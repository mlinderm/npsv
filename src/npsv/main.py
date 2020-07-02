#!/usr/bin/env python3
import argparse, io, logging, os, subprocess, shutil, sys, tempfile
from tqdm import tqdm
from shlex import quote
import pysam
import vcf
from npsv.npsv_options import *
from npsv.variant import variant_descriptor, write_record_to_indexed_vcf
from npsv.feature_extraction import (
    extract,
    header,
    Variant,
)
from npsv.random_variants import random_variants
from npsv.genotyper import genotype_vcf
from npsv.sample import Sample
import ray
from retry import retry


def make_argument_parser():
    parser = argparse.ArgumentParser(
        "npsv", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # "Core" options
    parser.add_argument("-i", "--input", help="Input VCF file", type=str, required=True)
    parser.add_argument("-b", "--bam", help="Input BAM file", type=str, required=True)
    parser.add_argument(
        "-o", "--output", help="Output directory", type=str, required=True
    )
    parser.add_argument(
        "--prefix",
        help="Prefix for feature files and genotypes",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-r",
        "--reference-sequence",
        help="Reference fasta file",
        type=str,
        dest="reference",
        required=True,
    )
    parser.add_argument(
        "--genome", help="BedTools Genome file", type=str, required=True
    )
    parser.add_argument(
        "-t",
        "--tempdir",
        help="Specify the temp directory",
        type=str,
        default=tempfile.gettempdir(),
    )
    parser.add_argument("--threads", help="Number of threads", type=int, default=1)
    logging_options = parser.add_mutually_exclusive_group()
    logging_options.add_argument(
        "-d",
        "--debug",
        help="Debug logging",
        action="store_const",
        dest="loglevel",
        const=logging.DEBUG,
        default=logging.WARNING,
    )
    logging_options.add_argument(
        "-v",
        "--verbose",
        help="Verbose logging",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
    )

    # Synthesizer options
    synth_options = parser.add_argument_group()
    synth_options.add_argument(
        "--n", help="Number of replicates", type=int, default=100
    )
    synth_options.add_argument(
        "--sim-ref",
        help="Simulate AC=0 genotype instead of randomly sample",
        dest="sim_ref",
        action="store_true",
        default=False,
    )
    synth_options.add_argument(
        "--reuse",
        help="Use existing synthetic BAM files if available in output directory",
        action="store_true",
        default=False,
    )
    synth_options.add_argument(
        "--profile", help="ART profile", type=str, default="HS25"
    )

    # Options used by "sub tools"
    data_options = parser.add_argument_group()
    add_data_options(data_options)

    feature_options = parser.add_argument_group()
    add_feature_options(feature_options)

    random_options = parser.add_argument_group()
    add_random_options(random_options)

    genotyping_options = parser.add_argument_group()
    add_genotyping_options(genotyping_options)

    return parser


def art_read_length(read_length, profile):
    """Make sure read length is compatible ART"""
    if profile in ("HS10", "HS20"):
        return min(read_length, 100)
    elif profile in ("HS25", "HSXn", "HSXt"):
        return min(read_length, 150)
    else:
        return read_length


def check_if_bwa_index_loaded(reference: str) -> bool:
    """Check if bwa index is loaded in shared memory
    
    Args:
        reference (str): Path to reference file
    
    Returns:
        bool: True if reference is already loaded
    """
    shared_name = os.path.basename(reference)
    indices = subprocess.check_output(
        "bwa shm -l", shell=True, universal_newlines=True, stderr=subprocess.DEVNULL
    )
    for index in indices.split("\n"):
        if index.startswith(shared_name):
            return True
    return False


def ray_iterator(obj_ids):
    while obj_ids:
        done, obj_ids = ray.wait(obj_ids)
        yield ray.get(done[0])


@ray.remote
def simulate_deletion(args, sample, record, variant_vcf_path, description):
    # logging.info(f"Extracting features for {description}")
    out_file = open(os.path.join(args.output, description + ".sim.tsv"), "w")

    # In correctly formatted VCF, POS is first base of event when zero-indexed, while
    # END is 1-indexed closed end or 0-indexed half-open end
    pos = record.POS
    end = int(record.sv_end)
    event_length = end - pos

    if args.sim_ref:
        # Generate the null model via simulation
        simulated_ACs = [0, 1, 2]
    else:
        # Generate random variants as the null model
        random_vcf_path = os.path.join(args.output, description + ".random.vcf")
        with open(random_vcf_path, "w") as random_vcf_file:
            random_variants(
                args.reference,
                args.genome,
                args.gaps,
                random_vcf_file,
                size=event_length,
                n=args.n,
                use_X=args.use_X,
                only_sex=args.only_sex,
                flank=args.flank,
            )

        # Extract features from the random data
        random_extract_args = argparse.Namespace(**vars(args))
        setattr(random_extract_args, "header", False)
        setattr(random_extract_args, "threads", 1)
        extract(
            random_extract_args,
            random_vcf_path,
            args.bam,
            sample=sample,
            ac=0,
            out_file=out_file,
            force_chrom=record.CHROM,
            force_pos=pos,
            force_end=end,
            insert_hist=True,
        )

        simulated_ACs = [1, 2]

    # Extract features from synthetic data.
    synth_extract_args = argparse.Namespace(**vars(args))
    setattr(synth_extract_args, "header", False)
    setattr(synth_extract_args, "threads", 1)

    # Generate FASTA path once (can be reused by realigner for every replicate)
    variant = Variant.from_pyvcf(record)
    fasta_path, ref_contig, alt_contig = variant.synth_fasta(args)

    for z in simulated_ACs:
        synthetic_bam_path = os.path.join(args.output, f"{description}_{z}.bam")
        if not args.reuse or not os.path.exists(synthetic_bam_path):
            # Check if shared reference is available
            shared_ref_arg = ""
            if check_if_bwa_index_loaded(args.reference):
                shared_ref_arg = f"-S {quote(os.path.basename(args.reference))}"

            # Synthetic data generation seems to fail randomly for some variants
            @retry(tries=2)
            def gen_synth_bam():
                # Generate synthetic BAM file
                synth_commandline = f"synthBAM \
                    -t {quote(args.tempdir)} \
                    -R {quote(args.reference)} \
                    {shared_ref_arg} \
                    -g {quote(args.genome)} \
                    -c {sample.mean_coverage / 2:0.1f} \
                    -m {sample.mean_insert_size} \
                    -s {sample.std_insert_size} \
                    -l {art_read_length(sample.read_length, args.profile)} \
                    -p {args.profile} \
                    -f {args.flank} \
                    -i {args.n} \
                    -z {z} \
                    -n {sample.name} \
                    {variant_vcf_path} {synthetic_bam_path}"
                synth_result = subprocess.run(
                    synth_commandline, shell=True, stderr=subprocess.PIPE,
                )
                if synth_result.returncode != 0 or not os.path.exists(synthetic_bam_path):
                    raise RuntimeError(f"Synthesis script failed to generate {synthetic_bam_path}")
            
            gen_synth_bam()

        for i in range(1, args.n + 1):
            try:
                single_sample_bam_file = tempfile.NamedTemporaryFile(
                    mode="wb", delete=False, suffix=".bam", dir=args.tempdir
                )
                single_sample_bam_file.close()

                pysam.view(  # pylint: disable=no-member
                    synthetic_bam_path,
                    "-b",
                    "-h",
                    "-r",
                    f"synth{i}",
                    "-o",
                    single_sample_bam_file.name,
                    catch_stdout=False,
                )
                pysam.index(  # pylint: disable=no-member
                    single_sample_bam_file.name, "-b"
                )

                extract(
                    synth_extract_args,
                    variant_vcf_path,
                    single_sample_bam_file.name,
                    ac=z,
                    sample=sample,
                    out_file=out_file,
                    input_fasta=fasta_path,
                    ref_contig=ref_contig,
                    alt_contig=alt_contig,
                    insert_hist=args.insert_hist,
                )
            finally:
                os.remove(single_sample_bam_file.name)
                os.remove(single_sample_bam_file.name + ".bai")

    out_file.close()
    return out_file.name


def main():
    parser = make_argument_parser()
    args = parser.parse_args()

    logging.basicConfig(level=args.loglevel)

    logging.info(
        f"Creating {args.output} output and {args.tempdir} temporary directories if they don't exist"
    )
    os.makedirs(args.output, exist_ok=True)
    os.makedirs(args.tempdir, exist_ok=True)

    # Initialize parallel computing setup
    ray.init(num_cpus=args.threads, temp_dir=args.tempdir, include_webui=False)

    # TODO: If library is not specified compute statistics, i.e. mean insert size, tec.
    if args.stats_path is not None:
        logging.info("Extracting BAM stats from NPSV stats file")
        sample = Sample.from_npsv(args.stats_path, bam_path=args.bam)
    elif None not in (
        args.fragment_mean,
        args.fragment_sd,
        args.read_length,
        args.depth,
    ):
        logging.info("Using Normal distribution for BAM stats")
        sample = Sample.from_distribution(
            args.bam,
            args.fragment_mean,
            args.fragment_sd,
            args.read_length,
            mean_coverage=args.depth,
        )
    else:
        raise parser.error(
            "Library information needed. Either provide distribution parameters or run `npsvg preprocess` to generate stats file."
        )

    # For each variant generate synthetic bam file(s) and extract relevant evidence
    record_results = []
    vcf_reader = vcf.Reader(filename=args.input)
    for record in tqdm(vcf_reader, desc="Preparing variants"):
        # npsv currently only supports deletions
        if not record.is_sv or record.var_subtype != "DEL":
            continue

        description = variant_descriptor(record)
        if record.ID is None:
            # Update ID to be more meaningful
            record.ID = description

        # Construct single variant VCF outside of worker so we don't need to pass the reader into the thread
        variant_vcf_path = os.path.join(args.output, description + ".vcf")
        if not args.reuse or not os.path.exists(variant_vcf_path + ".gz"):
            variant_vcf_path = write_record_to_indexed_vcf(
                record, vcf_reader, variant_vcf_path
            )
        else:
            # Variant file already exists, no need to recreate
            variant_vcf_path += ".gz"

        record_results.append(
            simulate_deletion.remote(
                args, sample, record, variant_vcf_path, description
            )
        )

    # Concatenate output files to create final result
    sim_tsv_path = os.path.join(args.output, args.prefix + ".sim.tsv")
    logging.info("Concatenating results in %s", sim_tsv_path)
    with open(sim_tsv_path, "w") as file:
        header(out_file=file, ac=True)
    with open(sim_tsv_path, "ab") as sink:
        for result in tqdm(
            ray_iterator(record_results),
            total=len(record_results),
            desc="Simulating variants",
        ):
            with open(result, "rb") as source:
                shutil.copyfileobj(source, sink)

    # Finished generating features from simulated data

    # Extract features from "real" data
    with open(
        os.path.join(args.output, args.prefix + ".real.tsv"), "w"
    ) as real_tsv_file:
        logging.info(
            "Extracting features from 'real' data (output in %s)", real_tsv_file.name
        )
        real_args = argparse.Namespace(**vars(args))
        setattr(real_args, "header", True)
        extract(real_args, args.input, args.bam, out_file=real_tsv_file, sample=sample, insert_hist=True)

    # Perform genotyping
    with open(os.path.join(args.output, args.prefix + ".npsv.vcf"), "w") as gt_vcf_file:
        logging.info("Determining genotypes (output in %s)", gt_vcf_file.name)
        genotyping_args = argparse.Namespace(**vars(args))
        genotype_vcf(
            genotyping_args,
            args.input,
            sim_tsv_path,
            real_tsv_file.name,
            gt_vcf_file,
            samples=[sample.name],
        )


if __name__ == "__main__":
    main()

