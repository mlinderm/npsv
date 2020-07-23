#!/usr/bin/env python3
import argparse, io, logging, os, subprocess, shutil, sys, tempfile
from tqdm import tqdm
from shlex import quote
import pysam
import vcf
from npsv.npsv_options import *
from npsv.variant import Variant, variant_descriptor, write_record_to_indexed_vcf
from npsv.feature_extraction import (
    extract,
    extract_variant_features,
    header,
)
from npsv.random_variants import random_deletions
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
    add_simulation_options(synth_options)

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
def simulate_and_extract(args, sample, variant, variant_vcf_path, description):
    extract_args = argparse.Namespace(**vars(args))
    setattr(extract_args, "header", False)
    setattr(extract_args, "threads", 1)

    sim_out_file = open(os.path.join(args.output, description + ".sim.tsv"), "w")
    real_out_file = open(os.path.join(args.output, description + ".real.tsv"), "w")

    # Extract features from randomly sampled variants as null model
    if args.sim_ref:
        # Generate the null model via simulation
        simulated_ACs = [0, 1, 2]
    else:
        simulated_ACs = [1, 2]

        # Generate random variants as the null model
        random_vcf_path = os.path.join(args.output, description + ".random.vcf")
        with open(random_vcf_path, "w") as random_vcf_file:
            random_deletions(
                args.reference,
                args.genome,
                args.gaps,
                random_vcf_file,
                size=variant.event_length,
                n=args.n,
                use_X=args.use_X,
                only_sex=args.only_sex,
                flank=args.flank,
            )

        # Extract features from the random data
        extract(
            extract_args,
            random_vcf_path,
            args.bam,
            sample=sample,
            ac=0,
            out_file=sim_out_file,
            force_variant=variant,
            insert_hist=True,
        )

    # Extract features from synthetic data.

    # Generate FASTA path once (can be reused by realigner for every replicate)
    fasta_path, ref_contig, alt_contig = variant.synth_fasta(args)

    # Check if shared reference is available
    shared_ref_arg = ""
    if check_if_bwa_index_loaded(args.reference):
        shared_ref_arg = f"-S {quote(os.path.basename(args.reference))}"

    stats_file_arg = ""
    hap_coverage = sample.mean_coverage / 2
    if args.covg_gc_bias and args.stats_path is not None:
        # Use the BAM stats to model GC bias in the simulation
        stats_file_arg = f"-j {quote(args.stats_path)}"
        hap_coverage *= sample.max_gc_normalized_coverage(limit=args.max_gc_norm_covg)

    for z in simulated_ACs:
        synthetic_bam_path = os.path.join(args.output, f"{description}_{z}.bam")
        if not args.reuse or not os.path.exists(synthetic_bam_path):
            # Synthetic data generation seems to fail randomly for some variants
            @retry(tries=2)
            def gen_synth_bam():
                # Generate synthetic BAM file
                synth_commandline = f"synthBAM \
                    -t {quote(args.tempdir)} \
                    -R {quote(args.reference)} \
                    {shared_ref_arg} \
                    -g {quote(args.genome)} \
                    -c {hap_coverage:0.1f} \
                    -m {sample.mean_insert_size} \
                    -s {sample.std_insert_size} \
                    -l {art_read_length(sample.read_length, args.profile)} \
                    -p {args.profile} \
                    -f {args.flank} \
                    -i {args.n} \
                    -z {z} \
                    -n {sample.name} \
                    {stats_file_arg} \
                    {variant_vcf_path} {synthetic_bam_path}"
                synth_result = subprocess.run(
                    synth_commandline, shell=True, stderr=subprocess.PIPE,
                )
                if synth_result.returncode != 0 or not os.path.exists(
                    synthetic_bam_path
                ):
                    raise RuntimeError(
                        f"Synthesis script failed to generate {synthetic_bam_path}"
                    )

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

                features = extract_variant_features(
                    extract_args,
                    variant,
                    single_sample_bam_file.name,
                    sample,
                    input_fasta=fasta_path,
                    ref_contig=ref_contig,
                    alt_contig=alt_contig,
                    insert_hist=args.insert_hist,
                )
                features.print_features(sim_out_file, ac=z)
            finally:
                os.remove(single_sample_bam_file.name)
                os.remove(single_sample_bam_file.name + ".bai")

    # Extract features from real data
    features = extract_variant_features(
        extract_args,
        variant,
        args.bam,
        sample,
        input_fasta=fasta_path,
        ref_contig=ref_contig,
        alt_contig=alt_contig,
        insert_hist=True,
    )
    features.print_features(real_out_file)

    sim_out_file.close()
    real_out_file.close()
    return sim_out_file.name, real_out_file.name


def main():
    parser = make_argument_parser()
    args = parser.parse_args()

    logging.basicConfig(level=args.loglevel)

    # Create any directories that are needed
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
        variant = Variant.from_pyvcf(record, args.reference)
        # npsv currently only supports deletions
        if variant is None:
            continue

        description = variant_descriptor(record)

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
            simulate_and_extract.remote(
                args, sample, variant, variant_vcf_path, description
            )
        )

    # Concatenate output files to create feature files
    sim_tsv_path = os.path.join(args.output, args.prefix + ".sim.tsv")
    real_tsv_path = os.path.join(args.output, args.prefix + ".real.tsv")
    logging.info("Extracting features (to %s and %s)", sim_tsv_path, real_tsv_path)

    with open(sim_tsv_path, "w") as file:
        header(out_file=file, ac=True)
    with open(real_tsv_path, "w") as file:
        header(out_file=file, ac=True)

    with open(sim_tsv_path, "ab") as sim_sink, open(real_tsv_path, "ab") as real_sink:
        for sim_result, real_result in tqdm(
            ray_iterator(record_results),
            total=len(record_results),
            desc="Extracting features",
        ):
            with open(sim_result, "rb") as source:
                shutil.copyfileobj(source, sim_sink)
            with open(real_result, "rb") as source:
                shutil.copyfileobj(source, real_sink)

    # Perform genotyping
    with open(os.path.join(args.output, args.prefix + ".npsv.vcf"), "w") as gt_vcf_file:
        logging.info("Determining genotypes (output in %s)", gt_vcf_file.name)
        genotyping_args = argparse.Namespace(**vars(args))
        genotype_vcf(
            genotyping_args,
            args.input,
            sim_tsv_path,
            real_tsv_path,
            gt_vcf_file,
            samples=[sample.name],
        )


if __name__ == "__main__":
    main()

