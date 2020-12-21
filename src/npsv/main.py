#!/usr/bin/env python3
import argparse, io, logging, os, re, subprocess, shutil, sys, tempfile
from tqdm import tqdm
from shlex import quote
import pysam
import vcf
from npsv.npsv_options import *
from npsv.variant import Variant, variant_descriptor, write_record_to_indexed_vcf
from npsv.feature_extraction import (
    extract,
    extract_features,
    coverage_over_region,
    Features
)
from npsv.random_variants import random_variants, CHROM_REGEX_SEX
from npsv.genotyper import genotype_vcf
from npsv.sample import Sample
from npsv.simulation import GNOMAD_MULT_COVG
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
        help="Prefix for output feature files and genotypes",
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
    for kind in ("DEL","INS"):
        synth_options.add_argument(
            f"--{kind}-n", dest=f"{kind}_n", help="Number of replicates", type=int, default=GENOTYPING_DEFAULTS[kind].get("n", 100),
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
    synth_options.add_argument(
        "--covg-gc-bias",
        dest="covg_gc_bias",
        help="Model GC bias in simulated coverage",
        action="store_true",
        default=False,
    )
    synth_options.add_argument(
        "--sim-depth",
        dest="sim_depth",
        help="Source of depth for simulated data",
        default="global",
        choices=("global", "chrom", "flank"),
    )
    synth_options.add_argument(
        "--max-flank-norm-covg",
        dest="max_flank_norm_covg",
        help="Max flank normalized coverage",
        type=float,
        default=2.0,
    )
    synth_options.add_argument(
        "--keep-synth-bams",
        dest="keep_synth_bams",
        help="Keep synthetic data files",
        action="store_true",
        default=False,
    )
    add_simulation_options(synth_options)
    add_hybrid_options(synth_options)

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

    # Create directory for synthetic files
    if args.keep_synth_bams:
        synth_dir = args.output
    else:
        synth_dir = tempfile.mkdtemp(dir=args.tempdir)
        setattr(extract_args, "tempdir", synth_dir)

    # If using hybrid mode, only simulate the number of replicates needed for single model
    # for variants larger than threshold
    if getattr(args, f"{variant.subtype}_gt_mode") == "hybrid" and variant.event_length >= getattr(args, f"{variant.subtype}_hybrid_threshold", GENOTYPING_DEFAULTS["ANY"]["hybrid_threshold"]):
        replicates = args.downsample
    else:
        replicates = getattr(args, f"{variant.subtype}_n", GENOTYPING_DEFAULTS["ANY"]["n"])

    sim_out_file = open(os.path.join(args.output, description + ".sim.tsv"), "w")
    real_out_file = open(os.path.join(args.output, description + ".real.tsv"), "w")

    # Extract features from randomly sampled variants as null model
    # TODO: Integrate duplications into random variant generation
    if args.sim_ref or variant.is_duplication or variant.event_length == 0:
        # Generate the null model via simulation
        simulated_ACs = [0, 1, 2]
    else:
        simulated_ACs = [1, 2]

        use_X = only_sex = False
        if sample.gender == 1 and re.match(CHROM_REGEX_SEX, variant.chrom):
            # For 46XY we only want to sample from haploid sex chromosomes for X, Y variants
            only_sex = True
        elif sample.gender == 2:
            # For 46XX we can include diploid X along with autosome
            use_X = True

        # Generate random variants as the null model
        random_vcf_path = os.path.join(synth_dir, description + ".random.vcf")
        with open(random_vcf_path, "w") as random_vcf_file:
            random_variants(
                variant,
                args.reference,
                args.genome,
                args.gaps,
                random_vcf_file,
                n=replicates,
                use_X=use_X,
                only_sex=only_sex,
                flank=args.flank,
            )

        # Extract features from the random data
        try:
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
        except Exception as e:
            print(description)
            raise


    # Extract features from synthetic data.

    # Generate FASTA path once (can be reused by realigner for every replicate)
    fasta_path, ref_contig, alt_contig = variant.synth_fasta(extract_args)

    # Check if shared reference is available
    shared_ref_arg = ""
    if check_if_bwa_index_loaded(extract_args.reference):
        shared_ref_arg = f"-S {quote(os.path.basename(extract_args.reference))}"

    if args.sim_depth == "global":
        hap_coverage = sample.mean_coverage / 2
    elif args.sim_depth == "chrom":
        hap_coverage = sample.chrom_mean_coverage(variant.chrom) / 2
    elif args.sim_depth == "flank":
        # Determine coverage from flanking regions around event
        _, left_flank_bases, left_flank_length = coverage_over_region(
            extract_args.bam,
            variant.left_flank_region_string(args.flank),
            min_mapq=extract_args.min_mapq, min_baseq=extract_args.min_baseq, min_anchor=extract_args.min_anchor
        )
        _, right_flank_bases, right_flank_length = coverage_over_region(
            extract_args.bam,
            variant.right_flank_region_string(args.flank),
            min_mapq=extract_args.min_mapq, min_baseq=extract_args.min_baseq, min_anchor=extract_args.min_anchor
        )
        total_flank_bases = left_flank_bases + right_flank_bases
        total_flank_length = left_flank_length + right_flank_length
        if total_flank_bases > 0 and total_flank_length > 0:
            flank_coverage = total_flank_bases / total_flank_length
            hap_coverage = min(flank_coverage, sample.mean_coverage * extract_args.max_flank_norm_covg) / 2
        else:
            hap_coverage = sample.mean_coverage / 2

    stats_file_arg = ""
    gnomad_covg_file_arg = ""
    if extract_args.gnomad_covg is not None:
        # Use gnomAD coverage file to model coverage bias in the simulation
        gnomad_covg_file_arg = f"-g {quote(extract_args.gnomad_covg)}"
        hap_coverage *= GNOMAD_MULT_COVG
    elif extract_args.covg_gc_bias and extract_args.stats_path is not None:
        # Use the BAM stats to model GC bias in the simulation
        stats_file_arg = f"-j {quote(extract_args.stats_path)}"
        hap_coverage *= sample.max_gc_normalized_coverage(limit=extract_args.max_gc_norm_covg)

    for z in simulated_ACs:
        synthetic_bam_path = os.path.join(synth_dir, f"{description}_{z}.bam")
        if not extract_args.reuse or not os.path.exists(synthetic_bam_path):
            # Synthetic data generation seems to fail randomly for some variants, so retry
            @retry(tries=2)
            def gen_synth_bam():
                # Generate synthetic BAM file
                synth_commandline = f"synthBAM \
                    -t {quote(extract_args.tempdir)} \
                    -R {quote(extract_args.reference)} \
                    {shared_ref_arg} \
                    -c {hap_coverage:0.1f} \
                    -m {sample.mean_insert_size} \
                    -s {sample.std_insert_size} \
                    -l {art_read_length(sample.read_length, extract_args.profile)} \
                    -p {extract_args.profile} \
                    -f {extract_args.flank} \
                    -i {replicates} \
                    -z {z} \
                    -n {sample.name} \
                    {stats_file_arg} \
                    {gnomad_covg_file_arg} \
                    {variant_vcf_path} \
                    {synthetic_bam_path}"
                synth_result = subprocess.run(
                    synth_commandline, shell=True, stderr=subprocess.PIPE,
                )
                if synth_result.returncode != 0 or not os.path.exists(
                    synthetic_bam_path
                ):
                    print(synth_result.stderr)
                    raise RuntimeError(
                        f"Synthesis script failed to generate {synthetic_bam_path}"
                    )

            gen_synth_bam()

        for i in range(1, replicates + 1):
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

                features = extract_features(
                    extract_args,
                    variant,
                    single_sample_bam_file.name,
                    sample,
                    input_fasta=fasta_path,
                    ref_contig=ref_contig,
                    alt_contig=alt_contig,
                    insert_hist=extract_args.insert_hist,
                )
                features.print_features(sim_out_file, ac=z)
            finally:
                os.remove(single_sample_bam_file.name)
                os.remove(single_sample_bam_file.name + ".bai")

    # Extract features from real data
    features = extract_features(
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

    # Cleanup temporary files created for this variant
    if synth_dir != args.output:
        shutil.rmtree(synth_dir)

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
    ray.init(num_cpus=args.threads, _temp_dir=args.tempdir, include_dashboard=False)

    # TODO: If library is not specified compute statistics, i.e. mean insert size, tec.
    if args.stats_path is not None:
        logging.info("Extracting BAM stats from NPSV stats file")
        sample = Sample.from_npsv(args.stats_path, bam_path=args.bam, ped_path=args.ped_path)
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
    observed_variants = {}
    record_results = []
    vcf_reader = vcf.Reader(filename=args.input)
    for i, record in enumerate(tqdm(vcf_reader, desc="Preparing variants")):
        variant = Variant.from_pyvcf(record, args.reference)
        # npsv currently only supports deletions
        if variant is None:
            continue

        # NPSV currently does not support variants with duplicate start and end coordinates
        description = variant_descriptor(record)
        if observed_variants.setdefault(description, i) != i:
            logging.warning("Skipping variant with duplicate description %s", description)
            continue

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
        Features.header(out_file=file, ac=True)
    with open(real_tsv_path, "w") as file:
        Features.header(out_file=file, ac=False)

    with open(sim_tsv_path, "ab") as sim_sink, open(real_tsv_path, "ab") as real_sink:
        for sim_result, real_result in tqdm(
            ray_iterator(record_results),
            total=len(record_results),
            desc="Extracting features",
        ):
            with open(sim_result, "rb") as source:
                shutil.copyfileobj(source, sim_sink)
            sim_sink.flush()
            with open(real_result, "rb") as source:
                shutil.copyfileobj(source, real_sink)
            real_sink.flush()

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

