import io, json, logging, math, os, re, subprocess, tempfile
from collections import Counter
from shlex import quote
from scipy.stats import norm
import numpy as np
import pandas as pd
import pysam
import pybedtools.bedtool as bed

NO_LIB_LABEL = "UNLABELED"
GENERIC_LIB_LABEL = "GENERIC"

GC_ROUNDING_DIGITS = 2


def sample_from_bam(bam):
    """Extract sample name from BAM"""
    read_groups = bam.header["RG"]
    samples = set([rg["SM"] for rg in read_groups])
    assert len(samples) == 1, "BAM file {} must contain a single sample".format(
        bam.filename
    )
    sample, = samples  # Extract single value from set
    return sample


def stats_from_dict(dict):
    lengths = np.array(dict.keys())
    counts = np.array(dict.values())
    num_reads = np.sum(counts)

    mean = np.sum(lengths * counts) / num_reads
    std = math.sqrt(np.sum(counts * np.square(lengths - mean)) / (num_reads - 1))

    return mean, std


class ParametricFragmentDensity:
    def __init__(self, mean, sd):
        self.mean = mean
        self.sd = sd

    def __getitem__(self, size):
        return norm.pdf(size, loc=self.mean, scale=self.sd)


class SparseFragmentDensity:
    def __init__(self, mean, sd, dens):
        self.mean = mean
        self.sd = sd
        self.dens = dens

    def __getitem__(self, size):
        density = self.dens.get(size)
        if density is None:
            density = norm.pdf(size, loc=self.mean, scale=self.sd)
        return density


class Library:
    """
    Implements the SVTyper Library interface
    """

    def __init__(self, mean, sd, read_length, dens):
        self.mean = mean
        self.sd = sd
        self.read_length = read_length
        self.dens = dens

    @property
    def mean_insert_size(self):
        return self.mean

    @property
    def std_insert_size(self):
        return self.sd

    @property
    def insert_size_density(self):
        return self.dens

    def gc_normalized_coverage(self, gc_fraction):
        return 1.0

    def chrom_normalized_coverage(self, chrom):
        return 1.0


class DistributionLibrary(Library):
    def __init__(self, mean, sd, read_length):
        Library.__init__(
            self, mean, sd, read_length, ParametricFragmentDensity(mean, sd)
        )


class HistogramLibrary(Library):
    def __init__(self, mean, sd, read_length, fragment_hist):
        # Convert histogram to density
        total_fragments = sum(fragment_hist.values())
        fragment_density = {}
        for length, count in fragment_hist.items():
            fragment_density[int(length)] = float(count) / total_fragments

        Library.__init__(
            self,
            mean,
            sd,
            read_length,
            SparseFragmentDensity(mean, sd, fragment_density),
        )


class NPSVStatsLibrary(DistributionLibrary):
    def __init__(
        self,
        mean,
        sd,
        read_length,
        gc_normalized_coverage={},
        chrom_normalized_coverage={},
    ):
        DistributionLibrary.__init__(self, mean, sd, read_length)
        self.gc_normalization_factors = gc_normalized_coverage
        self.chrom_normalization_factors = chrom_normalized_coverage

    def gc_normalized_coverage(self, gc_fraction):
        return self.gc_normalization_factors.get(
            round(gc_fraction, GC_ROUNDING_DIGITS), 1.0
        )

    def chrom_normalized_coverage(self, chrom):
        return self.chrom_normalization_factors.get(chrom, 1.0)


class Sample(object):
    def __init__(self, name, libraries, rg_to_lib, generic_lib, mean_coverage=None):
        self.name = name
        self.libraries = libraries
        self.rg_to_lib = rg_to_lib
        self.generic_lib = generic_lib
        self.mean_coverage = mean_coverage

    @classmethod
    def from_svtyper(cls, bam, json_path):
        """
        Construct libraries from SVTyper 'lib_info' JSON file
        :param bam: PySam AlignmentFile
        :param json_path: Path to SVTyper 'lib_info' JSON file
        :return: Library object
        """
        sample_name = sample_from_bam(bam)

        with open(json_path, "r") as file:
            lib_info = json.load(file)

            generic_hist = Counter()
            generic_read_length = []

            library_dict = {}
            rg_to_lib = {}

            libraries = lib_info[sample_name]["libraryArray"]
            for library in libraries:
                library_name = library.get("library_name", NO_LIB_LABEL)
                if library_name == "":
                    library_name = NO_LIB_LABEL

                library_obj = HistogramLibrary(
                    library["mean"],
                    library["sd"],
                    library["read_length"],
                    library["histogram"],
                )

                library_dict[library_name] = library_obj
                for read_group in library["readgroups"]:
                    rg_to_lib[read_group] = library_obj

                # Add entries to "generic" library
                for length, count in library["histogram"].iteritems():
                    generic_hist[int(length)] += int(count)
                generic_read_length.append(library["read_length"])

            # Construct generic library
            generic_mean, generic_std = stats_from_dict(generic_hist)
            generic_lib = HistogramLibrary(
                generic_mean,
                generic_std,
                int(np.mean(generic_read_length)),
                generic_hist,
            )

            return cls(sample_name, library_dict, rg_to_lib, generic_lib)

    @classmethod
    def from_distribution(
        cls,
        bam_path: str,
        mean: float,
        sd: float,
        read_length: int,
        mean_coverage: float = None,
    ) -> "Sample":
        """Create Sample object from distribution parameters
        
        Args:
            bam_path (str): Path to BAM file
            mean (float): Mean insert size
            sd (float): Standard deviation of the insert size
            read_length (int): Read length
            mean_coverage (float, optional): Sample mean coverage. Defaults to None.
        
        Returns:
            Sample: Sample object with library stats
        """
        with pysam.AlignmentFile(  # pylint: disable=no-member
            bam_path, mode="rb"
        ) as bam_reader:
            sample_name = sample_from_bam(bam_reader)
            library_object = DistributionLibrary(mean, sd, read_length)

            library_dict = {}
            rg_to_lib = {}
            for readgroup in bam_reader.header["RG"]:
                library_dict[readgroup.get("LB", NO_LIB_LABEL)] = library_object
                rg_to_lib[readgroup["ID"]] = library_object
            return cls(
                sample_name,
                library_dict,
                rg_to_lib,
                library_object,
                mean_coverage=mean_coverage,
            )

    @classmethod
    def from_npsv(
        cls, json_path: str, bam_path: str = None, min_gc_bin=100, max_gc_error=0.01
    ) -> "Sample":
        """Create Sample object from NPSV BAM stats JSON file
        
        Args:
            json_path (str): Path to JSON file
            bam_path (str, optional): Path to BAM file to override BAM path in JSON file. Defaults to None.
            min_obs_gc_bin (int, optional): Minimum observation of GC fraction to compute normalization factor. Defaults to 100.
        
        Returns:
            Sample: Sample object with library stats
        """
        with open(json_path, "r") as file:
            bam_info = json.load(file)

            # Filter GC entries with limited data
            gc_normalized_coverage = {}
            for gc, norm_covg in bam_info["gc_normalized_coverage"].items():
                if (
                    bam_info["gc_bin_count"].get(gc, 0) >= min_gc_bin
                    and bam_info["gc_normalized_coverage_error"].get(gc, 0)
                    <= max_gc_error
                ):
                    gc_normalized_coverage[float(gc)] = norm_covg

            library_object = NPSVStatsLibrary(
                bam_info["mean_insert_size"],
                bam_info["std_insert_size"],
                bam_info["read_length"],
                gc_normalized_coverage,
                bam_info.get("chrom_normalized_coverage", {}),
            )

            library_dict = {}
            rg_to_lib = {}

            # pylint: disable=no-member
            with pysam.AlignmentFile(bam_path or bam_info["bam"], "rb") as bam_file:
                for readgroup in bam_file.header["RG"]:
                    library_dict[readgroup.get("LB", NO_LIB_LABEL)] = library_object
                    rg_to_lib[readgroup["ID"]] = library_object

            return cls(
                bam_info["sample"],
                library_dict,
                rg_to_lib,
                library_object,
                mean_coverage=bam_info["mean_coverage"],
            )

    def has_read_group(self, read_group):
        return read_group in self.rg_to_lib

    def get_lib(self, read_group: str) -> Library:
        """Return Library object for specific read_group or generic library is read_group doesn't exist
        
        Args:
            read_group (str): Read group ID
        
        Returns:
            Library: Library object
        """
        return self.rg_to_lib.get(read_group, self.generic_lib)

    @property
    def read_length(self, read_group=None):
        return self.get_lib(read_group).read_length

    @property
    def mean_insert_size(self, read_group=None):
        return self.get_lib(read_group).mean

    @property
    def std_insert_size(self, read_group=None):
        return self.get_lib(read_group).sd

    def gc_mean_coverage(self, gc_fraction: float, read_group: str = None) -> float:
        """Return mean coverage for bins with this GC fraction
        
        If there is insufficient data for that GC fraction it will return the original mean coverage

        Args:
            gc_fraction (float): GC fraction in range 0-1
            read_group (str, optional): Read group ID. Defaults to None.
        
        Returns:
            float: Mean coverage
        """
        return (
            self.get_lib(read_group).gc_normalized_coverage(gc_fraction)
            * self.mean_coverage
        )

    def chrom_mean_coverage(self, chrom: str, read_group: str = None) -> float:
        """Return mean coverage for specific chromosome
        
        Args:
            chrom (str): Chromosome
            read_group (str, optional): Read group ID. Defaults to None.
        
        Returns:
            float: Mean coverage
        """
        return (
            self.get_lib(read_group).chrom_normalized_coverage(chrom)
            * self.mean_coverage
        )


def samtools_norm_coverage_group(table, mean_coverage):
    """Compute normalized coverage weighted by alignable bases"""
    return (np.sum(table.bases) / np.sum(table.align_len)) / mean_coverage


def compute_coverage_with_samtools(
    args, input_bam, mean_coverage, gc_window_size=20000
):
    # Generate GC and chromosome normalized coverage for the entire BAM file
    # Based on https://www.biostars.org/p/92744/
    logging.info(
        "Computing chromosome and GC normalized coverage for %s with samtools and bedtools",
        input_bam,
    )

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", dir=args.tempdir
    ) as window_bed_file:
        # pylint: disable=unexpected-keyword-arg
        windows_bed = bed.BedTool().window_maker(
            g=args.genome, w=gc_window_size, output=window_bed_file.name
        )

        windows_table = (
            # pylint: disable=no-member,unexpected-keyword-arg
            bed.BedTool(pysam.bedcov(window_bed_file.name, input_bam), from_string=True)
            .nucleotide_content(fi=args.reference)
            .to_dataframe(
                index_col=False,
                header=0,
                usecols=[0, 1, 2, 3, 7, 8, 10, 11, 12],
                names=[
                    "chrom",
                    "start",
                    "end",
                    "bases",
                    "num_C",
                    "num_G",
                    "num_N",
                    "num_oth",
                    "seq_len",
                ],
                dtype={"chrom": str},
            )
        )

    # Remove windows with no alignable data
    windows_table["align_len"] = (
        windows_table.seq_len - windows_table.num_N - windows_table.num_oth
    )
    windows_table = windows_table[windows_table.align_len != 0]

    # Compute normalized coverage by chromosome
    norm_coverage_by_chrom = (
        windows_table.groupby("chrom")
        .apply(samtools_norm_coverage_group, mean_coverage)
        .to_dict()
    )

    # Compute normalized coverage by GC bin
    gc_fraction = np.round(
        (windows_table.num_G + windows_table.num_C) / windows_table.align_len, 2
    )
    norm_coverage = (windows_table.bases / windows_table.align_len) / mean_coverage
    norm_coverage_by_gc = (
        norm_coverage.groupby(gc_fraction).agg(["count", "mean"]).to_dict()
    )

    return norm_coverage_by_chrom, norm_coverage_by_gc


def goleft_norm_coverage_group(table):
    """Compute normalized coverage weighted by alignable bases"""
    weights = table.align_len / np.sum(table.align_len)
    return np.sum(weights * table.norm_covg)


def compute_coverage_with_goleft(args, input_bam):
    # Generate GC and chromosome normalized coverage for the entire BAM file
    logging.info(
        "Computing chromosome and GC normalized coverage for %s with goleft and bedtools",
        input_bam,
    )
    with tempfile.TemporaryDirectory(dir=args.tempdir) as output_dir:
        prefix = os.path.basename(output_dir)

        commandline = "{exec} indexcov --directory {dir} --fai {fai} {bam}".format(
            exec=quote(args.goleft),
            dir=output_dir,
            fai=quote(args.reference + ".fai"),
            bam=quote(input_bam + ".crai" if input_bam.endswith("cram") else input_bam),
        )
        subprocess.check_call(
            commandline,
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # Compute chromosome and GC normalized coverage
        windows_table = (
            # pylint: disable=unexpected-keyword-arg
            bed.BedTool(fn=os.path.join(output_dir, prefix + "-indexcov.bed.gz"))
            .nucleotide_content(fi=args.reference)
            .to_dataframe(
                index_col=False,
                header=0,
                usecols=[0, 1, 2, 3, 7, 8, 10, 11, 12],
                names=[
                    "chrom",
                    "start",
                    "end",
                    "norm_covg",
                    "num_C",
                    "num_G",
                    "num_N",
                    "num_oth",
                    "seq_len",
                ],
                dtype={"chrom": str},
            )
        )
        # Remove windows with no alignable data
        windows_table["align_len"] = (
            windows_table.seq_len - windows_table.num_N - windows_table.num_oth
        )
        windows_table = windows_table[windows_table.align_len != 0]

        norm_coverage_by_chrom = (
            windows_table.groupby("chrom").apply(goleft_norm_coverage_group).to_dict()
        )

        return norm_coverage_by_chrom


def compute_bam_stats(args, input_bam: str, gc_window_size=20000):
    """Compute bam file stats, e.g. mean insert size, etc.
    
    Args:
        args ([type]): Arguments
        input_bam (str): Path to BAM file
        gc_window_size (int, optional): Window size for computing GC content. Defaults to 20000.
    
    Returns:
        dict: Stats (suitable for writing as JSON file)
    """
    # Generate stats for the entire BAM file
    logging.info("Computing coverage and insert size statistics with goleft")
    commandline = "{exec} covstats --fasta {ref} {bam}".format(
        exec=quote(args.goleft), ref=quote(args.reference), bam=quote(input_bam)
    )
    covstats = subprocess.check_output(commandline, shell=True, universal_newlines=True)
    covstats_table = pd.read_csv(io.StringIO(covstats), sep="\t")

    covstats_record = covstats_table[covstats_table.bam == input_bam].to_dict("records")
    assert len(covstats_record) == 1, "Couldn't find stats for BAM file"
    covstats_record = covstats_record[0]

    if args.picard_insert is not None:
        logging.info("Using Picard insert size statistics")
        insert_table = pd.read_csv(args.picard_insert, sep="\t", comment="#", nrows=1)
        insert_dict = insert_table.to_dict("records")[0]
        covstats_record["template_mean"] = insert_dict["MEAN_INSERT_SIZE"]
        covstats_record["template_sd"] = insert_dict["STANDARD_DEVIATION"]
    
    if args.picard_wgs is not None:
        logging.info("Using Picard mean coverage")
        wgs_table = pd.read_csv(args.picard_wgs, sep="\t", comment="#", nrows=1)
        wgs_dict = wgs_table.to_dict("records")[0]
        covstats_record["coverage"] = wgs_dict["MEAN_COVERAGE"]

    mean_coverage = covstats_record["coverage"]

    if args.picard_gc is not None:
        # "Fast" path. Use goleft to compute chromosomal coverage and extract GC bias from Picard metrics
        norm_coverage_by_chrom = compute_coverage_with_goleft(args, input_bam)

        # Extract GC bias from Picard metrics output
        gc_bias_table = pd.read_csv(args.picard_gc, sep="\t", comment="#")
        gc_bias_table.GC = gc_bias_table.GC / 100
        gc_bias_dict = (
            gc_bias_table[gc_bias_table.READS_USED == "ALL"]
            .set_index("GC")[["NORMALIZED_COVERAGE", "WINDOWS", "ERROR_BAR_WIDTH"]]
            .to_dict()
        )

        norm_coverage_by_gc = {
            "mean": gc_bias_dict["NORMALIZED_COVERAGE"],
            "count": gc_bias_dict["WINDOWS"],
            "error": gc_bias_dict["ERROR_BAR_WIDTH"],
        }
    else:
        # "Slow" path. Compute chromosomal coverage and GC bias with bedtools/samtools
        norm_coverage_by_chrom, norm_coverage_by_gc = compute_coverage_with_samtools(
            args, input_bam, mean_coverage, gc_window_size=gc_window_size
        )

    # Construct stats dictionary that can be written to JSON
    stats = {
        "sample": covstats_record["sample"],
        "bam": input_bam,
        "read_length": covstats_record["read_length"],
        "mean_insert_size": covstats_record["template_mean"],
        "std_insert_size": covstats_record["template_sd"],
        "mean_coverage": mean_coverage,
        "chrom_normalized_coverage": norm_coverage_by_chrom,
        "gc_normalized_coverage": norm_coverage_by_gc["mean"],
        "gc_bin_count": norm_coverage_by_gc["count"],
        "gc_normalized_coverage_error": norm_coverage_by_gc.get("error", {}),
    }
    return stats