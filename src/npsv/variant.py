import os, subprocess, tempfile
import vcf
import pysam
import pysam.bcftools as bcftools


def is_precise(record: vcf.model._Record):
    """Return true if SV has precise breakpoints"""
    # PyVCF approach for is_sv_precise is not flexible enough
    assert record.is_sv
    return not (
        record.INFO.get("IMPRECISE") is not None
        or record.INFO.get("CIPOS") is not None
        or record.INFO.get("CIEND") is not None
    )


def get_ci(record: vcf.model._Record, key: str, default_ci: int):
    """Get SV confidence interval or default for VCF record
    
    Arguments:
        record {vcf.model._Record} -- VCF SV record
        key {str} -- CI INFO key
        default_ci {int} -- Default value for CI if missing
    
    Returns:
        list -- Confidence interval
    """
    try:
        return [int(val) for val in record.INFO[key]]
    except KeyError:
        if is_precise(record):
            return [0, 0]
        else:
            return [-default_ci, default_ci]


def variant_descriptor(record: vcf.model._Record):
    return "{}_{}_{}_{}".format(
        record.CHROM, record.POS, record.sv_end, record.var_subtype
    )


def write_record_to_indexed_vcf(record, vcf_reader, vcf_path):
    vcf_writer = vcf.Writer(open(vcf_path, "w"), vcf_reader)
    vcf_writer.write_record(record)
    vcf_writer.close()
    # pylint: disable=no-member
    return pysam.tabix_index(vcf_path, preset="vcf", force=True)

class Variant(object):
    def __init__(self, record):
        self.record = record
        # Make sure every variant has a unique ID
        if self.record.ID is None:
            self.record.ID = variant_descriptor(record)
        
    @property
    def chrom(self):
        return self.record.CHROM

    @property
    def pos(self):
        return self.record.POS

    @property
    def id(self):
        return self.record.ID

    @property
    def end(self):
        """1-indexed closed end"""
        return int(self.record.sv_end)


    @property
    def is_precise(self):
        """Return true if SV has precise breakpoints"""
        # PyVCF approach for is_sv_precise is not flexible enough
        assert self.record.is_sv
        return not (
            self.record.INFO.get("IMPRECISE") is not None
            or self.record.INFO.get("CIPOS") is not None
            or self.record.INFO.get("CIEND") is not None
        )

    @property
    def is_deletion(self):
        return False

    def get_ci(self, key: str, default_ci: int):
        """Get SV confidence interval or default for VCF record
        
        Arguments:
            key {str} -- CI INFO key
            default_ci {int} -- Default value for CI if missing
        
        Returns:
            list -- Confidence interval
        """
        try:
            return [int(val) for val in self.record.INFO[key]]
        except KeyError:
            if self.is_precise:
                return [0, 0]
            else:
                return [-default_ci, default_ci]

    def region_string(self, flank=0):
        """Return 1-indexed fully closed region"""
        # In correctly formatted VCF, POS is first base of event when zero-indexed, while
        # END is 1-indexed closed end or 0-indexed half-open end
        return f"{self.record.CHROM}:{self.record.POS+1-flank}-{self.end+flank}"

    def right_flank_region_string(self, right_flank, left_flank=0):
        """Return 1-indexed fully closed region"""
        return f"{self.record.CHROM}:{self.end+1+left_flank}-{self.end+right_flank}"

    def reference_sequence(self, args, region=None, flank=0):
        try:
            if region is None:
                region = self.region_string(flank)

            # pylint: disable=no-member
            with pysam.FastaFile(args.reference) as ref_fasta:
                ref_seq = ref_fasta.fetch(region=region)
                return ref_seq
        except ValueError as err:
            print(region)
            raise err


    def to_minimal_vcf(self, args):
        raise NotImplementedError()

    def synth_fasta(self, args):
        # TODO: Normalize ref and alt contig lengths
        region = self.region_string(args.flank)
        ref_seq = self.reference_sequence(args, region=region)
        
        # Construct alternate alleles with bcftools consensus
        try:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", dir=args.tempdir
            ) as ref_allele_fasta:
                vcf_path = self.to_minimal_vcf(args)
                print(
                    f">{region}", ref_seq, sep="\n", file=ref_allele_fasta, flush=True
                )
                # pylint: disable=no-member
                alt_seq = bcftools.consensus("-f", ref_allele_fasta.name, vcf_path)
        finally:
            os.remove(vcf_path)
            os.remove(vcf_path + ".tbi")

        ref_contig = region.replace(":", "_").replace("-", "_")
        alt_contig = ref_contig + "_alt"
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta", dir=args.tempdir
        ) as allele_fasta:
            print(">", ref_contig, sep="", file=allele_fasta)
            print(ref_seq, file=allele_fasta)
            print(">", alt_contig, sep="", file=allele_fasta)
            allele_fasta.write(alt_seq[alt_seq.find("\n") + 1 :])

        return allele_fasta.name, ref_contig, alt_contig

    

    @classmethod
    def from_pyvcf(cls, record):
        if not record.is_sv:
            return None

        kind = record.var_subtype
        if kind.startswith("DEL"):
            return DeletionVariant(record)


class DeletionVariant(Variant):
    def __init__(self, record):
        Variant.__init__(self, record)

    @property
    def is_deletion(self):
        return True

    @property
    def event_length(self):
        # In correctly formatted VCF, POS is first base of event when zero-indexed, while
        # END is 1-indexed closed end or 0-indexed half-open end
        return int(self.record.sv_end) - self.record.POS

    @property
    def alt_length(self):
        assert len(self.record.ALT) == 1, "Multiple alternates are not supported"
        allele = self.record.ALT[0]
        if isinstance(allele, vcf.model._SV):
            # Symbolic allele
            return 1
        else:
            return len(allele)

    @property
    def as_bedtools_interval(self):
        # In correctly formatted VCF, POS is first base of event when zero-indexed, while
        # END is 1-indexed closed end or 0-indexed half-open end
        import pybedtools as bedtools
        return bedtools.Interval(self.record.CHROM, self.record.POS, self.record.sv_end)

    def to_minimal_vcf_record(self, info={}):
        # Mixin additional info fields
        base_info = { "SVTYPE": "DEL", "END": self.end, "SVLEN": -self.event_length }
        base_info.update(info)
        line = "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t.\t.\t{info_string}".format(
            chrom=self.record.CHROM,
            pos=self.record.POS,
            id=self.record.ID or ".",
            ref=self.record.REF,
            alt=self.record.ALT[0],
            info_string=";".join((f"{k}={v}" for (k, v) in base_info.items()))
        )
        return line

    def to_minimal_vcf(self, args, tempdir=None):
        # Create minimal VCF if needed for use with bcftools
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False, dir=(tempdir or args.tempdir)
        ) as vcf_file:
            print(
                """\
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">    
##ALT=<ID=DEL,Description="Deletion">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO""",
                file=vcf_file,
            )
            record = self.to_minimal_vcf_record()
            print(record, file=vcf_file)

        # Unfortunately tabix_index leaks file handles (in is_gzip_file function it appears)
        subprocess.check_call(
            "bgzip {0} && tabix {0}.gz".format(vcf_file.name), shell=True
        )
        return vcf_file.name + ".gz"

    