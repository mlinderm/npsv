import logging, os, subprocess, tempfile, textwrap
from backports.cached_property import cached_property
import vcf
import pysam


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


def overwrite_reader_samples(vcf_reader: vcf.Reader, samples):
    vcf_reader.samples = samples
    vcf_reader._sample_indexes = dict([(x,i) for (i,x) in enumerate(vcf_reader.samples)])


class Variant(object):
    def __init__(self, record, reference, padding=1):
        self.record = record
        # Make sure every variant has a unique ID
        if self.record.ID is None:
            self.record.ID = variant_descriptor(record)
        self.reference = reference
        self.padding = padding

        # Update padding for complex sequence-resolved variants
        ref_allele = record.REF
        alt_allele = record.ALT[0]
        if isinstance(alt_allele, vcf.model._Substitution):
            padding_string = os.path.commonprefix([ref_allele, str(alt_allele)])
            self.padding = len(padding_string)
        if self.padding > 1:
            logging.warning("%s has more than expected number of padding bases, if VCF normalized?", self.record.ID)    

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
    def event_length(self):
        """Length difference between reference and alternate alleles, i.e. absolute value of SVLEN field"""
        return abs(self.record.INFO["SVLEN"][0])

    @property
    def ref_length(self):
        """Length of reference allele including any padding bases"""
        raise NotImplementedError()

    @property
    def alt_length(self):
        """Length of alternate allele including any padding bases"""
        raise NotImplementedError()

    @property
    def subtype(self):
        return self.record.var_subtype

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

    @property
    def is_insertion(self):
        return False

    @property
    def is_duplication(self):
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
        """Return 1-indexed fully closed region describing inside the variant, e.g. the deleted region"""
        # In correctly formatted VCF, POS is first base of event when zero-indexed, while
        # END is 1-indexed closed end or 0-indexed half-open end
        return f"{self.chrom}:{self.pos+self.padding-flank}-{self.end+flank}"

    def left_flank_region_string(self, left_flank, right_flank=0):
        """Return 1-indexed fully closed region"""
        return f"{self.chrom}:{self.pos+self.padding-left_flank}-{self.pos+self.padding-1+right_flank}"
    
    def right_flank_region_string(self, right_flank, left_flank=0):
        """Return 1-indexed fully closed region"""
        return f"{self.chrom}:{self.end+1-left_flank}-{self.end+right_flank}"

    def ref_breakpoints(self, flank, contig=None):
        if contig is None:
            contig = self.chrom
        event_end = flank + self.ref_length - self.padding # ref_length includes padding base
        return (f"{contig}:{flank}-{flank + 1}", f"{contig}:{event_end}-{event_end+1}" if event_end > flank else None)

    def alt_breakpoints(self, flank, contig=None):
        if contig is None:
            contig = self.chrom
        event_end = flank + self.alt_length - self.padding # alt_length includes padding base
        return (f"{contig}:{flank}-{flank + 1}", f"{contig}:{event_end}-{event_end + 1}" if event_end > flank else None)

    def reference_sequence(self, region=None, flank=0):
        try:
            if region is None:
                region = self.region_string(flank)

            # pylint: disable=no-member
            with pysam.FastaFile(self.reference) as ref_fasta:
                ref_seq = ref_fasta.fetch(region=region)
                return ref_seq
        except ValueError as err:
            print(region)
            raise err

    def _alt_seq(self, flank, ref_seq):
        raise NotImplementedError()    

    def synth_fasta(
        self,
        args,
        allele_fasta=None,
        ac=1,
        ref_contig=None,
        alt_contig=None,
        line_width=60,
    ):
        region = self.region_string(args.flank)
        ref_seq = self.reference_sequence(region=region)
        if ac != 0:
            alt_seq = self._alt_seq(args.flank, ref_seq)
        
        if ref_contig is None:
            ref_contig = region.replace(":", "_").replace("-", "_")
        if alt_contig is None:
            alt_contig = ref_contig + "_alt"
        if allele_fasta is None:
            allele_fasta = tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".fasta", dir=args.tempdir
            )

        # Write out FASTA
        print(">", ref_contig, sep="", file=allele_fasta)
        if ac == 0 or ac == 1:
            for line in textwrap.wrap(ref_seq, width=line_width):
                print(line, file=allele_fasta)
        print(">", alt_contig, sep="", file=allele_fasta)
        if ac == 1 or ac == 2:
            for line in textwrap.wrap(alt_seq, width=line_width):
                print(line, file=allele_fasta)

        return allele_fasta.name, ref_contig, alt_contig

    @cached_property
    def ref_gc_fraction(self):
        gc = 0
        alignable_bases = 0
        for base in self.reference_sequence(flank=0):
            if base == "G" or base == "C" or base == "g" or base == "c":
                gc += 1
                alignable_bases += 1
            elif base != "N" and base != "n":
                alignable_bases += 1
        if alignable_bases > 0:
            return float(gc) / alignable_bases, alignable_bases
        else:
            return 0.0, 0

    def to_minimal_vcf(self, args):
        raise NotImplementedError()

    @classmethod
    def from_pyvcf(cls, record, reference=None):
        if not record.is_sv:
            return None
        if len(record.ALT) != 1:
            return None

        kind = record.var_subtype
        if kind.startswith("DEL"):
            return DeletionVariant(record, reference)
        elif kind.startswith("INS"):
            return InsertionVariant(record, reference)
        elif kind.startswith("DUP"):
            return DuplicationVariant(record, reference)


class DeletionVariant(Variant):
    def __init__(self, record, reference):
        Variant.__init__(self, record, reference)

    @property
    def subtype(self):
        return self.record.var_subtype

    @property
    def is_deletion(self):
        return True

    @property
    def ref_length(self):
        return self.end - self.pos + 1

    @property
    def alt_length(self):
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
        assert (
            "SVLEN" in self.record.INFO
        ), "SVLEN field not present in upstream variant"
        base_info = {
            "SVTYPE": "DEL",
            "END": self.end,
            "SVLEN": ",".join(map(str, self.record.INFO["SVLEN"])),
        }
        base_info.update(info)
        line = "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t.\t.\t{info_string}".format(
            chrom=self.record.CHROM,
            pos=self.record.POS,
            id=self.record.ID or ".",
            ref=self.record.REF,
            alt=self.record.ALT[0],
            info_string=";".join((f"{k}={v}" for (k, v) in base_info.items())),
        )
        return line

    def to_minimal_vcf(self, args, tempdir=None):
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

    def _alt_seq(self, flank, ref_seq):
        assert len(self.record.ALT) == 1, "Multiple alternates are not supported"
        allele = self.record.ALT[0]
        if isinstance(allele, vcf.model._SV):
            return ref_seq[: flank] + ref_seq[-flank :]
        elif isinstance(allele, vcf.model._Substitution):
            return ref_seq[: flank] + str(allele)[self.padding:] + ref_seq[-flank :]
        else:
            raise ValueError("Unsupported allele type")

    def gnomad_coverage_profile(
        self,
        args,
        gnomad_coverage: str,
        covg_file=None,
        ref_contig=None,
        alt_contig=None,
        line_width=60,
    ):
        region = self.region_string(args.flank)
        gnomad_coverage_tabix = pysam.TabixFile(gnomad_coverage)

        # Use a FASTQ-like +33 scheme for encoding depth
        ref_covg = "".join(
            map(
                lambda x: chr(min(round(float(x[2])) + 33, 126)),
                gnomad_coverage_tabix.fetch(region=region, parser=pysam.asTuple()),
            )
        )
        gnomad_coverage_tabix.close()

        assert len(self.record.ALT) == 1, "Multiple alternates are not supported"
        allele = self.record.ALT[0]
        if isinstance(allele, vcf.model._SV):
            alt_covg = ref_covg[: args.flank] + ref_covg[-args.flank :]
        elif isinstance(allele, vcf.model._Substitution):
            alt_covg = (
                ref_covg[: args.flank - self.padding]
                + ref_covg[args.flank - self.padding] * len(allele)
                + ref_covg[-args.flank :]
            )
        else:
            raise ValueError("Unsupported allele type")

        # Write out
        if ref_contig is None:
            ref_contig = region.replace(":", "_").replace("-", "_")
        if alt_contig is None:
            alt_contig = ref_contig + "_alt"
        if covg_file is None:
            covg_file = tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".fasta", dir=args.tempdir
            )

        print(ref_contig, ref_covg, sep="\t", file=covg_file)
        print(alt_contig, alt_covg, sep="\t", file=covg_file)

        return covg_file.name, ref_contig, alt_contig


class InsertionVariant(Variant):
    def __init__(self, record, reference):
        Variant.__init__(self, record, reference)

    @property
    def is_insertion(self):
        return True

    @property
    def ref_length(self):
        return 1

    @property
    def alt_length(self):
        allele = self.record.ALT[0]
        if isinstance(allele, vcf.model._SV):
            # Symbolic allele
            return self.event_length + 1
        else:
            return len(allele)

    def region_string(self, flank=0):
        if flank == 0:
            raise ValueError("Can't represent interior region of an insertion")
        return super().region_string(flank=flank)  

    def _alt_seq(self, flank, ref_seq):
        assert len(self.record.ALT) == 1, "Multiple alternates are not supported"
        allele = self.record.ALT[0]
        if isinstance(allele, vcf.model._SV):
            # Work with insertions as encoded by Manta
            if "SVINSSEQ" in self.record.INFO:
                allele = self.record.INFO["SVINSSEQ"][0]
                return ref_seq[: flank] + str(allele) + ref_seq[flank :]
            else:
                raise ValueError("Unsupported allele type")
        elif isinstance(allele, vcf.model._Substitution):
            return ref_seq[: flank - self.padding] + str(allele) + ref_seq[flank :]
        else:
            raise ValueError("Unsupported allele type")


class DuplicationVariant(Variant):
    def __init__(self, record, reference):
        Variant.__init__(self, record, reference)

    @property
    def is_duplication(self):
        return True

    @property
    def ref_length(self):
        return self.end - self.pos + 1

    def ref_breakpoints(self, flank, contig=None):
        if contig is None:
            contig = self.chrom
        duplication_junction = flank + self.ref_length - self.padding # ref_length includes padding base
        return (f"{contig}:{flank}-{flank + 1}", f"{contig}:{duplication_junction}-{duplication_junction + 1}")

    def alt_breakpoints(self, flank, contig=None):
        if contig is None:
            contig = self.chrom
        duplication_junction = flank + self.ref_length - self.padding # ref_length includes padding base
        return (f"{contig}:{duplication_junction}-{duplication_junction + 1}", None)

    def _alt_seq(self, flank, ref_seq):
        allele = self.record.ALT[0]
        if isinstance(allele, vcf.model._SV):
            event_length = self.event_length
            assert 2*flank + event_length == len(ref_seq), "Inconsistent length for reference sequence"
            return ref_seq[: flank] + 2 * ref_seq[flank:flank+event_length] + ref_seq[-flank :] 
        else:
            raise ValueError("Unsupported allele type")


def consensus_fasta(args, input_vcf: str, output_file):
    record = next(vcf.Reader(filename=input_vcf))
    variant = Variant.from_pyvcf(record, args.reference)
    variant.synth_fasta(
        args,
        allele_fasta=output_file,
        ac=args.ac,
        ref_contig=args.ref_contig,
        alt_contig=args.alt_contig,
    )

