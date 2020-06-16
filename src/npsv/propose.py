import vcf
import pybedtools as bedtools
from .variant import Variant, DeletionVariant
from npsv import npsva
from scipy.signal import peak_prominences, find_peaks


def propose_variants(args, input_path, output_file):
    """Generate alternate representations of the a variant"""

    if args.simple_repeats_bed:
        # TODO? Require certain overlap?
        simple_repeats_bed = bedtools.BedTool(args.simple_repeats_bed)
    else:
        simple_repeats_bed = None

    print(
        """\
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=HAS_PROPOSED,Number=1,Type=String,Description="Original variant with alternate proposed representations">      
##INFO=<ID=ORIGINAL,Number=1,Type=String,Description="This record is a proposed alternate representation for variant with this ID">      
##ALT=<ID=DEL,Description="Deletion">
##contig=<ID=1,length=249250621>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO""",
    )

    vcf_reader = vcf.Reader(filename=input_path)
    for record in vcf_reader:
        variant = Variant.from_pyvcf(record)
        repeats = simple_repeats_bed.tabix_intervals(variant.region_string()) if simple_repeats_bed else []
        if not repeats:
            print(variant.to_minimal_vcf_record(), file=output_file)
            continue

        var_seq = variant.reference_sequence(args)

        # Repeat properties to specify region for alignment
        #min_length = min([int(repeat.fields[3]) for repeat in repeats])
        min_start = min(repeat.start for repeat in repeats)
        max_end = max(repeat.end for repeat in repeats)

        ref_seq = variant.reference_sequence(
            args, region=f"{variant.chrom}:{min_start}-{max_end}"
        )
  
        # Slide allele down sequence looking for high scoring alignments
        scores = []
        for i in range(0, len(ref_seq)-len(var_seq)):
            matches = sum(c1 == c2 for c1, c2 in zip(var_seq, ref_seq[i:i+len(var_seq)]))
            scores.append(matches)

        # Require peak to be larger than random matches...
        # peaks + min_start is the 1-indexed start coordinate for first deleted base
        peaks, _ = find_peaks(scores, width=1, prominence=int(len(var_seq) * 0.25))
        peaks = peaks[peaks + min_start != variant.pos + 1] # Exclude exact match with original sequence
        
        # Generate original record
        print(variant.to_minimal_vcf_record({"HAS_PROPOSED": "True"}), file=output_file)
        
        # Generate alternate records
        for peak in peaks:
            # TODO: Realign allele sequence to get better end coordinate?
            record = "{chrom}\t{pos}\t.\t{ref}\t<DEL>\t.\t.\tSVTYPE=DEL;END={end};SVLEN={len};ORIGINAL={original}".format(
                chrom=variant.chrom,
                pos=peak + min_start - 1,
                ref=ref_seq[peak-1],
                end=peak + min_start + len(var_seq)-1,
                len=-len(var_seq),
                original=variant.id,
            )
            print(record, file=output_file)

       

        # # Get the reference sequence for the variant (with a padding base for constructing subsequent variants)
        # right_ref_seq = variant.reference_sequence(
        #     args, region=f"{variant.chrom}:{variant.pos+min_length}-{max_end}"
        # )
        # right_other_spans = npsva.alternate_variant_locations(
        #     right_ref_seq[1:], var_seq
        # )

        # # Generate original record
        # print(variant.to_minimal_vcf_record({"HAS_PROPOSED": "True"}), file=output_file)
        
        # # Generate alternate records
        # for (left_closed, right_closed) in right_other_spans:
        #     record = "{chrom}\t{pos}\t.\t{ref}\t<DEL>\t.\t.\tSVTYPE=DEL;END={end};SVLEN={len};ORIGINAL={original}".format(
        #         chrom=variant.chrom,
        #         pos=variant.end + left_closed,
        #         ref=right_ref_seq[left_closed],
        #         end=variant.end + right_closed,
        #         len=-(right_closed - left_closed),
        #         original=variant.id,
        #     )
        #     print(record, file=output_file)
