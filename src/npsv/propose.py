import numpy as np
from scipy.signal import peak_prominences, find_peaks
import vcf
import pysam
from .variant import Variant, DeletionVariant
from npsv import npsva

PEAK_FINDING_FLANK=5

def propose_variants(args, input_path, output_file):
    """Generate alternate representations of the a variant"""

    if args.simple_repeats_bed:
        # TODO? Require certain overlap?
        simple_repeats_bed = pysam.TabixFile(args.simple_repeats_bed)
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
        print(variant.to_minimal_vcf_record(), file=output_file)

        repeats = simple_repeats_bed.fetch(region=variant.region_string(), parser=pysam.asTuple()) if simple_repeats_bed else []
        if not repeats:
            continue
        
        proposed_variants = set()
        for repeat in repeats:
            consensus_length = int(repeat[3])
            repeat_length = consensus_length*float(repeat[4])
            
            if variant.event_length > repeat_length or (repeat_length - variant.event_length) / consensus_length < 1:
                continue
            
            event_repeat_count = round(variant.event_length / consensus_length)
            if event_repeat_count == 0:
                continue

            # Find event starts
            repeat_start = int(repeat[1]) - PEAK_FINDING_FLANK # Flank for peak finding
            repeat_end = int(repeat[2]) + PEAK_FINDING_FLANK

            # repeat is BED-file with half-open start
            ref_seq = variant.reference_sequence(
                args, region=f"{repeat[0]}:{repeat_start+1}-{repeat_end}"
            )
            
            consensus_seq = repeat[5]
            scores = []
            for i in range(0, len(ref_seq)-len(consensus_seq)):
                matches = sum(c1 == c2 for c1, c2 in zip(consensus_seq, ref_seq[i:i+len(consensus_seq)]))
                scores.append(matches)
        
            peaks, properties = find_peaks(scores, width=1, distance=consensus_length*0.8)
        
            # Enforce maximum number of potential alternate variants by selecting most prominent
            peaks = peaks[np.argsort(properties["prominences"])[:args.max_proposals]]

            # Generate alternate records
            for peak in peaks:
                alt_pos = peak + repeat_start  # 1-indexed base immediately before event
                alt_end = peak + repeat_start + event_repeat_count * consensus_length
                if alt_pos == variant.pos and alt_end == variant.end:
                    continue # Already have representation of this variant
                
                # TODO: Realign allele sequence to get better end coordinate?
                proposed_variants.add((alt_pos, alt_end, ref_seq[peak-1]))


        for pos, end, ref in sorted(proposed_variants):
            record = f"{variant.chrom}\t{pos}\t.\t{ref}\t<DEL>\t.\t.\tSVTYPE=DEL;END={end};SVLEN={-(end - pos)};ORIGINAL={variant.id}"
            print(record, file=output_file)

        
