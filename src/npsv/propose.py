import collections, logging, operator, sys
import numpy as np
from scipy.signal import peak_prominences, find_peaks
from scipy.special import logsumexp
from scipy.stats import binom, chi2
import vcf
import pysam
from intervaltree import Interval, IntervalTree
from .variant import Variant, variant_descriptor
from .genotyper import MAHAL_FEATURES

PEAK_FINDING_FLANK = 5
ORIGINAL_KEY = "ORIGINAL"
VCF_FORMAT = "GT:PL:DM:AD:CL:OGT:OPL:ODM:OAD"


def propose_variants(args, input_vcf: str, output_file):
    """Generate alternate representations for variants in a VCF file.

    Note that the resulting VCF file is not in sorted order.

    Args:
        args: Arguments from arg parseer
        input_vcf (str): Path to input VCF file
        output_file: Output file object
    """

    # Setup VCF reader and writer...
    vcf_reader = vcf.Reader(filename=input_vcf)
    assert (
        ORIGINAL_KEY not in vcf_reader.infos
    ), f"{ORIGINAL_KEY} already presented in VCF INFO field"
    vcf_reader.infos["ORIGINAL"] = vcf.parser._Info(
        "ORIGINAL",
        ".",
        "String",
        "This record is a proposed alternate representation for these variant IDs",
        None,
        None,
    )
    vcf_writer = vcf.Writer(output_file, vcf_reader)

    # Setup repeats file
    if args.simple_repeats_bed:
        simple_repeats_bed = pysam.TabixFile(args.simple_repeats_bed)
    else:
        logging.warn(
            "Proposing alternate variant representations requires BED file with simple repeats"
        )
        simple_repeats_bed = None

    proposed_variants = {}

    for record in vcf_reader:
        variant = Variant.from_pyvcf(record, args.reference)
        assert variant.id, "Variant proposal requires all variants to have a unique ID"
        vcf_writer.write_record(variant.record)
        
        if not variant.is_deletion:
            # Only deletions currently supported
            continue

        if variant.event_length > args.DEL_hybrid_threshold:
            continue

        # TODO? Require certain overlap?
        repeats = (
            simple_repeats_bed.fetch(
                region=variant.region_string(), parser=pysam.asTuple()
            )
            if simple_repeats_bed
            else []
        )
        if not repeats:
            continue

        for repeat in repeats:
            consensus_length = int(repeat[3])
            # Only propose variants for larger VNTRs
            if consensus_length < args.min_consensus_length:
                continue

            # Only propose variants if original variant is smaller than repeat region
            repeat_length = consensus_length * float(repeat[4])
            if (
                variant.event_length > repeat_length
                or (repeat_length - variant.event_length) / consensus_length < 1
            ):
                continue

            event_repeat_count = round(variant.event_length / consensus_length)
            if event_repeat_count == 0:
                continue

            # Find event starts (don't allow event to extend past end of repeat)
            repeat_start = int(repeat[1]) - PEAK_FINDING_FLANK  # Flank for peak finding
            repeat_end = int(repeat[2]) + PEAK_FINDING_FLANK - variant.event_length

            # repeat is BED-file with half-open start
            ref_seq = variant.reference_sequence(f"{repeat[0]}:{repeat_start+1}-{repeat_end}")

            consensus_seq = repeat[5]
            scores = []
            for i in range(0, len(ref_seq) - len(consensus_seq)):
                matches = sum(
                    c1 == c2
                    for c1, c2 in zip(
                        consensus_seq, ref_seq[i : i + len(consensus_seq)]
                    )
                )
                scores.append(matches)

            peaks, properties = find_peaks(
                scores, width=1, distance=consensus_length * 0.8
            )

            # Enforce maximum number of potential alternate variants by selecting most prominent
            peaks = peaks[np.argsort(properties["prominences"])[: args.max_proposals]]

            # Generate alternate records
            for peak in peaks:
                # TODO: Realign allele sequence to get better end coordinate?
                alt_pos = peak + repeat_start  # 1-indexed base immediately before event
                alt_end = peak + repeat_start + event_repeat_count * consensus_length
                if alt_pos == variant.pos and alt_end == variant.end:
                    continue  # Already have representation of this variant

                # TODO: Handle different variant types
                key = (variant.chrom, alt_pos, alt_end)
                if key in proposed_variants:
                    proposed_variants[key][4].append(variant.id)
                else:
                    proposed_variants[key] = (
                        variant.chrom,
                        alt_pos,
                        alt_end,
                        ref_seq[peak - 1],
                        [variant.id],
                    )

    for chrom, pos, end, ref, originals in proposed_variants.values():
        record = f"{chrom}\t{pos}\t.\t{ref}\t<DEL>\t.\t.\tSVTYPE=DEL;END={end};SVLEN={-(end - pos)};{ORIGINAL_KEY}={','.join(originals)}"
        if vcf_reader.samples:
            record += "\tGT"
        for _ in vcf_reader.samples:
            record += "\t."  # Define unknown genotypes for proposed variants
        print(record, file=output_file)


def dm2_to_prob(score, df=len(MAHAL_FEATURES)):
    with np.errstate(under="ignore"):
        prob = chi2.logsf(score, df)
        return np.exp(prob - logsumexp(prob))


def ad_to_prob(ad):
    ad = [int(reads) for reads in ad]
    assert len(ad) == 2, "Unexpected number of alleles in call"
    total_reads = sum(ad)
    return [
        binom.pmf(ad[1], total_reads, 0.05),
        binom.pmf(ad[1], total_reads, 0.5),
        binom.pmf(ad[1], total_reads, 0.95),
    ]


class RangeTree:
    def __init__(self):
        self._trees = collections.defaultdict(IntervalTree)

    def add(self, contig: str, start: int, end: int, data):
        return self._trees[contig].addi(start, end, data)

    def merge_overlaps(self, **kwargs):
        for tree in self._trees.values():
            tree.merge_overlaps(**kwargs)

    def values(self):
        for tree in self._trees.values():
            for _, _, data in tree.items():
                yield data


def refine_variants(args, input_vcf: str, output_file):
    """Identify the "best" representation for a variant

    Updates the genotypes for "original" variant with more similar alternate representation. Note that the 
    resulting VCF file is not in sorted order.

    Args:
        args: Arguments from arg parseer
        input_vcf (str): Path to input VCF file
        output_file: Output file object
    """

    orig_start_idx = 0 if args.include_orig_ref else 1

    # Setup VCF reader and writer...
    vcf_reader = vcf.Reader(filename=input_vcf)
    num_samples = len(vcf_reader.samples)

    # Add new format fields
    vcf_reader.formats["CL"] = vcf.parser._Format(
        "CL", "1", "String", "Call location used for genotype",
    )
    vcf_reader.formats["OGT"] = vcf.parser._Format(
        "OGT", "1", "String", "Genotype for the original variant",
    )
    vcf_reader.formats["OPL"] = vcf.parser._Format(
        "OPL", "G", "Integer", "Genotype likelihood for each genotype",
    )
    vcf_reader.formats["ODM"] = vcf.parser._Format(
        "ODM", "G", "Float", "Original Mahalanobis distance for each genotype",
    )
    vcf_reader.formats["OAD"] = vcf.parser._Format(
        "OAD", "R", "Integer", "Original read depth for each allele",
    )
    vcf_writer = vcf.Writer(args.output, vcf_reader)

    # TODO: Include data from original format
    AltCallData = collections.namedtuple(
        "AltCallData", ["GT", "PL", "DM", "AD", "CL", "OGT", "OPL", "ODM", "OAD"]
    )

    original_records = {}
    alternate_records = {}

    for record in vcf_reader:
        if ORIGINAL_KEY not in record.INFO:
            assert record.ID not in original_records, "Duplicate original variants"
            original_records[record.ID] = record
        else:
            originals = record.INFO[ORIGINAL_KEY]
            for original in originals:
                if original in alternate_records:
                    alternate_records[original].append(record)
                else:
                    alternate_records[original] = [record]

    # Determine variant groups 
    variant_ranges = RangeTree()
    for id, original_record in original_records.items():
        total_range = Interval(original_record.POS - 1, int(original_record.sv_end)) # Convert to 0-indexed half open
        for alternate_record in alternate_records.get(id, []):
            total_range = Interval(min(alternate_record.POS - 1, total_range.begin), max(int(alternate_record.sv_end), total_range.end))
        variant_ranges.add(original_record.CHROM, total_range.begin, total_range.end, [id])

    if args.merge_orig_blocks:    
        variant_ranges.merge_overlaps(data_reducer=operator.add, data_initializer=[])


    # Determine the best alternate representation(s) in each group
    closest_alts = {}
    for ids in variant_ranges.values():
        if args.select == "dm2":
            best_alts = [(float("Inf"), None, None)] * num_samples
        elif args.select == "prob":
            best_alts = [(0, None, None)] * num_samples
        
        for id in ids:
            possible_alternate_records = alternate_records.get(id, [])
            if args.include_orig_in_block:
                possible_alternate_records.append(original_records[id])
            for alternate_record in possible_alternate_records:
                for i, alternate_call in enumerate(alternate_record.samples):
                    best_score, *_ = best_alts[i]
                    
                    if args.select == "dm2" and alternate_call.data.DM:
                        min_alt_dist_idx = np.argmin(alternate_call.data.DM)
                        min_alt_dist = alternate_call.data.DM[min_alt_dist_idx]
                        if min_alt_dist_idx == 0 or min_alt_dist >= best_score:
                            continue
                        best_alts[i] = (min_alt_dist, id, alternate_record)
                    elif args.select == "prob" and alternate_call.data.AD:
                        alt_prob = ad_to_prob(alternate_call.data.AD)
                        max_alt_prob_idx = np.argmax(alt_prob)
                        max_alt_prob = alt_prob[max_alt_prob_idx]
                        if max_alt_prob_idx == 0 or max_alt_prob <= best_score:
                            continue
                        best_alts[i] = (max_alt_prob, id, alternate_record)
                    else:
                        continue

        for id in ids:
            closest_alts[id] = best_alts


    for id, record in original_records.items():
        if id not in closest_alts:
            # No alternate records present for this variant
            vcf_writer.write_record(record)
        else:
            # Update record with new FORMAT
            record.FORMAT = VCF_FORMAT
           
            # Identify the best alternate representation and genotype for each sample
            closest_alt = closest_alts[id]
            for i, call in enumerate(record.samples):
                orig_call = vcf.model._Call(
                    record,
                    call.sample,
                    AltCallData(
                        GT=call.data.GT,
                        PL=getattr(call.data, "PL", None),
                        DM=getattr(call.data, "DM", None),
                        AD=getattr(call.data, "AD", None),
                        CL=None,
                        OGT=None,
                        OPL=None,
                        ODM=None,
                        OAD=None,
                    ),
                )
                record.samples[i] = orig_call
                
                alt_score, alt_id, alt_record = closest_alt[i]
                if alt_record is None or alt_record is record:
                    # No alternate record to update with, or we are tying to update with ourselves (same original record)
                    continue
                
                if args.select == "dm2" and call.data.DM:
                    orig_dist = min(call.data.DM[orig_start_idx:])
                    if alt_score >= orig_dist:
                        continue 
                elif args.select == "prob" and best_call.data.AD:
                    orig_prob = max(ad_to_prob(best_call.data.AD)[orig_start_idx:])
                    if alt_score <= orig_dist:
                        continue 
                else:
                    continue

                alt_call = alt_record.samples[i]
                best_call_data = AltCallData(
                    GT=alt_call.data.GT if alt_id == id else "0/0",
                    PL=getattr(alt_call.data, "PL", None),
                    DM=getattr(alt_call.data, "DM", None),
                    AD=getattr(alt_call.data, "AD", None),
                    CL=variant_descriptor(alt_record),
                    OGT=call.data.GT,
                    OPL=getattr(call.data, "PL", None),
                    ODM=getattr(call.data, "DM", None),
                    OAD=getattr(call.data, "AD", None),
                )
                best_call = vcf.model._Call(record, call.sample, best_call_data)
                record.samples[i] = best_call

            vcf_writer.write_record(record)
