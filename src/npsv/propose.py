import collections, logging, sys
import numpy as np
from scipy.signal import peak_prominences, find_peaks
from scipy.special import logsumexp
from scipy.stats import binom, chi2
import vcf
import pysam
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

            # Find event starts
            repeat_start = int(repeat[1]) - PEAK_FINDING_FLANK  # Flank for peak finding
            repeat_end = int(repeat[2]) + PEAK_FINDING_FLANK

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


def refine_variants(args, input_vcf: str, output_file):
    """Identify the "best" representation for a variant

    Updates the genotypes for "original" variant with more similar alternate representation. Note that the 
    resulting VCF file is not in sorted order.

    Args:
        args: Arguments from arg parseer
        input_vcf (str): Path to input VCF file
        output_file: Output file object
    """

    # Setup VCF reader and writer...
    vcf_reader = vcf.Reader(filename=input_vcf)
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

    for id, record in original_records.items():
        if id not in alternate_records:
            # No alternate records present for this variant
            vcf_writer.write_record(record)
        else:
            # Update record with new FORMAT
            record.FORMAT = VCF_FORMAT

            # Identify best alternate representation and genotype for each sample
            for i, call in enumerate(record.samples):
                best_call = vcf.model._Call(
                    record,
                    call.sample,
                    AltCallData(
                        GT=call.data.GT,
                        PL=call.data.PL,
                        DM=call.data.DM,
                        AD=call.data.AD,
                        CL=None,
                        OGT=None,
                        OPL=None,
                        ODM=None,
                        OAD=None,
                    ),
                )
                if args.select == "dm2" and best_call.data.DM:
                    min_dist = min(best_call.data.DM[1:])
                elif args.select == "prob" and best_call.data.AD:
                    max_prob = max(ad_to_prob(best_call.data.AD)[1:])
                else:
                    record.samples[i] = best_call
                    continue

                # Identify other representations to see if we want to overwrite the genotype
                # Update if:
                # 1. Smallest DM^2 (largest GT prob) for alternate variant is for non-reference genotype, and
                # 2. Alternate variant's smallest DM^2 (largest GT prob) is less (greater) than original's non-reference DM^2 (GT prob)
                for alternate_record in alternate_records[id]:
                    alternate_call = alternate_record.samples[i]
                    if args.select == "dm2" and alternate_call.data.DM:
                        min_alt_dist_idx = np.argmin(alternate_call.data.DM)
                        min_alt_dist = alternate_call.data.DM[min_alt_dist_idx]
                        if min_alt_dist_idx == 0 or min_alt_dist >= min_dist:
                            continue
                        min_dist = min_alt_dist
                    elif args.select == "prob" and alternate_call.data.AD:
                        alt_prob = ad_to_prob(alternate_call.data.AD)
                        max_alt_prob_idx = np.argmax(alt_prob)
                        max_alt_prob = alt_prob[max_alt_prob_idx]
                        if max_alt_prob_idx == 0 or max_alt_prob <= max_prob:
                            continue
                        max_prob = max_alt_prob
                    else:
                        continue
                    
                    best_call_data = AltCallData(
                        GT=alternate_call.data.GT,
                        PL=alternate_call.data.PL,
                        DM=alternate_call.data.DM,
                        AD=alternate_call.data.AD,
                        CL=variant_descriptor(alternate_record),
                        OGT=call.data.GT,
                        OPL=call.data.PL,
                        ODM=call.data.DM,
                        OAD=call.data.AD,
                    )
                    best_call = vcf.model._Call(record, call.sample, best_call_data)

                record.samples[i] = best_call
            vcf_writer.write_record(record)
