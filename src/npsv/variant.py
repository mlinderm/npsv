import vcf
import pysam


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

