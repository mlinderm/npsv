import logging, math, sys, warnings
from typing import List
import vcf
import pybedtools.bedtool as bed
import pandas as pd
import numpy as np
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from scipy.special import logsumexp
from sklearn.ensemble import RandomForestClassifier
from sklearn.covariance import MinCovDet
import scipy.spatial.distance as distance
from scipy.stats import binom, chi2, zscore
from sklearn.metrics import accuracy_score
from scipy.special import comb

from .variant import variant_descriptor
from .feature_extraction import Features

# Suppress the future warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
if logging.getLogger(__name__).getEffectiveLevel() != logging.DEBUG:
    warnings.simplefilter(action="ignore", category=UserWarning)

VAR_COL = ["#CHROM", "START", "END", "TYPE", "SAMPLE"]
FEATURE_COL = Features.FEATURES
KLASS_COL = "AC"
SAMPLE_COL = "SAMPLE"
OUTPUT_COL = (
    VAR_COL
    + FEATURE_COL
    + [
        "AC",
        "REF_GENOTYPE_LIKELIHOOD",
        "HET_GENOTYPE_LIKELIHOOD",
        "HOM_GENOTYPE_LIKELIHOOD",
        "ALT_GENOTYPE_LIKELIHOOD",
        "REF_QUAL",
        "ALT_QUAL",
        "CLASSIFIER",
    ]
)

# Extracted and derived features
FEATURES = Features.FEATURES + [
    "REF_READ_REL",
    "ALT_READ_REL",
    "PROB_HOMREF",
    "PROB_HET",
    "PROB_HOMALT",
    "REF_SPAN_REL",
    "ALT_SPAN_REL",
]
ABSOLUTE_FEATURES = ["REF_SPLIT", "ALT_SPLIT", "REF_SPAN", "ALT_SPAN", "COVG"]

C_RANGE = np.logspace(-2, 10, 13)
GAMMA_RANGE = np.logspace(-9, 3, 13)
PARAM_GRID = dict(gamma=GAMMA_RANGE, C=C_RANGE)

VCF_COLUMN_HEADERS = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
]
VCF_FORMAT = "GT:GR"
AC_TO_GT = ("0/0", "0/1", "1/1")


def variant_to_string(variant):
    return "{0}:{1}-{2} {3} in {4}".format(*variant)


def phred(prob, max_val=99):
    if 1 - prob > 0:
        return abs(min(-10.0 * math.log10(1 - prob), max_val))
    else:
        return max_val


def get_tsv_columns(filename: str) -> List[str]:
    with open(filename, "r") as file:
        return file.readline().split()


def record_to_var_col(record: vcf.model._Record, sample):
    return (record.CHROM, record.POS, int(record.sv_end), record.var_subtype, sample)


def pred_to_vcf(real_data, pred) -> str:
    """Prediction to VCF call field for sample"""
    assert real_data.shape[0] == 1, "Real data should have just one variant"
    return "{gt}:{grr},{gra}".format(
        gt=AC_TO_GT[pred],
        grr=int(real_data["REF_SPLIT"].iloc[0]),
        gra=int(real_data["ALT_SPLIT"].iloc[0]),
    )


def multi_rf_classify(sim_data, real_data, features=FEATURE_COL, klass=KLASS_COL):
    x = sim_data[features]
    y = sim_data[klass]

    clf = RandomForestClassifier(max_depth=3, min_samples_leaf=5)
    clf = clf.fit(x, y)

    real_x = real_data[features]
    pred = clf.predict(real_x)

    return pred


def multi_svm_classify(sim_data, real_data, features=FEATURE_COL, klass=KLASS_COL):
    x = sim_data[features]
    y = sim_data[klass]

    # Normalize the training set
    scaler = StandardScaler()
    x = scaler.fit_transform(x)

    # Build the model
    clf = SVC(kernel="rbf")

    # Fit the model
    clf.fit(x, y)

    # Predict the real data
    real_x = scaler.transform(real_data[features])
    return clf.predict(real_x)


def single_mahalanobis(sim_data, real_data, features=FEATURE_COL, klass=KLASS_COL):
    """Classify single variant using Mahalanobis distance"""
    assert real_data.shape[0] == 1, "Real data should have just one variant"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        x = sim_data[features]
        y = sim_data[klass]
        real_x = real_data[features]

        def mahal_score(data):
            robust_cov = MinCovDet(assume_centered=False).fit(data)
            return robust_cov.mahalanobis(real_x)[0]

        score = x.groupby(y).apply(mahal_score)

        pred = score.idxmin()
        with np.errstate(under="ignore"):
            prob = chi2.logsf(score, len(features))
            prob = np.exp(prob - logsumexp(prob))

        return (pred, prob, score)


def single_svm_classify(sim_data, real_data, features=FEATURE_COL, klass=KLASS_COL):
    """Classify single variant using SVM"""
    assert real_data.shape[0] == 1, "Real data should have just one variant"

    x = sim_data[features]
    y = sim_data[klass]

    # Normalize the training set
    scaler = StandardScaler()
    x = scaler.fit_transform(x)

    # Grid-search to find the best gamma and C values for the rbf kernel
    cv = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=42)
    grid = GridSearchCV(
        SVC(probability=True, kernel="rbf"), param_grid=PARAM_GRID, cv=cv
    )

    # Fit the model
    grid.fit(x, y)

    # Predict the real data
    real_x = scaler.transform(real_data[features])

    pred = grid.predict(real_x)[0]
    prob = grid.predict_proba(real_x)[0]

    return (pred, prob)


def single_randomforest_classify(
    sim_data, real_data, features=FEATURE_COL, klass=KLASS_COL
):
    """Classify single variant using Random Forests"""
    assert real_data.shape[0] == 1, "Real data should have just one variant"

    x = sim_data[features]
    y = sim_data[klass]

    clf = RandomForestClassifier(max_depth=3, min_samples_leaf=5)
    clf = clf.fit(x, y)

    real_x = real_data[features]
    pred = clf.predict(real_x)[0]
    prob = clf.predict_proba(real_x)[0]

    return (pred, prob)


def add_derived_features(data):
    """Update data with derived features"""
    # Binomial probability features adapted from svviz2 and
    # https://www.sciencedirect.com/science/article/pii/S0092867418316337#sec4
    total_reads = data["REF_SPLIT"] + data["ALT_SPLIT"]
    data["REF_READ_REL"] = np.where(
        total_reads == 0, 0, data["REF_SPLIT"] / total_reads
    )
    data["ALT_READ_REL"] = np.where(
        total_reads == 0, 0, data["ALT_SPLIT"] / total_reads
    )

    data["PROB_HOMREF"] = binom.pmf(data["ALT_SPLIT"], total_reads, 0.05)
    data["PROB_HET"] = binom.pmf(data["ALT_SPLIT"], total_reads, 0.5)
    data["PROB_HOMALT"] = binom.pmf(data["ALT_SPLIT"], total_reads, 0.95)

    total_span = data["REF_SPAN"] + data["ALT_SPAN"]
    data["REF_SPAN_REL"] = np.where(total_span == 0, 0, data["REF_SPAN"] / total_span)
    data["ALT_SPAN_REL"] = np.where(total_span == 0, 0, data["ALT_SPAN"] / total_span)
    return data


def filter_by_zscore(data, features, remove_z):
    """Remove rows with |z scores| > remove_z"""
    return data[(np.abs(np.nan_to_num(zscore(data[features]))) < remove_z).all(axis=1)]


def genotype_vcf(
    args,
    input_vcf: str,
    input_sim: str,
    input_real: str,
    output_file=sys.stdout,
    remove_z=5.0,
):
    """Write new VCF with NSPV-determined genotypes

    Args:
        args (argparse.Namespace): Command arguments
        input_vcf (str): Path to input VCF file
        input_sim (str): Path to TSV of NPSV features for simulated data
        input_real (str): Path to TSV of NPSV features for real data
        output_file ([type], optional): File object for writing VCF. Defaults to sys.stdout.
    """
    # Load simulation and real data
    if not args.local and args.filter_bed is not None:
        # If provided, filter training data by BED file
        filter_bed = bed.BedTool(args.filter_bed)
        sim_bed = bed.BedTool(input_sim)
        sim_data = sim_bed.intersect(filter_bed, u=True, f=0.5).to_dataframe(
            header=None,
            na_values=".",
            names=get_tsv_columns(input_sim),
            dtype={"#CHROM": str, "SAMPLE": str, "AC": int},
        )

    else:
        sim_data = pd.read_table(
            input_sim, dtype={"#CHROM": str, "SAMPLE": str, "AC": int}
        )
    add_derived_features(sim_data)

    real_data = pd.read_table(input_real, dtype={"#CHROM": str, "SAMPLE": str})
    add_derived_features(real_data)

    # Filter to available features
    available_features = set(sim_data) & set(real_data)
    features = [feature for feature in FEATURES if feature in available_features]
    logging.info("Genotyping with features: %s", ", ".join(features))

    if not args.local:
        # Perform "global" genotyping by training on the entire simulated data
        downsample_sim_data = (
            sim_data.groupby(VAR_COL + [KLASS_COL])
            .head(args.downsample)
            .reset_index(drop=True)
        )

        # Drop rows with na, inf, etc. from training data
        with pd.option_context("mode.use_inf_as_na", True):
            downsample_sim_data = downsample_sim_data.dropna(axis=0, subset=features)

        # Expand here with additional classifiers
        logging.info(
            "Building global %s classifier based on simulated data", args.classifier
        )
        if args.classifier == "svm":
            pred = multi_svm_classify(downsample_sim_data, real_data, features=features)
        elif args.classifier == "rf":
            pred = multi_rf_classify(downsample_sim_data, real_data, features=features)
        else:
            raise NotImplementedError(f"Unknown classifier type: {args.classifier}")

        if logging.getLogger().isEnabledFor(logging.DEBUG) and KLASS_COL in real_data:
            logging.debug(
                "Accuracy compared to reported %s in 'real' data: %f",
                KLASS_COL,
                accuracy_score(real_data[KLASS_COL], pred),
            )

        # Report accuracy for Bayesian model as a comparison
        if logging.getLogger().isEnabledFor(logging.DEBUG) and KLASS_COL in real_data:
            accuracy = accuracy_score(
                real_data[KLASS_COL],
                np.argmax(
                    real_data[["PROB_HOMREF", "PROB_HET", "PROB_HOMALT"]].values, axis=1
                ),
            )
            logging.debug(
                "Accuracy for Bayesian model compared to reported %s in 'real' data: %f",
                KLASS_COL,
                accuracy,
            )
    else:
        # Prepare simulated data for per variant classifier
        grouped_sim_data = sim_data.groupby(VAR_COL)

    # Prepare real data to write out per-variant predictions
    grouped_real_data = real_data.groupby(VAR_COL)

    # Write new VCF
    # --------------------------------------
    vcf_reader = vcf.Reader(filename=input_vcf)  # Original VCF file

    # Add new fields to the header  (TODO: Add features as format field)
    vcf_reader.formats["GT"] = vcf.parser._Format("GT", 1, "String", "Genotype")
    vcf_reader.formats["GR"] = vcf.parser._Format(
        "GR",
        "R",
        "Integer",
        "Reads mapped by Paragraph to reference and alternate sequences",
    )

    # If original VCF is sites only...
    if len(vcf_reader._column_headers) < 9:
        vcf_reader._column_headers = VCF_COLUMN_HEADERS
    vcf_reader.samples = list(real_data[SAMPLE_COL].unique())  # Set sample names

    # Write new VCF entries, building local classifiers as needed
    vcf_writer = vcf.Writer(output_file, vcf_reader, lineterminator="")
    for record in vcf_reader:
        # Write sites-only and FORMAT columns (overwriting any original values)
        record.FORMAT = None
        record.samples = []
        vcf_writer.write_record(record)
        output_file.write(f"\t{VCF_FORMAT}")

        # Get prediction (we can't assume that the simulated data is in the same order as the VCF)
        for sample in vcf_reader.samples:
            group_vals = record_to_var_col(record, sample)

            try:
                real_group = grouped_real_data.get_group(group_vals)
                indices = grouped_real_data.indices[group_vals]
                if len(indices) > 1:
                    logging.warn(
                        "Multiple 'real' data rows for %s@%s. Skipping.",
                        sample,
                        variant_descriptor(record),
                    )
                    output_file.write("\t.")
                    continue

            except KeyError:
                logging.warn(
                    "No 'real' data found for %s@%s. Skipping.",
                    sample,
                    variant_descriptor(record),
                )
                output_file.write("\t.")
                continue

            assert len(indices) == 1, "Should only be only 'real data' entry"

            # Construct VCF call entry
            if not args.local:
                call = pred_to_vcf(real_group, pred[indices[0]])
            else:
                # Construct local classifier
                sim_group = grouped_sim_data.get_group(group_vals)

                # Drop rows with na, inf, etc. from training data
                with pd.option_context("mode.use_inf_as_na", True):
                    sim_group = sim_group.dropna(axis=0, subset=features)

                # Remove outliers from training data
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    sim_group = (
                        sim_group.groupby(KLASS_COL)
                        .apply(filter_by_zscore, features, remove_z)
                        .reset_index(drop=True)
                    )

                logging.debug("Classifying with %d observations", sim_group.shape[0])
                # If in debug mode, also print Mahalanobis distance
                if logging.getLogger().isEnabledFor(logging.DEBUG):
                    # Just use "absolute" features for Mahalanobis distance
                    _, _, mahal_score = single_mahalanobis(
                        sim_group, real_group, features=ABSOLUTE_FEATURES
                    )
                    logging.debug(
                        "%s@%s: Mahalanobis^2=%s",
                        sample,
                        variant_descriptor(record),
                        np.array2string(mahal_score.values, separator=","),
                    )

                # Classify variant and generate VCF genotype entry
                if args.classifier == "svm":
                    pred, _ = single_svm_classify(
                        sim_group, real_group, features=features
                    )
                elif args.classifier == "rf":
                    pred, _ = single_randomforest_classify(
                        sim_group, real_group, features=features
                    )
                else:
                    raise NotImplementedError(
                        f"Unknown classifier type: {args.classifier}"
                    )
                call = pred_to_vcf(real_group, pred)

            output_file.write("\t" + call)

        output_file.write("\n")  # Finish off the VCF line