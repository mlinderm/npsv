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
from sklearn.linear_model import LogisticRegression
from scipy.special import logsumexp
from sklearn.ensemble import RandomForestClassifier
from sklearn.covariance import MinCovDet
import scipy.spatial.distance as distance
from scipy.stats import binom, chi2, zscore
from sklearn.metrics import accuracy_score
from scipy.special import comb
from tqdm import tqdm
import xgboost as xgb

from .variant import Variant, variant_descriptor, overwrite_reader_samples
from .feature_extraction import Features

# Suppress the future warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
if logging.getLogger(__name__).getEffectiveLevel() != logging.DEBUG:
    warnings.simplefilter(action="ignore", category=UserWarning)

GENOTYPING_MODES = ["single", "variant", "hybrid"]
CLASSIFIERS = ["svm", "rf", "logistic", "xgboost"]

VAR_COL = ["#CHROM", "START", "END", "TYPE", "SAMPLE"]
FEATURE_COL = Features.FEATURES
KLASS_COL = "AC"
SAMPLE_COL = "SAMPLE"
TYPE_COL = "TYPE"
AD_COL = ["REF_READ", "ALT_READ"]

MAHAL_FEATURES = [
    "INSERT_LOWER",
    "INSERT_UPPER",
    "DHFFC",
    "REF_READ_REL",
    "ALT_READ_REL",
    "REF_WEIGHTED_SPAN_REL",
    "ALT_WEIGHTED_SPAN_REL",
]


CLASSIFIER_FEATURES = {
    "svm": [
        "SVLEN",
        "REF_READ_REL",
        "ALT_READ_REL",
        "REF_WEIGHTED_SPAN_REL",
        "ALT_WEIGHTED_SPAN_REL",
        "INSERT_LOWER",
        "INSERT_UPPER",
        "DHFC",
        "DHBFC",
        "DHFFC",
        "PROB_HOMREF",
        "PROB_HET",
        "PROB_HOMALT",
        "CLIP_PRIMARY",
    ],
    "rf": [
        "SVLEN",
        "ALT_READ_REL",
        "ALT_WEIGHTED_SPAN_REL",
        "INSERT_LOWER",
        "INSERT_UPPER",
        "DHFFC",
        "PROB_HOMREF",
        "PROB_HET",
        "PROB_HOMALT",
        "CLIP_PRIMARY",
    ],
    "logistic": [
        "SVLEN",
        "ALT_READ_REL",
        "ALT_WEIGHTED_SPAN_REL",
        "INSERT_LOWER",
        "INSERT_UPPER",
        "DHFFC",
        "PROB_HOMREF",
        "PROB_HET",
        "PROB_HOMALT",
        "CLIP_PRIMARY",
    ],
    "xgboost": [
        "SVLEN",
        "ALT_READ_REL",
        "ALT_WEIGHTED_SPAN_REL",
        "INSERT_LOWER",
        "INSERT_UPPER",
        "DHFFC",
        "PROB_HOMREF",
        "PROB_HET",
        "PROB_HOMALT",
        "CLIP_PRIMARY",
    ],
}


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
VCF_FORMAT = "GT:AD:PL:DM"
AC_TO_GT = ("0/0", "0/1", "1/1")

SVM_PARAM_GRID = [
    {"C": [1, 10, 100, 500, 1000, 5000, 10000], "gamma": ["scale", 0.001, 0.0055, .01, .055, 0.1, 0.55], "kernel": ['rbf']},
]

RF_PARAM_GRID = [
    {"n_estimators": [50, 100, 500, 1000], "max_depth": [None, 3, 5, 10], "min_samples_leaf": [1, 3, 5], "min_samples_split": [2, 5, 10]}
]

LOGISTIC_PARAM_GRID = [
    {"penalty": ["l2"], "C": [0.001, 0.01, 0.1, 1, 10, 100, 1000] }
]

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


def pred_to_vcf(real_data, pred, prob=None, dm2=None, ad=None) -> str:
    """Prediction to VCF call field for sample"""
    assert real_data.shape[0] == 1, "Real data should have just one variant"

    # Compute PL
    if prob is not None:
        pl = -10.0 * np.log10(prob + 1e-11)
        pl = np.clip(np.round(pl - np.min(pl)).astype(int), 0, 99)

    return "{gt}:{ad}:{pl}:{md}".format(
        gt=AC_TO_GT[pred] if pred in (0, 1, 2) else "./.",
        ad=",".join(map(str, real_data[AD_COL].to_numpy().squeeze().astype(int))),
        pl=",".join(map(str, pl)) if prob is not None else ".",
        md=",".join(map(lambda x: str(round(x, 1)), dm2)) if dm2 is not None else ".",
    )


def rf_classify(sim_data, real_data, features=FEATURE_COL, klass=KLASS_COL, param_search=False, threads=1, n_estimators=100, **kwargs):
    x = sim_data[features]
    y = sim_data[klass]

    # Fit the model
    clf = RandomForestClassifier(n_estimators=n_estimators, max_features="auto", **kwargs)
    
    if param_search:
        # Grid-search to find the best parameters
        clf = GridSearchCV(clf, param_grid=RF_PARAM_GRID, cv=5, n_jobs=threads)
    
    clf = clf.fit(x, y)

    if param_search:
        logging.info("Best RF params: %s", clf.best_params_)

    # Predict the real data
    real_x = real_data[features]
    pred = clf.predict(real_x)
    prob = clf.predict_proba(real_x)
    return (pred, prob)


def svm_classify(sim_data, real_data, features=FEATURE_COL, klass=KLASS_COL, param_search=False, gamma="scale", threads=1, **kwargs):
    x = sim_data[features]
    y = sim_data[klass]

    # Normalize the training set
    scaler = StandardScaler()
    x = scaler.fit_transform(x)

    # Build the model
    if gamma != "scale" and gamma != "auto":
        gamma = float(gamma)
    clf = SVC(kernel="rbf", probability=True, gamma=gamma, **kwargs)

    if param_search:
        # Grid-search to find the best parameters
        clf = GridSearchCV(clf, param_grid=SVM_PARAM_GRID, cv=5, n_jobs=threads)

    # Fit the model
    clf.fit(x, y)
    
    if param_search:
        logging.info("Best SVM params: %s", clf.best_params_)

    # Predict the real data
    real_x = scaler.transform(real_data[features])
    pred = clf.predict(real_x)
    prob = clf.predict_proba(real_x)
    return (pred, prob)


def single_mahalanobis(sim_data, real_data, features=FEATURE_COL, klass=KLASS_COL):
    """Classify single variant using Mahalanobis distance"""
    assert real_data.shape[0] == 1, "Real data should have just one variant"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        x = sim_data[features]
        y = sim_data[klass]
        real_x = real_data[features]

        def mahal_score(data):
            if data.shape[0] < 2:
                # Insufficient data for calculation
                return float("inf")
            robust_cov = MinCovDet(assume_centered=False).fit(data)
            return robust_cov.mahalanobis(real_x)[0]

        score = x.groupby(y).apply(mahal_score)
        assert len(score) == 3, "Missing classes in Mahalanobis distance calculation"
        pred = score.idxmin()
        with np.errstate(under="ignore"):
            prob = chi2.logsf(score, len(features))
            prob = np.exp(prob - logsumexp(prob))

        return (pred, prob, score)


def logistic_classify(sim_data, real_data, features=FEATURE_COL, klass=KLASS_COL, param_search=False, threads=1, **kwargs):
    x = sim_data[features]
    y = sim_data[klass]

    # Normalize the training set
    scaler = StandardScaler()
    x = scaler.fit_transform(x)
  
    # Build the model
    clf = LogisticRegression(**kwargs)

    if param_search:
        # Grid-search to find the best parameters
        clf = GridSearchCV(clf, param_grid=LOGISTIC_PARAM_GRID, cv=5, n_jobs=threads)

    # Fit the model
    clf.fit(x, y)
    
    if param_search:
        logging.info("Best Logistic params: %s", clf.best_params_)

    # Predict the real data
    real_x = scaler.transform(real_data[features])
    pred = clf.predict(real_x)
    prob = clf.predict_proba(real_x)
    return (pred, prob)


def xgboost_classify(sim_data, real_data, features=FEATURE_COL, klass=KLASS_COL, param_search=False, threads=1, **kwargs):
    x = sim_data[features]
    y = sim_data[klass]
  
    # Build the model
    clf = xgb.XGBClassifier(n_estimators=100, n_jobs=threads, **kwargs)

    # Fit the model
    clf.fit(x, y)

    # Predict the real data
    real_x = real_data[features]
    pred = clf.predict(real_x)
    prob = clf.predict_proba(real_x)
    return (pred, prob)


def add_derived_features(data):
    """Update data with derived features"""
    # Binomial probability features adapted from svviz2 and
    # https://www.sciencedirect.com/science/article/pii/S0092867418316337#sec4

    # Round since we can get fraction read counts when there are multiple breakpoints
    ref_reads = np.ceil(data["REF_READ"])
    alt_reads = np.ceil(data["ALT_READ"])
    total_reads = ref_reads + alt_reads
    data["REF_READ_REL"] = np.where(total_reads == 0, 0, ref_reads / total_reads)
    data["ALT_READ_REL"] = np.where(total_reads == 0, 0, alt_reads / total_reads)

    data["PROB_HOMREF"] = binom.pmf(alt_reads, total_reads, 0.05)
    data["PROB_HET"] = binom.pmf(alt_reads, total_reads, 0.5)
    data["PROB_HOMALT"] = binom.pmf(alt_reads, total_reads, 0.95)

    total_span = data["REF_WEIGHTED_SPAN"] + data["ALT_WEIGHTED_SPAN"]
    data["REF_WEIGHTED_SPAN_REL"] = np.where(total_span == 0, 0, data["REF_WEIGHTED_SPAN"] / total_span)
    data["ALT_WEIGHTED_SPAN_REL"] = np.where(total_span == 0, 0, data["ALT_WEIGHTED_SPAN"] / total_span)
    return data


def filter_by_zscore(data, features, remove_z):
    """Remove rows with |z scores| > remove_z"""
    return data[(np.abs(np.nan_to_num(zscore(data[features]), posinf=0.0, neginf=0.0)) < remove_z).all(axis=1)]


def genotype_vcf(
    args,
    input_vcf: str,
    input_sim: str,
    input_real: str,
    output_file=sys.stdout,
    remove_z=5.0,
    samples=[],
):
    """Write new VCF with NPSV-determined genotypes

    Args:
        args (argparse.Namespace): Command arguments
        input_vcf (str): Path to input VCF file
        input_sim (str): Path to TSV of NPSV features for simulated data
        input_real (str): Path to TSV of NPSV features for real data
        output_file ([type], optional): File object for writing VCF. Defaults to sys.stdout.
    """
    # Load simulation and real data
    if args.filter_bed is not None:
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
            input_sim, na_values=".", dtype={"#CHROM": str, "SAMPLE": str, "AC": int}
        )

    if sim_data.shape[0] == 0:
        # No data is available, copy input to output and exit
        vcf_reader = vcf.Reader(filename=input_vcf)
        if samples is not None:
            overwrite_reader_samples(vcf_reader, samples)
        vcf.Writer(output_file, vcf_reader)
        if not vcf_reader._reader.closed:
            vcf_reader._reader.close()
        return

    real_data = pd.read_table(
        input_real, na_values=".", dtype={"#CHROM": str, "SAMPLE": str}
    )

    # Add derived features
    add_derived_features(sim_data)
    add_derived_features(real_data)

    # Lines to add to VCF header
    classifier_vcf_metadata = {}

    if not ({args.DEL_gt_mode, args.INS_gt_mode}).isdisjoint({"single", "hybrid"}):
        single_pred = np.full(real_data.shape[0], -1, dtype=int)
        single_prob = np.full((real_data.shape[0], 3), 0, dtype=float)

        # Construct classifiers for each variant type
        typed_sim_data = sim_data.groupby(TYPE_COL)
        typed_real_data = real_data.groupby(TYPE_COL)
        for variant_type, sim_group in typed_sim_data:
            if getattr(args, f"{variant_type}_gt_mode", "single") == "single":
                single_classifier = getattr(args, f"{variant_type}_classifier", "rf")
            else:
                single_classifier = getattr(args, f"{variant_type}_hybrid_classifier")

            real_group = typed_real_data.get_group(variant_type)
            sim_group = (
                sim_group.groupby(VAR_COL + [KLASS_COL])
                .sample(args.downsample)
                .reset_index(drop=True)
            )

            with pd.option_context("mode.use_inf_as_na", True):
                # Drop any columns that are entirely NA
                sim_group = sim_group.dropna(axis=1, how="all")
                real_group = real_group.dropna(axis=1, how="all")

                # Filter to available features
                desired_features = CLASSIFIER_FEATURES[single_classifier]
                single_features = list(
                    set(desired_features) & set(sim_group) & set(real_group)
                )

                # Then drop rows with na, inf, etc. from training data
                sim_group = sim_group.dropna(axis=0, subset=single_features)

            # Expand here with additional classifiers
            logging.info(
                "Building 'single model' %s classifier for %s variants (%d observations) with features: %s",
                single_classifier,
                variant_type,
                sim_group.shape[0],
                ", ".join(single_features),
            )
            classifier_vcf_metadata[f"npsv_{variant_type}_single_classifier"] = [
                f"{single_classifier}({','.join(single_features)})"
            ]
            if single_classifier == "svm":
                pred, prob = svm_classify(
                    sim_group, real_group, features=single_features, param_search=args.param_search, threads=args.threads, gamma=args.svm_gamma, C=args.svm_C,
                )
            elif single_classifier == "rf":
                pred, prob = rf_classify(
                    sim_group, real_group, features=single_features, param_search=args.param_search, threads=args.threads, n_estimators=args.rf_n_estimators, max_depth=args.rf_max_depth,
                )
            elif single_classifier == "logistic":
                pred, prob = logistic_classify(
                    sim_group, real_group, features=single_features, param_search=args.param_search, threads=args.threads, penalty=args.lr_penalty, C=args.lr_C,
                )
            elif single_classifier == "xgboost":
                pred, prob = xgboost_classify(
                    sim_group, real_group, features=single_features, param_search=args.param_search, threads=args.threads,
                )

            # Reconstruct original vectors
            indices = typed_real_data.indices[variant_type]
            single_pred[indices] = pred
            if prob.shape[1] == 3:
                single_prob[indices] = prob

        # Report accuracy if actual class is defined
        if logging.getLogger().isEnabledFor(logging.DEBUG) and KLASS_COL in real_data:
            logging.debug(
                "Accuracy compared to reported %s in 'real' data: %f",
                KLASS_COL,
                accuracy_score(real_data[KLASS_COL], single_pred),
            )
            bayesian_accuracy = accuracy_score(
                real_data[KLASS_COL],
                np.argmax(
                    real_data[["PROB_HOMREF", "PROB_HET", "PROB_HOMALT"]].values,
                    axis=1,
                ),
            )
            logging.debug(
                "Accuracy for Bayesian model compared to reported %s in 'real' data: %f",
                KLASS_COL,
                bayesian_accuracy,
            )

    if not ({args.DEL_gt_mode, args.INS_gt_mode}).isdisjoint({"variant", "hybrid"}):
        if args.filter_bed is not None:
            # Clunky, but for the variant model we need unfiltered training data
            sim_data = pd.read_table(
                input_sim, na_values=".", dtype={"#CHROM": str, "SAMPLE": str, "AC": int}
            )
            add_derived_features(sim_data)
        
        # Prepare simulated data for per-variant classifier
        grouped_sim_data = sim_data.groupby(VAR_COL)

        # Filter to available features for per-variant classifier
        pervariant_features = {}
        for kind in ("DEL", "INS"):
            if getattr(args, f"{kind}_gt_mode") not in ("variant", "hybrid"):
                continue
            variant_classifier = getattr(args, f"{kind}_classifier")
            #desired_features = set(CLASSIFIER_FEATURES[variant_classifier]) - set(("SVLEN",))
            desired_features = set(CLASSIFIER_FEATURES[variant_classifier])
            pervariant_features[kind] = list(
                desired_features & set(sim_data) & set(real_data)
            )
            logging.info(
                "Building 'per-variant' %s classifiers for %s variants based on simulated data with features: %s",
                variant_classifier,
                kind,
                ", ".join(pervariant_features[kind]),
            )
            classifier_vcf_metadata[f"npsv_{kind}_variant_classifier"] = [
                f"{variant_classifier}({','.join(pervariant_features[kind])})"
            ]

    # Prepare real data to write out per-variant predictions
    grouped_real_data = real_data.groupby(VAR_COL)

    # Write new VCF
    # --------------------------------------
    vcf_reader = vcf.Reader(filename=input_vcf)  # Original VCF file

    # Add new fields to the header
    vcf_reader.metadata.update(classifier_vcf_metadata)
    vcf_reader.metadata["npsv_dm2"] = [f"mahal({','.join(MAHAL_FEATURES)})"]
    vcf_reader.formats["GT"] = vcf.parser._Format("GT", 1, "String", "Genotype")
    vcf_reader.formats["DM"] = vcf.parser._Format(
        "DM", "G", "Float", "Mahalanobis distance for each genotype",
    )
    vcf_reader.formats["PL"] = vcf.parser._Format(
        "PL", "G", "Integer", "Phred-scaled genotype likelihoods",
    )
    vcf_reader.formats["AD"] = vcf.parser._Format(
        "AD", "R", "Integer", "Read depth for each allele",
    )

    # If original VCF is sites only...
    if len(vcf_reader._column_headers) < 9:
        vcf_reader._column_headers = VCF_COLUMN_HEADERS
    # Set sample names
    overwrite_reader_samples(
        vcf_reader, list(set(real_data[SAMPLE_COL]) | set(samples))
    )

    # Write new VCF entries, building local classifiers as needed
    vcf_writer = vcf.Writer(output_file, vcf_reader, lineterminator="")
    for record in tqdm(vcf_reader, desc="Genotyping variants"):
        variant = Variant.from_pyvcf(record)

        # Write sites-only and FORMAT columns (overwriting any original or potentially invalidated values)
        record.INFO.pop("AC", None)
        record.INFO.pop("AN", None)
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
            indices = indices[0]

            # Construct VCF call entry
            variant_type = record.var_subtype
            gt_mode = getattr(args, f"{variant_type}_gt_mode", "single")  
            if gt_mode == "single" or (gt_mode == "hybrid"
                and variant.event_length >= getattr(args, f"{variant_type}_hybrid_threshold", 1000)
            ):
                call = pred_to_vcf(
                    real_group, single_pred[indices], single_prob[indices], ad=real_group[AD_COL].to_numpy().squeeze()
                )
            else:
                # Construct local classifier
                sim_group = grouped_sim_data.get_group(group_vals)
                if args.variant_downsample:
                    sim_group = (
                        sim_group.groupby(VAR_COL + [KLASS_COL])
                        .head(args.variant_downsample)
                        .reset_index(drop=True)
                    )
                
                with pd.option_context("mode.use_inf_as_na", True), warnings.catch_warnings():
                    warnings.simplefilter("ignore")

                    # Drop any columns that are entirely NA
                    sim_group = sim_group.dropna(axis=1, how="all")
                    real_group = real_group.dropna(axis=1, how="all")

                    # Update the available features
                    avail_features = list(
                        set(pervariant_features[variant_type]) & set(sim_group) & set(real_group)
                    )

                    # Then drop rows with na, inf, etc. from training data
                    sim_group = sim_group.dropna(axis=0, subset=avail_features)

                    # Remove outliers from training data
                    rows_before = sim_group.shape[0]
                    sim_group = (
                        sim_group.groupby(KLASS_COL)
                        .apply(filter_by_zscore, avail_features, remove_z)
                        .reset_index(drop=True)
                    )
                    # TODO: Check that all three classes are still present
                    logging.debug(
                        "Classifying with %d observations after removing %d outliers",
                        sim_group.shape[0], rows_before - sim_group.shape[0],
                    )

                mahal_score = None
                if args.dm2:
                    try:
                        _, _, mahal_score = single_mahalanobis(
                            sim_group,
                            real_group,
                            features=list(set(avail_features) & set(MAHAL_FEATURES)),
                        )
                    except:
                        pass

                try:
                    variant_classifier = getattr(args, f"{variant_type}_classifier")
                    if variant_classifier == "svm":
                        pred, prob = svm_classify(
                            sim_group, real_group, features=avail_features, gamma=args.svm_gamma, C=args.svm_C,
                        )
                    elif variant_classifier == "rf":
                        pred, prob = rf_classify(
                            sim_group, real_group, features=avail_features, n_estimators=args.rf_n_estimators, max_depth=args.rf_max_depth,
                        )
                    elif variant_classifier == "logistic":
                        pred, prob = logistic_classify(
                            sim_group, real_group, features=avail_features, penalty=args.lr_penalty, C=args.lr_C,
                        )
                    elif variant_classifier == "xgboost":
                        pred, prob = xgboost_classify(
                            sim_group, real_group, features=avail_features,
                        )
                    call = pred_to_vcf(
                        real_group, pred.item(0), prob[0,], dm2=mahal_score, ad=real_group[AD_COL].to_numpy().squeeze()
                    )
                except ValueError as e:
                    logging.error(
                        "Genotyping error for %s: %s", variant_descriptor(record), e
                    )
                    call = "./."

            output_file.write("\t" + call)

        output_file.write("\n")  # Finish off the VCF line
