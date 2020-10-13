import logging, os, re, warnings
import vcf
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from npsv.feature_extraction import Features
from npsv.variant import variant_descriptor
from npsv.genotyper import add_derived_features, filter_by_zscore

FEATURE_COL = Features.FEATURES
VARIANT_COL = ["#CHROM", "START", "END", "TYPE"]

def plot_hist(data, col, colors, ax):
    sns.distplot(data.loc[data["AC"] == "REF",col], kde=False, color=colors[0], ax=ax)
    sns.distplot(data.loc[data["AC"] == "HET",col], kde=False, color=colors[1], ax=ax)
    sns.distplot(data.loc[data["AC"] == "HOM",col], kde=False, color=colors[2], ax=ax)
    ax.axvline(data.loc[data["AC"] == "Act",col].values[0], 0, 1, color=colors[3])

def plot_features(
    args, sim_path: str, real_path: str, vcf_path: str, out_dir_path: str
):
    """Generate pairwise plot of simulated and 'real' features
    
    Args:
        args (argparse.Namespace): Additional command line arguments
        sim_path (str): Path to NPSV features from 'simulated' data
        real_path (str): Path to NPSV features from 'real' data
        vcf_path (str): Path to input VCF file
        out_dir_path (str): Directory for plot files
    """
    # Create output directory if it doesn't exist
    os.makedirs(out_dir_path, exist_ok=True)
    logging.info("Generating plots in %s", out_dir_path)

    # Group the data to prepare for querying variants
    sim_data = pd.read_table(sim_path, dtype={"#CHROM": str, "AC": int})
    add_derived_features(sim_data)
    sim_data = sim_data.groupby(VARIANT_COL)

    real_data = pd.read_table(real_path, dtype={"#CHROM": str})
    add_derived_features(real_data)
    real_data = real_data.groupby(VARIANT_COL)

    # Depending on feature extractor, not all features may be available
    available_features = set(sim_data.obj) & set(real_data.obj)
    features = [feature for feature in FEATURE_COL if feature in available_features]

    vcf_reader = vcf.Reader(filename=vcf_path)
    for record in vcf_reader:
        variant = (
            record.CHROM,
            int(record.POS),
            int(record.sv_end),
            record.var_subtype,
        )

        try:
            current_sim = sim_data.get_group(variant)
            current_real = real_data.get_group(variant)
        except KeyError:
            # No data available for this variant, skipping
            logging.debug(
                "No simulated or real data found for %s. Skipping.",
                variant_descriptor(record),
            )
            continue
        current_real["AC"] = [-1]

        # Remove outliers with Z score above threshold
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            current_sim = (
                current_sim.groupby("AC")
                .apply(filter_by_zscore, features, 5)
                .reset_index(drop=True)
            )

        plot_data = current_sim.append(current_real)
        # Don't yet know how to encode AC directly (need strings for plotting)
        plot_data["AC"] = pd.Categorical(
            plot_data["AC"], categories=[0, 1, 2, -1]
        ).rename_categories(["REF", "HET", "HOM", "Act"])

        colors = sns.mpl_palette("Set1", 3) + [(0, 0, 0)]  # Actual data is black
        markers = { "REF": "o", "HET": "o", "HOM": "o", "Act": "s"}
        
        fig, ((ax11, ax12, ax13, ax14), (ax21, ax22, ax23, ax24)) = plt.subplots(2, 4, figsize=(14, 8))
        
        sns.scatterplot(ax=ax11, x="REF_READ", y="ALT_READ", data=plot_data, hue="AC", style="AC", markers=markers, palette=colors)
        sns.scatterplot(ax=ax12, x="REF_WEIGHTED_SPAN", y="ALT_WEIGHTED_SPAN", data=plot_data, hue="AC", style="AC", markers=markers, palette=colors)
        sns.scatterplot(ax=ax13, x="INSERT_LOWER", y="INSERT_UPPER", data=plot_data, hue="AC", style="AC", markers=markers, palette=colors)
        plot_hist(ax=ax14, col="CLIP_PRIMARY", data=plot_data, colors=colors)

        plot_hist(ax=ax21, col="COVERAGE", data=plot_data, colors=colors)
        plot_hist(ax=ax22, col="DHFC", data=plot_data, colors=colors)
        plot_hist(ax=ax23, col="DHBFC", data=plot_data, colors=colors)
        plot_hist(ax=ax24, col="DHFFC", data=plot_data, colors=colors)

        fig.suptitle("{}:{}-{}".format(*variant), size=16)
        fig.subplots_adjust(top=0.95, wspace=0.3, hspace=0.3)

        # Save plot to file name based on variant descriptor
        description = variant_descriptor(record)
        logging.info("Plotting variant into %s.png", description)
        plt.savefig(os.path.join(out_dir_path, f"{description}.png"))
