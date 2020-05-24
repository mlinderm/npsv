import logging, os, re, warnings
import vcf
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from npsv.variant import variant_descriptor

FEATURE_COL = ["REF_SPAN", "ALT_SPAN", "REF_SPLIT", "ALT_SPLIT", "COVG"]
VARIANT_COL = ["#CHROM", "START", "END", "TYPE"]


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
    sim_data = pd.read_table(sim_path, dtype={"#CHROM": str, "AC": int}).groupby(
        VARIANT_COL
    )
    real_data = pd.read_table(real_path, dtype={"#CHROM": str}).groupby(VARIANT_COL)

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
                .apply(
                    lambda g: g[
                        (np.abs(np.nan_to_num(stats.zscore(g[features]))) < 5).all(
                            axis=1
                        )
                    ]
                )
                .reset_index(drop=True)
            )

        plot_data = current_sim.append(current_real)
        # Don't yet know how to encode AC directly (needs to strings for plotting)
        plot_data["AC"] = pd.Categorical(
            plot_data["AC"], categories=[0, 1, 2, -1]
        ).rename_categories(["REF", "HET", "HOM", "Act"])

        sns.set(style="ticks", color_codes=True)
        colors = sns.mpl_palette("Set1", 3) + [(0, 0, 0)]  # Actual data is black
        fig = sns.pairplot(
            plot_data,
            vars=features,
            hue="AC",
            diag_kind="hist",
            markers=["o", "o", "o", "x"],
            palette=colors,
        )
        fig.add_legend()
        fig.fig.suptitle("{}:{}-{}".format(*variant), size=16)
        fig.fig.subplots_adjust(top=0.95)

        # Save plot to file name based on variant descriptor
        description = variant_descriptor(record)
        logging.info("Plotting variant into %s.png", description)
        plt.savefig(os.path.join(out_dir_path, variant_descriptor(record) + ".png"))