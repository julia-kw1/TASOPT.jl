"""Plot SAD payload-range results."""

from pathlib import Path

import matplotlib
import pandas as pd


matplotlib.use("Agg")
import matplotlib.pyplot as plt

pt = 1.0 / 72.27
width = 345 * pt
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=plt.get_cmap("Dark2").colors)

def set_aiaa_style():
    plt.rcParams.update({
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 8,
        "lines.markersize": 3,
        "figure.autolayout": True,
        "savefig.dpi": 350,
        "axes.grid": True,
        "grid.linestyle": ":",
        "grid.linewidth": 0.4,
    })
    return plt


set_aiaa_style()


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "results" / "outputs"
IMAGE_DIR = ROOT / "results" / "images"

PAYLOAD_RANGE_SUFFIX = "_payload_range_comparison.csv"


def plot_configuration(csv_file):
    image_file = IMAGE_DIR / f"{csv_file.stem}.png"
    data = pd.read_csv(csv_file)
    title = data["title"].dropna().iloc[0]
    fig, ax = plt.subplots(figsize=(width, 0.72 * width))

    model_data = data[data["series_type"] == "model"]
    for (altitude_ft, mach), series in model_data.groupby(["altitude_ft", "mach"], sort=False):
        ax.plot(
            series["range_nmi"],
            series["oew_plus_payload_lb"],
            label=f"FL{altitude_ft / 100:.0f} / M{mach:.2f}",
            marker="o",
        )

    reference = data[data["series_type"] == "published_reference"]
    ax.plot(
        reference["range_nmi"],
        reference["oew_plus_payload_lb"],
        linestyle="--",
        marker="D",
        label="Published reference",
    )
    ax.set(
        xlabel="Range (nmi)",
        ylabel="OEW + Payload (lb)",
        title=f"{title} Payload-Range Comparison",
        xlim=(0, None),
    )
    ax.legend(loc="upper right")
    ax.yaxis.set_major_formatter("{x:,.0f}")
    image_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(image_file)
    plt.close(fig)
    print(f"Saved: {image_file.resolve()}")


def main():
    for csv_file in sorted(OUTPUT_DIR.glob(f"*{PAYLOAD_RANGE_SUFFIX}")):
        plot_configuration(csv_file)


if __name__ == "__main__":
    main()
