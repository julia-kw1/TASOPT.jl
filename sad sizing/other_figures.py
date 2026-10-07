"""Plot SAD payload-Mach and fuselage-radius sweep results."""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd


matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

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


# Presentation color palette
COLOR_ISO = "#EBF8FF"          # Soft blue for Isothermal Layer
COLOR_UPPER_STRATO = "#F7FAFC"
COLOR_CRUISE = "#BEE3F8"      # Mid-tone blue for cruise band
COLOR_BASELINE = "#4A5568"    # Executive slate grey for E175
COLOR_TARGET = "#FF6C1D"      # Accent red for target & gap arrow
COLOR_GAP_BG = "#FFF5F5"      # Light tint for gap shading
TEXT_MAIN = "#121D33"          # Main text color (dark slate)
COLOR_LINE_BOUNDARY = "#FFFFFF"   # Crisp white boundary lines like reference image
COLOR_TEXT_ATM = "#FFFFFF"        # White text for readability over dark upper gradient


ISOTHERMAL_LAYER = [36089, 65617] # Isothermal layer bounds in feet (11 km to 20 km)

def plot_altitude_visual():
    y_min, y_max = 29000, 72000
    e175_ceiling = 41000
    injection_altitude = 65000

    fig, ax = plt.subplots(figsize=(width, 0.85 * width))

    # 1. Continuous Blue Sky Gradient (Light troposphere -> Deep space blue upper stratosphere)
    gradient = np.linspace(0, 1, 256).reshape(256, 1)
    cmap = LinearSegmentedColormap.from_list("sky_gradient", ["#FFFFFF", "#7498BE", "#223D62"])
    ax.imshow(gradient, aspect='auto', cmap=cmap, extent=[0, 1.05, y_min, y_max], origin='lower', zorder=0)

    # Element Styling
    COLOR_GAP_BG = "#FF6B6B"

    # Isothermal / Upper Stratosphere Boundary (65,617 ft)
    ax.hlines(y=65617, xmin=0.0, xmax=1.0, color=COLOR_LINE_BOUNDARY, linestyle="--", linewidth=0.8, alpha=0.9, zorder=2)
    ax.hlines(y=36089, xmin=0.0, xmax=1.0, color=COLOR_BASELINE, linestyle="--", linewidth=0.8, alpha=0.5, zorder=2)

    # Layer Labels inside gradient regions
    ax.text(0.83, 37600, "Isothermal\n     Layer", ha='left', va='center', alpha=0.5, fontsize=7.5, fontstyle='italic', color=COLOR_BASELINE, fontweight='bold', zorder=3)
    ax.text(0.83, 67500, "      Upper\nStratosphere", ha='left', va='center', fontsize=7.5, fontstyle='italic', color=COLOR_TEXT_ATM, fontweight='bold', zorder=3)

    # 3. Airline Cruise Altitude Band (35,000 ft to 39,000 ft)
    ax.axhspan(35000, 39000, color="#ACB5BF", alpha=0.45, zorder=2)
    ax.text(0.42, 37000, "Airline Cruise Altitude (35k-39k ft)", ha='center', va='center', 
            fontweight='bold', fontsize=7, color="#1A202C", zorder=4)

    # 4. Capability Gap Tint Highlight
    ax.axhspan(e175_ceiling, injection_altitude, xmin=0.1, xmax=0.69, color=COLOR_GAP_BG, alpha=0.20, zorder=2)

    # 5. Baseline E175 Ceiling Line
    ax.hlines(y=e175_ceiling, xmin=0.1, xmax=0.7, color=COLOR_BASELINE, linestyle="-", linewidth=2, zorder=5)
    ax.text(0.425, e175_ceiling + 1000, f"E175 Ceiling ({e175_ceiling:,} ft)", ha='center', 
            va='bottom', fontweight='bold', fontsize=7.5, color=COLOR_BASELINE, zorder=6)

    # 6. Injection Target Altitude Line
    ax.hlines(y=injection_altitude, xmin=0.1, xmax=0.7, color=COLOR_TARGET, linestyle="--", linewidth=2, zorder=5)
    ax.text(0.425, injection_altitude - 1200, f"Injection Target ({injection_altitude:,} ft)", ha='center',
            va='top', fontweight='bold', fontsize=7.5, color=COLOR_TARGET, zorder=6)

    # 7. Gap Bracket & Rotated Label
    x_bracket = 0.72
    gap_ft = (injection_altitude - e175_ceiling) / 1000
    ax.annotate(
        "",
        xy=(x_bracket, injection_altitude),
        xytext=(x_bracket, e175_ceiling),
        arrowprops=dict(arrowstyle="<->", color=COLOR_TARGET, lw=1.5, mutation_scale=10),
        zorder=6,
    )
    ax.text(
        x_bracket + 0.03,
        (e175_ceiling + injection_altitude) / 2,
        f"~{gap_ft:,.0f}k ft Gap",
        ha="left",
        va="center",
        fontweight="bold",
        fontsize=8,
        color=COLOR_TARGET,
        zorder=6,
    )

    # Primary Y-axis formatting
    ax.set_ylim(y_min, y_max)
    ax.set_ylabel("Altitude (ft)", fontweight="bold", color="#1A202C")
    ax.set_yticks(np.arange(30000, 80000, 10000))
    ax.yaxis.set_major_formatter("{x:,.0f}")
    ax.set_xlim(0, 1.05)
    ax.get_xaxis().set_visible(False)

    # Clean Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    output_path = IMAGE_DIR / "altitude_visual_ladder.png"
    fig.savefig(output_path, dpi=350, bbox_inches="tight")
    plt.close(fig)


def plot_mission_profile():
    # takeoff, climb, cruise, descent, landing
    x_vals = np.array([0, 0.05, 0.22, 0.78, 0.95, 1.0])
    y_vals = np.array([0, 0, 65000, 65000, 0, 0])
    y_min, y_max = 10000, 72000

    fig, ax = plt.subplots(figsize=(0.85*width, 0.42 * width))

    # fill for stratsophere
    ax.axhspan(ISOTHERMAL_LAYER[1], 80000, color=COLOR_CRUISE, alpha=0.45, zorder=1)
    ax.text(0.01, 69000, "Upper Stratosphere", fontsize=6.5, fontstyle="italic", color="#828F9F", zorder=2)

    # fill under mission
    ax.fill_between(x_vals, y_vals, color="#E2E8F0", alpha=0.35, zorder=2)

    # plot mission profile line
    ax.plot(x_vals, y_vals, linewidth=1.8, color=COLOR_BASELINE, zorder=4)

    # highlight cruise segment
    ax.plot([0.22, 0.78], [65000, 65000], color=COLOR_TARGET, linewidth=2.4, zorder=5)

    # hollow waypoint markers
    ax.scatter(x_vals[1:-1], y_vals[1:-1], color="white", s=14, zorder=6, edgecolor=COLOR_BASELINE, lw=1.2)

    ax.text(
        0.50,
        49000,
        "Cruise + Dispersion \nTarget: 65,000 ft",
        ha="center",
        va="center",
        fontsize=7.5,
        fontweight="bold",
        color=COLOR_TARGET,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFF5F5", edgecolor=COLOR_TARGET, lw=0.8),
        zorder=7,
    )

    # Phase text labels along profile legs
    ax.text(0.135, 34000, "Climb", ha="right", va="center", fontsize=7, fontweight="bold", color=TEXT_MAIN, rotation=65)
    ax.text(0.865, 34000, "Descent", ha="left", va="center", fontsize=7, fontweight="bold", color=TEXT_MAIN, rotation=-65)

    # Takeoff & Landing labels
    ax.text(0.00, -7000, "Takeoff", ha="left", va="top", fontsize=7, fontweight="bold", color=TEXT_MAIN)
    ax.text(1.00, -7000, "Landing", ha="right", va="top", fontsize=7, fontweight="bold", color=TEXT_MAIN)

    # Y-Axis & Frame Formatting
    ax.set_ylim(-10000, 78000)
    ax.set_ylabel("Altitude (ft)", fontweight="bold", fontsize=8)
    ax.set_yticks([0, 30000, 65000])
    ax.yaxis.set_major_formatter("{x:,.0f}")

    # Remove X-axis line & ticks
    ax.set_xlim(-0.02, 1.02)
    ax.get_xaxis().set_visible(False)

    # Clean spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    
    output_path = IMAGE_DIR / "mission_profile.png"
    fig.savefig(output_path, dpi=350, bbox_inches="tight")
    plt.close(fig)


def main():
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    plot_altitude_visual()
    plot_mission_profile()

if __name__ == "__main__":
    main()