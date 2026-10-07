"""Plot SAD payload-Mach and fuselage-radius sweep results."""

from pathlib import Path

import matplotlib
import numpy as np
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
RADIUS_CSV = OUTPUT_DIR / "e175_cfm56_radius_payload_ceiling_sweep.csv"
PAYLOAD_MACH_SUFFIX = "_payload_mach_sweep.csv"


def contour_data(data, x_column, y_column, value_column):
    x = np.sort(data[x_column].unique())
    y = np.sort(data[y_column].unique())
    values = data.pivot(index=y_column, columns=x_column, values=value_column)
    return x, y, values.loc[y, x].to_numpy()


def rounded_limits(values, step):
    return (np.floor(np.min(values) / step) * step, np.ceil(np.max(values) / step) * step)


def rounded_ticks(values, step):
    lower, upper = rounded_limits(values, step)
    return np.arange(lower, upper + step, step)


def endpoint_ticks(values, step):
    """Rounded interior ticks with the sampled endpoints shown exactly."""
    lower = float(np.min(values))
    upper = float(np.max(values))
    interior = rounded_ticks(values, step)
    interior = interior[(interior > lower) & (interior < upper)]
    return np.unique(np.concatenate(([lower], interior, [upper])))


def save_figure(fig, filename):
    image_file = IMAGE_DIR / filename
    fig.savefig(image_file)
    print(f"Saved: {image_file.resolve()}")
    plt.close(fig)


def plot_title(configuration_title, metric):
    if configuration_title.startswith("R_f = "):
        radius = configuration_title.split()[2]
        configuration_title = rf"$R_f = {radius}\ \mathrm{{in.}}$"
    return f"{configuration_title} {metric}"


def plot_payload_mach_contour(data, value_column, metric, colorbar_label, filename):
    payload_column = "so2_payload_lb" if "so2_payload_lb" in data else "payload_lb"
    payload, mach, values = contour_data(data, payload_column, "mach", value_column)
    title = plot_title(data["title"].iloc[0], metric)

    fig, ax = plt.subplots(figsize=(width, 0.72 * width))
    contours = ax.contourf(payload, mach, values, levels=12)
    colorbar = fig.colorbar(contours, ax=ax)
    colorbar.set_label(colorbar_label)
    ax.set(
        xlabel="SO$_2$ payload (lb)" if payload_column == "so2_payload_lb" else "Payload (lb)",
        ylabel="Cruise Mach number",
        title=title,
        ylim=(np.min(mach), np.max(mach)),
    )
    save_figure(fig, filename)


def plot_payload_mach_sweep(csv_file):
    data = pd.read_csv(csv_file)
    configuration = csv_file.stem.removesuffix("_payload_mach_sweep")
    plots = (
        ("ceiling_ft", "Service Ceiling", "Service ceiling (ft)", "absolute_ceiling"),
        ("cruise_cl", "$C_L$ at Service Ceiling", "$C_L$", "cl_at_ceiling"),
        ("specific_excess_power_fpm", "Specific Excess Power", "Specific excess power (ft/min)",
         "ps_35000ft"),
    )
    for value_column, metric, colorbar_label, output_suffix in plots:
        plot_payload_mach_contour(
            data,
            value_column,
            metric,
            colorbar_label,
            f"{output_suffix}__{configuration}.png",
        )


def plot_radius_ceiling(data):
    payload_column = "so2_payload_lb" if "so2_payload_lb" in data else "payload_lb"
    capacity_column = (
        "so2_payload_capacity_lb"
        if "so2_payload_capacity_lb" in data
        else "payload_capacity_lb"
    )
    radius, payload, ceiling = contour_data(data, "radius_in", payload_column, "ceiling_ft")
    radius_data = data.groupby("radius_in", sort=True).first()
    capacity = radius_data[capacity_column]
    if "selected_radius" in radius_data:
        selected = radius_data[radius_data["selected_radius"]]
        selected_radius = selected.index[0]
    else:
        eligible = (
            radius_data[radius_data["design_payload_feasible"]]
            if "design_payload_feasible" in radius_data
            else radius_data
        )
        selected_radius = eligible["mtow_lb"].idxmin()
    reference_payload = (
        data["design_so2_payload_lb"].iloc[0]
        if "design_so2_payload_lb" in data
        else data["excess_thrust_reference_payload_lb"].iloc[0]
    )

    fig, ax = plt.subplots(figsize=(width, 0.72 * width))
    contours = ax.contourf(radius, payload, ceiling, levels=12)
    colorbar = fig.colorbar(contours, ax=ax)
    colorbar.set_label("Service ceiling (ft)")
    ax.plot(radius, capacity.loc[radius], color="black",
            label="Maximum SO$_2$ capacity")
    ax.scatter(
        selected_radius,
        reference_payload,
        marker="o",
        s=70,
        edgecolor="black",
        zorder=3,
        label=f"Selected radius ({selected_radius:.0f} in)",
    )
    ax.set(
        xlabel="Fuselage radius (in)",
        ylabel="SO$_2$ payload (lb)",
        title="E175 with CFM56 Service Ceiling",
    )
    ax.set_xticks(endpoint_ticks(radius, 5.0))
    ax.set_xlim(np.min(radius), np.max(radius))
    ax.set_yticks(rounded_ticks(payload, 2_500.0))
    ax.set_ylim(rounded_limits(payload, 2_500.0)[0], np.max(payload))
    ax.legend(loc="upper left")
    save_figure(fig, "absolute_ceiling__e175_cfm56_radius_payload.png")


def plot_mass_trade(data):
    mass_data = data.groupby("radius_in", sort=True).first()
    radius = mass_data.index.to_numpy()
    capacity_column = (
        "so2_payload_capacity_lb"
        if "so2_payload_capacity_lb" in mass_data
        else "payload_capacity_lb"
    )

    fig, (mtow_ax, capacity_ax) = plt.subplots(2, 1, figsize=(width, 0.72 * width))
    mtow_ax.plot(radius, mass_data["mtow_lb"])
    capacity_ax.plot(radius, mass_data[capacity_column])
    mtow_ax.set_title("Maximum Takeoff Weight")
    capacity_ax.set_title("Maximum SO$_2$ Payload Capacity")

    for ax, column, ylabel, step in (
        (mtow_ax, "mtow_lb", "MTOW (lb)", 2_500.0),
        (capacity_ax, capacity_column, "SO$_2$ capacity (lb)", 2_500.0),
    ):
        ax.set(xlabel="Fuselage radius (in)", ylabel=ylabel,
               xlim=(np.min(radius), np.max(radius)))
        ax.set_xticks(endpoint_ticks(radius, 5.0))
        ax.set_yticks(rounded_ticks(mass_data[column], step))
        lower_limit = rounded_limits(mass_data[column], step)[0]
        upper_limit = mass_data[column].max()
        headroom = 0.05 * (upper_limit - lower_limit)
        ax.set_ylim(lower_limit, upper_limit + headroom)

    save_figure(fig, "mass_trade__e175_cfm56_radius.png")


def main():
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    for csv_file in sorted(OUTPUT_DIR.glob(f"*{PAYLOAD_MACH_SUFFIX}")):
        name = csv_file.name.removesuffix(PAYLOAD_MACH_SUFFIX)
        plot_payload_mach_sweep(csv_file)
        print(f"Plotted {name}: {csv_file.name}")

    if RADIUS_CSV.is_file():
        radius_data = pd.read_csv(RADIUS_CSV)
        plot_radius_ceiling(radius_data)
        plot_mass_trade(radius_data)
        print(f"Plotted radius sweep: {RADIUS_CSV.name}")


if __name__ == "__main__":
    main()
