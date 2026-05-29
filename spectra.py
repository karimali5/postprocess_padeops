from utilities import plot_spectra, plot_spectrum
import numpy as np
import os

S800 = {"name":"S800", "path": "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH800/slices_t092440_n014462", "timestamp": "092440", "nstamp": "014462"}
S250 = {"name":"S250", "path": "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH250/slices_t092897_n014397", "timestamp": "092897", "nstamp": "014397"}
stations = [10,20,30,40,50,60,70,80,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540]
xfarm = 120
yfarm = (95.61+63.11)/2
zfarm = (95-126/2)/126
Lfarm = 40
Wfarm = 95.61-63.11
Zfarm = 1
hubheight = 95/126
nx = 2200
ny = 1000
nz = 800
Lx = 698.4126984
Ly = 158.7301587
Lz = 63.49206349
dx = Lx/nx
dy = Ly/ny
dz = Lz/nz
sim = S800
field = "ddy_tauxy"

label_latex = {
    "ddz_tauxz":r"$\partial_z \Delta \tau_{xz}$", 
    "ddy_tauxy":r"$\partial_y \Delta \tau_{xy}$", 
    "ddx_dp":r"$\partial_x \Delta p$", 
    "ddz_dp":r"$\partial_z \Delta p$", 
    }

saved_settings={
    "S800":{
        "ddz_tauxz":{
            "smooth_bins":200,
            "xright": 2*np.pi/(3.5*dy),
            "ybot": 5*10**(-4),
            "wavenumber_marker_ymax": 2*10**(1),
            "wavenumber_marker_label_y": 2*10**(0),
        },
        "ddy_tauxy":{
            "smooth_bins":200,
            "xright": 2*np.pi/(3.5*dy),
            "ybot": 1*10**(-5),
            "wavenumber_marker_ymax": 1*10**(2),
            "wavenumber_marker_label_y": 9*10**(0),
        },
        "ddx_dp":{
            "smooth_bins":200,
            "xright": 2*np.pi/(3.5*dy),
            "ybot": 1*10**(-4),
            "wavenumber_marker_ymax": 8*10**(0),
            "wavenumber_marker_label_y": 1*10**(0),
        },
        "ddz_dp":{
            "smooth_bins":200,
            "xright": 2*np.pi/(3.5*dy),
            "ybot": 1*10**(-4),
            "ytop": 1*10**(1),
            "wavenumber_marker_ymax": 2*10**(0),
            "wavenumber_marker_label_y": 6*10**(-1),
        },
    },

    "S250":{
        "ddz_tauxz":{
            "smooth_bins":200,
            "xright": 2*np.pi/(3.5*dy),
            "ybot": 5*10**(-4),
            "wavenumber_marker_ymax": 2*10**(1),
            "wavenumber_marker_label_y": 2*10**(0),
        },
        "ddy_tauxy":{
            "smooth_bins":200,
            "xright": 2*np.pi/(3.5*dy),
            "ybot": 1*10**(-5),
            "wavenumber_marker_ymax": 1*10**(2),
            "wavenumber_marker_label_y": 9*10**(0),
        },
        "ddx_dp":{
            "smooth_bins":200,
            "xright": 2*np.pi/(3.5*dy),
            "ybot": 1*10**(-4),
            "wavenumber_marker_ymax": 8*10**(0),
            "wavenumber_marker_label_y": 1*10**(0),
        },
        "ddz_dp":{
            "smooth_bins":200,
            "xright": 2*np.pi/(3.5*dy),
            "ybot": 1*10**(-4),
            "ytop": 1*10**(1),
            "wavenumber_marker_ymax": 2*10**(0),
            "wavenumber_marker_label_y": 6*10**(-1),
        }
    },
}

def saved_setting(key):
    return saved_settings.get(sim["name"], {}).get(field, {}).get(key, None)

# Spectrum file selection
# kind: Spectrum.F90 output family. Options include horizontal, x/streamwise,
#       y/spanwise, zsummary, and height.
# field: Field name used by Spectrum.F90 after sanitizing spaces/symbols.
# directory: Folder containing Spectrum.F90 CSV output.
spectrum_selection = {
    "kind": "horizontal",
    "field": None,
    "directory": None,
}

# Data transforms
# fieldgain: Multiplier applied to spectral energy E.
# lengthgain: Multiplier applied to z, z_centroid, and z_spread.
# kgain: Multiplier applied to k, kx, and ky.
# background: Optional spectrum CSV/dict to subtract.
# percentage: If background is set, plot percent change relative to background.
# normalize: None, "parseval"/"variance", "max", "integral"/"sum", True, or a numeric divisor.
# normalize_reference: CSV used when normalize="parseval"; use one file to preserve relative curve magnitudes.
# premultiply: Plot kE instead of E for line spectra.
# density: Label helper for spectra written with write_density=.true.
# skip_zero: Remove k=0 before log-axis plotting.
# k_position: Plot CSV wavenumbers as bin "center", "lower", or "upper" edges.
data_style = {
    "fieldgain": 1.0,
    "lengthgain": 1.0,
    "kgain": 1.0,
    "background": None,
    "percentage": False,
    "normalize": "parseval",
    "normalize_reference": os.path.join(sim["path"], f"spectrum_{field}.csv"),
    "premultiply": True,
    "density": True,
    "skip_zero": True,
    "k_position": None,
}

# Line plot style
# color: Matplotlib line color. None uses the default color cycle.
# linewidth: Line width.
# linestyle: Line style.
# marker: Marker style. None disables markers.
# markersize: Marker size. None uses Matplotlib default.
# alpha: Line transparency.
line_style = {
    "color": "k",
    "linewidth": 1.2,
    "linestyle": "-",
    "marker": None,
    "markersize": None,
    "alpha": None,
}

# Height-resolved spectrum style
# cmap: Colormap for spectrum_height_<field>.csv.
# norm: Color normalization. "log" is useful for spectra; use "linear" for raw colors.
# vmin/vmax: Color limits.
# colorbar: Draw local colorbars for height-resolved panels.
# colorbar_label: Colorbar label. None uses default spectrum label.
# colorbar_pad: Spacing between axes and colorbar.
# colorbar_size: Width of matched vertical colorbar.
height_style = {
    "cmap": "viridis",
    "norm": "log",
    "vmin": None,
    "vmax": None,
    "colorbar": True,
    "colorbar_label": None,
    "colorbar_pad": 0.08,
    "colorbar_size": "2%",
}

# Y-z plane spectra reduction
# yzplane_reduce: How to collapse spectrum_yzplane_<field>.csv over x.
#                 Use "sum" for total active wake/farm energy, "mean", "median", or "none".
# x_bounds: Optional (xmin, xmax) before reduction, useful for farm/wake subregions.
# x_normalize: If True with yzplane_reduce="sum", divide by the number of x planes.
yzplane_style = {
    "yzplane_reduce": "sum",
    "x_bounds": None,
    "x_normalize": False,
}

# Spectrum smoothing
# smooth: Enable post-processing smoothing for display/comparison.
# smooth_method: "logbin" groups high-resolution k data into logarithmic bins;
#                "moving" applies a centered moving window on the existing k grid.
# smooth_bins: Number of log bins for smooth_method="logbin".
# smooth_window: Moving-window width for smooth_method="moving"; even values are rounded up.
# smooth_stat: Bin/window statistic. Use "mean", "median", "sum", or "weighted_mean".
# smooth_height: Apply the same smoothing independently at each z for height spectra.
smoothing_style = {
    "smooth": True,
    "smooth_method": "logbin",
    "smooth_bins": saved_setting("smooth_bins"),
    "smooth_window": 5,
    "smooth_stat": "mean",
    "smooth_height": True,
}

# Axes and quantity selection
# xscale/yscale: Axis scales, usually "log" for line spectra.
# xlim/ylim: Direct axis limits.
# klim/zlim: Aliases for spectrum k-limits and height z-limits.
# xleft/xright/ybot/ytop: One-sided axis limits.
# xlabel/ylabel/title: Manual labels and title. None uses defaults where possible.
# grid/grid_which: Toggle grid and choose "major", "minor", or "both".
# quantity: Quantity to plot from CSV. Use E, z_centroid, or z_spread when available.
axis_style = {
    "xscale": "log",
    "yscale": "log",
    "xlim": None,
    "ylim": None,
    "klim": None,
    "zlim": None,
    "xleft": None,
    "xright": saved_setting("xright"),
    "ybot": saved_setting("ybot"),
    "ytop": saved_setting("ytop"),
    "xlabel": r"$k_\alpha$",
    "ylabel": r"$k_\alpha E_\alpha/\langle q^2 \rangle$",
    "title": label_latex[field],
    "grid": False,
    "grid_which": "both",
    "quantity": "E",
}

# Reference slopes
# references: List of power-law guide lines for line spectra. Each entry needs
#             {"slope": value, "anchor": (k0, E0)} and can include label/color/style.
reference_style = {
    "references": [
        # Example:
        # {"slope": -5/3, "anchor": (0.1, 1.0), "label": r"$k^{-5/3}$"},
    ],
}

# Wavenumber scale markers
# turbine_diameter: Physical/nondimensional turbine diameter D. Plots k=2*pi/D.
# farm_length: Physical/nondimensional farm length L_f. Plots k=2*pi/L_f.
# turbine_diameter_label/farm_length_label: Text labels beside the marker lines.
# wavenumber_markers: Optional custom marker list. Each entry can use {"length": value}
#                     for k=2*pi/length or {"k": value} for direct placement.
# wavenumber_marker_label_y: Label height near the x-axis in axes coordinates by default.
# wavenumber_marker_label_x_nudge: Rightward text offset in points from the marker line.
# wavenumber_marker_ymin/ymax: Optional data-coordinate vertical extent for marker lines.
# wavenumber_marker_*: Default style for all marker lines/text.
wavenumber_marker_style = {
    "turbine_diameter": 1,
    "farm_length": Lfarm,
    "turbine_diameter_label": r"$2\pi$",
    "farm_length_label": r"$2\pi D/L_f$",
    "wavenumber_markers": 
    [
        {"length": 2.5, "label":r"$4\pi D/s_y$"},
        {"length": 5, "label":r"$2\pi D/s_x$"}, 
        {"length": Lx, "label":r"$2\pi D/L_x$"}, 
        {"length": Ly, "label":r"$2\pi D/L_y$"},
    ],
    "wavenumber_marker_color": "k",
    "wavenumber_marker_linestyle": "--",
    "wavenumber_marker_alpha": 0.7,
    "wavenumber_marker_linewidth": 0.8,
    "wavenumber_marker_ymin": None,
    "wavenumber_marker_ymax": saved_setting("wavenumber_marker_ymax"),
    "wavenumber_marker_label": True,
    "wavenumber_marker_label_y": saved_setting("wavenumber_marker_label_y"),
    "wavenumber_marker_label_x_nudge": -8.0,
    "wavenumber_marker_fontsize": 6,
}

# Figure layout and export
# figsize: Figure size in inches.
# dpi: Saved image resolution.
# panel: If True, make one subplot per simulation. Height spectra force panel mode.
# ncols: Number of subplot columns in panel mode.
# annotation: Panel labels, using the same options as plot_slices.
# stamp and stamp_*: Optional figure-level annotation.
# tight_layout: Apply Matplotlib tight layout before saving.
# bbox_inches/pad_inches: Savefig bounding-box controls.
# export: Output image filename.
layout_style = {
    "figsize": (5.2, 3.4),
    "dpi": 300,
    "panel": False,
    "ncols": 1,
    "annotation": "letter_name",
    "stamp": None,
    "stamp_location": (0.98, 0.98),
    "stamp_ha": "right",
    "stamp_va": "top",
    "tight_layout": True,
    "bbox_inches": "tight",
    "pad_inches": 0.02,
    "export": f"spectra_{sim["name"]}_{field}.png",
}

# Legend options
# legend: Draw legend for line spectra.
# legend_loc: Matplotlib legend location. None lets Matplotlib choose.
# legend_ncols: Number of legend columns.
# legend_frameon: Toggle legend frame.
# legend_kwargs: Extra keyword arguments passed to ax.legend.
legend_style = {
    "legend": True,
    "legend_loc": "upper left",
    "legend_ncols": 3,
    "legend_frameon": False,
    "legend_kwargs": {},
}

# Font sizes
# tick_fontsize: Tick labels.
# label_fontsize: Axis labels.
# title_fontsize: Axis title.
# legend_fontsize: Legend text.
# annotation_fontsize: Panel annotation text.
# stamp_fontsize: Figure stamp text.
# colorbar_tick_fontsize: Colorbar tick labels.
# colorbar_label_fontsize: Colorbar label.
font_style = {
    "tick_fontsize": 7,
    "label_fontsize": 8,
    "title_fontsize": 9,
    "legend_fontsize": 7,
    "annotation_fontsize": 7,
    "stamp_fontsize": 9,
    "colorbar_tick_fontsize": 7,
    "colorbar_label_fontsize": 7,
}

all_kwargs = {
    **spectrum_selection,
    **data_style,
    **line_style,
    **height_style,
    **yzplane_style,
    **smoothing_style,
    **axis_style,
    **reference_style,
    **wavenumber_marker_style,
    **layout_style,
    **legend_style,
    **font_style,
}


def plot_all_spectra(simulations, **overrides):
    kwargs = all_kwargs.copy()
    kwargs.update(overrides)
    return plot_spectra(simulations, **kwargs)


def plot_one_spectrum(filename, **overrides):
    kwargs = all_kwargs.copy()
    kwargs.update(overrides)
    return plot_spectrum(filename, **kwargs)


example_simulations = [
    # Use explicit filenames:
    # {"filename": "/path/to/spectrum_u.csv", "name": "case A"},
    # {"filename": "/path/to/spectrum_u.csv", "name": "case B"},
    #
    # Or let utilities build names from Spectrum.F90 conventions:
    # {"directory": "/path/to/output", "field": "u", "kind": "horizontal", "name": "case A"},
    # {"directory": "/path/to/output", "field": "u", "kind": "x", "name": "case A"},
    {
        "filename": os.path.join(sim["path"],f"spectrum_{field}.csv"),
        "name": r"Horizontal $(\alpha=h)$",
        "plot_kwargs": {
              "color": "k",
              "linestyle": "-",
              "linewidth": 1.4,
              "k_position": "lower",
          },
    },

    {
        "filename": os.path.join(sim["path"], f"spectrum_x_{field}.csv"),
        "name": r"Streamwise $(\alpha=x)$",
        "plot_kwargs": {
              "color": "tab:blue",
              "linestyle": "-",
              "linewidth": 1.4,
              "k_position": "center",
          },
    },

    {
        "filename": os.path.join(sim["path"],f"spectrum_y_{field}.csv"),
        "name": r"Lateral $(\alpha=y)$",
        "plot_kwargs": {
              "color": "tab:red",
              "linestyle": "-",
              "linewidth": 1.4,
              "k_position": "center",
          },
    },

]

if __name__ == "__main__":
    if example_simulations:
        plot_all_spectra(example_simulations)
