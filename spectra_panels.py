import os

import matplotlib.pyplot as plt
import numpy as np

from utilities import (
    _configure_latex_fonts,
    _plot_fontsize,
    _plot_wavenumber_markers,
    _prepare_plot_spectrum_data,
)


S800 = {
    "name": "S800",
    "path": "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH800/slices_t092440_n014462",
    "timestamp": "092440",
    "nstamp": "014462",
}
S250 = {
    "name": "S250",
    "path": "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH250/slices_t092897_n014397",
    "timestamp": "092897",
    "nstamp": "014397",
}

hubheight = 95 / 126
nx = 2200
ny = 1000
nz = 800
Lx = 698.4126984
Ly = 158.7301587
Lz = 63.49206349
dx = Lx / nx
dy = Ly / ny
dz = Lz / nz
Lfarm = 40

field_labels = {
    "ddx_dp": r"$\partial_x \Delta p$",
    "ddy_dp": r"$\partial_y \Delta p$",
    "ddz_dp": r"$\partial_z \Delta p$",
    "ddx_tauxx": r"$\partial_x \Delta \tau_{xx}$",
    "ddy_tauxy": r"$\partial_y \Delta \tau_{xy}$",
    "ddz_tauxz": r"$\partial_z \Delta \tau_{xz}$",
}

panel_defs = [
    {
        "kind": "horizontal",
        "title": r"Horizontal $(\alpha=h)$",
        "coord": "k",
        "prefix": "spectrum_",
        "k_position": "lower",
    },
    {"kind": "x", "title": r"Streamwise $(\alpha=x)$", "coord": "kx", "prefix": "spectrum_x_", "k_position": "center"},
    {"kind": "y", "title": r"Lateral $(\alpha=y)$", "coord": "ky", "prefix": "spectrum_y_", "k_position": "center"},
]

fields = [
    {"field": "ddx_dp", "label": field_labels["ddx_dp"], "color": "tab:blue"},
    {"field": "ddy_tauxy", "label": field_labels["ddy_tauxy"], "color": "tab:red"},
    {"field": "ddz_tauxz", "label": field_labels["ddz_tauxz"], "color": "k"},
]

plot_style = {
    "sim": S800,
    "figsize": (8, 4),
    "dpi": 300,
    "export": "spectra_panels_S800.png",
    "xscale": "log",
    "yscale": "log",
    "xleft": 1 * np.pi / Lx,
    "xright": 2 * np.pi / (3.5 * dy),
    "ylim": None,
    "ybot": 1e-5,
    "ytop": 1e2,
    "y_margin": 0.06,
    "xlabel": r"$k_\alpha$",
    "ylabel": r"$k_\alpha E_\alpha/\langle q^2\rangle$",
    "normalize": "parseval",
    "premultiply": True,
    "density": True,
    "skip_zero": True,
    "smooth": True,
    "smooth_method": "logbin",
    "smooth_bins": 300,
    "smooth_stat": "mean",
    "linewidth": 1.5,
    "alpha": 1.0,
    "legend_ncols": 3,
    "legend_loc": "lower left",
    "legend_bbox_to_anchor": (0.27, 1.02),
    "legend_frameon": False,
    "panel_label_x": 0.015,
    "panel_label_y": 1.02,
    "tight_layout": False,
    "bbox_inches": "tight",
    "pad_inches": 0.03,
    "tick_fontsize": 7,
    "label_fontsize": 9,
    "legend_fontsize": 9,
    "annotation_fontsize": 8,
    "wavenumber_markers": [
        {"length": 1, "label": r"$2\pi$"},
        {"length": 2.5, "label": r"$4\pi D/s_y$"},
        {"length": 5, "label": r"$2\pi D/s_y$"},
        {"length": Lfarm, "label": r"$2\pi D/L_f$"},
        {"length": Lx, "label": r"$2\pi D/L_x$"},
        {"length": Ly, "label": r"$2\pi D/L_y$"},
    ],
    "wavenumber_marker_color": "k",
    "wavenumber_marker_linestyle": "--",
    "wavenumber_marker_alpha": 0.55,
    "wavenumber_marker_linewidth": 0.8,
    "wavenumber_marker_ymin": None,
    "wavenumber_marker_ymax": 3e1,
    "wavenumber_marker_label": False,
    "wavenumber_marker_fontsize": 6,
}


def spectrum_filename(directory, field, panel):
    return os.path.join(directory, f"{panel['prefix']}{field}.csv")


def prepare_panel_data(fields, panels, **kwargs):
    sim = kwargs["sim"]
    prepared = []
    printed_refs = set()

    for field_spec in fields:
        field_name = field_spec["field"]
        reference = field_spec.get(
            "normalize_reference",
            kwargs.get("normalize_reference", os.path.join(sim["path"], f"spectrum_{field_name}.csv")),
        )
        for panel in panels:
            panel_kwargs = kwargs.copy()
            panel_kwargs.update(field_spec.get("plot_kwargs", {}))
            panel_kwargs["_printed_normalization_references"] = printed_refs
            panel_kwargs["normalize_reference"] = reference
            panel_kwargs["k_position"] = panel.get("k_position", panel_kwargs.get("k_position", None))
            filename = spectrum_filename(sim["path"], field_name, panel)
            data = _prepare_plot_spectrum_data(filename, **panel_kwargs)
            prepared.append(
                {
                    "field": field_spec,
                    "panel": panel,
                    "data": data,
                }
            )
    return prepared


def spectrum_xy(data, coord_key):
    return data[coord_key], data["E"]


def apply_shared_limits(axes, prepared, panels, **kwargs):
    x_values = []
    for item in prepared:
        coord_key = item["panel"]["coord"]
        x, y = spectrum_xy(item["data"], coord_key)
        finite = np.isfinite(x) & np.isfinite(y)
        if np.any(finite):
            x_values.append(x[finite])

    if x_values:
        all_x = np.concatenate(x_values)
        x_left = kwargs.get("xleft", None)
        x_right = kwargs.get("xright", None)
        if x_left is None:
            positive = all_x[all_x > 0]
            if positive.size:
                x_left = np.nanmin(positive)
        if x_right is None:
            x_right = np.nanmax(all_x)
        for ax in axes:
            ax.set_xlim(left=x_left, right=x_right)

    if kwargs.get("ylim", None) is not None:
        for ax in axes:
            ax.set_ylim(kwargs["ylim"])
    else:
        y_bottom = kwargs.get("ybot", kwargs.get("ymin", 0.0))
        y_top = kwargs.get("ytop", kwargs.get("ymax", None))
        for ax in axes:
            ax.set_ylim(bottom=y_bottom, top=y_top)


def plot_spectra_panels(fields, panels=panel_defs, **kwargs):
    _configure_latex_fonts()
    prepared = prepare_panel_data(fields, panels, **kwargs)

    fig, axes = plt.subplots(
        len(panels),
        1,
        figsize=kwargs.get("figsize", (5.4, 6.4)),
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]

    panel_lookup = {panel["kind"]: ax for panel, ax in zip(panels, axes)}
    legend_handles = {}

    for item in prepared:
        field_spec = item["field"]
        panel = item["panel"]
        data = item["data"]
        ax = panel_lookup[panel["kind"]]
        x, y = spectrum_xy(data, panel["coord"])
        line, = ax.plot(
            x,
            y,
            color=field_spec.get("color", None),
            linestyle=field_spec.get("linestyle", "-"),
            linewidth=field_spec.get("linewidth", kwargs.get("linewidth", 1.5)),
            alpha=field_spec.get("alpha", kwargs.get("alpha", 1.0)),
            label=field_spec.get("label", field_spec["field"]),
        )
        legend_handles.setdefault(field_spec["field"], line)

    apply_shared_limits(axes, prepared, panels, **kwargs)

    letters = "abcdefghijklmnopqrstuvwxyz"
    for i, (panel, ax) in enumerate(zip(panels, axes)):
        ax.set_xscale(kwargs.get("xscale", "log"))
        ax.set_yscale(kwargs.get("yscale", "linear"))
        ax.tick_params(axis="both", which="major", labelsize=_plot_fontsize(kwargs, "tick_fontsize", 7))
        if kwargs.get("grid", False):
            ax.grid(True, which=kwargs.get("grid_which", "both"), alpha=0.25, linewidth=0.5)
        _plot_wavenumber_markers(ax, **kwargs)
        ax.set_ylabel(kwargs.get("ylabel", r"$k_\alpha E_\alpha/\langle q^2\rangle$"), fontsize=_plot_fontsize(kwargs, "label_fontsize", 8))
        ax.text(
            kwargs.get("panel_label_x", 0.015),
            kwargs.get("panel_label_y", 0.92),
            f"({letters[i]}) {panel['title']}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=_plot_fontsize(kwargs, "annotation_fontsize", 8),
            clip_on=False,
        )

    axes[-1].set_xlabel(kwargs.get("xlabel", r"$k_\alpha$"), fontsize=_plot_fontsize(kwargs, "label_fontsize", 8))

    handles = list(legend_handles.values())
    labels = [handle.get_label() for handle in handles]
    axes[0].legend(
        handles,
        labels,
        loc=kwargs.get("legend_loc", "lower left"),
        bbox_to_anchor=kwargs.get("legend_bbox_to_anchor", (0.0, 1.02)),
        ncols=kwargs.get("legend_ncols", len(handles)),
        frameon=kwargs.get("legend_frameon", False),
        fontsize=_plot_fontsize(kwargs, "legend_fontsize", 7),
    )

    fig.align_ylabels(axes)
    if kwargs.get("tight_layout", False):
        fig.tight_layout()
    else:
        fig.subplots_adjust(left=0.15, right=0.98, bottom=0.08, top=0.9, hspace=0.20)

    export = kwargs.get("export", "spectra_panels.png")
    print(export)
    fig.savefig(
        export,
        dpi=kwargs.get("dpi", 300),
        bbox_inches=kwargs.get("bbox_inches", "tight"),
        pad_inches=kwargs.get("pad_inches", 0.03),
    )
    if kwargs.get("close", True):
        plt.close(fig)
    return prepared


if __name__ == "__main__":
    plot_spectra_panels(fields, **plot_style)
