import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import csv
import netCDF4 as nc
import numpy as np
import warnings

lw = 1.
hubheight = 95 / 126


def load_profile(file_path):
    """
    Load a vertical profile from a CSV file.

    Returns
    -------
    z : list
        List of vertical coordinates.
    profile : list
        List of horizontally averaged values.
    """
    z = []
    profile = []
    try:
        csvfile = open(file_path, 'r')
    except FileNotFoundError:
        warnings.warn(f"Profile CSV file not found: {file_path}", stacklevel=2)
        raise

    with csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip the header row
        for row in reader:
            z.append(float(row[0]))
            profile.append(float(row[1]))
    return z, profile


def load_horizontal_average_slice(file_path, scale=1.0):
    with nc.Dataset(file_path, "r") as ds:
        dimensions = set(ds.dimensions.keys())
        for name, variable in ds.variables.items():
            if name in dimensions or getattr(variable, "ndim", 0) != 2:
                continue
            return float(np.nanmean(np.asarray(variable[:]) * scale))
    raise ValueError(f"No 2D slice field found in {file_path}.")


def plot_profile_panels(profile_panels, **kwargs):
    """
    Plot side-by-side vertical profile panels.
    """
    export = kwargs.get('export', None)
    ylabel = kwargs.get('ylabel', None)
    figsize = kwargs.get('figsize', (6, 2.4))
    dpi = kwargs.get('dpi', 300)
    ylim = kwargs.get('ylim', None)
    ybot = kwargs.get('ybot', kwargs.get('ymin', None))
    ytop = kwargs.get('ytop', kwargs.get('ymax', None))
    grid = kwargs.get('grid', False)
    plotfarm = kwargs.get('plotfarm', False)
    farm_bounds = kwargs.get('farm_bounds', None)
    legend = kwargs.get('legend', True)
    legend_loc = kwargs.get('legend_loc', None)
    legend_ncols = kwargs.get('legend_ncols', 1)
    legend_frameon = kwargs.get('legend_frameon', False)
    legend_kwargs = kwargs.get('legend_kwargs', {})
    tight_layout = kwargs.get('tight_layout', True)
    close = kwargs.get('close', False)
    hide_top_right_spines = kwargs.get('hide_top_right_spines', False)
    tick_fontsize = kwargs.get('tick_fontsize', 7)
    label_fontsize = kwargs.get('label_fontsize', 9)
    title_fontsize = kwargs.get('title_fontsize', 9)
    legend_fontsize = kwargs.get('legend_fontsize', 9)
    annotation = kwargs.get('annotation', "letter_name")
    linewidth = kwargs.get('linewidth', 1.2)
    linestyle = kwargs.get('linestyle', '-')
    marker = kwargs.get('marker', None)
    alpha = kwargs.get('alpha', None)
    color = kwargs.get('color', None)
    zorder = kwargs.get('zorder', None)
    wspace = kwargs.get('wspace', 0.28)
    tighten_whitespace = kwargs.get('tighten_whitespace', False)
    subplot_adjust = kwargs.get('subplot_adjust', {})

    with plt.rc_context({"text.usetex": True, "font.family": "serif"}):
        fig, axes = plt.subplots(
            1,
            len(profile_panels),
            figsize=figsize,
            sharey=True,
            squeeze=False,
        )
        axes = axes[0, :]

        for panel_index, (ax, panel) in enumerate(zip(axes, profile_panels)):
            panel_profiles = panel.get("profiles", panel) if isinstance(panel, dict) else panel
            panel_xlabel = panel.get("xlabel", None) if isinstance(panel, dict) else None
            panel_name = panel.get("name", None) if isinstance(panel, dict) else None
            panel_title = panel.get("title", None) if isinstance(panel, dict) else None
            panel_legend = panel.get("legend", legend) if isinstance(panel, dict) else legend
            panel_legend_loc = panel.get("legend_loc", legend_loc) if isinstance(panel, dict) else legend_loc
            panel_legend_kwargs = legend_kwargs.copy()
            if isinstance(panel, dict):
                panel_legend_kwargs.update(panel.get("legend_kwargs", {}))
            panel_zeroline = panel.get("zeroline", False) if isinstance(panel, dict) else False
            panel_shaded_regions = panel.get("shaded_regions", []) if isinstance(panel, dict) else []
            panel_horizontal_lines = panel.get("horizontal_lines", []) if isinstance(panel, dict) else []
            panel_auto_xlim_visible = panel.get("auto_xlim_visible", False) if isinstance(panel, dict) else False
            panel_xtick_spacing = panel.get("xtick_spacing", None) if isinstance(panel, dict) else None
            panel_xlim = panel.get("xlim", None) if isinstance(panel, dict) else None
            panel_xlft = panel.get("xlft", panel.get("xleft", None)) if isinstance(panel, dict) else None
            panel_xrght = panel.get("xrght", panel.get("xright", None)) if isinstance(panel, dict) else None
            panel_ylim = panel.get("ylim", ylim) if isinstance(panel, dict) else ylim
            panel_ybot = panel.get("ybot", panel.get("ymin", ybot)) if isinstance(panel, dict) else ybot
            panel_ytop = panel.get("ytop", panel.get("ymax", ytop)) if isinstance(panel, dict) else ytop
            visible_values = []

            for profile in panel_profiles:
                z, values = load_profile(profile['file'])
                z = np.asarray(z) / hubheight
                profile_gain = profile.get('gain', profile.get('xgain', 1.0))
                values = np.asarray(values) * profile_gain
                transform = profile.get("transform", None)
                if transform is not None:
                    values = transform(values)
                if panel_auto_xlim_visible:
                    visible = np.ones_like(z, dtype=bool)
                    if panel_ylim is not None:
                        visible &= z >= min(panel_ylim)
                        visible &= z <= max(panel_ylim)
                    if panel_ybot is not None:
                        visible &= z >= panel_ybot
                    if panel_ytop is not None:
                        visible &= z <= panel_ytop
                    visible_values.extend(values[visible])
                ax.plot(
                    values,
                    z,
                    label=profile.get('label', profile.get('name', 'p')),
                    color=profile.get('color', color),
                    linestyle=profile.get('linestyle', profile.get('style', linestyle)),
                    linewidth=profile.get('linewidth', linewidth),
                    marker=profile.get('marker', marker),
                    alpha=profile.get('alpha', alpha),
                    zorder=profile.get('zorder', zorder),
                )

            if panel_xlim is not None:
                ax.set_xlim(panel_xlim)
            if panel_xlft is not None:
                ax.set_xlim(left=panel_xlft)
            if panel_xrght is not None:
                ax.set_xlim(right=panel_xrght)
            if panel_ylim is not None:
                ax.set_ylim(panel_ylim)
            if panel_ybot is not None:
                ax.set_ylim(bottom=panel_ybot)
            if panel_ytop is not None:
                ax.set_ylim(top=panel_ytop)

            if plotfarm:
                if farm_bounds is None:
                    rotor_bottom = (hubheight - 0.5) / hubheight
                    rotor_top = (hubheight + 0.5) / hubheight
                else:
                    rotor_bottom, rotor_top = farm_bounds
                ax.axhspan(
                    rotor_bottom,
                    rotor_top,
                    facecolor=kwargs.get('farm_facecolor', 'tab:red'),
                    edgecolor=kwargs.get('farm_edgecolor', 'none'),
                    alpha=kwargs.get('farm_alpha', 0.12),
                    linewidth=kwargs.get('farm_linewidth', 0),
                    zorder=kwargs.get('farm_zorder', 0),
                )
            for region in panel_shaded_regions:
                ax.axhspan(
                    region["bottom"],
                    region["top"],
                    facecolor=region.get("facecolor", "grey"),
                    edgecolor=region.get("edgecolor", "none"),
                    alpha=region.get("alpha", 0.12),
                    linewidth=region.get("linewidth", 0),
                    zorder=region.get("zorder", 0),
                )
            for line in panel_horizontal_lines:
                ax.axhline(
                    line["y"],
                    color=line.get("color", "0.45"),
                    linewidth=line.get("linewidth", 0.7),
                    linestyle=line.get("linestyle", "--"),
                    alpha=line.get("alpha", 1.0),
                    zorder=line.get("zorder", 1),
                )
            if panel_zeroline:
                ax.axvline(
                    0,
                    color=kwargs.get('zeroline_color', '0.4'),
                    linewidth=kwargs.get('zeroline_width', 0.7),
                    linestyle=kwargs.get('zeroline_style', ':'),
                    zorder=kwargs.get('zeroline_zorder', 0),
                )

            if panel_xlabel is not None:
                ax.set_xlabel(panel_xlabel, fontsize=label_fontsize)
            if panel_index == 0 and ylabel is not None:
                ax.set_ylabel(ylabel, fontsize=label_fontsize)
            if panel_title is None and annotation == "letter_name":
                letter = chr(ord("a") + panel_index)
                panel_title = f"({letter}) {panel_name}" if panel_name else f"({letter})"
            if panel_title is not None:
                ax.set_title(panel_title, fontsize=title_fontsize, loc='left')

            ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
            if grid:
                ax.grid(True, alpha=kwargs.get('grid_alpha', 0.25), linewidth=kwargs.get('grid_linewidth', 0.5))
            if hide_top_right_spines:
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.tick_params(top=False, right=False)
            if panel_legend:
                ax.legend(
                    loc=panel_legend_loc,
                    ncol=legend_ncols,
                    frameon=legend_frameon,
                    fontsize=legend_fontsize,
                    **panel_legend_kwargs,
                )
            if panel_auto_xlim_visible and visible_values:
                xmin = min(visible_values)
                xmax = max(visible_values)
                dx = xmax - xmin
                pad = 0.05 * dx if dx > 0 else 0.5
                ax.set_xlim(xmin - pad, xmax + pad)
            if panel_xtick_spacing is not None:
                ax.xaxis.set_major_locator(MultipleLocator(panel_xtick_spacing))

        if tight_layout:
            fig.tight_layout()
            fig.subplots_adjust(wspace=wspace)
        if tighten_whitespace:
            default_adjust = {
                "left": 0.09,
                "right": 0.995,
                "bottom": 0.19,
                "top": 0.88,
                "wspace": wspace,
            }
            default_adjust.update(subplot_adjust)
            fig.subplots_adjust(**default_adjust)
        if export is None:
            export = "baseflow.png"
        print(export)
        fig.savefig(export, dpi=dpi, bbox_inches=kwargs.get('bbox_inches', None), pad_inches=kwargs.get('pad_inches', 0.1))
        if close:
            plt.close(fig)
        return fig, axes


profiles_dir = "/anvil/scratch/x-kali/PadeOpsSims/INV800-5K/profiles"
slices_dir = "/anvil/scratch/x-kali/PadeOpsSims/INV800-5K/slices"

inversion_bottom = load_horizontal_average_slice(f"{slices_dir}/Run08_t101691_INVH0.nc", scale=1/126) / hubheight
inversion_top = load_horizontal_average_slice(f"{slices_dir}/Run08_t101691_INVH2.nc", scale=1/126) / hubheight
inversion_lines = [
    {"y": inversion_bottom, "color": "0.45", "linestyle": "--", "linewidth": 0.7, "zorder": 1},
    {"y": inversion_top, "color": "0.45", "linestyle": "--", "linewidth": 0.7, "zorder": 1},
]

panel_a_legend_loc = "upper left"
panel_a_legend_bbox_to_anchor = (0.3, 0.38)
panel_a_legend_borderaxespad = 0.0

profile_panels = [
    {
        "profiles": [
            {
                "file": f"{profiles_dir}/u_HA_z.csv",
                "label": r"$\langle {\overline{u}}^b \rangle_{xy}$",
                "color": "k",
                "style": "-",
                "linewidth": lw,
            },
            {
                "file": f"{profiles_dir}/v_HA_z.csv",
                "label": r"$\langle {\overline{v}}^b \rangle_{xy}$",
                "color": "k",
                "style": "--",
                "linewidth": lw,
            },
        ],
        "xlabel": r"$u,\,v$",
        "name": r"Wind speed components",
        "legend_loc": panel_a_legend_loc,
        "legend_kwargs": {
            "bbox_to_anchor": panel_a_legend_bbox_to_anchor,
            "borderaxespad": panel_a_legend_borderaxespad,
        },
        "zeroline": False,
        "horizontal_lines": inversion_lines,
        "xtick_spacing": 0.25,
    },
    {
        "profiles": [
            {
                "file": f"{profiles_dir}/theta_HA_z.csv",
                "label": r"$\langle {\overline{\theta}}^b \rangle_{xy}$",
                "color": "k",
                "style": "-",
                "linewidth": lw,
                "transform": lambda theta: theta - 273.15,
            },
        ],
        "xlabel": r"$\theta~[^\circ\mathrm{C}]$",
        "name": r"Potential Temperature $\langle {\overline{\theta}}^b \rangle_{xy}$",
        "legend": False,
        "auto_xlim_visible": True,
        "horizontal_lines": inversion_lines,
        "xtick_spacing": 1,
        # "shaded_regions": [
        #     {
        #         "bottom": inversion_bottom,
        #         "top": inversion_top,
        #         "facecolor": "0.7",
        #         "alpha": 0.35,
        #         "zorder": 0,
        #     },
        # ],
    },
    {
        "profiles": [
            {
                "file": f"{profiles_dir}/tke_HA_z.csv",
                "label": r"$\langle {\overline{k}}^b \rangle_{xy}$",
                "color": "k",
                "style": "-",
                "linewidth": lw,
            },
        ],
        "xlabel": r"$k$",
        "name": r"TKE $\langle {\overline{k}}^b \rangle_{xy}$",
        "legend": False,
        "horizontal_lines": inversion_lines,
        "xtick_spacing": 0.001,
    },
]

plot_style = {
    "figsize": (5.2, 1.8),
    "dpi": 200,
    "hide_top_right_spines": True,
    "legend": True,
    "legend_loc": "lower left",
    "legend_frameon": False,
    "legend_fontsize": 7,
    "tick_fontsize": 6,
    "label_fontsize": 7,
    "title_fontsize": 7,
    "tight_layout": True,
    "wspace": 0.2,
    "tighten_whitespace": True,
    "plotfarm": True,
    "farm_alpha": 0.12,
    "farm_zorder": 0,
    "annotation": "letter_name",
    "export": "baseflow.png",
}

axis_style = {
    "ylabel": r"$z/z_h$",
    "ybot": 0,
    "ytop": 14,
}

if __name__ == "__main__":
    plot_profile_panels(
        profile_panels,
        **plot_style,
        **axis_style,
    )
