import matplotlib.pyplot as plt
import csv
import numpy as np
import warnings

xfarm = 160
Lfarm = 40
Ly =  158.7301587
Lz =  63.49206349
rmsgain = np.sqrt(1/(Ly*Lz))
lw = 1.

def load_profile(file_path):
    """
    Load a profile from one CSV file or sum profiles from a list of CSV files.

    Parameters
    ----------
    file_path : str or list of str
        Path to the CSV file containing the profile data. If a list is passed,
        the y-values from all files are added together.

    Returns
    -------
    x : list
        List of x-coordinates.
    y : list
        List of y-coordinates.
    """
    if isinstance(file_path, (list, tuple)):
        x_sum = None
        y_sum = None
        for path in file_path:
            x, y = load_profile(path)
            x = np.asarray(x)
            y = np.asarray(y)
            if x_sum is None:
                x_sum = x
                y_sum = y
            else:
                if x.shape != x_sum.shape or not np.allclose(x, x_sum):
                    raise ValueError(f"Profile x-coordinates do not match for {path}.")
                y_sum = y_sum + y
        if x_sum is None:
            raise ValueError("file_path list must include at least one file.")
        return x_sum, y_sum

    x = []
    y = []
    try:
        csvfile = open(file_path, 'r')
    except FileNotFoundError:
        warnings.warn(f"Profile CSV file not found: {file_path}", stacklevel=2)
        raise

    with csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip the header row
        for row in reader:
            x.append(float(row[0]))
            y.append(float(row[1]))
    return x, y

def _axis_transform_function(axis_name, transformer):
    if transformer is None:
        return None
    if not isinstance(transformer, dict):
        raise TypeError(f"{axis_name}axis_transformer must be a dictionary or None.")

    transform = (
        transformer.get("function")
        or transformer.get("transform")
        or transformer.get("callable")
        or transformer.get("func")
    )
    if not callable(transform):
        raise TypeError(
            f"{axis_name}axis_transformer must include a callable under "
            "'function', 'transform', 'callable', or 'func'."
        )
    return transform

def _axis_inverse_function(axis_name, transformer):
    if transformer is None:
        return None

    inverse = (
        transformer.get("inverse")
        or transformer.get("inverse_function")
        or transformer.get("inverse_transform")
    )
    if inverse is not None and not callable(inverse):
        raise TypeError(f"{axis_name}axis_transformer inverse must be callable.")
    return inverse

def _transform_values(values, transform):
    values = np.asarray(values)
    if transform is None:
        return values
    return np.asarray(transform(values))

def _transform_scalar(value, transform):
    if value is None or transform is None:
        return value
    return float(np.asarray(transform(np.asarray(value))))

def _apply_profile_offset(x, y, offset_x):
    if offset_x is None or offset_x == 0:
        return y

    x = np.asarray(x)
    y = np.asarray(y)
    sort_order = np.argsort(x)
    x_sorted = x[sort_order]
    y_sorted = y[sort_order]

    if offset_x < x_sorted[0] or offset_x > x_sorted[-1]:
        raise ValueError(
            f"Profile offset x-coordinate {offset_x} is outside the profile "
            f"x-range [{x_sorted[0]}, {x_sorted[-1]}]."
        )

    return y - np.interp(offset_x, x_sorted, y_sorted)

def _apply_axis_transformer(ax, axis_name, transformer, transform_axes_data, label_fontsize):
    if transformer is None:
        return

    transform = _axis_transform_function(axis_name, transformer)
    inverse = _axis_inverse_function(axis_name, transformer)
    label = transformer.get("label", transformer.get(f"{axis_name}label", None))
    ticks = transformer.get("ticks", None)

    if ticks is not None:
        ticklabels = np.asarray(ticks)
        if transform_axes_data:
            tick_positions = ticklabels
        elif inverse is None:
            raise TypeError(
                f"{axis_name}axis_transformer with 'ticks' must include a callable "
                "'inverse', 'inverse_function', or 'inverse_transform'."
            )
        else:
            tick_positions = [inverse(tick) for tick in ticklabels]
    else:
        if axis_name == "x":
            tick_positions = ax.get_xticks()
        else:
            tick_positions = ax.get_yticks()
        ticklabels = tick_positions if transform_axes_data else [transform(tick) for tick in tick_positions]

    ticklabels = [f"{tick:g}" for tick in ticklabels]
    if axis_name == "x":
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(ticklabels)
        if label is not None:
            ax.set_xlabel(label, fontsize=label_fontsize)
    else:
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(ticklabels)
        if label is not None:
            ax.set_ylabel(label, fontsize=label_fontsize)

def plot_profiles(profiles, **kwargs):
    """
    Plot the profiles of the simulations.

    Parameters
    ----------
    profiles : list of dict
        List of profile dictionaries to plot.
    **kwargs : dict
        Additional keyword arguments to pass to the plotting function.

    Returns
    -------
    None
    """
    export = kwargs.get('export', None)
    xlabel = kwargs.get('xlabel', None)
    ylabel = kwargs.get('ylabel', None)
    title = kwargs.get('title', None)
    figsize = kwargs.get('figsize', (6, 4))
    dpi = kwargs.get('dpi', 300)
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)
    xlft = kwargs.get('xlft', kwargs.get('xleft', None))
    xrght = kwargs.get('xrght', kwargs.get('xright', None))
    ybot = kwargs.get('ybot', kwargs.get('ymin', None))
    ytop = kwargs.get('ytop', kwargs.get('ymax', None))
    xscale = kwargs.get('xscale', None)
    yscale = kwargs.get('yscale', None)
    xaxis_transformer = kwargs.get('xaxis_transformer', kwargs.get('xaxis_transform', None))
    yaxis_transformer = kwargs.get('yaxis_transformer', kwargs.get('yaxis_transform', None))
    transform_axes_data = kwargs.get('transform_axes_data', False)
    plotfarm = kwargs.get('plotfarm', False)
    farm_bounds = kwargs.get('farm_bounds', None)
    xfarm = kwargs.get('xfarm', None)
    Lfarm = kwargs.get('Lfarm', None)
    zeroline = kwargs.get('zeroline', False)
    grid = kwargs.get('grid', False)
    legend = kwargs.get('legend', True)
    legend_loc = kwargs.get('legend_loc', None)
    legend_ncols = kwargs.get('legend_ncols', 1)
    legend_frameon = kwargs.get('legend_frameon', False)
    legend_kwargs = kwargs.get('legend_kwargs', {})
    tight_layout = kwargs.get('tight_layout', True)
    close = kwargs.get('close', False)
    hide_top_right_spines = kwargs.get('hide_top_right_spines', False)
    tick_fontsize = kwargs.get('tick_fontsize', 10)
    label_fontsize = kwargs.get('label_fontsize', 11)
    title_fontsize = kwargs.get('title_fontsize', 9)
    legend_fontsize = kwargs.get('legend_fontsize', 10)
    linewidth = kwargs.get('linewidth', 1.2)
    linestyle = kwargs.get('linestyle', '-')
    marker = kwargs.get('marker', None)
    alpha = kwargs.get('alpha', None)
    color = kwargs.get('color', None)
    zorder = kwargs.get('zorder', None)
    ygain = kwargs.get('ygain', kwargs.get('yscale_factor', 1.0))

    x_transform = _axis_transform_function("x", xaxis_transformer) if transform_axes_data else None
    y_transform = _axis_transform_function("y", yaxis_transformer) if transform_axes_data else None

    with plt.rc_context({"text.usetex": True, "font.family": "serif"}):
        fig = plt.figure(figsize=figsize)
        ax = plt.axes()

        for profile in profiles:
            x,y = load_profile(profile['file'])
            profile_ygain = profile.get('ygain', profile.get('yscale_factor', ygain))
            y = np.asarray(y) * profile_ygain
            y = _apply_profile_offset(x, y, profile.get('offset', None))
            ax.plot(
                _transform_values(x, x_transform),
                _transform_values(y, y_transform),
                label=profile.get('label', profile.get('name', 'p')),
                color=profile.get('color', color),
                linestyle=profile.get('linestyle', profile.get('style', linestyle)),
                linewidth=profile.get('linewidth', linewidth),
                marker=profile.get('marker', marker),
                alpha=profile.get('alpha', alpha),
                zorder=profile.get('zorder', zorder),
            )

        if plotfarm:
            if farm_bounds is None:
                if xfarm is None or Lfarm is None:
                    raise ValueError("plotfarm=True requires either farm_bounds or both xfarm and Lfarm.")
                farm_bounds = (xfarm, xfarm + Lfarm)
            farm_left = _transform_scalar(farm_bounds[0], x_transform)
            farm_right = _transform_scalar(farm_bounds[1], x_transform)
            ax.axvspan(
                min(farm_left, farm_right),
                max(farm_left, farm_right),
                facecolor=kwargs.get('farm_facecolor', 'tab:red'),
                edgecolor=kwargs.get('farm_edgecolor', 'none'),
                alpha=kwargs.get('farm_alpha', 0.12),
                linewidth=kwargs.get('farm_linewidth', 0),
                zorder=kwargs.get('farm_zorder', 0),
            )

        if zeroline:
            ax.axhline(0, color=kwargs.get('zeroline_color', '0.4'), linewidth=kwargs.get('zeroline_width', 0.7), linestyle=kwargs.get('zeroline_style', '--'), zorder=0)

        if xscale is not None:
            ax.set_xscale(xscale)
        if yscale is not None:
            ax.set_yscale(yscale)

        if xlim is not None:
            ax.set_xlim((_transform_scalar(xlim[0], x_transform), _transform_scalar(xlim[1], x_transform)))
        if ylim is not None:
            ax.set_ylim((_transform_scalar(ylim[0], y_transform), _transform_scalar(ylim[1], y_transform)))
        if xlft is not None:
            ax.set_xlim(left=_transform_scalar(xlft, x_transform))
        if xrght is not None:
            ax.set_xlim(right=_transform_scalar(xrght, x_transform))
        if ybot is not None:
            ax.set_ylim(bottom=_transform_scalar(ybot, y_transform))
        if ytop is not None:
            ax.set_ylim(top=_transform_scalar(ytop, y_transform))

        if xlabel is not None:
            ax.set_xlabel(xlabel, fontsize=label_fontsize)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=label_fontsize)
        if title is not None:
            ax.set_title(title, fontsize=title_fontsize)

        _apply_axis_transformer(ax, "x", xaxis_transformer, transform_axes_data, label_fontsize)
        _apply_axis_transformer(ax, "y", yaxis_transformer, transform_axes_data, label_fontsize)

        ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
        if grid:
            ax.grid(True, alpha=kwargs.get('grid_alpha', 0.25), linewidth=kwargs.get('grid_linewidth', 0.5))
        if hide_top_right_spines:
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(top=False, right=False)
        if legend:
            ax.legend(
                loc=legend_loc,
                ncol=legend_ncols,
                frameon=legend_frameon,
                fontsize=legend_fontsize,
                **legend_kwargs,
            )
        if tight_layout:
            plt.tight_layout()
        if export is None:
            export = "streamwise_profiles.png"
        print(export)
        fig.savefig(export, dpi=dpi, bbox_inches=kwargs.get('bbox_inches', None), pad_inches=kwargs.get('pad_inches', 0.1))
        if close:
            plt.close(fig)
        return fig, ax

def plot_profile_panels(profile_panels, **kwargs):
    """
    Plot stacked profile panels.

    Each panel dictionary can override the global y-axis limits with
    ``ylim``, ``ybot``/``ymin``, or ``ytop``/``ymax``.
    """
    export = kwargs.get('export', None)
    xlabel = kwargs.get('xlabel', None)
    ylabel = kwargs.get('ylabel', None)
    figsize = kwargs.get('figsize', (5, 5))
    dpi = kwargs.get('dpi', 300)
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)
    xlft = kwargs.get('xlft', kwargs.get('xleft', None))
    xrght = kwargs.get('xrght', kwargs.get('xright', None))
    ybot = kwargs.get('ybot', kwargs.get('ymin', None))
    ytop = kwargs.get('ytop', kwargs.get('ymax', None))
    xscale = kwargs.get('xscale', None)
    yscale = kwargs.get('yscale', None)
    xaxis_transformer = kwargs.get('xaxis_transformer', kwargs.get('xaxis_transform', None))
    yaxis_transformer = kwargs.get('yaxis_transformer', kwargs.get('yaxis_transform', None))
    transform_axes_data = kwargs.get('transform_axes_data', False)
    plotfarm = kwargs.get('plotfarm', False)
    farm_bounds = kwargs.get('farm_bounds', None)
    xfarm = kwargs.get('xfarm', None)
    Lfarm = kwargs.get('Lfarm', None)
    zeroline = kwargs.get('zeroline', False)
    grid = kwargs.get('grid', False)
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
    title_fontsize = kwargs.get('title_fontsize', 10)
    title_loc = kwargs.get('title_loc', 'left')
    legend_fontsize = kwargs.get('legend_fontsize', 9)
    linewidth = kwargs.get('linewidth', 1.2)
    linestyle = kwargs.get('linestyle', '-')
    marker = kwargs.get('marker', None)
    alpha = kwargs.get('alpha', None)
    color = kwargs.get('color', None)
    zorder = kwargs.get('zorder', None)
    ygain = kwargs.get('ygain', kwargs.get('yscale_factor', 1.0))
    hspace = kwargs.get('hspace', 0.12)

    x_transform = _axis_transform_function("x", xaxis_transformer) if transform_axes_data else None
    y_transform = _axis_transform_function("y", yaxis_transformer) if transform_axes_data else None

    with plt.rc_context({"text.usetex": True, "font.family": "serif"}):
        fig, axes = plt.subplots(
            len(profile_panels),
            1,
            figsize=figsize,
            sharex=True,
            squeeze=False,
        )
        axes = axes[:, 0]

        for panel_index, (ax, panel) in enumerate(zip(axes, profile_panels)):
            panel_profiles = panel.get("profiles", panel) if isinstance(panel, dict) else panel
            panel_ylabel = panel.get("ylabel", ylabel) if isinstance(panel, dict) else ylabel
            panel_title = panel.get("title", None) if isinstance(panel, dict) else None
            #panel_label = panel.get("panel_label", f"({chr(ord('a') + panel_index)})") if isinstance(panel, dict) else f"({chr(ord('a') + panel_index)})"
            panel_label=None
            panel_ylim = panel.get("ylim", ylim) if isinstance(panel, dict) else ylim
            panel_ybot = panel.get("ybot", panel.get("ymin", ybot)) if isinstance(panel, dict) else ybot
            panel_ytop = panel.get("ytop", panel.get("ymax", ytop)) if isinstance(panel, dict) else ytop
            panel_legend_ncols = panel.get("legend_ncols", legend_ncols) if isinstance(panel, dict) else legend_ncols

            for profile in panel_profiles:
                x, y = load_profile(profile['file'])
                profile_ygain = profile.get('ygain', profile.get('yscale_factor', ygain))
                y = np.asarray(y) * profile_ygain
                y = _apply_profile_offset(x, y, profile.get('offset', None))
                ax.plot(
                    _transform_values(x, x_transform),
                    _transform_values(y, y_transform),
                    label=profile.get('label', profile.get('name', 'p')),
                    color=profile.get('color', color),
                    linestyle=profile.get('linestyle', profile.get('style', linestyle)),
                    linewidth=profile.get('linewidth', linewidth),
                    marker=profile.get('marker', marker),
                    alpha=profile.get('alpha', alpha),
                    zorder=profile.get('zorder', zorder),
                )

            if plotfarm:
                if farm_bounds is None:
                    if xfarm is None or Lfarm is None:
                        raise ValueError("plotfarm=True requires either farm_bounds or both xfarm and Lfarm.")
                    farm_bounds = (xfarm, xfarm + Lfarm)
                farm_left = _transform_scalar(farm_bounds[0], x_transform)
                farm_right = _transform_scalar(farm_bounds[1], x_transform)
                ax.axvspan(
                    min(farm_left, farm_right),
                    max(farm_left, farm_right),
                    facecolor=kwargs.get('farm_facecolor', 'tab:red'),
                    edgecolor=kwargs.get('farm_edgecolor', 'none'),
                    alpha=kwargs.get('farm_alpha', 0.12),
                    linewidth=kwargs.get('farm_linewidth', 0),
                    zorder=kwargs.get('farm_zorder', 0),
                )

            if zeroline:
                ax.axhline(0, color=kwargs.get('zeroline_color', '0.4'), linewidth=kwargs.get('zeroline_width', 0.7), linestyle=kwargs.get('zeroline_style', '--'), zorder=0)
            if xscale is not None:
                ax.set_xscale(xscale)
            if yscale is not None:
                ax.set_yscale(yscale)
            if xlim is not None:
                ax.set_xlim((_transform_scalar(xlim[0], x_transform), _transform_scalar(xlim[1], x_transform)))
            if panel_ylim is not None:
                ax.set_ylim((_transform_scalar(panel_ylim[0], y_transform), _transform_scalar(panel_ylim[1], y_transform)))
            if xlft is not None:
                ax.set_xlim(left=_transform_scalar(xlft, x_transform))
            if xrght is not None:
                ax.set_xlim(right=_transform_scalar(xrght, x_transform))
            if panel_ybot is not None:
                ax.set_ylim(bottom=_transform_scalar(panel_ybot, y_transform))
            if panel_ytop is not None:
                ax.set_ylim(top=_transform_scalar(panel_ytop, y_transform))

            if panel_ylabel is not None:
                ax.set_ylabel(panel_ylabel, fontsize=label_fontsize)
            if panel_label is not None or panel_title is not None:
                title_parts = [part for part in (panel_label, panel_title) if part]
                ax.set_title(" ".join(title_parts), fontsize=title_fontsize, loc=title_loc)

            _apply_axis_transformer(ax, "y", yaxis_transformer, transform_axes_data, label_fontsize)
            ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
            if panel_index == len(axes) - 1:
                if xlabel is not None:
                    ax.set_xlabel(xlabel, fontsize=label_fontsize)
                _apply_axis_transformer(ax, "x", xaxis_transformer, transform_axes_data, label_fontsize)
            else:
                ax.tick_params(axis='x', labelbottom=False)

            if grid:
                ax.grid(True, alpha=kwargs.get('grid_alpha', 0.25), linewidth=kwargs.get('grid_linewidth', 0.5))
            if hide_top_right_spines:
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.tick_params(top=False, right=False)
            if legend:
                ax.legend(
                    loc=legend_loc,
                    ncol=panel_legend_ncols,
                    frameon=legend_frameon,
                    fontsize=legend_fontsize,
                    **legend_kwargs,
                )

        if tight_layout:
            fig.tight_layout()
            fig.subplots_adjust(hspace=hspace)
        if export is None:
            export = "streamwise_profile_panels.png"
        print(export)
        fig.savefig(export, dpi=dpi, bbox_inches=kwargs.get('bbox_inches', None), pad_inches=kwargs.get('pad_inches', 0.1))
        if close:
            plt.close(fig)
        return fig, axes

profiles = [

    {
        'file': "/anvil/scratch/x-kali/PadeOpsSims/INV800-5K/profiles/Net_tkeadvect_rms_x.csv",
        'label': r'$\Delta \mathcal{A}$',
        'color': 'tab:red',
        'style': '-',
        'ygain': rmsgain*100000,
        "linewidth": lw,
        #'offset': 140,
    },

    #######

    {
        'file': "/anvil/scratch/x-kali/PadeOpsSims/INV800-5K/profiles/Net_PrsRedist_rms_x.csv",
        'label': r'$\Delta \mathcal{R}$',
        'color': 'tab:brown',
        'style': '-',
        'ygain': rmsgain*100000,
        "linewidth": lw,
        #'offset': 140,
    },

    #######

    {
        'file': "/anvil/scratch/x-kali/PadeOpsSims/INV800-5K/profiles/Net_Prod_rms_x.csv",
        'label': r'$\Delta \mathcal{P}$',
        'color': 'k',
        'style': '-',
        'ygain': rmsgain*100000,
        "linewidth": lw,
        #'offset': 140,
    },

     #######
    
    {
        'file': "/anvil/scratch/x-kali/PadeOpsSims/INV800-5K/profiles/Net_Bouyancy_rms_x.csv",
        'label': r'$\Delta \mathcal{B}$',
        'color': 'cyan',
        'style': '-',
        'ygain': rmsgain*100000,
        "linewidth": lw,
        #'offset': 140,
    },

    ####
    {
        'file': "/anvil/scratch/x-kali/PadeOpsSims/INV800-5K/profiles/Net_SGSTransport_rms_x.csv",
        'label': r'$\Delta \mathcal{D}$',
        'color': 'tab:orange',
        'style': '-',
        'ygain': rmsgain*100000,
        "linewidth": lw,
        #'offset': 120,
    },

    ####
    {
        'file': "/anvil/scratch/x-kali/PadeOpsSims/INV800-5K/profiles/Net_Transport_rms_x.csv",
        'label': r'$\Delta \mathcal{T}$',
        'color': 'tab:green',
        'style': '-',
        'ygain': rmsgain*100000,
        "linewidth": lw,
        'offset': 120,
    },

    ####

    {
        'file': "/anvil/scratch/x-kali/PadeOpsSims/INV800-5K/profiles/Net_Dissipation_rms_x.csv",
        'label': r'$\Delta \varepsilon$',
        'color': 'tab:blue',
        'style': '-',
        'ygain': rmsgain*100000,
        "linewidth": lw,
        #'offset': 120,
    },

]


axis_style = {
    "xlft": -0.5*Lfarm + xfarm,
    "xrght": 3.*Lfarm + xfarm,
    "xlabel": r"$x$",
    "ylabel": r"$\mathrm{RMS}~(\times 10^{-5})$",
    "xaxis_transformer": {
        "function": lambda x: (x - xfarm) / Lfarm,
        "inverse": lambda xt: xt * Lfarm + xfarm,
        "ticks": np.arange(-0.5, 3.1, 0.5),
        "label": r"$(x-x_0)/L_p$",
    },
    "transform_axes_data": False,
    "ytop": 6,
    "ybot": 0,
}

overlay_style = {
    "plotfarm": True,
    "xfarm": xfarm,
    "Lfarm": Lfarm,
}

plot_style = {
    "figsize": (5, 2.7),
    "dpi": 200,
    "hide_top_right_spines": True,
    "legend": True,
    "legend_loc": "best",
    "legend_frameon": False,
    "tight_layout": True,
    "legend_fontsize": 8,
    "hspace": 0.3,
    "export": "nettke.png",
}

if __name__ == "__main__":
    profile_panels = [
        {"profiles": profiles[:], "ytop":8, "title":"", "ylabel": r"$\mathrm{RMS}~(\times 10^{-5})$"},
    ]
    plot_profile_panels(
        profile_panels,
        **plot_style,
        **axis_style,
        **overlay_style,
    )
