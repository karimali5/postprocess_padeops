import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import warnings
from matplotlib.colors import TwoSlopeNorm, LogNorm
from matplotlib.ticker import FuncFormatter, MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.transforms import ScaledTranslation
import csv
from typing import Dict, Optional, Any, Tuple
from pathlib import Path
import netCDF4 as nc
import imageio.v2 as imageio
from scipy.ndimage import gaussian_filter
from matplotlib.patches import Rectangle

def _configure_latex_fonts():
    plt.rcdefaults()
    texbin = "/apps/spack/anvil/apps/texlive/20200406-gcc-11.2.0-eeavxnm/bin/x86_64-linux"
    if os.path.isdir(texbin):
        path = os.environ.get("PATH", "")
        path_parts = path.split(os.pathsep) if path else []
        if texbin not in path_parts:
            os.environ["PATH"] = os.pathsep.join([texbin, *path_parts])
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "text.latex.preamble": r"\usepackage{lmodern}",  # modern CM
    })

def _format_tick_max_two_decimals(value, _position=None):
    if abs(value) < 0.005:
        value = 0.0
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _plot_fontsize(kwargs, key, default):
    return kwargs.get(key, kwargs.get("fontsize", kwargs.get("font_size", default)))

def _add_matched_vertical_colorbar(fig, cf, ax, pad=0.08, size="2%"):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=size, pad=pad)
    return fig.colorbar(cf, cax=cax, orientation="vertical")

def _axis_transform_function(transformer):
    if transformer is None:
        return None

    transform = (
        transformer.get("function")
        or transformer.get("transform")
        or transformer.get("callable")
        or transformer.get("func")
    )

    if not callable(transform):
        raise TypeError(
            "axis transformer must include a callable under "
            "'function', 'transform', 'callable', or 'func'."
        )

    return transform

def streamwise_bulk_profiles(fields, **kwargs):
    _configure_latex_fonts()

    export = kwargs.get("export", None)
    plot_farm = kwargs.get("plotfarm", kwargs.get("plot_farm", False))
    farm_bounds = kwargs.get("farm_bounds", kwargs.get("farmsize", None))
    zeroline = kwargs.get("zeroline", False)
    scale = kwargs.get("scale", 1.0)
    xlim = kwargs.get("xlim", None)
    xleft = kwargs.get("xlft", kwargs.get("xleft", None))
    xright = kwargs.get("xrght", kwargs.get("xright", None))
    ylim = kwargs.get("ylim", None)
    ybot = kwargs.get("ybot", kwargs.get("ymin", None))
    ytop = kwargs.get("ytop", kwargs.get("ymax", None))
    xaxis_transformer = kwargs.get("xaxis_transformer", kwargs.get("xaxis_transform", None))
    x_transform = _axis_transform_function(xaxis_transformer)
    xstep = kwargs.get("xstep", None)
    xlabel = kwargs.get("xlabel", None)
    ylabel = kwargs.get("ylabel", None)
    title = kwargs.get("title", None)
    data_key = kwargs.get("data_key", "data")
    skip_nan = kwargs.get("skip_nan", True)
    inclusive = kwargs.get("inclusive", True)
    profile_value = kwargs.get("profile_value", kwargs.get("profile_stat", "sum"))
    bounds_default = kwargs.get("bounds", None)
    fieldgain_default = kwargs.get("fieldgain", 1)
    lengthgain_default = kwargs.get("lengthgain", 1)
    smooth_default = kwargs.get("smooth", False)
    l_smooth_default = kwargs.get("l_smooth", 5)
    linewidth_default = kwargs.get("linewidth", kwargs.get("lw", 1.2))
    linestyle_default = kwargs.get("linestyle", kwargs.get("ls", "-"))
    marker_default = kwargs.get("marker", None)
    alpha_default = kwargs.get("alpha", None)
    color_default = kwargs.get("color", None)
    zorder_default = kwargs.get("zorder", None)
    figsize = kwargs.get("figsize", (6, 4))
    dpi = kwargs.get("dpi", 300)
    legend = kwargs.get("legend", True)
    legend_ncols = kwargs.get("legend_ncols", kwargs.get("ncols", 1))
    legend_loc = kwargs.get("legend_loc", None)
    legend_frameon = kwargs.get("legend_frameon", False)
    legend_kwargs = kwargs.get("legend_kwargs", {})
    grid = kwargs.get("grid", False)
    tight_layout = kwargs.get("tight_layout", True)
    close = kwargs.get("close", True)
    tick_fontsize = _plot_fontsize(kwargs, "tick_fontsize", 6)
    label_fontsize = _plot_fontsize(kwargs, "label_fontsize", 7)
    title_fontsize = _plot_fontsize(kwargs, "title_fontsize", 8)
    legend_fontsize = _plot_fontsize(kwargs, "legend_fontsize", 6)

    if not fields:
        raise ValueError("fields must contain at least one field dictionary.")

    if profile_value not in ("mean", "sum"):
        raise ValueError("profile_value must be 'mean' or 'sum'.")

    fig, ax = plt.subplots(figsize=figsize)
    profiles = []

    for ifield, field in enumerate(fields):
        if not isinstance(field, dict):
            raise TypeError(f"Field {ifield} must be a dictionary.")

        bounds = field.get("bounds", bounds_default)
        label = field.get("label", field.get("name", None))
        files = field.get("files", field.get("filename", field.get("file", None)))

        if files is None:
            raise KeyError(f"Field {ifield} is missing required key 'files'.")
        if bounds is None:
            raise KeyError(f"Field {ifield} is missing required key 'bounds'.")

        data = read_netcdf_slice(
            files,
            varname=field.get("varname", kwargs.get("varname", None)),
            indexing=field.get("indexing", kwargs.get("indexing", "ij")),
        )

        fieldgain = field.get("fieldgain", fieldgain_default)
        lengthgain = field.get("lengthgain", lengthgain_default)
        if fieldgain != 1:
            data[data_key] *= fieldgain
        if lengthgain != 1:
            data["x1"] *= lengthgain
            data["x2"] *= lengthgain
            data["X1"] *= lengthgain
            data["X2"] *= lengthgain

        if field.get("smooth", smooth_default):
            data = gaussian_filter_slice(data, field.get("l_smooth", l_smooth_default))

        results = average_slice_over_bounds(
            data,
            bounds,
            data_key=data_key,
            skip_nan=field.get("skip_nan", skip_nan),
            inclusive=field.get("inclusive", inclusive),
        )
        profiles.append(results)

        if results["output_type"] == "profile":
            if x_transform is not None:
                xaxis = [x_transform(x_) for x_ in results["coord"]]
            else:
                xaxis = results["coord"]

            ax.plot(
                xaxis,
                results[profile_value] * scale,
                label=label,
                color=field.get("color", color_default),
                linewidth=field.get("linewidth", field.get("lw", linewidth_default)),
                linestyle=field.get("linestyle", field.get("ls", linestyle_default)),
                marker=field.get("marker", marker_default),
                alpha=field.get("alpha", alpha_default),
                zorder=field.get("zorder", zorder_default),
            )

    if plot_farm and farm_bounds is not None:
        if len(farm_bounds) >= 4:
            farm_bounds = (farm_bounds[0], farm_bounds[0] + farm_bounds[2])
        if x_transform is not None:
            farm_bounds_ = (x_transform(farm_bounds[0]), x_transform(farm_bounds[1]))
        else:
            farm_bounds_ = farm_bounds
        ax.axvspan(farm_bounds_[0], farm_bounds_[1], color='grey', alpha=0.2, zorder=-2)
    
    if zeroline:
        ax.axhline(y=0, color='k', alpha=0.2, linestyle='--', zorder=-1)

    if xlim is not None:
        if x_transform is not None:
            xlim_ = (x_transform(xlim[0]), x_transform(xlim[1]))
        else:
            xlim_ = xlim
        ax.set_xlim(xlim_)
    if xleft is not None:
        ax.set_xlim(left=x_transform(xleft) if x_transform is not None else xleft)
    if xright is not None:
        ax.set_xlim(right=x_transform(xright) if x_transform is not None else xright)

    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    elif xaxis_transformer is not None and xaxis_transformer.get("label", None) is not None:
        ax.set_xlabel(xaxis_transformer["label"], fontsize=label_fontsize)
    else:
        ax.set_xlabel('x', fontsize=label_fontsize)

    if xstep is not None:
        ax.xaxis.set_major_locator(MultipleLocator(xstep))
    elif xaxis_transformer is not None and xaxis_transformer.get("ticks", None) is not None:
        ticks = xaxis_transformer["ticks"]
        tick_positions = [x_transform(tick) for tick in ticks] if x_transform is not None else ticks
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([f"{tick:g}" for tick in ticks])

    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=label_fontsize)

    if ylim is not None:
        ax.set_ylim(ylim)
    if ybot is not None:
        ax.set_ylim(bottom=ybot)
    if ytop is not None:
        ax.set_ylim(top=ytop)

    if title is not None:
        ax.set_title(title, fontsize=title_fontsize)
    if grid:
        ax.grid(True, alpha=0.25, linewidth=0.5)
    if legend:
        legend_options = {
            "frameon": legend_frameon,
            "ncols": legend_ncols,
            "fontsize": legend_fontsize,
        }
        if legend_loc is not None:
            legend_options["loc"] = legend_loc
        legend_options.update(legend_kwargs)
        ax.legend(**legend_options)

    ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)

    if tight_layout:
        plt.tight_layout()

    if export is None:
        export = "streamwise_profiles.png"
    print(export)
    fig.savefig(
        export,
        dpi=dpi,
        bbox_inches=kwargs.get("bbox_inches", "tight"),
        pad_inches=kwargs.get("pad_inches", 0.02),
    )
    if close:
        plt.close(fig)

    return profiles

def average_slice_over_bounds(
    d: Dict[str, Any],
    bounds: Dict[str, Tuple[float, float]],
    data_key: str = "data",
    skip_nan: bool = True,
    inclusive: bool = True,
) -> Dict[str, Any]:
    """
    Average a 2D slice over bounds applied to one or both dimensions.

    If bounds are provided for both dimensions, returns a scalar average.

    If bounds are provided for only one dimension, returns a 1D average profile
    along the unbounded dimension.

    Parameters
    ----------
    d : dict
        Dictionary returned by read_slice_netcdf_for_contour().
        Must contain "x1", "x2", "X1", "X2", data_key, and "dims".

    bounds : dict
        Bounds keyed by actual dimension names stored in d["dims"].

        Examples
        --------
        {"y": (0.0, 1000.0), "z": (50.0, 200.0)}
        {"z": (50.0, 200.0)}

        The order of keys does not matter.

    data_key : str
        Key for the 2D field. Default is "data".

    skip_nan : bool
        If True, ignore NaNs in the average.
        If False, any genuine NaN inside the selected region propagates.

    inclusive : bool
        If True, use inclusive bounds.
        If False, use strict bounds.

    Returns
    -------
    dict
        For two bounded dimensions:
            output_type = "scalar"

        For one bounded dimension:
            output_type = "profile"
    """

    if "dims" not in d:
        raise KeyError("Input dictionary must contain key 'dims'.")

    if len(d["dims"]) != 2:
        raise ValueError(f"d['dims'] must contain exactly two names. Got {d['dims']}.")

    dim1, dim2 = tuple(d["dims"])
    valid_dims = {dim1, dim2}

    if len(bounds) == 0:
        raise ValueError("At least one bounded dimension must be provided.")

    if len(bounds) > 2:
        raise ValueError("This function only supports 2D slices.")

    unknown_dims = [dim for dim in bounds if dim not in valid_dims]
    if unknown_dims:
        raise KeyError(
            f"Bounds contain unknown dimension(s): {unknown_dims}. "
            f"Valid dimensions are {tuple(d['dims'])}."
        )

    x1 = np.asarray(d["x1"])
    x2 = np.asarray(d["x2"])

    X1 = np.asarray(d["X1"])
    X2 = np.asarray(d["X2"])
    data = np.asarray(d[data_key])

    if X1.shape != data.shape or X2.shape != data.shape:
        raise ValueError(
            f"Shape mismatch: X1{X1.shape}, X2{X2.shape}, data{data.shape}. "
            "Expected X1, X2, and data to have the same shape."
        )

    # Infer which dimension corresponds to each data axis.
    #
    # For indexing="ij":
    #     data.shape = (len(x1), len(x2))
    #     axis_dims = (dim1, dim2)
    #
    # For indexing="xy":
    #     data.shape = (len(x2), len(x1))
    #     axis_dims = (dim2, dim1)
    if data.shape == (x1.size, x2.size):
        axis_dims = (dim1, dim2)
    elif data.shape == (x2.size, x1.size):
        axis_dims = (dim2, dim1)
    else:
        raise ValueError(
            "Could not infer data-axis order from coordinate lengths. "
            f"data.shape={data.shape}, len(x1)={x1.size}, len(x2)={x2.size}."
        )

    axis_of_dim = {
        axis_dims[0]: 0,
        axis_dims[1]: 1,
    }

    coord_1d = {
        dim1: x1,
        dim2: x2,
    }

    coord_2d = {
        dim1: X1,
        dim2: X2,
    }

    mask = np.ones(data.shape, dtype=bool)
    normalized_bounds = {}

    for dim, bound_pair in bounds.items():
        lower, upper = sorted(bound_pair)
        C = coord_2d[dim]

        if inclusive:
            mask &= (C >= lower) & (C <= upper)
        else:
            mask &= (C > lower) & (C < upper)

        normalized_bounds[dim] = (lower, upper)

    reduced_dims = tuple(dim for dim in (dim1, dim2) if dim in bounds)
    remaining_dims = tuple(dim for dim in (dim1, dim2) if dim not in bounds)

    # ------------------------------------------------------------------
    # Case 1: both dimensions bounded -> scalar average
    # ------------------------------------------------------------------
    if len(remaining_dims) == 0:
        values = data[mask]

        if values.size == 0:
            raise ValueError(
                "No grid points found inside the requested rectangle. "
                f"Requested bounds: {normalized_bounds}"
            )

        if skip_nan:
            values = values[~np.isnan(values)]
            if values.size == 0:
                raise ValueError("All selected values inside the bounds are NaN.")

            total = np.sum(values)
            count = values.size
            mean = total / count

        else:
            total = np.sum(values)
            count = values.size
            mean = total / count

        return {
            "output_type": "scalar",
            "mean": mean,
            "sum": total,
            "count": count,
            "mask": mask,
            "bounds": normalized_bounds,
            "reduced_dims": reduced_dims,
            "remaining_dims": remaining_dims,
            "dims": (dim1, dim2),
            "axis_dims": axis_dims,
        }

    # ------------------------------------------------------------------
    # Case 2: one dimension bounded -> profile along remaining dimension
    # ------------------------------------------------------------------
    bounded_dim = reduced_dims[0]
    remaining_dim = remaining_dims[0]

    reduce_axis = axis_of_dim[bounded_dim]

    if skip_nan:
        selected = np.where(mask, data, np.nan)

        total = np.nansum(selected, axis=reduce_axis)
        count = np.sum(~np.isnan(selected), axis=reduce_axis)

        with np.errstate(invalid="ignore", divide="ignore"):
            mean = total / count

    else:
        # Sum only in-bounds cells without allowing out-of-bounds cells
        # to inject NaNs. Genuine NaNs inside the selected slab should propagate.
        data_zero_outside = np.where(mask, data, 0.0)

        total_no_nan = np.sum(
            np.where(np.isnan(data_zero_outside), 0.0, data_zero_outside),
            axis=reduce_axis,
        )

        has_nan_in_bounds = np.any(mask & np.isnan(data), axis=reduce_axis)

        total = np.where(has_nan_in_bounds, np.nan, total_no_nan)
        count = np.sum(mask, axis=reduce_axis)

        with np.errstate(invalid="ignore", divide="ignore"):
            mean = total / count

    if np.all(count == 0):
        raise ValueError(
            "No grid points found inside the requested slab. "
            f"Requested bounds: {normalized_bounds}"
        )

    coord = coord_1d[remaining_dim]

    # Sanity check: profile length should match remaining coordinate length.
    if mean.shape[0] != coord.size:
        raise RuntimeError(
            "Internal shape mismatch in profile output. "
            f"Profile length is {mean.shape[0]}, but len({remaining_dim})={coord.size}. "
            f"axis_dims={axis_dims}, reduce_axis={reduce_axis}."
        )

    return {
        "output_type": "profile",
        "mean": mean,
        "sum": total,
        "count": count,
        "coord": coord,
        "coord_name": remaining_dim,
        "mask": mask,
        "bounds": normalized_bounds,
        "reduced_dims": reduced_dims,
        "remaining_dims": remaining_dims,
        "dims": (dim1, dim2),
        "axis_dims": axis_dims,
    }

def inversion_height(sim, **kwargs):
    """
    Return an inversion-height slice and a matching colorbar label.

    Parameters
    ----------
    sim : dict
        Simulation dictionary with keys "path" and "timestamp".
    RID : int, optional
        Run ID to read. Default is 9.
    BRID : int, optional
        Baseline run ID used for deficits and percentages. Default is 8.
    mode : {"h0", "h1", "h2", "dh"}, optional
        Inversion-height field to construct. "h1" is the midpoint between
        INVH0 and INVH2, and "dh" is INVH2 - INVH0.

    Returns
    -------
    tuple
        (slice_dict, label)
    """

    def _inversion_height_filename(sim, rid, suffix):
        return os.path.join(
            sim["path"],
            f"Run{rid:02d}_t{sim['timestamp']}_INVH{suffix}.nc",
        )

    def _read_inversion_height(sim, rid, mode, scale):
        h0 = read_netcdf_slice(_inversion_height_filename(sim, rid, "0"), scale=scale)

        if mode == "h0":
            return h0

        h2 = read_netcdf_slice(_inversion_height_filename(sim, rid, "2"), scale=scale)
        _validate_compatible_slices(h0, h2, "combine inversion-height")

        full = _copy_slice_dict(h0)
        if mode == "h1":
            full["data"] = 0.5 * (h0["data"] + h2["data"])
        elif mode == "h2":
            full["data"] = h2["data"]
        elif mode == "dh":
            full["data"] = h2["data"] - h0["data"]
        return full

    rid = kwargs.get("RID", kwargs.get("rid", 9))
    brid = kwargs.get("BRID", kwargs.get("brid", 8))
    mode = kwargs.get("mode", "h1")
    deficit = kwargs.get("deficit", True)
    percentage = kwargs.get("percentage", True)
    smooth = kwargs.get("smooth", True)
    l_smooth = kwargs.get("l_smooth", 3.5)
    scale = kwargs.get("scale", 1)

    valid_modes = ("h0", "h1", "h2", "dh")
    if mode not in valid_modes:
        raise ValueError(f"mode must be one of {valid_modes}. Got {mode!r}.")

    if mode == "dh":
        label = r"\eta"
        reference_label = r"\eta_0"
    else:
        label = rf"h_{mode[1]}"
        reference_label = rf"h_{{{mode[1]}, 0}}"

    full = _read_inversion_height(sim, rid, mode, scale)
    base = None

    if deficit or percentage:
        base = _read_inversion_height(sim, brid, mode, scale)
        _validate_compatible_slices(full, base, "compare inversion-height")

    if deficit:
        full["data"] = full["data"] - base["data"]
        label = r"\Delta " + label

    if percentage:
        with np.errstate(invalid="ignore", divide="ignore"):
            full["data"] = full["data"] / base["data"] * 100.0
        label = label + rf"/ {reference_label} (\%)"

    if smooth:
        full = gaussian_filter_slice(full, l_smooth=l_smooth)

    return full, "$" + label + "$"

def transect_netcdf_slice(file, axis_name, value, **kwargs):
    """
    Return a 1D transect from a 2D slice by fixing one axis at a specified value.

    Parameters
    ----------
    slice_dict : dict
        Dictionary returned by read_slice_netcdf_for_contour(..., indexing="ij").
        Must contain "x1", "x2", "data", and "dims".
    axis_name : str
        Name of the axis to fix, e.g. "x", "y", or "z".
    value : float
        Coordinate value at which to take the transect.
    bounds_error : bool, optional
        If True, raise an error when value is outside the coordinate range.
        If False, clamp to the nearest boundary.

    Returns
    -------
    dict
        Dictionary containing:
          - "coord": 1D coordinate along the transect
          - "values": interpolated 1D field values
          - "fixed_axis": axis_name
          - "fixed_value": requested value
          - "transect_axis": remaining axis name
          - "dims": tuple describing returned transect
    """

    bounds_error = kwargs.get("bounds_error", True)
    plot_farm = kwargs.get("plot_farm", False)
    farm_bounds = kwargs.get("farm_bounds", None)
    export = kwargs.get("export", None)
    ylabel = kwargs.get("ylabel", None)
    zeroline = kwargs.get("zeroline", False)
    xleft = kwargs.get("xleft", None)
    xright = kwargs.get("xright", None)
    plot = kwargs.get("plot", True)
    scale = kwargs.get("scale", 1.0)
    smooth = kwargs.get("smooth", False)
    l_smooth = kwargs.get("l_smooth", 1.0)

    _configure_latex_fonts()

    slice_dict = read_netcdf_slice(file)
    if smooth:
        slice_dict = gaussian_filter_slice(slice_dict, l_smooth)

    dims = tuple(slice_dict["dims"])
    arr = np.asarray(slice_dict["data"], dtype=float)

    if len(dims) != 2:
        raise ValueError(f"Expected a 2D slice with two dims, got dims={dims}")

    if axis_name not in dims:
        raise ValueError(
            f"Axis {axis_name!r} not found in slice dims {dims}. "
            f"Available axes are {dims}."
        )

    if arr.shape != (len(slice_dict["x1"]), len(slice_dict["x2"])):
        raise ValueError(
            "Expected data shape to match (len(x1), len(x2)). "
            f"Got data.shape={arr.shape}, len(x1)={len(slice_dict['x1'])}, "
            f"len(x2)={len(slice_dict['x2'])}. "
            "Use read_slice_netcdf_for_contour(..., indexing='ij')."
        )

    fixed_dim_index = dims.index(axis_name)

    if fixed_dim_index == 0:
        fixed_coord = np.asarray(slice_dict["x1"], dtype=float)
        transect_coord = np.asarray(slice_dict["x2"], dtype=float)
        transect_axis = dims[1]
    else:
        fixed_coord = np.asarray(slice_dict["x2"], dtype=float)
        transect_coord = np.asarray(slice_dict["x1"], dtype=float)
        transect_axis = dims[0]

    vmin = np.nanmin(fixed_coord)
    vmax = np.nanmax(fixed_coord)

    if bounds_error and not (vmin <= value <= vmax):
        raise ValueError(
            f"Requested {axis_name}={value} is outside coordinate range "
            f"[{vmin}, {vmax}]."
        )

    value_clamped = np.clip(value, vmin, vmax)

    # Handle decreasing coordinate arrays by reversing interpolation direction.
    increasing = fixed_coord[0] <= fixed_coord[-1]

    if increasing:
        coord_interp = fixed_coord
        arr_interp = arr
    else:
        coord_interp = fixed_coord[::-1]
        if fixed_dim_index == 0:
            arr_interp = arr[::-1, :]
        else:
            arr_interp = arr[:, ::-1]

    # Exact or interpolated index location.
    hi = np.searchsorted(coord_interp, value_clamped)

    if hi == 0:
        lo = hi = 0
        weight = 0.0
    elif hi >= len(coord_interp):
        lo = hi = len(coord_interp) - 1
        weight = 0.0
    else:
        lo = hi - 1
        xlo = coord_interp[lo]
        xhi = coord_interp[hi]
        weight = (value_clamped - xlo) / (xhi - xlo)

    if fixed_dim_index == 0:
        values_lo = arr_interp[lo, :]
        values_hi = arr_interp[hi, :]
    else:
        values_lo = arr_interp[:, lo]
        values_hi = arr_interp[:, hi]

    values = (1.0 - weight) * values_lo + weight * values_hi

    if plot:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(transect_coord, values, c = 'k', linewidth=1.0)
        ax.set_xlabel(f"${transect_axis}$")
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        else:
            ax.set_ylabel(f"${slice_dict.get('varname', 'F')}$")
        
        if zeroline:
            ax.axhline(y=0, color='k', alpha=0.2, linestyle='--', zorder=-1)

        if xleft is not None:
            ax.set_xlim(left=xleft)
        if xright is not None:
            ax.set_xlim(right=xright)

        if plot_farm and farm_bounds is not None:
            ax.axvspan(farm_bounds[0], farm_bounds[1], color='grey', alpha=0.2, zorder=-2)
        if export is None:
            export = f"transect_{axis_name}.png"
        print(export)
        fig.savefig(export, dpi=300)
        plt.close(fig)

    return {
        "coord": transect_coord,
        "values": values*scale,
        "fixed_axis": axis_name,
        "fixed_value": value,
        "fixed_value_used": value_clamped,
        "transect_axis": transect_axis,
        "dims": (transect_axis,),
    }

def interpolate_streamtube(tube, x_target):

    def _interpolate_streamtube_npz(trajectories, x_query):
        """
        Interpolate all streamtube points at a requested x-location.

        Parameters
        ----------
        trajectories : ndarray
            Array of shape (npoints, 3, nxstations).
        x_query : float
            Streamwise location where y and z are requested.

        Returns
        -------
        yq : ndarray
            Interpolated y-coordinates of shape (npoints,).
        zq : ndarray
            Interpolated z-coordinates of shape (npoints,).
        """

        trajectories = np.asarray(trajectories["trajectories"])

        if trajectories.ndim != 3:
            raise ValueError("trajectories must have shape (npoints, 3, nxstations).")

        if trajectories.shape[1] != 3:
            raise ValueError("Second dimension of trajectories must have size 3 for x, y, z.")

        x = trajectories[:, 0, :]
        y = trajectories[:, 1, :]
        z = trajectories[:, 2, :]

        npoints = trajectories.shape[0]

        yq = np.empty(npoints, dtype=float)
        zq = np.empty(npoints, dtype=float)

        for p in range(npoints):
            xp = x[p, :]
            yp = y[p, :]
            zp = z[p, :]

            # Remove invalid / non-finite trajectory entries if present
            valid = np.isfinite(xp) & np.isfinite(yp) & np.isfinite(zp)

            if np.count_nonzero(valid) < 2:
                yq[p] = np.nan
                zq[p] = np.nan
                continue

            xp = xp[valid]
            yp = yp[valid]
            zp = zp[valid]

            # np.interp requires increasing x
            order = np.argsort(xp)
            xp = xp[order]
            yp = yp[order]
            zp = zp[order]

            if x_query < xp[0] or x_query > xp[-1]:
                yq[p] = np.nan
                zq[p] = np.nan
                continue

            yq[p] = np.interp(x_query, xp, yp)
            zq[p] = np.interp(x_query, xp, zp)

        return yq, zq

    def _interpolate_streamtube_csv(data, x_target):
        """
        Interpolate streamtube coordinates at a target x station.

        Parameters
        ----------
        data : dict
            Dictionary of the form:

                data[x]["d1"] -> y coordinates
                data[x]["d2"] -> z coordinates

            Each x station must have the same number of contour points,
            ordered consistently around the contour.

        x_target : float
            Target x station where the contour is interpolated.

        Returns
        -------
        y_interp : numpy.ndarray
            Interpolated y coordinates.

        z_interp : numpy.ndarray
            Interpolated z coordinates.
        """

        keys = [k for k in data.keys() if isinstance(k, (int, float))]
        x_stations = np.array(sorted(keys), dtype=float)

        if x_target < x_stations[0] or x_target > x_stations[-1]:
            warnings.warn(
                f"x_target={x_target} is outside the available range "
                f"[{x_stations[0]}, {x_stations[-1]}]. Returning None.",
                RuntimeWarning,
            )
            return None, None

        # Exact match
        if x_target in data:
            return (
                np.array(data[x_target]["d1"], dtype=float),
                np.array(data[x_target]["d2"], dtype=float),
            )

        # Find neighboring stations
        i_right = np.searchsorted(x_stations, x_target)
        i_left = i_right - 1

        x_left = x_stations[i_left]
        x_right = x_stations[i_right]

        y_left = np.array(data[x_left]["d1"], dtype=float)
        z_left = np.array(data[x_left]["d2"], dtype=float)

        y_right = np.array(data[x_right]["d1"], dtype=float)
        z_right = np.array(data[x_right]["d2"], dtype=float)

        if y_left.shape != y_right.shape or z_left.shape != z_right.shape:
            raise ValueError(
                f"Contour point count mismatch between x={x_left} and x={x_right}."
            )

        weight = (x_target - x_left) / (x_right - x_left)

        y_interp = (1.0 - weight) * y_left + weight * y_right
        z_interp = (1.0 - weight) * z_left + weight * z_right

        return y_interp, z_interp

    if not isinstance(tube, dict):
        raise TypeError("tube must be a streamtube dictionary returned by read_streamtube().")

    source = tube.get("source")
    if source == "csv":
        return _interpolate_streamtube_csv(tube, x_target)
    elif source == "npz":
        return _interpolate_streamtube_npz(tube, x_target)
    else:
        raise ValueError(f"Unsupported streamtube source '{source}'.")

def gaussian_filter_slice(slice_data, l_smooth = 3.5):
    smoothed = _copy_slice_dict(slice_data)
    dx = np.mean(np.diff(slice_data["x1"]))
    dy = np.mean(np.diff(slice_data["x2"]))
    sigma_x = l_smooth / dx
    sigma_y = l_smooth / dy
    sigma = (sigma_x, sigma_y)  # tune in grid-cell units
    d = gaussian_filter(slice_data['data'], sigma=sigma, mode="nearest")
    smoothed['data'] = d
    return smoothed

def _copy_slice_dict(slice_data):
    # Keep caller-owned slice dictionaries immutable from utility operations.
    copied = slice_data.copy()
    for key in ("x1", "x2", "X1", "X2", "data"):
        if key in copied:
            copied[key] = np.array(copied[key], copy=True)
    return copied

def read_netcdf_slice(filename, **kwargs):
    read_kwargs_base = kwargs.copy()
    smooth = read_kwargs_base.pop("smooth", False)
    l_smooth = read_kwargs_base.pop("l_smooth", 3.5)

    def _read_netcdf(
        path: str,
        varname: Optional[str] = None,
        indexing: str = "ij",
        scale = 1,
    ) -> Dict[str, Any]:
        """
        Read a 2D slice NetCDF written by export_slice_to_netcdf() and return
        a dict with X1, X2, data ready for contour/contourf.

        Parameters
        ----------
        path : str
            NetCDF file path.
        varname : str | None
            Name of the 2D field variable. If None, auto-detects the first 2D
            variable whose dimensions match (x1_name, x2_name) or (x2_name, x1_name).
        x1_name, x2_name : str
            Coordinate variable names (defaults match Fortran defaults).
        indexing : {"ij","xy"}
            Passed to numpy.meshgrid. Use "ij" to match Fortran array ordering
            (slice(nx, ny) written as dims [x1, x2]). Use "xy" if you prefer the
            plotting convention.

        Returns
        -------
        dict with keys:
        - "x1", "x2": 1D coordinate arrays
        - "X1", "X2": 2D meshgrid arrays
        - "data": 2D array aligned with X1/X2
        - "varname": resolved variable name
        - "dims": tuple of dimension names for the returned data
        """
        with nc.Dataset(path, "r") as ds:
            dims = list(ds.dimensions.keys())

            x1_name = dims[0]
            x2_name = dims[1]

            x1 = np.asarray(ds.variables[x1_name][:])
            x2 = np.asarray(ds.variables[x2_name][:])

            # Resolve data variable
            if varname is None:
                cand = []
                for name, v in ds.variables.items():
                    if name in (x1_name, x2_name):
                        continue
                    if getattr(v, "ndim", 0) != 2:
                        continue
                    dims = tuple(getattr(v, "dimensions", ()))
                    # Prefer exact match to Fortran dim order [x1, x2]
                    if dims == (x1_name, x2_name):
                        cand.insert(0, name)
                    elif dims == (x2_name, x1_name):
                        cand.append(name)

                if not cand:
                    raise KeyError(
                        "Could not auto-detect a 2D field variable. "
                        f"Please pass varname explicitly. Variables: {list(ds.variables.keys())}"
                    )
                varname = cand[0]

            if varname not in ds.variables:
                raise KeyError(
                    f"Variable '{varname}' not found in {path}. "
                    f"Available: {list(ds.variables.keys())}"
                )

            v = ds.variables[varname]
            data = np.asarray(v[:, :])
            dims = tuple(v.dimensions)

            # Ensure data is aligned as (x1, x2) to match meshgrid('ij') and Fortran write
            if dims == (x1_name, x2_name):
                data_x1x2 = data
            elif dims == (x2_name, x1_name):
                data_x1x2 = data.T
                dims = (x1_name, x2_name)
            else:
                raise ValueError(
                    f"Field '{varname}' has dims {dims}, expected ({x1_name},{x2_name}) "
                    f"or ({x2_name},{x1_name})."
                )

            # Build meshgrid compatible with contour/contourf
            X1, X2 = np.meshgrid(x1, x2, indexing=indexing)
            # For indexing="ij": X1,X2 shape = (len(x1), len(x2)) -> matches data_x1x2
            # For indexing="xy": X1,X2 shape = (len(x2), len(x1)) -> need transpose
            if indexing == "xy":
                data_out = data_x1x2.T
            elif indexing == "ij":
                data_out = data_x1x2
            else:
                raise ValueError("indexing must be 'ij' or 'xy'")

            return {
                "x1": x1,
                "x2": x2,
                "X1": X1,
                "X2": X2,
                "data": data_out * scale,
                "varname": varname,
                "dims": dims,
            }
    
    if isinstance(filename, list):
        if not filename:
            raise ValueError("filename list must contain at least one slice.")
        scales = read_kwargs_base.get("scale", 1)
        if np.isscalar(scales):
            scales = [scales] * len(filename)
        elif len(scales) != len(filename):
            raise ValueError(
                "When filename is a list, scale must be a scalar or a list "
                "with the same length as filename."
            )
        first = True
        for file_, scale_ in zip(filename, scales):
            if isinstance(file_, dict):
                d_ = _copy_slice_dict(file_)
                d_["data"] *= scale_
            else:
                read_kwargs = read_kwargs_base.copy()
                read_kwargs["scale"] = scale_
                d_ = _read_netcdf(file_, **read_kwargs)
            if first:
                d = _copy_slice_dict(d_)
                first = False
            else:
                if d['data'].shape != d_['data'].shape:
                    raise ValueError(
                        f"Cannot add slices with shapes {d['data'].shape} and {d_['data'].shape}."
                    )
                if not (np.array_equal(d['x1'], d_['x1']) and np.array_equal(d['x2'], d_['x2'])):
                    raise ValueError("Cannot add slices with different x1/x2 coordinates.")
                d['data'] += d_['data']            
    elif isinstance(filename, str):
        d = _read_netcdf(filename, **read_kwargs_base)
    elif isinstance(filename, dict):
        d = _copy_slice_dict(filename)
    else:
        raise TypeError("filename must be a path, list of paths/dicts, or slice dictionary")
    if smooth:
        d = gaussian_filter_slice(d, l_smooth=l_smooth)
    return d

def read_streamtube(file):
    file = os.fspath(file)

    def _read_streamtube_csv(filename):
        """
        Read a CSV file with header:
            x,y1,z1,y2,z2,...,yN,zN

        Returns
        -------
        data : dict
            Dictionary keyed by x station.

            data[x]["d1"] -> numpy array of y coordinates
            data[x]["d2"] -> numpy array of z coordinates
        """

        data = {"source": "csv"}

        with open(filename, "r", newline="") as f:
            reader = csv.reader(f)

            header = next(reader)

            if header[0].strip() != "x":
                raise ValueError("First column must be 'x'.")

            coord_headers = header[1:]

            if len(coord_headers) % 2 != 0:
                raise ValueError("Coordinate columns must come in y,z pairs.")

            npoints = len(coord_headers) // 2

            for row in reader:
                if not row:
                    continue

                x = float(row[0])

                if len(row) != 1 + 2 * npoints:
                    raise ValueError(
                        f"Row at x={x} has {len(row)} columns, "
                        f"but expected {1 + 2 * npoints}."
                    )

                y = []
                z = []

                for p in range(npoints):
                    y_val = float(row[1 + 2*p])
                    z_val = float(row[1 + 2*p + 1])

                    y.append(y_val)
                    z.append(z_val)

                data[x] = {
                    "d1": np.array(y),
                    "d2": np.array(z),
                }

        return data

    def _read_streamtube_npz(file):
        """
        Read streamtube trajectories from an NPZ file.

        Expected NPZ content
        --------------------
        trajectories : ndarray
            Array of shape (npoints, 3, nxstations), where

                trajectories[:, 0, :] = x
                trajectories[:, 1, :] = y
                trajectories[:, 2, :] = z

        Returns
        -------
        data : dict
            Dictionary containing:

            data["source"]       -> "npz"
            data["trajectories"] -> ndarray of shape (npoints, 3, nxstations)
        """

        npz_data = np.load(file)
        trajectories = npz_data["trajectories"]

        trajectories = np.asarray(trajectories, dtype=float)

        if trajectories.ndim != 3:
            raise ValueError("Stored trajectories must have shape (npoints, 3, nxstations).")

        if trajectories.shape[1] != 3:
            raise ValueError("Second dimension of trajectories must have size 3 for x, y, z.")

        data = {
            "source": "npz",
            "trajectories": trajectories,
        }

        return data

    if file.endswith(".csv"):
        return _read_streamtube_csv(file)
    elif file.endswith(".npz"):
        return _read_streamtube_npz(file)
    else:
        raise ValueError("Streamtube file must have extension '.csv' or '.npz'.")

def make_mp4(png_files, output_mp4, fps=10, quality=8):
    """
    Create an MP4 animation from an ordered list of PNG files.

    Parameters
    ----------
    png_files : list[str] or list[pathlib.Path]
        Ordered list of PNG files.
    output_mp4 : str or pathlib.Path
        Output MP4 filename.
    fps : int or float, optional
        Frames per second.
    quality : int, optional
        Encoding quality for imageio/ffmpeg. Higher is better.
    """
    png_files = [Path(f) for f in png_files]
    output_mp4 = Path(output_mp4)

    if not png_files:
        raise ValueError("png_files is empty.")

    for file in png_files:
        if not file.exists():
            raise FileNotFoundError(f"File not found: {file}")

    output_mp4.parent.mkdir(parents=True, exist_ok=True)

    with imageio.get_writer(
        output_mp4,
        fps=fps,
        codec="libx264",
        quality=quality,
        macro_block_size=1,
    ) as writer:
        for file in png_files:
            frame = np.asarray(imageio.imread(file))
            if frame.ndim == 2:
                frame = np.stack([frame, frame, frame], axis=-1)
            if frame.shape[-1] == 4:
                frame = frame[:, :, :3]
            height, width = frame.shape[:2]
            pad_height = height % 2
            pad_width = width % 2
            if pad_height or pad_width:
                frame = np.pad(
                    frame,
                    ((0, pad_height), (0, pad_width), (0, 0)),
                    mode="edge",
                )
            writer.append_data(frame)

    print(f"Saved MP4 to: {output_mp4}")

def make_gif(png_files, output_gif, fps=5):
    """
    Create a GIF animation from an ordered list of PNG files.

    Parameters
    ----------
    png_files : list[str] or list[pathlib.Path]
        Ordered list of PNG files.
    output_gif : str or pathlib.Path
        Output GIF filename.
    fps : int or float, optional
        Frames per second.
    """
    png_files = [Path(f) for f in png_files]
    output_gif = Path(output_gif)

    if not png_files:
        raise ValueError("png_files is empty.")
    if fps <= 0:
        raise ValueError("fps must be positive.")

    for file in png_files:
        if not file.exists():
            raise FileNotFoundError(f"File not found: {file}")

    output_gif.parent.mkdir(parents=True, exist_ok=True)
    duration = 1.0 / fps

    with imageio.get_writer(output_gif, mode="I", duration=duration) as writer:
        for file in png_files:
            writer.append_data(imageio.imread(file))

    print(f"Saved GIF to: {output_gif}")

def _validate_compatible_slices(primary, secondary, operation):
    if primary["data"].shape != secondary["data"].shape:
        raise ValueError(
            f"Cannot {operation} slices with shapes "
            f"{primary['data'].shape} and {secondary['data'].shape}."
        )
    if not (
        np.array_equal(primary["x1"], secondary["x1"])
        and np.array_equal(primary["x2"], secondary["x2"])
    ):
        raise ValueError(f"Cannot {operation} slices with different x1/x2 coordinates.")

def _prepare_plot_slice_data(filename, **kwargs):
    fieldgain = kwargs.get('fieldgain', 1)
    lengthgain = kwargs.get('lengthgain', 1)
    background = kwargs.get('background', None)
    percentage = kwargs.get('percentage', False)
    smooth = kwargs.get('smooth', False)
    l_smooth = kwargs.get('l_smooth', 5)

    d = read_netcdf_slice(filename)
    d["data"] *= fieldgain

    if background is not None:
        bg = read_netcdf_slice(background)
        bg["data"] *= fieldgain
        _validate_compatible_slices(d, bg, "subtract")

        d["data"] -= bg["data"]
        if percentage:
            d["data"] = d["data"] / bg["data"] * 100

    d["X1"] *= lengthgain
    d["X2"] *= lengthgain

    if smooth:
        d = gaussian_filter_slice(d, l_smooth = l_smooth)

    return d

def _make_shared_norm_and_levels(slices, **kwargs):
    deficit = kwargs.get('deficit', True)
    s = kwargs.get('s', None)
    cmap_reg = kwargs.get('cmap_reg', 'viridis')
    cmap_def = kwargs.get('cmap_def', 'PRGn')
    levels = kwargs.get('levels', 200)
    vmin = kwargs.get('vmin', None)
    vmax = kwargs.get('vmax', None)
    colorbar_range = kwargs.get('colorbar_range', 'norm')

    data_min = min(np.nanmin(d["data"]) for d in slices)
    data_max = max(np.nanmax(d["data"]) for d in slices)

    if deficit:
        scale = s
        if scale is None:
            lower = data_min if vmin is None else vmin
            upper = data_max if vmax is None else vmax
            scale = max(abs(lower), abs(upper))
        if scale == 0:
            scale = 1.0
        print(scale)
        norm = TwoSlopeNorm(vmin=-scale, vcenter=0.0, vmax=scale)
        level_count = levels if np.isscalar(levels) else len(levels)
        if colorbar_range == 'data':
            lower = data_min if vmin is None else vmin
            upper = data_max if vmax is None else vmax
            if lower == upper:
                lower -= 0.5
                upper += 0.5
            plot_levels = np.linspace(lower, upper, int(level_count))
        elif colorbar_range == 'norm':
            plot_levels = np.linspace(-scale, scale, int(level_count))
        else:
            raise ValueError("colorbar_range must be 'norm' or 'data'.")
        return cmap_def, norm, plot_levels, scale

    lower = data_min if vmin is None else vmin
    upper = data_max if vmax is None else vmax
    if lower == upper:
        lower -= 0.5
        upper += 0.5
    norm = TwoSlopeNorm(vmin=lower, vcenter=(lower + upper)/2, vmax=upper)
    return cmap_reg, norm, levels, None

def _format_simulation_label(index, name, annotation):
    if annotation is None or annotation is False:
        return None

    letter = chr(ord("a") + index)
    if annotation == "letter":
        return f"({letter})"
    if annotation == "name":
        return name
    if annotation == "letter_name":
        return f"({letter}) {name}" if name else f"({letter})"
    if callable(annotation):
        return annotation(index, name)
    return f"({letter}) {name}" if name else f"({letter})"

def _simulation_plot_kwargs(global_kwargs, sim):
    panel_kwargs = global_kwargs.copy()
    for key in ("fieldgain", "s"):
        if key in sim:
            panel_kwargs[key] = sim[key]
    panel_kwargs.update(sim.get("plot_kwargs", {}))
    panel_kwargs["blh_file"] = sim.get("blh_file", panel_kwargs.get("blh_file", None))
    return panel_kwargs

def _simulation_has_local_kw(sim, key):
    return key in sim or key in sim.get("plot_kwargs", {})

def _plot_slice_on_axis(ax, d, cmap, norm, plot_levels, simulation_label=None, **kwargs):
    lengthgain = kwargs.get('lengthgain', 1)
    equal_aspect = kwargs.get('equal_aspect', True)
    ytop = kwargs.get('ytop', None)
    ybot = kwargs.get('ybot', None)
    xlft = kwargs.get('xlft', None)
    xrght = kwargs.get('xrght', None)
    plotfarm = kwargs.get('plotfarm', False)
    farmsize = kwargs.get('farmsize', None)
    stamp = kwargs.get('stamp', None)
    plot_blh = kwargs.get('plot_blh', False)
    blh_file = kwargs.get('blh_file', None)
    xaxis_transformer = kwargs.get('xaxis_transformer', None)
    yaxis_transformer = kwargs.get('yaxis_transformer', None)
    xlabel = kwargs.get('xlabel', None)
    ylabel = kwargs.get('ylabel', None)
    label_fontsize = _plot_fontsize(kwargs, 'label_fontsize', 7)
    tick_fontsize = _plot_fontsize(kwargs, 'tick_fontsize', 6)
    stamp_fontsize = _plot_fontsize(kwargs, 'stamp_fontsize', 8)
    annotation_fontsize = _plot_fontsize(kwargs, 'annotation_fontsize', 6)
    stamp_location = kwargs.get('stamp_location', (0.78, 0.85))
    stamp_ha = kwargs.get('stamp_ha', 'left')
    stamp_va = kwargs.get('stamp_va', 'baseline')
    aspect = kwargs.get('aspect', None)
    contours = kwargs.get('contours', False)
    bounding_shell = kwargs.get('bounding_shell', None)
    bounding_shell_threshold = kwargs.get('bounding_shell_threshold', -0.01)
    bounding_shell_local = kwargs.get('bounding_shell_local', False)
    fixed_colorbar_ticks = kwargs.get('fixed_colorbar_ticks', False)
    plot_streamtube = kwargs.get('plot_streamtube', False)
    streamtube = kwargs.get('streamtube', None)
    station = kwargs.get('station', None)
    transform_axes_data = kwargs.get('transform_axes_data', False)

    def _get_transform(axis_name, transformer):
        if transformer is None:
            return None
        if not isinstance(transformer, dict):
            raise TypeError(f"{axis_name}axis_transformer must be a dictionary or None.")
        transform = (
            transformer.get("function")
            or transformer.get("transform")
            or transformer.get("callable")
        )
        if not callable(transform):
            raise TypeError(
                f"{axis_name}axis_transformer must include a callable under "
                "'function', 'transform', or 'callable'."
            )
        return transform

    x_transform = _get_transform("x", xaxis_transformer) if transform_axes_data else None
    y_transform = _get_transform("y", yaxis_transformer) if transform_axes_data else None

    def _transform_values(values, transform):
        if transform is None:
            return values
        transformed = transform(np.asarray(values))
        return np.asarray(transformed)

    def _transform_scalar(value, transform):
        if transform is None or value is None:
            return value
        return float(np.asarray(transform(np.asarray(value))))

    def _transform_rect(bounds):
        x0, y0, width, height = bounds
        x_left = _transform_scalar(x0, x_transform)
        x_right = _transform_scalar(x0 + width, x_transform)
        y_bottom = _transform_scalar(y0, y_transform)
        y_top = _transform_scalar(y0 + height, y_transform)
        return (
            min(x_left, x_right),
            min(y_bottom, y_top),
            abs(x_right - x_left),
            abs(y_top - y_bottom),
        )

    X1_plot = _transform_values(d["X1"], x_transform)
    X2_plot = _transform_values(d["X2"], y_transform)

    def _draw_field():
        if contours:
            return ax.contourf(
                X1_plot, X2_plot, d["data"],
                levels=plot_levels, cmap=cmap, norm=norm,
                extend=("both" if fixed_colorbar_ticks else None))
        return ax.pcolormesh(
            X1_plot, X2_plot, d["data"],
            shading='auto', cmap=cmap, norm=norm
        )

    def _add_bounding_shell():
        if not bounding_shell:
            return

        # Local threshold scales the contour level by the current field magnitude.
        if bounding_shell_local:
            threshold = bounding_shell_threshold * np.nanmax(np.abs(np.asarray(d["data"])))
        else:
            threshold = bounding_shell_threshold
        ax.contour(X1_plot, X2_plot, d["data"], levels=[threshold], colors='k', linewidths=0.5, alpha=0.2)

    def _add_streamtube():
        if not (plot_streamtube and streamtube is not None and station is not None):
            return

        tube = read_streamtube(streamtube["file"])
        ytube, ztube = interpolate_streamtube(tube, station)
        ax.plot(
            _transform_values(ytube*lengthgain, x_transform),
            _transform_values(ztube*lengthgain, y_transform),
            color=streamtube["color"],
            lw=streamtube["linewidth"],
        )

    def _add_farm():
        if not (plotfarm and farmsize is not None):
            return

        farm_x, farm_y, farm_width, farm_height = _transform_rect(farmsize)
        rect = Rectangle(
            (farm_x, farm_y),
            farm_width,
            farm_height,
            facecolor='none',
            edgecolor='tab:red',
            linestyle='--',
            linewidth=0.75,
            alpha=0.4,
        )
        ax.add_patch(rect)

    def _add_blh():
        if not (plot_blh and blh_file is not None):
            return

        for blh_ in blh_file:
            data = blh_['data']
            ax.plot(
                _transform_values(np.array(data['coord']), x_transform),
                _transform_values(np.array(data['values']), y_transform),
                color=blh_['color'],
                linestyle=blh_.get('style', '-'),
                linewidth=0.5,
                alpha=0.5,
            )

    def _apply_axes_limits_and_aspect():
        if(xlft is not None): ax.set_xlim(left=_transform_scalar(xlft, x_transform))
        if(xrght is not None): ax.set_xlim(right=_transform_scalar(xrght, x_transform))
        if(ybot is not None): ax.set_ylim(bottom=_transform_scalar(ybot, y_transform))
        if(ytop is not None): ax.set_ylim(top=_transform_scalar(ytop, y_transform))
        if(equal_aspect): ax.set_aspect("equal")
        if aspect is not None: ax.set_aspect(aspect)

    def _apply_axis_transformer(axis_name, transformer):
        if transformer is None:
            return

        if not isinstance(transformer, dict):
            raise TypeError(f"{axis_name}axis_transformer must be a dictionary or None.")

        transform = (
            transformer.get("function")
            or transformer.get("transform")
            or transformer.get("callable")
        )
        label_key = f"{axis_name}label"
        label = transformer.get("label", transformer.get(label_key, None))
        ticks = transformer.get("ticks", None)
        inverse = (
            transformer.get("inverse")
            or transformer.get("inverse_function")
            or transformer.get("inverse_transform")
        )

        if not callable(transform):
            raise TypeError(
                f"{axis_name}axis_transformer must include a callable under "
                "'function', 'transform', or 'callable'."
            )

        if ticks is not None:
            ticklabels = np.asarray(ticks)
            if transform_axes_data:
                tick_positions = ticklabels
            elif not callable(inverse):
                raise TypeError(
                    f"{axis_name}axis_transformer with 'ticks' must include a callable "
                    "'inverse', 'inverse_function', or 'inverse_transform'."
                )
            else:
                tick_positions = [inverse(x) for x in ticklabels]
        else:
            if axis_name == "x":
                tick_positions = ax.get_xticks()
            else:
                tick_positions = ax.get_yticks()
            ticklabels = tick_positions if transform_axes_data else [transform(x) for x in tick_positions]

        ticklabels = [f"{v:g}" for v in ticklabels]
        if axis_name == "x":
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(ticklabels)
        else:
            ax.set_yticks(tick_positions)
            ax.set_yticklabels(ticklabels)
        if label is not None:
            if axis_name == "x":
                ax.set_xlabel(label, fontsize=label_fontsize)
            else:
                ax.set_ylabel(label, fontsize=label_fontsize)

    
    cf = _draw_field()
    ax.set_xlabel(xlabel, fontsize=label_fontsize)
    ax.set_ylabel(ylabel, fontsize=label_fontsize)

    _add_bounding_shell()
    _add_streamtube()
    _apply_axes_limits_and_aspect()
    _add_farm()
    _add_blh()

    if stamp is not None:
        ax.text(
            stamp_location[0], stamp_location[1], stamp,
            transform=ax.transAxes,
            fontsize=stamp_fontsize,
            ha=stamp_ha,
            va=stamp_va,
        )

    if simulation_label is not None:
        ax.text(
            0.02, 1.02, simulation_label,
            transform=ax.transAxes,
            va="bottom", ha="left",
            fontsize=annotation_fontsize,
            clip_on=False,
        )

    _apply_axis_transformer("x", xaxis_transformer)
    _apply_axis_transformer("y", yaxis_transformer)
    ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)

    return cf

def _add_slice_colorbar(fig, axes, cf, scale, **kwargs):
    deficit = kwargs.get('deficit', True)
    fixed_colorbar_ticks = kwargs.get('fixed_colorbar_ticks', False)
    colorbar_tick_count = kwargs.get('colorbar_tick_count', 9)
    colorbar_orient = kwargs.get('colorbar_orient', 'vertical')
    colorbar_label = kwargs.get('colorbar_label', None)
    colorbar_pad = kwargs.get('colorbar_pad', 0.08)
    colorbar_fraction = kwargs.get('colorbar_fraction', 0.025)
    colorbar_size = kwargs.get('colorbar_size', kwargs.get('colorbar_thickness', '2%'))
    tick_fontsize = kwargs.get(
        'colorbar_tick_fontsize',
        _plot_fontsize(kwargs, 'tick_fontsize', 6),
    )
    label_fontsize = _plot_fontsize(
        kwargs,
        'colorbar_label_fontsize',
        6 if colorbar_orient == 'vertical' else 8,
    )

    if colorbar_orient == 'vertical':
        if isinstance(axes, (list, tuple, np.ndarray)):
            cbar = fig.colorbar(
                cf, ax=axes, orientation=colorbar_orient, pad=colorbar_pad, fraction=colorbar_fraction
            )
        else:
            cbar = _add_matched_vertical_colorbar(
                fig, cf, axes, pad=colorbar_pad, size=colorbar_size
            )
        if colorbar_label is not None:
            cbar.set_label(colorbar_label, fontsize=label_fontsize, labelpad=6)
    else:
        cax = fig.add_axes([0.27, 0.065, 0.46, 0.032])
        cbar = fig.colorbar(cf, cax=cax, orientation=colorbar_orient)
        if colorbar_label is not None:
            cbar.ax.text(
                1.07, 0.5, colorbar_label,
                transform=cbar.ax.transAxes,
                va="center", ha="left", fontsize=label_fontsize,
            )

    if deficit and fixed_colorbar_ticks:
        ticks = np.linspace(-scale, scale, colorbar_tick_count)
        cbar.set_ticks(ticks)

    cbar.formatter = FuncFormatter(_format_tick_max_two_decimals)
    cbar.update_ticks()
    cbar.ax.tick_params(labelsize=tick_fontsize)

    return cbar

def _add_local_slice_colorbars(fig, axes, cfs, color_specs, local_kwargs_list=None, **kwargs):
    if local_kwargs_list is None:
        local_kwargs_list = [kwargs] * len(cfs)

    for ax, cf, color_spec, local_kwargs in zip(axes, cfs, color_specs, local_kwargs_list):
        cbar_kwargs = kwargs.copy()
        cbar_kwargs.update(local_kwargs)
        deficit = cbar_kwargs.get('deficit', True)
        fixed_colorbar_ticks = cbar_kwargs.get('fixed_colorbar_ticks', False)
        colorbar_tick_count = cbar_kwargs.get('colorbar_tick_count', 9)
        colorbar_orient = cbar_kwargs.get('colorbar_orient', 'vertical')
        colorbar_label = cbar_kwargs.get('colorbar_label', None)
        colorbar_pad = cbar_kwargs.get('colorbar_pad', 0.03)
        colorbar_size = cbar_kwargs.get(
            'colorbar_size',
            cbar_kwargs.get('colorbar_thickness', '2%'),
        )
        tick_fontsize = cbar_kwargs.get(
            'colorbar_tick_fontsize',
            _plot_fontsize(cbar_kwargs, 'tick_fontsize', 6),
        )
        label_fontsize = _plot_fontsize(
            cbar_kwargs,
            'colorbar_label_fontsize',
            6 if colorbar_orient == 'vertical' else 8,
        )
        scale = color_spec[3]
        if colorbar_orient == 'vertical':
            cbar = _add_matched_vertical_colorbar(
                fig, cf, ax, pad=colorbar_pad, size=colorbar_size
            )
            if colorbar_label is not None:
                cbar.set_label(colorbar_label, fontsize=label_fontsize, labelpad=6)
        else:
            cbar = fig.colorbar(cf, ax=ax, orientation=colorbar_orient, pad=0.16, fraction=0.08)
            if colorbar_label is not None:
                cbar.ax.text(
                    1.03, 0.5, colorbar_label,
                    transform=cbar.ax.transAxes,
                    va="center", ha="left", fontsize=label_fontsize,
                )

        if deficit and fixed_colorbar_ticks and scale is not None:
            ticks = np.linspace(-scale, scale, colorbar_tick_count)
            cbar.set_ticks(ticks)

        cbar.formatter = FuncFormatter(_format_tick_max_two_decimals)
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=tick_fontsize)

def _adjust_slice_layout(fig, **kwargs):
    colorbar_orient = kwargs.get('colorbar_orient', 'vertical')
    hspace = kwargs.get('hspace', 0.18)
    wspace = kwargs.get('wspace', None)
    top = 0.91 if kwargs.get('_has_figure_title', False) else 0.96

    if colorbar_orient == 'vertical':
        fig.subplots_adjust(
            left=0.12, right=0.92, bottom=0.14, top=top,
            hspace=hspace, wspace=wspace,
        )
    else:
        fig.subplots_adjust(
            left=0.12, right=0.98, bottom=0.30, top=top,
            hspace=hspace, wspace=wspace,
        )

def _slice_figure_size(nrows=1, ncols=1, **kwargs):
    if "figsize" in kwargs and kwargs["figsize"] is not None:
        return kwargs["figsize"]

    colorbar_orient = kwargs.get('colorbar_orient', 'vertical')
    equal_aspect = kwargs.get('equal_aspect', True)

    # Single panels keep a compact publication-friendly default. Multi-panel
    # figures use slightly smaller panels so two-column layouts do not sprawl.
    if nrows == 1 and ncols == 1:
        panel_width = 5.2
        panel_height = 3.0 if equal_aspect else 2.8
    else:
        panel_width = 6.0 if ncols == 1 else 3.7
        panel_height = 0.72 if equal_aspect and ncols == 1 else (2.8 if equal_aspect else 2.45)

    width = panel_width * ncols
    height = panel_height * nrows

    if colorbar_orient == 'vertical':
        width += 0.55
    else:
        height += 0.45

    return (width, height)

def plot_slice(filename, **kwargs):
    _configure_latex_fonts()

    fig, ax = plt.subplots(figsize=_slice_figure_size(**kwargs))
    export = kwargs.get('export', None)

    d = _prepare_plot_slice_data(filename, **kwargs)
    cmap, norm, plot_levels, scale = _make_shared_norm_and_levels([d], **kwargs)
    cf = _plot_slice_on_axis(ax, d, cmap, norm, plot_levels, **kwargs)
    _add_slice_colorbar(fig, ax, cf, scale, **kwargs)
    _adjust_slice_layout(fig, **kwargs)
    
    if export is None:
        export = "slice.png"
    print(export)
    fig.savefig(export, dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)

def plot_slices(simulations, **kwargs):
    _configure_latex_fonts()

    if not simulations:
        raise ValueError("simulations must contain at least one simulation dictionary.")

    for i, sim in enumerate(simulations):
        if not isinstance(sim, dict):
            raise TypeError(f"Simulation {i} must be a dictionary.")
        if "filename" not in sim:
            raise KeyError(f"Simulation {i} is missing required key 'filename'.")

    export = kwargs.get('export', None)
    annotation = kwargs.get('annotation', 'letter_name')
    figure_stamp = kwargs.get('stamp', None)
    local_colorbars = kwargs.get('local_colorbars', False)
    stamp_fontsize = _plot_fontsize(kwargs, 'stamp_fontsize', 6)
    stamp_location = kwargs.get('stamp_location', (0.98, 0.985))
    stamp_ha = kwargs.get('stamp_ha', 'right')
    stamp_va = kwargs.get('stamp_va', 'top')
    ncols = kwargs.pop('ncols', 1)
    nrows = int(np.ceil(len(simulations) / ncols))
    panel_order = kwargs.get('panel_order', 'row')
    if panel_order in ('row', 'row_first', 'row-major', 'row_major'):
        panel_indices = list(range(nrows * ncols))
    elif panel_order in ('column', 'column_first', 'column-major', 'column_major'):
        panel_indices = [
            row * ncols + col
            for col in range(ncols)
            for row in range(nrows)
        ]
    else:
        raise ValueError("panel_order must be 'row' or 'column'.")
    panel_indices = panel_indices[:len(simulations)]
    figsize = _slice_figure_size(nrows=nrows, ncols=ncols, **kwargs)

    prepared = []
    for sim in simulations:
        panel_kwargs = _simulation_plot_kwargs(kwargs, sim)
        prepared.append(_prepare_plot_slice_data(sim["filename"], **panel_kwargs))

    use_local_colorbars = local_colorbars
    if use_local_colorbars:
        color_specs = [
            _make_shared_norm_and_levels(
                [d],
                **_simulation_plot_kwargs(kwargs, sim),
            )
            for sim, d in zip(simulations, prepared)
        ]
    else:
        color_specs = [
            _make_shared_norm_and_levels(prepared, **kwargs)
        ] * len(prepared)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes_flat = axes.ravel()
    panel_axes = [axes_flat[index] for index in panel_indices]
    cfs = []

    for i, (sim, d, ax, ax_index) in enumerate(zip(simulations, prepared, panel_axes, panel_indices)):
        panel_kwargs = _simulation_plot_kwargs(kwargs, sim)
        if figure_stamp is not None and "stamp" not in sim.get("plot_kwargs", {}):
            panel_kwargs["stamp"] = None
        label = _format_simulation_label(i, sim.get("name", None), annotation)
        cmap, norm, plot_levels, scale = color_specs[i]
        cf = _plot_slice_on_axis(ax, d, cmap, norm, plot_levels, simulation_label=label, **panel_kwargs)
        cfs.append(cf)
        if ax_index // ncols < nrows - 1:
            ax.set_xlabel(None)
            ax.tick_params(axis='x', labelbottom=False)

    # Remove unused axes when the simulation count is odd.
    used_indices = set(panel_indices)
    for index, ax in enumerate(axes_flat):
        if index in used_indices:
            continue
        fig.delaxes(ax)

    if not use_local_colorbars:
        _add_slice_colorbar(
            fig,
            panel_axes,
            cfs[-1],
            color_specs[-1][3],
            **kwargs,
        )
    if figure_stamp is not None:
        fig.text(
            stamp_location[0], stamp_location[1], figure_stamp,
            ha=stamp_ha, va=stamp_va,
            fontsize=stamp_fontsize,
        )
    layout_kwargs = kwargs.copy()
    layout_kwargs["_has_figure_title"] = figure_stamp is not None
    if use_local_colorbars and kwargs.get('colorbar_orient', 'vertical') == 'horizontal':
        layout_kwargs["hspace"] = max(layout_kwargs.get("hspace", 0.18), 0.62)
    _adjust_slice_layout(fig, **layout_kwargs)
    if use_local_colorbars:
        _add_local_slice_colorbars(
            fig,
            panel_axes,
            cfs,
            color_specs,
            local_kwargs_list=[
                _simulation_plot_kwargs(kwargs, sim)
                for sim in simulations
            ],
            **kwargs,
        )

    if export is None:
        export = "slices.png"
    print(export)
    fig.savefig(export, dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)

def _sanitize_spectrum_field_name(field):
    clean = str(field).strip()
    for char in (" ", "/", "\\", ":", ",", ";", "{", "}", "(", ")", "[", "]"):
        clean = clean.replace(char, "_")
    return clean

def spectrum_filename(directory, field, kind="horizontal"):
    kind_map = {
        "horizontal": "spectrum_{field}.csv",
        "isotropic": "spectrum_{field}.csv",
        "k": "spectrum_{field}.csv",
        "x": "spectrum_x_{field}.csv",
        "streamwise": "spectrum_x_{field}.csv",
        "kx": "spectrum_x_{field}.csv",
        "y": "spectrum_y_{field}.csv",
        "spanwise": "spectrum_y_{field}.csv",
        "ky": "spectrum_y_{field}.csv",
        "zsummary": "spectrum_zsummary_{field}.csv",
        "summary": "spectrum_zsummary_{field}.csv",
        "vertical_summary": "spectrum_zsummary_{field}.csv",
        "height": "spectrum_height_{field}.csv",
        "height_resolved": "spectrum_height_{field}.csv",
        "yzplane": "spectrum_yzplane_{field}.csv",
        "yz_plane": "spectrum_yzplane_{field}.csv",
        "x_ky": "spectrum_yzplane_{field}.csv",
    }
    normalized = str(kind).lower()
    if normalized not in kind_map:
        raise ValueError(f"Unknown spectrum kind '{kind}'.")
    clean_field = _sanitize_spectrum_field_name(field)
    return os.path.join(os.fspath(directory), kind_map[normalized].format(field=clean_field))

def _spectrum_kind_from_columns(columns):
    cols = tuple(columns)
    colset = set(cols)
    if colset == {"k", "E"}:
        return "horizontal"
    if colset == {"kx", "E"}:
        return "x"
    if colset == {"ky", "E"}:
        return "y"
    if colset == {"k", "E", "z_centroid", "z_spread"}:
        return "zsummary"
    if colset == {"z", "k", "E"}:
        return "height"
    if colset == {"x", "ky", "E"}:
        return "yzplane"
    raise ValueError(f"Unrecognized spectrum CSV columns: {cols}")

def read_spectrum_csv(filename, **kwargs):
    path = os.fspath(filename)
    fieldgain = kwargs.get("fieldgain", kwargs.get("scale", 1.0))
    kgain = kwargs.get("kgain", 1.0)
    lengthgain = kwargs.get("lengthgain", 1.0)

    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding=None)
    if data.dtype.names is None:
        raise ValueError(f"Spectrum CSV '{path}' must have a header row.")

    columns = tuple(name.strip() for name in data.dtype.names)
    kind = _spectrum_kind_from_columns(columns)

    def col(name):
        if name not in data.dtype.names:
            return None
        return np.atleast_1d(np.asarray(data[name], dtype=float))

    spectrum = {
        "source": "csv",
        "filename": path,
        "kind": kind,
        "columns": columns,
        "E": col("E") * fieldgain,
    }

    if kind in ("horizontal", "zsummary", "height"):
        spectrum["k"] = col("k") * kgain
        spectrum["coord_name"] = "k"
    if kind == "x":
        spectrum["kx"] = col("kx") * kgain
        spectrum["coord_name"] = "kx"
    if kind == "y":
        spectrum["ky"] = col("ky") * kgain
        spectrum["coord_name"] = "ky"
    if kind == "yzplane":
        spectrum["x"] = col("x") * lengthgain
        spectrum["ky"] = col("ky") * kgain
        spectrum["coord_name"] = ("x", "ky")
    if kind == "zsummary":
        spectrum["z_centroid"] = col("z_centroid") * lengthgain
        spectrum["z_spread"] = col("z_spread") * lengthgain
    if kind == "height":
        spectrum["z"] = col("z") * lengthgain
        spectrum["coord_name"] = ("z", "k")

    return spectrum

def _resolve_spectrum_filename(spec, **kwargs):
    if "filename" in spec:
        return spec["filename"]
    directory = spec.get("directory", spec.get("dir", kwargs.get("directory", None)))
    field = spec.get("field", kwargs.get("field", None))
    kind = spec.get("kind", kwargs.get("kind", "horizontal"))
    if directory is None or field is None:
        raise KeyError("Spectrum specification must include 'filename' or both 'directory' and 'field'.")
    return spectrum_filename(directory, field, kind=kind)

def _copy_spectrum_dict(data):
    copied = data.copy()
    for key in ("k", "kx", "ky", "x", "z", "E", "z_centroid", "z_spread"):
        if key in copied:
            copied[key] = np.array(copied[key], copy=True)
    return copied

def _reduce_yzplane_spectrum(data, **kwargs):
    if data.get("kind") != "yzplane":
        return data

    reduce_method = kwargs.get("yzplane_reduce", kwargs.get("x_reduce", "sum"))
    x_bounds = kwargs.get("x_bounds", kwargs.get("xbounds", None))
    if reduce_method in (None, False, "none", "raw"):
        return data

    x = np.asarray(data["x"], dtype=float)
    ky = np.asarray(data["ky"], dtype=float)
    energy = np.asarray(data["E"], dtype=float)
    mask = np.isfinite(x) & np.isfinite(ky) & np.isfinite(energy)
    if x_bounds is not None:
        mask &= (x >= x_bounds[0]) & (x <= x_bounds[1])
    if not np.any(mask):
        raise ValueError("No yz-plane spectrum rows remain after applying x_bounds.")

    x = x[mask]
    ky = ky[mask]
    energy = energy[mask]
    ky_values = np.unique(ky)
    reduced_energy = np.empty(len(ky_values), dtype=float)
    counts = np.empty(len(ky_values), dtype=int)

    for i, ky_value in enumerate(ky_values):
        ky_mask = ky == ky_value
        values = energy[ky_mask]
        counts[i] = np.count_nonzero(ky_mask)
        if reduce_method in ("sum", "integral"):
            reduced_energy[i] = np.nansum(values)
        elif reduce_method in ("mean", "average", "avg"):
            reduced_energy[i] = np.nanmean(values)
        elif reduce_method == "median":
            reduced_energy[i] = np.nanmedian(values)
        else:
            raise ValueError("yzplane_reduce must be 'sum', 'mean', 'median', or 'none'.")

    if kwargs.get("x_normalize", False) and reduce_method in ("sum", "integral"):
        reduced_energy = reduced_energy / np.maximum(counts, 1)

    return {
        "source": data.get("source", "csv"),
        "filename": data.get("filename", None),
        "kind": "y",
        "original_kind": "yzplane",
        "columns": ("ky", "E"),
        "ky": ky_values,
        "E": reduced_energy,
        "coord_name": "ky",
        "yzplane_reduce": reduce_method,
        "x_bounds": x_bounds,
        "x_count": counts,
    }

def _weighted_spectrum_stat(values, weights, stat):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(values)
    if not np.any(valid):
        return np.nan
    values = values[valid]
    weights = weights[valid]
    if stat in ("sum", "integral"):
        return np.nansum(values)
    if stat in ("median",):
        return np.nanmedian(values)
    if stat in ("weighted", "weighted_mean", "energy_weighted"):
        total_weight = np.nansum(weights)
        if total_weight == 0:
            return np.nanmean(values)
        return np.nansum(values * weights) / total_weight
    return np.nanmean(values)

def _logbin_spectrum(data, coord_key, bins=80, stat="mean"):
    coord = np.asarray(data[coord_key], dtype=float)
    positive = coord > 0
    if np.count_nonzero(positive) < 2:
        return data

    coord_pos = coord[positive]
    edges = np.logspace(np.log10(np.nanmin(coord_pos)), np.log10(np.nanmax(coord_pos)), int(bins) + 1)
    bin_index = np.digitize(coord, edges) - 1
    smooth = _copy_spectrum_dict(data)
    keep = []

    for ibin in range(int(bins)):
        mask = bin_index == ibin
        if not np.any(mask):
            continue
        keep.append(ibin)

    nbin = len(keep)
    if nbin == 0:
        return data

    smooth[coord_key] = np.empty(nbin, dtype=float)
    for key in ("E", "z_centroid", "z_spread"):
        if key in data and len(data[key]) == len(coord):
            smooth[key] = np.empty(nbin, dtype=float)

    for out_i, ibin in enumerate(keep):
        mask = bin_index == ibin
        weights = data["E"][mask] if "E" in data else np.ones(np.count_nonzero(mask))
        smooth[coord_key][out_i] = _weighted_spectrum_stat(coord[mask], weights, "weighted_mean")
        for key in ("E", "z_centroid", "z_spread"):
            if key in data and len(data[key]) == len(coord):
                key_stat = "weighted_mean" if key in ("z_centroid", "z_spread") and stat == "mean" else stat
                smooth[key][out_i] = _weighted_spectrum_stat(data[key][mask], weights, key_stat)

    smooth["smoothed"] = True
    smooth["smooth_method"] = "logbin"
    smooth["smooth_bins"] = bins
    smooth["smooth_stat"] = stat
    return smooth

def _moving_spectrum(data, coord_key, window=5, stat="mean"):
    window = int(window)
    if window <= 1:
        return data
    if window % 2 == 0:
        window += 1
    half_window = window // 2
    coord = np.asarray(data[coord_key], dtype=float)
    smooth = _copy_spectrum_dict(data)
    smooth_keys = [key for key in ("E", "z_centroid", "z_spread") if key in data and len(data[key]) == len(coord)]

    for i in range(len(coord)):
        left = max(0, i - half_window)
        right = min(len(coord), i + half_window + 1)
        weights = data["E"][left:right] if "E" in data else np.ones(right - left)
        for key in smooth_keys:
            key_stat = "weighted_mean" if key in ("z_centroid", "z_spread") and stat == "mean" else stat
            smooth[key][i] = _weighted_spectrum_stat(data[key][left:right], weights, key_stat)

    smooth["smoothed"] = True
    smooth["smooth_method"] = "moving"
    smooth["smooth_window"] = window
    smooth["smooth_stat"] = stat
    return smooth

def _smooth_height_spectrum(data, **kwargs):
    method = kwargs.get("smooth_method", "logbin")
    bins = kwargs.get("smooth_bins", 80)
    window = kwargs.get("smooth_window", kwargs.get("smooth_bins", 5))
    stat = kwargs.get("smooth_stat", "mean")
    z_values = np.unique(data["z"])
    pieces = []

    for z_value in z_values:
        mask = data["z"] == z_value
        row = {
            "kind": "horizontal",
            "k": data["k"][mask],
            "E": data["E"][mask],
        }
        if method == "logbin":
            row = _logbin_spectrum(row, "k", bins=bins, stat=stat)
        elif method in ("moving", "rolling"):
            row = _moving_spectrum(row, "k", window=window, stat=stat)
        else:
            raise ValueError("smooth_method must be 'logbin' or 'moving'.")
        pieces.append((z_value, row))

    smooth = _copy_spectrum_dict(data)
    smooth["z"] = np.concatenate([np.full(len(row["k"]), z_value) for z_value, row in pieces])
    smooth["k"] = np.concatenate([row["k"] for _, row in pieces])
    smooth["E"] = np.concatenate([row["E"] for _, row in pieces])
    smooth["smoothed"] = True
    smooth["smooth_method"] = method
    smooth["smooth_stat"] = stat
    return smooth

def _smooth_spectrum(data, **kwargs):
    if not kwargs.get("smooth", False):
        data["smoothed"] = False
        return data

    method = kwargs.get("smooth_method", "logbin")
    stat = kwargs.get("smooth_stat", "mean")
    if data.get("kind") == "height":
        if not kwargs.get("smooth_height", True):
            data["smoothed"] = False
            return data
        return _smooth_height_spectrum(data, **kwargs)

    coord_key = "k" if "k" in data else ("kx" if "kx" in data else ("ky" if "ky" in data else None))
    if coord_key is None:
        data["smoothed"] = False
        return data
    if method == "logbin":
        return _logbin_spectrum(data, coord_key, bins=kwargs.get("smooth_bins", 80), stat=stat)
    if method in ("moving", "rolling"):
        return _moving_spectrum(data, coord_key, window=kwargs.get("smooth_window", kwargs.get("smooth_bins", 5)), stat=stat)
    raise ValueError("smooth_method must be 'logbin' or 'moving'.")

def _shift_spectrum_wavenumber_position(data, **kwargs):
    position = kwargs.get("k_position", kwargs.get("wavenumber_position", "center"))
    if position in (None, "center", "centroid", "bin_center"):
        data["k_position"] = "center"
        return data
    if position not in ("lower", "left", "edge", "upper", "right"):
        raise ValueError("k_position must be 'center', 'lower', or 'upper'.")

    shifted = _copy_spectrum_dict(data)
    for coord_key in ("k", "kx", "ky"):
        if coord_key not in shifted:
            continue
        coord = np.asarray(shifted[coord_key], dtype=float)
        if coord.size < 2:
            continue
        dk = np.empty_like(coord)
        dk[1:] = np.diff(coord)
        dk[0] = dk[1]
        if position in ("lower", "left", "edge"):
            shifted[coord_key] = coord - 0.5 * dk
        else:
            shifted[coord_key] = coord + 0.5 * dk
        shifted[coord_key] = np.maximum(shifted[coord_key], 0.0)
    shifted["k_position"] = "lower" if position in ("lower", "left", "edge") else "upper"
    return shifted

def _filter_nonpositive_spectrum_wavenumbers(data):
    coord_key = "k" if "k" in data else ("kx" if "kx" in data else ("ky" if "ky" in data else None))
    if coord_key is None:
        return data
    coord = np.asarray(data[coord_key], dtype=float)
    finite = np.isfinite(coord)
    scale = np.nanmax(np.abs(coord[finite])) if np.any(finite) else 1.0
    tolerance = np.finfo(float).eps * max(scale, 1.0) * 100.0
    mask = coord > tolerance
    filtered = _copy_spectrum_dict(data)
    for key in ("k", "kx", "ky", "z", "E", "z_centroid", "z_spread"):
        if key in filtered and len(filtered[key]) == len(mask):
            filtered[key] = filtered[key][mask]
    return filtered

def _spectrum_coordinate_key(data):
    if "k" in data:
        return "k"
    if "kx" in data:
        return "kx"
    if "ky" in data:
        return "ky"
    return None

def _spectrum_bin_width(coord):
    coord = np.asarray(coord, dtype=float)
    finite = np.isfinite(coord)
    values = coord[finite]
    if values.size < 2:
        return 1.0
    diffs = np.diff(np.sort(values))
    diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
    if diffs.size == 0:
        return 1.0
    return np.nanmedian(diffs)

def _spectrum_parseval_energy(data, density=False):
    energy = np.asarray(data["E"], dtype=float)
    total = np.nansum(energy)
    if data.get("kind") == "height" and "z" in data:
        z_count = len(np.unique(data["z"]))
        if z_count > 0:
            total = total / z_count
    if data.get("original_kind") == "yzplane" and data.get("yzplane_reduce") in ("sum", "integral"):
        x_count = np.asarray(data.get("x_count", []), dtype=float)
        x_count = x_count[np.isfinite(x_count) & (x_count > 0)]
        if x_count.size > 0:
            total = total / np.nanmedian(x_count)

    if not density:
        return total

    coord_key = _spectrum_coordinate_key(data)
    if coord_key is None:
        return total
    return total * _spectrum_bin_width(data[coord_key])

def _spectrum_parseval_reference_energy(reference, **kwargs):
    if isinstance(reference, dict):
        data = _copy_spectrum_dict(reference)
    else:
        data = read_spectrum_csv(reference, **kwargs)
    data = _reduce_yzplane_spectrum(data, **kwargs)
    return _spectrum_parseval_energy(data, density=kwargs.get("density", kwargs.get("write_density", False)))

def _normalization_reference_key(reference):
    if isinstance(reference, dict):
        return id(reference)
    return os.fspath(reference)

def _prepare_plot_spectrum_data(filename, **kwargs):
    fieldgain = kwargs.get("fieldgain", kwargs.get("scale", 1.0))
    background = kwargs.get("background", None)
    percentage = kwargs.get("percentage", False)
    normalize = kwargs.get("normalize", None)
    premultiply = kwargs.get("premultiply", False)
    skip_zero = kwargs.get("skip_zero", kwargs.get("skip_k0", False))
    xscale = kwargs.get("xscale", "log")

    if isinstance(filename, dict):
        d = _copy_spectrum_dict(filename)
        if fieldgain != 1:
            d["E"] *= fieldgain
    else:
        d = read_spectrum_csv(filename, **kwargs)
    d = _reduce_yzplane_spectrum(d, **kwargs)

    if background is not None:
        bg = read_spectrum_csv(background, **kwargs) if not isinstance(background, dict) else _copy_spectrum_dict(background)
        bg = _reduce_yzplane_spectrum(bg, **kwargs)
        if d["kind"] != bg["kind"]:
            raise ValueError("Cannot compare spectrum background with a different kind.")
        coord_key = "k" if "k" in d else ("kx" if "kx" in d else "ky")
        if not np.array_equal(d[coord_key], bg[coord_key]):
            raise ValueError("Cannot subtract spectra with different coordinates.")
        if percentage:
            d["E"] = (d["E"] - bg["E"]) / bg["E"] * 100
        else:
            d["E"] = d["E"] - bg["E"]

    if normalize in ("parseval", "variance", "energy", "physical"):
        reference = kwargs.get("normalize_reference", kwargs.get("normalization_reference", None))
        if reference is None:
            normalize = _spectrum_parseval_energy(d, density=kwargs.get("density", kwargs.get("write_density", False)))
        else:
            normalize = _spectrum_parseval_reference_energy(reference, **kwargs)
            printed = kwargs.setdefault("_printed_normalization_references", set())
            reference_key = _normalization_reference_key(reference)
            if reference_key not in printed:
                print(f"Parseval normalization from {reference_key}: {normalize:.16e}")
                printed.add(reference_key)

    d = _shift_spectrum_wavenumber_position(d, **kwargs)

    coord_key = _spectrum_coordinate_key(d)
    if skip_zero or xscale == "log":
        d = _filter_nonpositive_spectrum_wavenumbers(d)

    d = _smooth_spectrum(d, **kwargs)
    coord_key = _spectrum_coordinate_key(d)

    if premultiply and coord_key is not None:
        d["E"] = d["E"] * d[coord_key]
        d["premultiplied"] = True
    else:
        d["premultiplied"] = False

    if normalize is not None and normalize is not False:
        if normalize in ("max", "maximum", True):
            denom = np.nanmax(np.abs(d["E"]))
        elif normalize in ("integral", "sum"):
            denom = np.nansum(d["E"])
        else:
            denom = float(normalize)
        if denom != 0:
            d["E"] = d["E"] / denom
        d["normalized"] = normalize
    else:
        d["normalized"] = None

    return d

def _spectrum_line_xy(data, quantity="E"):
    if quantity not in data:
        raise KeyError(f"Spectrum quantity '{quantity}' is not available for kind '{data['kind']}'.")
    if "k" in data:
        return data["k"], data[quantity], "k"
    if "kx" in data:
        return data["kx"], data[quantity], "kx"
    if "ky" in data:
        return data["ky"], data[quantity], "ky"
    raise ValueError("Height-resolved spectra are not line spectra.")

def _default_spectrum_xlabel(coord_key):
    labels = {
        "k": r"$k$",
        "kx": r"$k_x$",
        "ky": r"$k_y$",
    }
    return labels.get(coord_key, coord_key)

def _default_spectrum_ylabel(data, quantity="E", density=False):
    if quantity == "z_centroid":
        return r"$z_c$"
    if quantity == "z_spread":
        return r"$\sigma_z$"
    if data.get("premultiplied", False):
        return r"$kE(k)$"
    return r"$E(k)$" if not density else r"$E(k)\,/\,\Delta k$"

def _plot_spectrum_references(ax, references, **kwargs):
    if not references:
        return
    for ref in references:
        slope = ref["slope"]
        k0, e0 = ref["anchor"]
        xlim = ref.get("xlim", ax.get_xlim())
        x = np.asarray(ref.get("x", np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 80)))
        x = x[x > 0]
        y = e0 * (x / k0) ** slope
        ax.plot(
            x,
            y,
            color=ref.get("color", "0.35"),
            linestyle=ref.get("linestyle", ref.get("ls", "--")),
            linewidth=ref.get("linewidth", ref.get("lw", 0.8)),
            alpha=ref.get("alpha", 0.75),
            label=ref.get("label", None),
            zorder=ref.get("zorder", 1),
        )

def _collect_wavenumber_markers(**kwargs):
    markers = []
    marker_defs = kwargs.get("wavenumber_markers", None)
    if marker_defs:
        if isinstance(marker_defs, dict):
            marker_defs = [marker_defs]
        markers.extend(marker_defs)

    for length_key, label_default in (
        ("turbine_diameter", r"$2\pi/D$"),
        ("farm_length", r"$2\pi/L_f$"),
    ):
        length = kwargs.get(length_key, None)
        if length is None:
            continue
        markers.append({
            "length": length,
            "label": kwargs.get(f"{length_key}_label", label_default),
        })
    return markers

def _plot_wavenumber_markers(ax, **kwargs):
    markers = _collect_wavenumber_markers(**kwargs)
    if not markers:
        return

    default_color = kwargs.get("wavenumber_marker_color", "k")
    default_linestyle = kwargs.get("wavenumber_marker_linestyle", "--")
    default_alpha = kwargs.get("wavenumber_marker_alpha", 0.5)
    default_linewidth = kwargs.get("wavenumber_marker_linewidth", kwargs.get("wavenumber_marker_lw", 0.8))
    default_label = kwargs.get("wavenumber_marker_label", True)
    default_ymin = kwargs.get("wavenumber_marker_ymin", None)
    default_ymax = kwargs.get("wavenumber_marker_ymax", None)

    ymin, ymax = ax.get_ylim()
    for marker in markers:
        if np.isscalar(marker):
            marker = {"length": marker}
        if "k" in marker:
            k_value = float(marker["k"])
        else:
            length = float(marker["length"])
            if length <= 0:
                continue
            k_value = 2.0 * np.pi / length
        marker_ymin = marker.get("ymin", default_ymin)
        marker_ymax = marker.get("ymax", default_ymax)
        line_kwargs = {
            "color": marker.get("color", default_color),
            "linestyle": marker.get("linestyle", marker.get("ls", default_linestyle)),
            "alpha": marker.get("alpha", default_alpha),
            "linewidth": marker.get("linewidth", marker.get("lw", default_linewidth)),
            "zorder": marker.get("zorder", 0),
        }
        if marker_ymin is None and marker_ymax is None:
            ax.axvline(k_value, **line_kwargs)
        else:
            if marker_ymin is None:
                marker_ymin = ymin
            if marker_ymax is None:
                marker_ymax = ymax
            ax.vlines(k_value, marker_ymin, marker_ymax, **line_kwargs)
        label = marker.get("label", None)
        if label is not None and marker.get("show_label", default_label):
            label_y = marker.get("label_y", kwargs.get("wavenumber_marker_label_y", None))
            x_nudge = marker.get(
                "label_x_nudge",
                kwargs.get("wavenumber_marker_label_x_nudge", 3.0),
            )
            if label_y is None:
                label_y = 0.04
                label_transform = ax.get_xaxis_transform()
            else:
                label_transform = ax.transData
            label_transform = label_transform + ScaledTranslation(
                x_nudge / 72.0,
                0.0,
                ax.figure.dpi_scale_trans,
            )
            ax.text(
                k_value,
                label_y,
                label,
                rotation=90,
                ha=marker.get("ha", "left"),
                va=marker.get("va", "bottom"),
                fontsize=marker.get("fontsize", _plot_fontsize(kwargs, "wavenumber_marker_fontsize", 6)),
                color=marker.get("color", default_color),
                alpha=marker.get("text_alpha", marker.get("alpha", default_alpha)),
                transform=label_transform,
            )
    ax.set_ylim(ymin, ymax)

def _plot_spectrum_line_on_axis(ax, data, label=None, **kwargs):
    quantity = kwargs.get("quantity", "E")
    x, y, coord_key = _spectrum_line_xy(data, quantity=quantity)
    ax.plot(
        x,
        y,
        label=label,
        color=kwargs.get("color", None),
        linewidth=kwargs.get("linewidth", kwargs.get("lw", 1.2)),
        linestyle=kwargs.get("linestyle", kwargs.get("ls", "-")),
        marker=kwargs.get("marker", None),
        markersize=kwargs.get("markersize", kwargs.get("ms", None)),
        alpha=kwargs.get("alpha", None),
        zorder=kwargs.get("zorder", None),
    )
    return coord_key

def _height_spectrum_grid(data):
    z_values = np.unique(data["z"])
    k_values = np.unique(data["k"])
    grid = np.full((len(z_values), len(k_values)), np.nan)
    z_index = {value: i for i, value in enumerate(z_values)}
    k_index = {value: i for i, value in enumerate(k_values)}
    for z, k, energy in zip(data["z"], data["k"], data["E"]):
        grid[z_index[z], k_index[k]] = energy
    return k_values, z_values, grid

def _plot_height_spectrum_on_axis(ax, data, **kwargs):
    cmap = kwargs.get("cmap", "viridis")
    norm_type = kwargs.get("norm", kwargs.get("color_norm", "log"))
    vmin = kwargs.get("vmin", None)
    vmax = kwargs.get("vmax", None)
    k, z, energy = _height_spectrum_grid(data)
    if norm_type == "log":
        positive = energy[np.isfinite(energy) & (energy > 0)]
        if positive.size == 0:
            norm = None
        else:
            norm = LogNorm(
                vmin=np.nanmin(positive) if vmin is None else vmin,
                vmax=np.nanmax(positive) if vmax is None else vmax,
            )
    elif norm_type in (None, "linear"):
        norm = None
    else:
        norm = norm_type
    K, Z = np.meshgrid(k, z)
    return ax.pcolormesh(K, Z, energy, shading="auto", cmap=cmap, norm=norm, vmin=None if norm else vmin, vmax=None if norm else vmax)

def _apply_spectrum_axes(ax, coord_key=None, data=None, **kwargs):
    xscale = kwargs.get("xscale", "log")
    yscale = kwargs.get("yscale", "log")
    if data is not None and data.get("kind") == "height":
        yscale = kwargs.get("yscale", "linear")
    if xscale is not None:
        ax.set_xscale(xscale)
    if yscale is not None:
        ax.set_yscale(yscale)

    if kwargs.get("xlim", None) is not None:
        ax.set_xlim(kwargs["xlim"])
    if kwargs.get("ylim", None) is not None:
        ax.set_ylim(kwargs["ylim"])
    if kwargs.get("klim", None) is not None:
        ax.set_xlim(kwargs["klim"])
    if kwargs.get("zlim", None) is not None:
        ax.set_ylim(kwargs["zlim"])
    if kwargs.get("xleft", kwargs.get("xlft", None)) is not None:
        ax.set_xlim(left=kwargs.get("xleft", kwargs.get("xlft", None)))
    if kwargs.get("xright", kwargs.get("xrght", None)) is not None:
        ax.set_xlim(right=kwargs.get("xright", kwargs.get("xrght", None)))
    if kwargs.get("ybot", kwargs.get("ymin", None)) is not None:
        ax.set_ylim(bottom=kwargs.get("ybot", kwargs.get("ymin", None)))
    if kwargs.get("ytop", kwargs.get("ymax", None)) is not None:
        ax.set_ylim(top=kwargs.get("ytop", kwargs.get("ymax", None)))

    xlabel = kwargs.get("xlabel", None)
    ylabel = kwargs.get("ylabel", None)
    quantity = kwargs.get("quantity", "E")
    density = kwargs.get("density", kwargs.get("write_density", False))
    label_fontsize = _plot_fontsize(kwargs, "label_fontsize", 8)
    tick_fontsize = _plot_fontsize(kwargs, "tick_fontsize", 7)
    if xlabel is None and coord_key is not None:
        xlabel = _default_spectrum_xlabel(coord_key)
    if ylabel is None and data is not None:
        ylabel = r"$z$" if data.get("kind") == "height" else _default_spectrum_ylabel(data, quantity, density)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
    if kwargs.get("title", None) is not None:
        ax.set_title(kwargs["title"], fontsize=_plot_fontsize(kwargs, "title_fontsize", 9))
    if kwargs.get("grid", False):
        ax.grid(True, which=kwargs.get("grid_which", "both"), alpha=0.25, linewidth=0.5)
    ax.tick_params(axis="both", which="major", labelsize=tick_fontsize)

def plot_spectrum(filename, **kwargs):
    _configure_latex_fonts()
    export = kwargs.get("export", None)
    dpi = kwargs.get("dpi", 300)
    close = kwargs.get("close", True)
    fig, ax = plt.subplots(figsize=kwargs.get("figsize", (5.2, 3.4)))
    data = _prepare_plot_spectrum_data(filename, **kwargs)

    if data["kind"] == "height":
        cf = _plot_height_spectrum_on_axis(ax, data, **kwargs)
        _apply_spectrum_axes(ax, coord_key="k", data=data, **kwargs)
        _plot_wavenumber_markers(ax, **kwargs)
        cbar = _add_matched_vertical_colorbar(
            fig,
            cf,
            ax,
            pad=kwargs.get("colorbar_pad", 0.08),
            size=kwargs.get("colorbar_size", kwargs.get("colorbar_thickness", "2%")),
        )
        cbar.set_label(
            kwargs.get("colorbar_label", _default_spectrum_ylabel(data, "E", kwargs.get("density", False))),
            fontsize=_plot_fontsize(kwargs, "colorbar_label_fontsize", 7),
        )
        cbar.ax.tick_params(labelsize=kwargs.get("colorbar_tick_fontsize", _plot_fontsize(kwargs, "tick_fontsize", 7)))
    else:
        label = kwargs.get("label", None)
        coord_key = _plot_spectrum_line_on_axis(ax, data, label=label, **kwargs)
        _apply_spectrum_axes(ax, coord_key=coord_key, data=data, **kwargs)
        _plot_wavenumber_markers(ax, **kwargs)
        _plot_spectrum_references(ax, kwargs.get("references", None))
        if kwargs.get("legend", label is not None or bool(kwargs.get("references", None))):
            ax.legend(
                frameon=kwargs.get("legend_frameon", False),
                fontsize=_plot_fontsize(kwargs, "legend_fontsize", 7),
                loc=kwargs.get("legend_loc", None),
                ncols=kwargs.get("legend_ncols", kwargs.get("ncols", 1)),
            )

    if kwargs.get("tight_layout", True):
        plt.tight_layout()
    if export is None:
        export = "spectrum.png"
    print(export)
    fig.savefig(
        export,
        dpi=dpi,
        bbox_inches=kwargs.get("bbox_inches", "tight"),
        pad_inches=kwargs.get("pad_inches", 0.02),
    )
    if close:
        plt.close(fig)
    return data

def plot_spectra(simulations, **kwargs):
    _configure_latex_fonts()
    if not simulations:
        raise ValueError("simulations must contain at least one simulation dictionary.")
    kwargs.setdefault("_printed_normalization_references", set())

    prepared = []
    for i, sim in enumerate(simulations):
        if not isinstance(sim, dict):
            raise TypeError(f"Simulation {i} must be a dictionary.")
        panel_kwargs = kwargs.copy()
        panel_kwargs.update(sim.get("plot_kwargs", {}))
        filename = _resolve_spectrum_filename(sim, **kwargs)
        prepared.append(_prepare_plot_spectrum_data(filename, **panel_kwargs))

    height_mode = any(d["kind"] == "height" for d in prepared)
    panel = kwargs.get("panel", height_mode)
    export = kwargs.get("export", None)
    dpi = kwargs.get("dpi", 300)
    close = kwargs.get("close", True)

    if panel:
        ncols = kwargs.get("ncols", 1)
        nrows = int(np.ceil(len(simulations) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=kwargs.get("figsize", (5.2*ncols, 3.2*nrows)), squeeze=False)
        axes_flat = axes.ravel()
        cfs = []
        for i, (sim, data, ax) in enumerate(zip(simulations, prepared, axes_flat)):
            panel_kwargs = kwargs.copy()
            panel_kwargs.update(sim.get("plot_kwargs", {}))
            label = _format_simulation_label(i, sim.get("name", None), kwargs.get("annotation", "letter_name"))
            if data["kind"] == "height":
                cf = _plot_height_spectrum_on_axis(ax, data, **panel_kwargs)
                cfs.append((ax, cf, panel_kwargs, data))
                _apply_spectrum_axes(ax, coord_key="k", data=data, **panel_kwargs)
                _plot_wavenumber_markers(ax, **panel_kwargs)
            else:
                coord_key = _plot_spectrum_line_on_axis(ax, data, label=sim.get("name", None), **panel_kwargs)
                _apply_spectrum_axes(ax, coord_key=coord_key, data=data, **panel_kwargs)
                _plot_wavenumber_markers(ax, **panel_kwargs)
                _plot_spectrum_references(ax, panel_kwargs.get("references", None))
                if panel_kwargs.get("legend", False):
                    ax.legend(frameon=panel_kwargs.get("legend_frameon", False), fontsize=_plot_fontsize(panel_kwargs, "legend_fontsize", 7))
            if label is not None:
                ax.text(0.02, 1.02, label, transform=ax.transAxes, va="bottom", ha="left", fontsize=_plot_fontsize(panel_kwargs, "annotation_fontsize", 7), clip_on=False)

        for ax in axes_flat[len(simulations):]:
            fig.delaxes(ax)
        if cfs and kwargs.get("colorbar", True):
            for ax, cf, panel_kwargs, data in cfs:
                cbar = _add_matched_vertical_colorbar(
                    fig,
                    cf,
                    ax,
                    pad=panel_kwargs.get("colorbar_pad", 0.08),
                    size=panel_kwargs.get("colorbar_size", panel_kwargs.get("colorbar_thickness", "2%")),
                )
                cbar.set_label(
                    panel_kwargs.get("colorbar_label", _default_spectrum_ylabel(data, "E", panel_kwargs.get("density", False))),
                    fontsize=_plot_fontsize(panel_kwargs, "colorbar_label_fontsize", 7),
                )
                cbar.ax.tick_params(labelsize=panel_kwargs.get("colorbar_tick_fontsize", _plot_fontsize(panel_kwargs, "tick_fontsize", 7)))
    else:
        fig, ax = plt.subplots(figsize=kwargs.get("figsize", (5.2, 3.4)))
        coord_key = None
        for sim, data in zip(simulations, prepared):
            panel_kwargs = kwargs.copy()
            panel_kwargs.update(sim.get("plot_kwargs", {}))
            coord_key = _plot_spectrum_line_on_axis(ax, data, label=sim.get("name", None), **panel_kwargs)
        _apply_spectrum_axes(ax, coord_key=coord_key, data=prepared[0], **kwargs)
        _plot_wavenumber_markers(ax, **kwargs)
        _plot_spectrum_references(ax, kwargs.get("references", None))
        if kwargs.get("legend", True):
            legend_options = {
                "frameon": kwargs.get("legend_frameon", False),
                "fontsize": _plot_fontsize(kwargs, "legend_fontsize", 7),
                "ncols": kwargs.get("legend_ncols", kwargs.get("ncols", 1)),
            }
            if kwargs.get("legend_loc", None) is not None:
                legend_options["loc"] = kwargs["legend_loc"]
            legend_options.update(kwargs.get("legend_kwargs", {}))
            ax.legend(**legend_options)

    if kwargs.get("stamp", None) is not None:
        fig.text(
            kwargs.get("stamp_location", (0.98, 0.98))[0],
            kwargs.get("stamp_location", (0.98, 0.98))[1],
            kwargs["stamp"],
            ha=kwargs.get("stamp_ha", "right"),
            va=kwargs.get("stamp_va", "top"),
            fontsize=_plot_fontsize(kwargs, "stamp_fontsize", 9),
        )
    if kwargs.get("tight_layout", True):
        plt.tight_layout()
    if export is None:
        export = "spectra.png"
    print(export)
    fig.savefig(
        export,
        dpi=dpi,
        bbox_inches=kwargs.get("bbox_inches", "tight"),
        pad_inches=kwargs.get("pad_inches", 0.02),
    )
    if close:
        plt.close(fig)
    return prepared
