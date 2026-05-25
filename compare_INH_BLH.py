from utilities import inversion_height, plot_slices, read_netcdf_slice
import numpy as np
import os

xfarm = 120
yfarm = (95.61+63.11)/2
zfarm = (95-126/2)/126
Lfarm = 40
Wfarm = 95.61-63.11
Zfarm = 1
hubheight = 95/126

S800 = {"path": "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH800/slices_t092440_n014462", "timestamp": "092440", "nstamp": "014462"}
S250 = {"path": "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH250/slices_t092897_n014397", "timestamp": "092897", "nstamp": "014397"}

h1_prcntg_S800, h1l = inversion_height(S800, RID=9, BRID=8, mode='h1', deficit=True, percentage=True, smooth=True, l_smooth=5)
dh_prcntg_S800, dhl = inversion_height(S800, RID=9, BRID=8, mode='dh', deficit=True, percentage=True, smooth=True, l_smooth=5)
h1_prcntg_S250, _ = inversion_height(S250, RID=9, BRID=8, mode='h1', deficit=True, percentage=True, smooth=True, l_smooth=5)
dh_prcntg_S250, _ = inversion_height(S250, RID=9, BRID=8, mode='dh', deficit=True, percentage=True, smooth=True, l_smooth=5)

def blh_fractional_change(sim, run=9, base_run=8, scale=1/126, smooth=True, l_smooth=2):
    smooth_kwargs = {"smooth": smooth, "l_smooth": l_smooth}
    run_file = os.path.join(sim["path"], f'Run{run:02d}_t{sim["timestamp"]}_SL_BLH.nc')
    base_file = os.path.join(sim["path"], f'Run{base_run:02d}_t{sim["timestamp"]}_SL_BLH.nc')

    delta = read_netcdf_slice(
        [run_file, base_file],
        scale=[scale, -scale],
        **smooth_kwargs,
    )
    base = read_netcdf_slice(base_file, scale=scale, **smooth_kwargs)
    delta["data"] = delta["data"] / base["data"] * 100
    return delta

BLH_S800 = blh_fractional_change(S800)
BLH_S250 = blh_fractional_change(S250)

plot_style = {
    "ncols": 2,
    "panel_order": "column",
    "figsize": (10, 5),
    "wspace": 0.3,
    "annotation": "letter_name",
}

field_style = {
    "fieldgain": 1,
    "lengthgain": 1.0,
    "deficit": True,
    "local_colorbars": True,
    "cmap_reg": "viridis",
    "cmap_def": "PRGn",
    "vmin": None,
    "vmax": None,
    "levels": 200,
    "colorbar_range": "data",
    "contours": True,
    "smooth": False,
    "l_smooth": 5,
    "percentage": False,
}

axis_style = {
    "equal_aspect": True,
    "ytop": 140,
    "xlft": 0,
    "xrght": 540,
    "aspect": None,
    "xaxis_transformer": {
        "function": lambda x: (x - xfarm) / Lfarm,
        "inverse": lambda xt: xt * Lfarm + xfarm,
        "ticks": np.arange(-2, 10.1, 1),
        "label": r"$(x-x_0)/L_f$",
    },
    "yaxis_transformer": {
        "function": lambda x: (x - yfarm) / Wfarm,
        "inverse": lambda xt: xt * Wfarm + yfarm,
        "ticks": np.arange(-2, 2.1, 1),
        "label": r"$(y-y_0)/W_f$",
    },
    "transform_axes_data": False,
}

overlay_style = {
    "plotfarm": True,
    "farmsize": [xfarm, yfarm-Wfarm/2, Lfarm, Wfarm],
    "plot_blh": False,
    "background": None,
    "bounding_shell": False,
    "bounding_shell_threshold": -0.02,
    "bounding_shell_local": True,
    "plot_streamtube": False,
    "streamtube": None,
    "station": None,
}

colorbar_style = {
    "colorbar_orient": "vertical",
    "fixed_colorbar_ticks": False,
}

font_style = {
    "tick_fontsize": 7,
    "label_fontsize": 7,
    "colorbar_tick_fontsize": 6,
    "colorbar_label_fontsize": 7,
    "annotation_fontsize": 8,
    "stamp_fontsize": 9,
    "stamp_location": (0.96, 0.94),
    "stamp_ha": "right",
    "stamp_va": "top",
}

plot_slices(
    [
        {
            "filename": h1_prcntg_S800,
            "name": r"$H_{\textrm{inv}}=800\,m$",
            "plot_kwargs":{
                "colorbar_label": h1l
            }
        },
        {
            "filename": dh_prcntg_S800,
            "name": r"$H_{\textrm{inv}}=800\,m$",
            "plot_kwargs":{
                "colorbar_label": dhl
            }
        },
        {
            "filename": BLH_S800,
            "name": r"$H_{\textrm{inv}}=800\,m$",
            "plot_kwargs":{
                "colorbar_label": r"$\Delta H/H_0\,(\%)$"
            }
        },
        
        {
            "filename": h1_prcntg_S250,
            "name": r"$H_{\textrm{inv}}=250\,m$",
            "plot_kwargs":{
                "colorbar_label": h1l
            }
        },
        {
            "filename": dh_prcntg_S250,
            "name": r"$H_{\textrm{inv}}=250\,m$",
            "plot_kwargs":{
                "colorbar_label": dhl
            }
        },
        {
            "filename": BLH_S250,
            "name": r"$H_{\textrm{inv}}=250\,m$",
            "plot_kwargs":{
                "colorbar_label": r"$\Delta H/H_0\,(\%)$"
            }
        },
    
    ],
    **plot_style,
    **field_style,
    **axis_style,
    **overlay_style,
    **colorbar_style,
    **font_style,
    export=f"INVH.png",
)
