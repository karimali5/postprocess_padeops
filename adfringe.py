from utilities import plot_slices, transect_netcdf_slice, read_netcdf_slice
import os
import numpy as np

plot_style = {
    "figsize": (4, 4),
    "panel_order": "row",
    "annotation": "letter_name",
}

field_style = {
    "fieldgain": 1,
    "lengthgain": 1.0,
    "deficit": True,
    "s": 1.1,
    "cmap_reg": "viridis",
    "cmap_def": "PRGn",
    "vmin": None,
    "vmax": None,
    "levels": 200,
    "contours": True,
    "smooth": False,
    "l_smooth": 5,
    "percentage": False,
}

axis_style = {
    "equal_aspect": False,
    "ytop": None,
    "ybot": None,
    "xlft": None,
    "xrght": None,
    "xlabel": r"$y$",
    "ylabel": r"$z$",
    "aspect": None,
    # "xaxis_transformer": {
    #     "function": lambda x: (x - yfarm) / Wfarm,
    #     "inverse": lambda xt: xt * Wfarm + yfarm,
    #     "ticks": np.arange(-2.25, 1.01, 0.25),
    #     "label": r"$(y-y_0)/W_f$",
    # },
    # "yaxis_transformer": {
    #     "function": lambda z: z / hubheight,
    #     "inverse": lambda zt: zt * hubheight,
    #     "ticks": np.arange(0, 16, 5),
    #     "label": r"$z/z_h$",
    # },
    "transform_axes_data": False,
}

overlay_style = {
    "plotfarm": False,
    "farmsize": None,
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
    "colorbar_orient": "horizontal",
    "colorbar_label": r"$\Delta w$",
    "fixed_colorbar_ticks": True,
}

font_style = {
    "tick_fontsize": 7,
    "label_fontsize": 8,
    "colorbar_tick_fontsize": 7,
    "colorbar_label_fontsize": 8,
    "annotation_fontsize": 8,
    "stamp_fontsize": 9,
    "stamp_location": (0.96, 0.94),
    "stamp_ha": "right",
    "stamp_va": "top",
}

def decimal_to_text(x):
    s = str(x)
    return s.replace(".", "p")

maindir="/anvil/scratch/x-kali/PadeOpsSims/test_fringeAD"
stations = [2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5]

for x in stations:
    plot_slices(
        [
            {
                "filename": [
                    os.path.join(maindir, "no_fringeAD", "Run06_comp_deficit_budget0_term01_t084852_n014769_SL_x="+decimal_to_text(x)+".nc"),
                ],
                "name": r"No AD Fringe",
            },

            {
                "filename": [
                    os.path.join(maindir, "orig", "Run06_comp_deficit_budget0_term01_t088000_n014278_SL_x="+decimal_to_text(x)+".nc"),
                ],
                "name": r"AD Fringe + No y taper",
            },

            {
                "filename": [
                    os.path.join(maindir, "smoothed_in_y", "Run06_comp_deficit_budget0_term01_t094000_n020277_SL_x="+decimal_to_text(x)+".nc"),
                ],
                "name": r"AD Fringe + y taper",
            },
        ],
        stamp=rf"$x={x:.2f}$",
        **plot_style,
        **field_style,
        **axis_style,
        **overlay_style,
        **colorbar_style,
        **font_style,
        export=f"delta_u_x={x}.png",
    )
