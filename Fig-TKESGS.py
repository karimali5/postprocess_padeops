from utilities import plot_slices, transect_netcdf_slice, read_netcdf_slice
import os
import numpy as np

S800 = {"path": "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH800/slices_t092440_n014462", "timestamp": "092440", "nstamp": "014462"}
S250 = {"path": "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH250/slices_t092897_n014397", "timestamp": "092897", "nstamp": "014397"}

S800_2K = {"path": "/anvil/scratch/x-kali/PadeOpsSims/INV800-2K/slices", "timestamp": "101643", "nstamp": "019053"}
S800_5K = {"path": "/anvil/scratch/x-kali/PadeOpsSims/INV800-5K/slices", "timestamp": "101691", "nstamp": "019053"}

stations = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540,560,580,600,620,640,660]
xfarm = 160
yfarm = (95.61+63.11)/2
zfarm = (95 -126/2)/126
Lfarm = 40
Wfarm = 95.61-63.11
Zfarm = 1
hubheight = 95/126

# S800_h0 = np.mean(read_netcdf_slice(os.path.join(S800['path'], f'Run09_t{S800["timestamp"]}_INVH0.nc'), scale=1/126)["data"])
# S250_h0 = np.mean(read_netcdf_slice(os.path.join(S250['path'], f'Run09_t{S250["timestamp"]}_INVH0.nc'), scale=1/126)["data"])

S800_2K_h0 = np.mean(read_netcdf_slice(os.path.join(S800_2K['path'], f'Run09_t{S800_2K["timestamp"]}_INVH0.nc'), scale=1/126)["data"])
S800_5K_h0 = np.mean(read_netcdf_slice(os.path.join(S800_5K['path'], f'Run09_t{S800_5K["timestamp"]}_INVH0.nc'), scale=1/126)["data"])

def blh_files(sim, x):
    return [
        {'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_INVH0.nc'), 'x', x, plot=False, scale=1/126, smooth=True, l_smooth=5), 'color':"k", "style":"--"},
        {'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_INVH2.nc'), 'x', x, plot=False, scale=1/126, smooth=True, l_smooth=5), 'color':"k", "style":"--"},
        #{'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_SL_BLH.nc'), 'x', x, plot=False, smooth=True, l_smooth=2), 'color':"tab:red", "style":"-"},
    ]

def blh_files_y(sim, y):
    return [
        {'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_INVH0.nc'), 'y', y, plot=False, scale=1/126, smooth=True, l_smooth=5), 'color':"k", "style":"--"},
        {'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_INVH2.nc'), 'y', y, plot=False, scale=1/126, smooth=True, l_smooth=5), 'color':"k", "style":"--"},
        #{'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_SL_BLH.nc'), 'y', y, plot=False, smooth=True, l_smooth=2), 'color':"tab:red", "style":"-"},
    ]

plot_style = {
    "figsize": (7, 3),
    "panel_order": "row",
    "annotation": True,
    "ncols":1,
    "wspace": 0.08,
    "hspace": 0.34,
    "hide_top_right_spines": True,
    "local_colorbars": True,
}

field_style = {
    "fieldgain": 1,
    "lengthgain": 1.0,
    "deficit": True,
    # "s": 14,
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
    "equal_aspect": True,
    "ytop": 4*hubheight,
    "ybot": None,
    "xrght": 0.47*Wfarm+yfarm,
    "xlft": -0.35*Wfarm+yfarm,
    "xlabel": r"$y$",
    "ylabel": r"$z$",
    "aspect": None,
    "xaxis_transformer": {
        "function": lambda x: (x - yfarm) / Wfarm,
        "inverse": lambda xt: xt * Wfarm + yfarm,
        "ticks": [-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4],
        "label": r"$(y-y_0)/W_p$",
    },
    # "xaxis_transformer": {
    #         "function": lambda x: (x - xfarm) / Lfarm,
    #         "inverse": lambda xt: xt * Lfarm + xfarm,
    #         "ticks": np.arange(-3.75, 10.51, 0.5),
    #         "label": r"$(x-x_0)/L_f$",
    #     },
    "yaxis_transformer": {
        "function": lambda z: z / hubheight,
        "inverse": lambda zt: zt * hubheight,
        "ticks": np.arange(0, 4.1, 2),
        "label": r"$z/z_h$",
    },
    "transform_axes_data": False,
}

overlay_style = {
    "plotfarm": False,
    #"farmsize": [yfarm-Wfarm/2, zfarm, Wfarm, Zfarm],
    "farmsize": [xfarm, zfarm, Lfarm, Zfarm],
    "plot_blh": True,
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
    "colorbar_label_fontsize": 6,
    "fixed_colorbar_ticks": True,
    "colorbar_tick_count": 3,
    "colorbar_tick_fontsize": 6,
    "colorbar_label": r"$\times 10^{-4}$",
    "colorbar_label_location": "bottom",
    "colorbar_size": "1.6%",
    "colorbar_pad": 0.05,
    "colorbar_label_x": 3,
    "colorbar_label_y": 0.9,
    # "colorbar_groups": [
    #     {"indices": range(0, 1), "location": "right", "s": 14, "label": r"$\times 10^{-3}$"},
    # ],
}

font_style = {
    "tick_fontsize": 7,
    "label_fontsize": 7,
    "annotation_fontsize": 7,
    "stamp_fontsize": 9,
    "stamp_location": (0.96, 0.94),
    "stamp_ha": "right",
    "stamp_va": "top",
}

#for x in stations:
plot_slices(
        [
            {
                "filename": [
                    os.path.join(S800_5K["path"], f"Run09_comp_deficit_budget3_term04_t"+S800_5K["timestamp"]+"_n"+S800_5K["nstamp"]+f"_SL_x={165}.nc"),
                ],
                "name": r"$\mathcal{D}_{b \Delta}$",
                #"blh_file": blh_files_y(S800_5K, yfarm),
                'fieldgain': 10000,
                's': 1,
            },

            {
                "filename": [
                    os.path.join(S800_5K["path"], f"Run09_comp_deficit_budget3_term05_t"+S800_5K["timestamp"]+"_n"+S800_5K["nstamp"]+f"_SL_x={165}.nc"),
                ],
                "name": r"$\mathcal{D}_{\Delta b}$",
                #"blh_file": blh_files_y(S800_5K, yfarm),
                'fieldgain': 10000,
                's': 1,
            },

            {
                "filename": [
                    os.path.join(S800_5K["path"], f"Run09_comp_deficit_budget3_term06_t"+S800_5K["timestamp"]+"_n"+S800_5K["nstamp"]+f"_SL_x={165}.nc"),
                ],
                "name": r"$\mathcal{D}_{\Delta \Delta}$",
                #"blh_file": blh_files_y(S800_5K, yfarm),
                'fieldgain': 10000,
                's': 1,
            },

             {
                "filename": [
                    os.path.join(S800_5K["path"], f"Run09_comp_deficit_budget3_term04_t"+S800_5K["timestamp"]+"_n"+S800_5K["nstamp"]+f"_SL_x={165}.nc"),
                    os.path.join(S800_5K["path"], f"Run09_comp_deficit_budget3_term05_t"+S800_5K["timestamp"]+"_n"+S800_5K["nstamp"]+f"_SL_x={165}.nc"),
                    os.path.join(S800_5K["path"], f"Run09_comp_deficit_budget3_term06_t"+S800_5K["timestamp"]+"_n"+S800_5K["nstamp"]+f"_SL_x={165}.nc"),
                ],
                "name": r"$\Delta \mathcal{D}$",
                #"blh_file": blh_files_y(S800_5K, yfarm),
                'fieldgain': 10000,
                's': 1,
            },
    
        ],
        #stamp=rf"$x-x_0={(x-xfarm)/Lfarm:.2f}L_f$",
        **plot_style,
        **field_style,
        **axis_style,
        **overlay_style,
        **colorbar_style,
        **font_style,
        export=f"TKESGS_x=165.png",
    )
