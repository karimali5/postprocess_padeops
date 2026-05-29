from utilities import plot_slices, transect_netcdf_slice, read_netcdf_slice
import os
import numpy as np

S800 = {"path": "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH800/slices_t092440_n014462", "timestamp": "092440", "nstamp": "014462"}
S250 = {"path": "/anvil/scratch/x-kali/PadeOpsSims/EXT-BLH250/slices_t092897_n014397", "timestamp": "092897", "nstamp": "014397"}
stations = [10,20,30,40,50,60,70,80,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540]
xfarm = 120
yfarm = (95.61+63.11)/2
zfarm = (95-126/2)/126
Lfarm = 40
Wfarm = 95.61-63.11
Zfarm = 1
hubheight = 95/126

S800_h0 = np.mean(read_netcdf_slice(os.path.join(S800['path'], f'Run09_t{S800["timestamp"]}_INVH0.nc'), scale=1/126)["data"])
S250_h0 = np.mean(read_netcdf_slice(os.path.join(S250['path'], f'Run09_t{S250["timestamp"]}_INVH0.nc'), scale=1/126)["data"])
sim = S800

def blh_files(sim, cutaxis):
    return [
        {'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_INVH0.nc'), cutaxis['axis'], cutaxis['value'], plot=False, scale=1/126, smooth=True, l_smooth=5), 'color':"k", "style":"--"},
        {'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_INVH2.nc'), cutaxis['axis'], cutaxis['value'], plot=False, scale=1/126, smooth=True, l_smooth=5), 'color':"k", "style":"--"},
        {'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_SL_BLH.nc'), cutaxis['axis'], cutaxis['value'], plot=False, smooth=True, l_smooth=2), 'color':"tab:red", "style":"-"},
    ]

plot_style = {
    "figsize": (8, 6),
    "panel_order": "row",
    "annotation": "letter_name",
    "hspace": 0.35,
}

field_style = {
    "lengthgain": 1.0,
    "deficit": True,
    "s": None,
    "local_colorbars": True,
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
    "ytop": 25*hubheight,
    "ybot": None,
    # "xlft": xfarm-2.5*Lfarm,
    # "xrght": xfarm+10*Lfarm,
    #"xlft": -2.2,
    #"xrght": 1.2,
    "xlabel": r"$y$",
    "ylabel": r"$z$",
    "aspect": None,
    "xaxis_transformer": {
        "function": lambda x: (x - yfarm) / Wfarm,
        "inverse": lambda xt: xt * Wfarm + yfarm,
        "ticks": np.arange(-2.25, 2.26, 0.25),
        "label": r"$(y-y_0)/W_f$",
    },
    # "xaxis_transformer": {
    #     "function": lambda x: (x - xfarm) / Lfarm,
    #     "inverse": lambda xt: xt * Lfarm + xfarm,
    #     "ticks": np.arange(-2.5, 10.01, 0.5),
    #     "label": r"$(x-x_0)/L_f$",
    # },
    "yaxis_transformer": {
        "function": lambda z: z / hubheight,
        "inverse": lambda zt: zt * hubheight,
        "ticks": np.arange(0, 26, 5),
        "label": r"$z/z_h$",
    },
    "transform_axes_data": False,
}

overlay_style = {
    "plotfarm": True,
    "farmsize": [yfarm-Wfarm/2, zfarm, Wfarm, Zfarm],
    #"farmsize": [xfarm, zfarm, Lfarm, Zfarm],
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
    "colorbar_label": None,
    "fixed_colorbar_ticks": True,
    "colorbar_size": "1.0%",
    "colorbar_tick_count": 5,
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

for x in stations:
    station_blh_files = blh_files(sim, {"axis":'x', 'value':x})
    plot_slices(
        [
            {
                "filename": [
                    os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term01_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
                ],
                "name": r"$\Delta u$",
                "blh_file": station_blh_files,
                "fieldgain": 1,
                "s":0.14,
            },

            {
                "filename": [
                    os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term03_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
                ],
                "name": r"$\Delta w$",
                "blh_file": station_blh_files,
                "fieldgain": 1,
                "s":0.05,
            },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term18_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term21_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term24_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\partial_z \tau_{xz}$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     "s":13,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term04_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta p\,\left(\times 10^{-3}\right)$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":20,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term18_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$-\partial_x \Delta p\,\left(\times 10^{-3}\right)$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     "s":14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term19_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$-\partial_y \Delta p\,\left(\times 10^{-3}\right)$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     "s":14,
            # },

            {
                "filename": [
                    os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term20_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
                ],
                "name": r"$-\partial_z \Delta p\,\left(\times 10^{-3}\right)$",
                "blh_file": station_blh_files,
                "fieldgain": -1000,
                "s":14,
            },

            {
                "filename": [
                    os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term17_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
                ],
                "name": r"$\Delta w_b\,\left(\times 10^{-3}\right)$",
                "blh_file": station_blh_files,
                "fieldgain": 1000,
                "s":6,
            },
            
        ],
        stamp=rf"$x-x_0={(x-xfarm)/Lfarm:.2f}L_f$",
        **plot_style,
        **field_style,
        **axis_style,
        **overlay_style,
        **colorbar_style,
        **font_style,
        export=f"S800_buoyancy_x={x}.png",
    )

# station_blh_files = blh_files(sim, {"axis":'y', 'value':yfarm})
# plot_slices(
#     [
#         {
#             "filename": [
#                 os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term01_t"+sim["timestamp"]+"_n"+sim["nstamp"]+"_SL_y=79p3625.nc"),
#             ],
#             "name": r"$\Delta u$",
#             "blh_file": station_blh_files,
#             "fieldgain": 1,
#             "s":0.14,
#         },

#         {
#             "filename": [
#                 os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term18_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_y=79p3625.nc"),
#                 os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term21_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_y=79p3625.nc"),
#                 os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term24_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_y=79p3625.nc"),
#             ],
#             "name": r"$-\partial_z \Delta \tau_{xz} (\times 10^{-3})$",
#             "blh_file": station_blh_files,
#             "fieldgain": -1000,
#             #"s":13,
#         },

#         {
#             "filename": [
#                 os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term17_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_y=79p3625.nc"),
#                 os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term20_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_y=79p3625.nc"),
#                 os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term23_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_y=79p3625.nc"),
#             ],
#             "name": r"$-\partial_y \Delta \tau_{xy} (\times 10^{-3})$",
#             "blh_file": station_blh_files,
#             "fieldgain": -1000,
#             #"s":13,
#         },

#         {
#             "filename": [
#                 os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term18_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_y=79p3625.nc"),
#             ],
#             "name": r"$-\partial_x \Delta p\,\left(\times 10^{-3}\right)$",
#             "blh_file": station_blh_files,
#             "fieldgain": -1000,
#             #"s":14,
#         },

#         {
#             "filename": [
#                 os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term20_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_y=79p3625.nc"),
#             ],
#             "name": r"$-\partial_z \Delta p\,\left(\times 10^{-3}\right)$",
#             "blh_file": station_blh_files,
#             "fieldgain": -1000,
#             #"s":14,
#         },
        
#     ],

#     #stamp=rf"$x-x_0={(x-xfarm)/Lfarm:.2f}L_f$",
#     **plot_style,
#     **field_style,
#     **axis_style,
#     **overlay_style,
#     **colorbar_style,
#     **font_style,
#     export=f"S250_midYplane.png",
# )
