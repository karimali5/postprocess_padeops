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
zfarm = (95-126/2)/126
Lfarm = 40
Wfarm = 95.61-63.11
Zfarm = 1
hubheight = 95/126

# S800_h0 = np.mean(read_netcdf_slice(os.path.join(S800['path'], f'Run09_t{S800["timestamp"]}_INVH0.nc'), scale=1/126)["data"])
# S250_h0 = np.mean(read_netcdf_slice(os.path.join(S250['path'], f'Run09_t{S250["timestamp"]}_INVH0.nc'), scale=1/126)["data"])

sim = S800_5K

def blh_files(sim, cutaxis):
    return [
        {'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_INVH0.nc'), cutaxis['axis'], cutaxis['value'], plot=False, scale=1/126, smooth=True, l_smooth=5), 'color':"k", "style":"--"},
        {'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_INVH2.nc'), cutaxis['axis'], cutaxis['value'], plot=False, scale=1/126, smooth=True, l_smooth=5), 'color':"k", "style":"--"},
        {'data': transect_netcdf_slice(os.path.join(sim['path'], f'Run09_t{sim["timestamp"]}_SL_BLH.nc'), cutaxis['axis'], cutaxis['value'], plot=False, smooth=True, l_smooth=2), 'color':"tab:red", "style":"-"},
    ]

plot_style = {
    "figsize": (5.5, 3),
    #"figsize": (8, 2.8),
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
    "zero_ground_first_row": False,
}

axis_style = {
    "equal_aspect": True,
    # "ytop": 3*hubheight,
    # "ybot": 15/126,
    # "xrght": 0.65*Wfarm+yfarm,
    # "xlft": -0.5*Wfarm+yfarm,
    "xlabel": r"$y$",
    "ylabel": r"$z$",
    "aspect": None,
    "xaxis_transformer": {
        "function": lambda x: (x - yfarm) / Wfarm,
        "inverse": lambda xt: xt * Wfarm + yfarm,
        #"ticks": np.arange(-0.4, 0.51, 0.1),
        "ticks": [-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4],
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
        "ticks": np.arange(0, 3.1, 1),
        "label": r"$z/z_h$",
    },
    "transform_axes_data": False,
}

overlay_style = {
    "plotfarm": False,
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
    #station_blh_files = blh_files(sim, {"axis":'x', 'value':x})
    plot_slices(
        [
            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget1_term01_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget1_term07_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget1_term07_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),

            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget1_term04_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget1_term12_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget1_term12_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),

            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget1_term06_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget1_term15_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget1_term15_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta k$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 0.5,
            #     #"s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term01_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term02_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term03_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term04_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term05_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta (u_j \partial_j k)\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     #"s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term06_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term07_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term08_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term09_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term10_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term11_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term12_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta \mathcal{P}_k\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":1.5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term07_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term08_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term09_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term10_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term11_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term12_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta \mathcal{P}_k - \Delta \mathcal{P}_{6}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     #"s":0.14,
            # },

            {
                "filename": [
                    os.path.join(sim["path"], f"Run09_comp_deficit_budget6_term07_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
                ],
                #"scale": [1,-100],
                "name": r"$\overline{u}^b \partial_x \Delta \overline{v}$",
                #"blh_file": station_blh_files,
                "fieldgain": 1000,
                #"s":0.03,
            },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget2_term06_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     #"scale": [1,-100],
            #     "name": r"$-\partial_j \overline{{u_3^\prime}^b \Delta u_j^\prime}~(\times10^{-3})$",
            #     #"blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget2_term09_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     #"scale": [1,-100],
            #     "name": r"$-\partial_j \overline{\Delta {u_3^\prime}  {u_j^\prime}^b}~(\times10^{-3})$",
            #     #"blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget2_term03_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget2_term06_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget2_term09_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     #"scale": [1,-100],
            #     "name": r"$-\Delta \partial_j \overline{{u_3^\prime}  {u_j^\prime}}~(\times10^{-3})$",
            #     #"blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term07_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\overline{\Delta u_i' u_j'^b} \partial_j \overline{ \Delta u_i}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":1.5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term08_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\overline{u_i'^b \Delta u_j'} \partial_j \overline{ \Delta u_i}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":1.5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term09_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\overline{u_i'^b u_j'^b} \partial_j \overline{ \Delta u_i}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":1.5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term10_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\overline{\Delta u_i' \Delta u_j'} \partial_j \overline{u_i^b}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":1.5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term11_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\overline{\Delta u_i' u_j'^b} \partial_j \overline{u_i^b}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":1.5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget4_term12_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\overline{u_i'^b \Delta u_j'} \partial_j \overline{u_i^b}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":1.5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term13_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term14_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term15_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term16_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term17_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term18_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term19_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta \mathcal{T}_k\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     #"s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term13_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta \mathcal{T}_{13}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     #"s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term14_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta \mathcal{T}_{14}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     #"s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term15_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta \mathcal{T}_{15}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     #"s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term16_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta \mathcal{T}_{16}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     #"s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term17_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta \mathcal{T}_{17}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     #"s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term18_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta \mathcal{T}_{18}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     #"s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget3_term19_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta \mathcal{T}_{19}\,(\times10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     #"s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term01_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta u$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1,
            #     "s":0.14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term07_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$U_0 \partial_x \Delta u\,\left(\times 10^{-3}\right)$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":100,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term01_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta u \partial_x \Delta u\,\left(\times 10^{-3}\right)$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":30,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term08_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$V_0 \partial_y \Delta u\,\left(\times 10^{-3}\right)$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term02_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta v \partial_y \Delta u\,\left(\times 10^{-3}\right)$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":15,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget5_term06_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$U_0^\prime \Delta w\,\left(\times 10^{-3}\right)$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":10,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term03_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta w$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1,
            #     "s":0.05,
            # },

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

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term20_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$-\partial_z \Delta p\,\left(\times 10^{-3}\right)$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     "s":14,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget0_term17_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta w_b\,\left(\times 10^{-3}\right)$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":6,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget8_term01_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta u_i u_j^b \partial_j u_i^b\,(\times 10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget8_term02_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$u_i^b u_j^b \partial_j \Delta u_i\,(\times 10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget8_term03_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta u_i u_j^b \partial_j \Delta u_i\,(\times 10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget8_term04_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$u_i^b \Delta u_j \partial_j u_i^b\,(\times 10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget8_term05_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta u_i \Delta u_j \partial_j u_i^b\,(\times 10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget8_term06_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$u_i^b \Delta u_j \partial_j \Delta u_i\,(\times 10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget8_term07_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$\Delta u_i \Delta u_j \partial_j \Delta u_i\,(\times 10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": 1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget8_term11_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$- \Delta u_i \partial_i \Delta p\,(\times 10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget8_term12_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$- u_i^b \partial_i \Delta p\,(\times 10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     "s":5,
            # },

            # {
            #     "filename": [
            #         os.path.join(sim["path"], f"Run09_comp_deficit_budget8_term13_t"+sim["timestamp"]+"_n"+sim["nstamp"]+f"_SL_x={x}.nc"),
            #     ],
            #     "name": r"$-\Delta u_i \partial_i p^b\,(\times 10^{-3})$",
            #     "blh_file": station_blh_files,
            #     "fieldgain": -1000,
            #     "s":5,
            # },
            
        ],
        #stamp=rf"$x-x_0={(x-xfarm)/Lfarm:.2f}L_f$",
        **plot_style,
        **field_style,
        **axis_style,
        **overlay_style,
        **colorbar_style,
        **font_style,
        export=f"U_dv={x}.png",
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
