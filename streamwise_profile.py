from utilities import streamwise_bulk_profiles, read_netcdf_slice
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
dz = 10/126

S800_h2 = np.mean(read_netcdf_slice(os.path.join(S800['path'], f'Run09_t{S800["timestamp"]}_INVH2.nc'), scale=1/126)["data"])
S250_h2 = np.mean(read_netcdf_slice(os.path.join(S250['path'], f'Run09_t{S250["timestamp"]}_INVH2.nc'), scale=1/126)["data"])

S800_h0 = np.mean(read_netcdf_slice(os.path.join(S800['path'], f'Run09_t{S800["timestamp"]}_INVH0.nc'), scale=1/126)["data"])

profiles = streamwise_bulk_profiles(
      [
          {
              "files":
                [
                    os.path.join(S800["path"], f"Run09_comp_deficit_budget5_term01_t"+S800["timestamp"]+"_n"+S800["nstamp"]+f"_integ_y.nc"),
                    os.path.join(S800["path"], f"Run09_comp_deficit_budget5_term02_t"+S800["timestamp"]+"_n"+S800["nstamp"]+f"_integ_y.nc"),
                    os.path.join(S800["path"], f"Run09_comp_deficit_budget5_term03_t"+S800["timestamp"]+"_n"+S800["nstamp"]+f"_integ_y.nc"),
                    os.path.join(S800["path"], f"Run09_comp_deficit_budget5_term04_t"+S800["timestamp"]+"_n"+S800["nstamp"]+f"_integ_y.nc"),
                    os.path.join(S800["path"], f"Run09_comp_deficit_budget5_term05_t"+S800["timestamp"]+"_n"+S800["nstamp"]+f"_integ_y.nc"),
                    os.path.join(S800["path"], f"Run09_comp_deficit_budget5_term06_t"+S800["timestamp"]+"_n"+S800["nstamp"]+f"_integ_y.nc"),
                    os.path.join(S800["path"], f"Run09_comp_deficit_budget5_term07_t"+S800["timestamp"]+"_n"+S800["nstamp"]+f"_integ_y.nc"),
                    os.path.join(S800["path"], f"Run09_comp_deficit_budget5_term08_t"+S800["timestamp"]+"_n"+S800["nstamp"]+f"_integ_y.nc"),
                    os.path.join(S800["path"], f"Run09_comp_deficit_budget5_term09_t"+S800["timestamp"]+"_n"+S800["nstamp"]+f"_integ_y.nc"),
                  ],
              "bounds": {"z": (0, 1000*S800_h2)},
              #"label": r"$\mathcal{A}\left(\Delta u; 0,L_y; 0, 15 z_h\right)$",
              "label": r"S800",
              "name": "A",
              "color": "k",
              "zorder": 2,
              "linewidth": 1.2,
              "linestyle": "-",
              "marker": None,
              "alpha": 1.0,
              "fieldgain": dz/(Wfarm*Zfarm),
              "lengthgain": 1.0,
              "smooth": False,
              "l_smooth": 5,
              "skip_nan": True,
              "inclusive": True,
              "varname": None,
              "indexing": "ij",
          },

          {
              "files": [
                  os.path.join(S250["path"], f"Run09_comp_deficit_budget5_term01_t"+S250["timestamp"]+"_n"+S250["nstamp"]+f"_integ_y.nc"),
                  os.path.join(S250["path"], f"Run09_comp_deficit_budget5_term02_t"+S250["timestamp"]+"_n"+S250["nstamp"]+f"_integ_y.nc"),
                  os.path.join(S250["path"], f"Run09_comp_deficit_budget5_term03_t"+S250["timestamp"]+"_n"+S250["nstamp"]+f"_integ_y.nc"),
                  os.path.join(S250["path"], f"Run09_comp_deficit_budget5_term04_t"+S250["timestamp"]+"_n"+S250["nstamp"]+f"_integ_y.nc"),
                  os.path.join(S250["path"], f"Run09_comp_deficit_budget5_term05_t"+S250["timestamp"]+"_n"+S250["nstamp"]+f"_integ_y.nc"),
                  os.path.join(S250["path"], f"Run09_comp_deficit_budget5_term06_t"+S250["timestamp"]+"_n"+S250["nstamp"]+f"_integ_y.nc"),
                  os.path.join(S250["path"], f"Run09_comp_deficit_budget5_term07_t"+S250["timestamp"]+"_n"+S250["nstamp"]+f"_integ_y.nc"),
                  os.path.join(S250["path"], f"Run09_comp_deficit_budget5_term08_t"+S250["timestamp"]+"_n"+S250["nstamp"]+f"_integ_y.nc"),
                  os.path.join(S250["path"], f"Run09_comp_deficit_budget5_term09_t"+S250["timestamp"]+"_n"+S250["nstamp"]+f"_integ_y.nc"),
                  ],
              "bounds": {"z": (0, 1000*S250_h2)},
              "label": r"S250",
              "name": "A",
              "color": "tab:red",
              "zorder": 3,
              "linewidth": 1.2,
              "linestyle": "-",
              "marker": None,
              "alpha": 1.0,
              "fieldgain": dz/(Wfarm*Zfarm),
              "lengthgain": 1.0,
              "smooth": False,
              "l_smooth": 5,
              "skip_nan": True,
              "inclusive": True,
              "varname": None,
              "indexing": "ij",
          },
      ],
      bounds=None,
      export="Adu.png",
      plotfarm=True,
      farm_bounds=(xfarm, xfarm+Lfarm),
      zeroline=True,
      scale=1.0,
      xlim=None,
      xlft=10,
      xrght=540,
      ylim=None,
      ybot=None,
      ytop=None,
      xaxis_transformer={
          "function": lambda x: (x-xfarm)/Lfarm,
          "ticks": None,
          "label": r"$(x-x_0)/L_f$",
      },
      xstep=None,
      xlabel=None,
      ylabel=r"$\left\langle D_t\Delta u \right\rangle_{yz}$",
      title=None,
      data_key="data",
      skip_nan=True,
      inclusive=True,
      profile_value="sum",
      fieldgain=1.0,
      lengthgain=1.0,
      smooth=False,
      l_smooth=5,
      linewidth=1.2,
      linestyle="-",
      marker=None,
      alpha=None,
      color=None,
      zorder=None,
      figsize=(6, 4),
      dpi=300,
      legend=True,
      legend_ncols=1,
      legend_loc=None,
      legend_frameon=False,
      legend_kwargs={},
      grid=False,
      tight_layout=True,
      close=True,
      tick_fontsize=7,
      label_fontsize=9,
      title_fontsize=10,
      legend_fontsize=10,
      bbox_inches="tight",
      pad_inches=0.02,
      varname=None,
      indexing="ij",
  )