#todo: add functionality for plotting diagnostics from simulation- add options so that the worms plot can be done alone as a check
from configure import FileExplorer
from compile_releases import StoreReleases
from plot_data import (PlotData, CatchData, plot_recruit_stat, plot_retain, plot_SOIN_recruit_dom_paths,
                       plot_BSSI_recruit_dom_paths, plot_temp_SG, plot_catch_points, plot_linreg, plot_arrivals,
                       plot_SG_rec_area, plot_worms, plot_poster_dom_paths)
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
compile_folder = 'A:/Cian_sinmod/opendrift/' + 'compile/'
analysis_folder = 'A:/Cian_sinmod/opendrift/' + 'analysis/'
trajectory_folder = 'A:/Cian_sinmod/opendrift/' + 'trajectory/'

# s2.1: results for catch data:
cdata = CatchData()
# catch facts
#cdata.catch_facts()
# figure 1:
#cdata.plot_fishing_season()
# figure 2:
#cdata.plot_lat_lon()
# figure 3:
#plot_catch_points(compile_folder, analysis_folder)

# s2.2: results for recruitment to SG:
# figure 1:
#plot_recruit_stat(compile_folder, analysis_folder)
#plot_linreg(compile_folder, analysis_folder)
# figure 2:
#plot_arrivals(compile_folder, analysis_folder)
#plot_retain(compile_folder, analysis_folder)
# figure 3:
#plot_SOIN_recruit_dom_paths(compile_folder, analysis_folder)
#plot_BSSI_recruit_dom_paths(compile_folder, analysis_folder)
# figure 4:
#plot_temp_SG(compile_folder, analysis_folder)
#plot_temp_month_SG(compile_folder, analysis_folder)
# figure 5: worm plots:
plot_worms(compile_folder, analysis_folder, trajectory_folder)
# SG area plot- conceptual figure;
#plot_SG_rec_area(compile_folder, analysis_folder)

# poster plot:
#plot_poster_dom_paths(compile_folder, analysis_folder)
breakpoint()
