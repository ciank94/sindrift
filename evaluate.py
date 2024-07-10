#todo: add functionality for plotting diagnostics from simulation- add options so that the worms plot can be done alone as a check
from configure import FileExplorer
from compile_releases import StoreReleases
from plot_data import PlotData
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
compile_folder = 'A:/Cian_sinmod/opendrift/' + 'compile/'
analysis_folder = 'A:/Cian_sinmod/opendrift/' + 'analysis/'
years = ['2016', '2017', '2018', '2019']
release_number = 10
# switches for different analysis
case_dom = False
case_srec = False
case_env = True

for y in years:
    p_plot = PlotData(key='SOIN', year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
    if case_dom:
        p_plot.plot_dom_paths(release_n=release_number)
    if case_srec:
        p_plot.plot_site_recruit_t(release_n=release_number)
    if case_env:
        p_plot.plot_hist_environment()






    #filename = file_prefix + 'CG_lat.npy'

    #CG_lat = np.load(filename)

    #filename = file_prefix + 'temp_exp.npy'
    #temp_exp = np.load(filename)


# Ready the file paths for analysis


# pseudocode:
# load files using numpy- get data needed


# pld.plot_CG_paths(df)
#
#     counter_r = counter_r + 1
#     print('Analysing file: ' + file_i)
#     pld = PlotData(fpath, file_i)
#
#     # count recruits to area:
#     recruit_number, recruit_time = pld.count_recruits()
#     df.recruit_values(recruit_number, recruit_time, counter_r, counter_y)
#
#     # pld.polygon_boundaries()
#     pld.init_plot()
#     pld.plot_background()
#     CG = pld.analysis_df.variables['CG']
#     plt.plot(CG[:, 0], CG[:, 1])
#     breakpoint()
#
#     pld.analysis_df.close()  # close file after use- due to file locking in local python script;
#
#     # times
#     pld.read_trajectory_df()
#     pld.plot_init()
#     if counter_r == 0:
#         pld.plot_recruits()
#     pld.trajectory_df.close()  # close file after use- due to file locking in local python script;
#     df.time_bounds(pld, counter_r, counter_y)
#
#     if counter_r == 0:
#         dom_paths = pld.dom_paths
#     else:
#         dom_paths = dom_paths + pld.dom_paths
#
#     #
#     # pld.plot_trajectory()
#
# pld.plot_dom_paths(dom_paths)
# del pld



        #breakpoint()
