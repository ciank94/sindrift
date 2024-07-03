#todo: add functionality for plotting diagnostics from simulation- add options so that the worms plot can be done alone as a check
from configure import FileExplorer
from results import PlotData, ReadAnalysis, StoreRelease
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc

# Ready the file paths for analysis
node_name = 'local'
file_loc = 'remote'
fpath = FileExplorer(node=node_name, model_name='cmems', key="SOIN")
year_range = [2016, 2016]
release_range = [1, 3]
if node_name == 'local':
    if file_loc == 'remote':
        fpath.mounted_paths()
    else:
        fpath.local_phys_states()


for year in range(year_range[0], year_range[1] + 1):
    fpath.search_path(year=year, release_start = release_range[0], release_end = release_range[1])  # select files which should be plotted
    df = StoreRelease(fpath)  # store information about simulations;
    for file_i in fpath.file_list:
        # read analysis file for storage:
        rd = ReadAnalysis(fpath, analysis_file=file_i)
        rd.read_analysis_df()

        # if it is the first release in the list, initialise variables to be stored
        df.counter()
        if df.counter_r == 0:
            df.init_variables(rd)

        df.store_variables(rd)
        rd.analysis_df.close()

        if df.counter_r == 0:
            rd.read_trajectory_df()
            df.store_trajectory_variables(rd)
            rd.trajectory_df.close()

    df.write_data()  # write data on recruits to file;

    pld = PlotData(df)  # initialise a plotting function
    pld.plot_time_recruits(df)
    pld.plot_site_recruits(df)  # plot recruits from sites
    pld.plot_dom_paths(df)  # plot the dominant pathways for year
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
