#todo: add functionality for plotting diagnostics from simulation- add options so that the worms plot can be done alone as a check
from configure import FileExplorer
from compile_releases import StoreReleases
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc

# Ready the file paths for analysis
node_name = 'local'
file_loc = 'remote'
fpath = FileExplorer(node=node_name, model_name='cmems', key="SOIN")
if node_name == 'local':
    if file_loc == 'remote':
        fpath.mounted_paths()
    else:
        fpath.local_phys_states()

fpath.search_path(year=2016, release_start=1, release_end=3)  # select files which should be analysed


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
