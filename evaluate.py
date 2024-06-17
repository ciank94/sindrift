#todo: add functionality for plotting diagnostics from simulation- add options so that the worms plot can be done alone as a check
from configure import FileExplorer
from compare import PlotData, CompareData
import numpy as np

# Ready the file paths for analysis
fpath = FileExplorer(node='remote', model_name='cmems', key="APSO")
fpath.mounted_paths()

# setup analysis over several years
y1 = 2017
y2 = 2020
r1 = 1 # release start
r2 = 7 # release end
df = CompareData(y1, y2, r1, r2)

counter_y = -1
for y in range(y1, y2+1, 1):
    counter_y = counter_y + 1
    fpath.search_path(year=y, release_start=r1, release_end=r2)  # select files which should be plotted

    counter_r = -1
    dom_paths = 0
    for file_i in fpath.file_list:
        counter_r = counter_r + 1
        print('Analysing file: ' + file_i)
        pld = PlotData(fpath, file_i)
        breakpoint()

        # count recruits to area:
        recruit_number, recruit_time = pld.count_recruits()
        df.recruit_values(recruit_number, recruit_time, counter_r, counter_y)

        #pld.polygon_boundaries()

        pld.analysis_df.close()  # close file after use- due to file locking in local python script;

        # times
        pld.read_trajectory_df()
        if counter_r == 0:
            pld.plot_recruits()
        pld.trajectory_df.close()  # close file after use- due to file locking in local python script;
        df.time_bounds(pld, counter_r, counter_y)

        if counter_r == 0:
            dom_paths = pld.dom_paths
        else:
            dom_paths = dom_paths + pld.dom_paths


        #
        #pld.plot_trajectory()

    pld.plot_dom_paths(dom_paths)
    del pld

        #breakpoint()
breakpoint()
