#todo: add functionality for plotting diagnostics from simulation- add options so that the worms plot can be done alone as a check
from configure import FileExplorer
from plot_data import PlotData
import numpy as np

# Ready the file paths for analysis
fpath = FileExplorer(node='local', model_name='cmems', key="APSO")
#fpath.mounted_paths()
fpath.search_path()  # select files which should be plotted
n_files = np.shape(fpath.file_list)[0]
recruit_number = np.zeros([n_files])
counter = -1
for file_v in fpath.file_list:
    counter = counter + 1
    analysis_file = fpath.analysis_path + file_v
    print('Analysing file: ' + analysis_file)
    pld = PlotData(fpath, analysis_file)
    pld.analysis_df.close()  # close file after use- due to file locking in local python script;
    pld.read_trajectory_df()
    recruit_number[counter] = pld.plot_recruits()
    pld.plot_trajectory()
    pld.plot_dom_paths()
    pld.trajectory_df.close()  # close file after use- due to file locking in local python script;
    breakpoint()
