#todo: add functionality for plotting diagnostics from simulation- add options so that the worms plot can be done alone as a check
from configure import FileExplorer

# Ready the file paths for analysis
fpath = FileExplorer(node='local', model_name='sinmod')
fpath.search_path()  # select files which should be plotted

breakpoint()



import sys
#from post_process import Process, locate_files
#from plotting import FuseData, PlotData
#import numpy as np

# sys.path.insert(0, 'C:/Users/ciank/PycharmProjects/sinmod/opendrift')  # add opendrift local path
node = 'local'
model_name = 'sinmod'
infile_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/results/'
outfile_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/results/'
#infile_path = outfile_path

#trajectory_file = infile_path + "APSO_2001_R1_trajectory.nc"
#tp_file = outfile_path + "APSO_2001_R1_trajectory_tp.nc"


# Parameters for analysis
key_id = "SG800"
release_f = 1
release_e = 4
release_ids = np.arange(release_f, release_e+1, 1)
years = np.arange(2020, 2021, 1)

for year_id in years:


    # Fuse datasets and get statistics
    fuse = FuseData(infile_path, outfile_path, year_id, release_ids, key_id)
    # Plot relevant analysis
    p = PlotData(fuse)
    p.plot_dom_paths(fuse)