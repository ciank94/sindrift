# Firstly, specify location of netcdf input files
import sys
from post_process import Process
from plotting import FuseData, PlotData
import numpy as np

# sys.path.insert(0, 'C:/Users/ciank/PycharmProjects/sinmod/opendrift')  # add opendrift local path
#infile_path = 'A:/Cian_sinmod/copernicus_client/results/'
outfile_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/results/'
infile_path = outfile_path



# Parameters for analysis
key_ids = ["APSO"]
year_f = 2001
year_e = 2001
year_ids = np.arange(year_f, year_e+1, 1)
release_f = 11
release_e = 11
release_ids = np.arange(release_f, release_e+1, 1)

# Fuse datasets and get statistics
fuse = FuseData(infile_path, outfile_path, year_ids, release_ids, key_ids)

# Plot relevant analysis
p = PlotData(fuse)
p.plot_dom_paths(fuse)

#y_end = 2011

#trajectory_file = infile_path + "APSO_2001_R10_trajectory.nc"
#y_i  = 2001
#r_i=10
#key="APSO"
#test=False
#Process(trajectory_file, outfile_path, y_i, r_i, key, test)
breakpoint()