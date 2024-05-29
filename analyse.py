# Firstly, specify location of netcdf input files
import sys
from post_process import Process
from plotting import FuseData, PlotData
import numpy as np

# sys.path.insert(0, 'C:/Users/ciank/PycharmProjects/sinmod/opendrift')  # add opendrift local path
infile_path = 'A:/Cian_sinmod/copernicus_client/results/'
outfile_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/results/'
#infile_path = outfile_path

#trajectory_file = infile_path + "APSO_2001_R1_trajectory.nc"
#tp_file = outfile_path + "APSO_2001_R1_trajectory_tp.nc"


# Parameters for analysis
key_id = "APSO"

release_f = 1
release_e = 25
release_ids = np.arange(release_f, release_e+1, 1)
years = np.arange(2000, 2004, 1)

for year_id in years:


    # Fuse datasets and get statistics
    fuse = FuseData(infile_path, outfile_path, year_id, release_ids, key_id)
    # Plot relevant analysis
    p = PlotData(fuse)
    p.plot_dom_paths(fuse)



breakpoint()
