# Firstly, specify location of netcdf input files
node = 'local'
if node == 'local':
    import sys
    sys.path.insert(0, 'C:/Users/ciank/PycharmProjects/sinmod/opendrift')  # add opendrift local path
    infile_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/results/'
    outfile_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/results/'
else:  # assume job is on server
    infile_path = '/cluster/projects/nn9828k/Cian_sinmod/copernicus_client/results/'
    outfile_path = infile_path

from plotting import FuseData, PlotData
import numpy as np

# Parameters for analysis
key_ids = ["APSO"]
year_f = 2000
year_e = 2000
year_ids = np.arange(year_f, year_e+1, 1)
release_f = 1
release_e = 1
release_ids = np.arange(release_f, release_e+1, 1)

# Fuse datasets and get statistics
fuse = FuseData(infile_path, year_ids, release_ids, key_ids)

# Plot relevant analysis
p = PlotData(fuse)
p.plot_dom_paths(fuse)
breakpoint()
#y_end = 2011