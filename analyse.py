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
key_ids = ["APSO"]
year_f = 2000
year_e = 2001
year_ids = np.arange(year_f, year_e+1, 1)
release_f = 1
release_e = 3
release_ids = np.arange(release_f, release_e+1, 1)

# Fuse datasets and get statistics
fuse = FuseData(infile_path, outfile_path, year_ids, release_ids, key_ids)

# Plot relevant analysis
p = PlotData(fuse)
p.plot_dom_paths(fuse)
breakpoint()

#y_end = 2011
import netCDF4 as nc
f = nc.Dataset(tp_file)
lat = f.variables['lat']
lon = f.variables['lon']
lat_m = np.zeros([np.shape(lat)[0]])
lon_m = np.zeros([np.shape(lat)[0]])
for i in range(0, np.shape(lat)[0]):
    lat_m[i] = np.nanmedian(lat[i, :])
    lon_m[i] = np.nanmedian(lon[i, :])
breakpoint()
#breakpoint()
y_i  = 2001
r_i=0
key="APSO"
test=False
#Process(trajectory_file, outfile_path, y_i, r_i, key, test)
breakpoint()
