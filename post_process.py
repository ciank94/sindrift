import os
import numpy as np
import netCDF4 as nc
from netCDF4 import num2date
import geopandas as gpd
from shapely.geometry import Polygon
import xarray as xr
import matplotlib.pyplot as plt


class PostProcess:

    def __init__(self, fpath, file_i, test):
        # files for analysis
        self.fpath = fpath
        self.analysis_file = self.fpath.analysis_path + file_i
        self.analysis_df = nc.Dataset(self.analysis_file, mode='r+')
        self.tr_file = fpath.trajectory_path + self.analysis_df.trajectory_file
        self.read_trajectory()

        # add from analysis file
        self.lon_bin_vals = self.analysis_df.variables['lon_bin_vals']
        self.lat_bin_vals = self.analysis_df.variables['lat_bin_vals']

        if test:
            self.p_stride = np.ceil(self.shp_p/20).astype(int)
            self.t_stride = np.ceil(self.shp_t/20).astype(int)
        else:
            self.p_stride = 1
            self.t_stride = 1
        return

    def read_trajectory(self):
        # extract data from trajectory file
        self.trajectory_df = nc.Dataset(self.tr_file, mode='r')
        times = self.trajectory_df.variables['time']
        self.dates = num2date(times, times.units)
        self.lat = np.array(self.trajectory_df['lat'][:])
        self.lat[self.lat > 1e06] = np.nan
        self.lon = np.array(self.trajectory_df['lon'][:])
        self.lon[self.lon > 1e06] = np.nan
        self.shp_p = np.shape(self.lat)[0]
        self.shp_t = np.shape(self.lat)[1]
        return
        # self.z = self.nc_file['z']

    def trajectory_analysis(self):
        # Set up a netcdf file for storing output from the analysis:
        self.analysis_df.variables['dom_paths'][:] = 0  # Initialize dom_path matrix
        self.analysis_df.variables['recruit_SG_north'][:] = 0

        # find bins for each particle
        id_lat = np.digitize(self.lat[:,::self.t_stride], self.lat_bin_vals)
        id_lon = np.digitize(self.lon[:,::self.t_stride], self.lon_bin_vals)

        dom_paths = np.zeros([np.shape(self.lon_bin_vals)[0], np.shape(self.lat_bin_vals)[0]])
        self.transit_hours = np.zeros(self.shp_p)
        self.visit_index = np.zeros(self.shp_p)

        # loop through particles and check their positions
        for p_i in range(0, self.shp_p, self.p_stride):
            self.p_i = p_i
            print("Particle analysis : " + str(self.p_i + 1) + " of " + str(self.shp_p))  # print particle id number in analysis
            idx_p = id_lat[p_i, :] < np.shape(self.lat_bin_vals)[0]
            dom_paths[id_lon[p_i, idx_p], id_lat[p_i, idx_p]] = dom_paths[id_lon[p_i, idx_p], id_lat[p_i, idx_p]] + 1
            # analyse retention/ recruitment in the following function
            self.retention(self.lon[p_i, :], self.lat[p_i, :])

        self.analysis_df.variables['dom_paths'][:] = dom_paths
        self.analysis_df.variables['recruit_SG_north'][:, 0] = self.transit_hours
        self.analysis_df.variables['recruit_SG_north'][:, 1] = self.visit_index
        return

    def temporal_analysis(self):
        CG_lon_lat = np.zeros([self.shp_t, 2])
        CG_lon_lat[:, 0] = np.nanmedian(self.lon, 0)
        CG_lon_lat[:, 1] = np.nanmedian(self.lat, 0)
        try:
            self.analysis_df.variables['CG'][:] = CG_lon_lat
        except:
            self.analysis_df.variables['CG'][:] = np.nan
        return

    def retention(self, lon, lat):
        lon_lim = [-39.5, -35]
        lat_lim = [-54, -53]
        idx_recruit = (lon > lon_lim[0]) & (lon < lon_lim[1]) & (lat > lat_lim[0]) & (lat < lat_lim[1])
        if not np.sum(idx_recruit) == 0:
            visits = np.where(idx_recruit>0)
            self.visit_index[self.p_i] = visits[0][0]
            t_to_poly = self.dates[visits[0][0]] - self.dates[0]
            self.transit_hours[self.p_i] = ((t_to_poly.days * 24) +
                                            np.floor(t_to_poly.seconds * 1 / (60 * 60)).astype(int))
        return

    def init_ncfile(self):
        # As I am appending to a file, I should check if dimensions or variables already exist
        dimension_key_dict = {'transit_info': 2, 'CG' : 2, 'n_parts': self.shp_p, 'obs': self.shp_t}

        for dimension in dimension_key_dict:
            if dimension not in self.analysis_df.dimensions.keys():
                self.analysis_df.createDimension(dimension, dimension_key_dict[dimension])

        variable_key_dict = {'dom_paths': {'datatype':'i4','dimensions':('n_lon_bins', 'n_lat_bins'), 'description': 'unique particle visits'},
                             'recruit_SG_north': {'datatype':'i4','dimensions':('n_parts', 'transit_info'), 'description': 'transit hours to polygon'},
                             'CG': {'datatype':'f4','dimensions':('obs', 'CG'), 'description': 'average longitude and latitude over time'}}

        for variable in variable_key_dict:
            if variable not in self.analysis_df.variables.keys():
                self.analysis_df.createVariable(variable, variable_key_dict[variable]['datatype'], variable_key_dict[variable]['dimensions'])
                self.analysis_df[variable].description = variable_key_dict[variable]['description']
        return












