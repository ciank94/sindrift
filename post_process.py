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

    def trajectory_analysis(self):
        # set up a netcdf file for storing output from the analysis:
        self.initialize_variables()

        # find bins for each particle
        self.id_lat = np.digitize(self.lat[:,::self.t_stride], self.lat_bin_vals)
        self.id_lon = np.digitize(self.lon[:,::self.t_stride], self.lon_bin_vals)

        # loop through particles and check their positions
        self.particle_analysis()
        self.environment_analysis()
        return

    def initialize_variables(self):
        # initialize variables to be stored
        self.dom_paths = np.zeros([np.shape(self.lon_bin_vals)[0], np.shape(self.lat_bin_vals)[0]])
        self.recruit_dom_paths = np.zeros([np.shape(self.lon_bin_vals)[0], np.shape(self.lat_bin_vals)[0]])
        self.transit_hours = np.zeros(self.shp_p)
        self.visit_index = np.zeros(self.shp_p)
        self.CG_lon = np.zeros([self.shp_t])
        self.CG_lat = np.zeros([self.shp_t])
        self.chl_exp = np.zeros([self.shp_t])
        self.o2_exp = np.zeros([self.shp_t])
        self.temp_exp = np.zeros([self.shp_t])
        self.analysis_df.variables['dom_paths'][:] = 0  # initialize dom_path matrix
        if self.analysis_df.bounds == 'NEMO':
            self.analysis_df.variables['recruit_SG_north'][:] = 0  # initialize recruitment to SG north
        elif self.analysis_df.bounds == 'SGret':
            self.analysis_df.variables['retain_SG_north'][:] = 0  # initialize recruitment to SG north
        return

    def particle_analysis(self):
        if self.analysis_df.bounds == 'NEMO':
            for p_i in range(0, self.shp_p, self.p_stride):
                self.p_i = p_i
                print("Particle analysis : " + str(self.p_i + 1) + " of " + str(self.shp_p))  # print particle id number in analysis
                idx_p = self.id_lat[p_i, :] < np.shape(self.lat_bin_vals)[0]
                self.dom_paths[self.id_lon[p_i, idx_p], self.id_lat[p_i, idx_p]] = (
                        self.dom_paths[self.id_lon[p_i, idx_p], self.id_lat[p_i, idx_p]] + 1)
                # analyse retention/ recruitment in the following function
                self.recruit_SG_north(self.lon[p_i, :], self.lat[p_i, :])
                if self.visit_index[self.p_i] > 0:
                    idx_v = self.visit_index[self.p_i].astype(int)
                    range_v = range(idx_v-100, idx_v)
                    self.recruit_dom_paths[self.id_lon[p_i, range_v], self.id_lat[p_i, range_v]] = (
                            self.recruit_dom_paths[self.id_lon[p_i, range_v], self.id_lat[p_i, range_v]] +1)

            self.analysis_df.variables['recruit_dom_paths'][:] = self.recruit_dom_paths
            self.analysis_df.variables['dom_paths'][:] = self.dom_paths
            self.analysis_df.variables['recruit_SG_north'][:, 0] = self.transit_hours
            self.analysis_df.variables['recruit_SG_north'][:, 1] = self.visit_index
            self.analysis_df['time'][:] = self.times[:]
            self.analysis_df['time'].units = self.times.units

        if self.analysis_df.bounds == 'SGret':
            for p_i in range(0, self.shp_p, self.p_stride):
                self.p_i = p_i
                print("Particle analysis : " + str(self.p_i + 1) + " of " + str(self.shp_p))  # print particle id number in analysis
                idx_p = self.id_lat[p_i, :] < np.shape(self.lat_bin_vals)[0]
                self.dom_paths[self.id_lon[p_i, idx_p], self.id_lat[p_i, idx_p]] = (
                        self.dom_paths[self.id_lon[p_i, idx_p], self.id_lat[p_i, idx_p]] + 1)
                # analyse retention/ recruitment in the following function
                self.retain_SG_north(self.lon[p_i, :], self.lat[p_i, :])

            self.analysis_df.variables['dom_paths'][:] = self.dom_paths
            self.analysis_df.variables['retain_SG_north'][:, 0] = self.transit_hours
            self.analysis_df.variables['retain_SG_north'][:, 1] = self.visit_index
        return




    def environment_analysis(self):
        self.read_bio()   # read NEMO ecosystem model data
        self.read_phys()  # read NEMO physics model data

        for i in range(0, self.shp_t):
            self.t_i = i
            self.CG_lon[i] = np.nanmean(self.lon[:, i])
            self.CG_lat[i] = np.nanmean(self.lat[:, i])
            self.subset_bio()
            self.subset_phys()

        self.analysis_df.variables['CG_lon'][:] = self.CG_lon
        self.analysis_df.variables['CG_lat'][:] = self.CG_lat
        self.analysis_df.variables['chl_exp'][:] = self.chl_exp
        self.analysis_df.variables['o2_exp'][:] = self.o2_exp
        self.analysis_df.variables['temp_exp'][:] = self.temp_exp
        self.analysis_df.variables['lat_init'][:] = self.lat[:, 0]
        self.analysis_df.variables['lon_init'][:] = self.lon[:, 0]

        # close files:
        self.bio_states.close()
        self.phys_states.close()
        return

    def retain_SG_north(self, lon1, lat1):
        lon2 = self.lon[0, :]
        lat2 = self.lat[0, :]
        phi_1 = np.radians(lat1)
        phi_2 = np.radians(lat2)
        R = 6371000
        delta_phi = np.radians(lat2 - lat1)
        delta_lambda = np.radians(lon2 - lon1)

        a = np.sin(delta_phi / 2.0) ** 2 + \
            np.cos(phi_1) * np.cos(phi_2) * \
            np.sin(delta_lambda / 2.0) ** 2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        dist_km = R*c/1000
        idx_retain = dist_km > 200
        if not np.sum(idx_retain) == 0:
             visits = np.where(idx_retain > 0)
             self.visit_index[self.p_i] = visits[0][0]
             t_to_poly = self.dates[visits[0][0]] - self.dates[0]
             self.transit_hours[self.p_i] = ((t_to_poly.days * 24) +
                                             np.floor(t_to_poly.seconds * 1 / (60 * 60)).astype(int))

        return


    def recruit_SG_north(self, lon, lat):
        lon_lim = [-39.4, -35]
        lat_lim = [-55, -53.2]
        idx_recruit = (lon > lon_lim[0]) & (lon < lon_lim[1]) & (lat > lat_lim[0]) & (lat < lat_lim[1])
        if not np.sum(idx_recruit) == 0:
            visits = np.where(idx_recruit > 0)
            self.visit_index[self.p_i] = visits[0][0]
            t_to_poly = self.dates[visits[0][0]] - self.dates[0]
            self.transit_hours[self.p_i] = ((t_to_poly.days * 24) +
                                            np.floor(t_to_poly.seconds * 1 / (60 * 60)).astype(int))
        return


    def subset_phys(self):
        t_diff = abs(self.dates[self.t_i] - self.phys_times)
        id_t = t_diff.argmin(0)
        lon_diff = abs(self.CG_lon[self.t_i] - self.phys_lons)
        id_lon = lon_diff.argmin(0)
        lat_diff = abs(self.CG_lat[self.t_i] - self.phys_lats)
        id_lat = lat_diff.argmin(0)
        id_depth = 0
        self.temp_exp[self.t_i] = self.temp[id_t, id_depth, id_lat, id_lon]
        return

    def subset_bio(self):
        t_diff = abs(self.dates[self.t_i] - self.bio_times)
        id_t = t_diff.argmin(0)
        lon_diff = abs(self.CG_lon[self.t_i] - self.bio_lons)
        id_lon = lon_diff.argmin(0)
        lat_diff = abs(self.CG_lat[self.t_i] - self.bio_lats)
        id_lat = lat_diff.argmin(0)
        id_depth = 0
        self.chl_exp[self.t_i] = self.chl[id_t, id_depth, id_lat, id_lon]
        self.o2_exp[self.t_i] = self.o2[id_t, id_depth, id_lat, id_lon]
        return


    def read_bio(self):
        self.bio_file = self.fpath.phys_states_path + 'CMEMS_GLBIO_D_full_' + str(
        self.analysis_df.sim_start_year) + '.nc'
        self.bio_states = nc.Dataset(self.bio_file, mode='r')
        self.bio_times = num2date(self.bio_states['time'], self.bio_states['time'].units)
        self.chl = self.bio_states['chl']
        self.o2 = self.bio_states['o2']
        self.bio_lats = self.bio_states['latitude'][:]
        self.bio_lons = self.bio_states['longitude'][:]
        return

    def read_phys(self):
        self.phys_file = (self.fpath.phys_states_path + 'CMEMS_GLPHYS_D_full_' +
                          str(self.analysis_df.sim_start_year) + '.nc')
        self.phys_states = nc.Dataset(self.phys_file, mode='r')
        self.phys_times = num2date(self.bio_states['time'], self.bio_states['time'].units)
        self.temp = self.phys_states['thetao']
        self.phys_lats = self.phys_states['latitude'][:]
        self.phys_lons = self.phys_states['longitude'][:]
        return

    def read_trajectory(self):
        # extract data from trajectory file
        self.trajectory_df = nc.Dataset(self.tr_file, mode='r')
        self.times = self.trajectory_df.variables['time']
        self.dates = num2date(self.times, self.times.units)
        self.lat = np.array(self.trajectory_df['lat'][:])
        self.lat[self.lat > 1e06] = np.nan
        self.lat[self.lat == 0] = np.nan
        self.lon = np.array(self.trajectory_df['lon'][:])
        self.lon[self.lon > 1e06] = np.nan
        self.lon[self.lon == 0] = np.nan
        self.shp_p = np.shape(self.lat)[0]
        self.shp_t = np.shape(self.lat)[1]
        return

    def init_ncfile(self):
        # As I am appending to a file, I should check if dimensions or variables already exist
        dimension_key_dict = {'transit_info': 2, 'trajectory': self.shp_p, 'obs': self.shp_t}

        for dimension in dimension_key_dict:
            if dimension not in self.analysis_df.dimensions.keys():
                self.analysis_df.createDimension(dimension, dimension_key_dict[dimension])

        if self.analysis_df.bounds == 'NEMO':
            variable_key_dict = {'dom_paths': {'datatype': 'i4', 'dimensions': ('n_lon_bins', 'n_lat_bins'),
                                               'description': 'unique particle visits into SG'},
                                 'recruit_dom_paths':{'datatype': 'i4', 'dimensions': ('n_lon_bins', 'n_lat_bins'),
                                               'description': 'unique particle visits'},
                                 'recruit_SG_north': {'datatype': 'i4', 'dimensions': ('trajectory', 'transit_info'),
                                                      'description': 'transit hours to polygon'},
                                 'CG_lon': {'datatype': 'f4', 'dimensions': ('obs',),
                                            'description': 'average longitude over time'},
                                 'CG_lat': {'datatype': 'f4', 'dimensions': ('obs',),
                                            'description': 'average latitude over time'},
                                 'o2_exp': {'datatype': 'f4', 'dimensions': ('obs',),
                                            'description': 'experienced oxygen concentration along trajectory'},
                                 'chl_exp': {'datatype': 'f4', 'dimensions': ('obs',),
                                             'description': 'experienced chlorophyll concentration along trajectory'},
                                 'temp_exp': {'datatype': 'f4', 'dimensions': ('obs',),
                                              'description': 'experienced temperature along trajectory'},
                                 'lat_init': {'datatype': 'f4', 'dimensions': ('trajectory',),
                                              'description': 'initial latitude'},
                                 'lon_init': {'datatype': 'f4', 'dimensions': ('trajectory',),
                                              'description': 'initial longitude'},
                                 'time': {'datatype': 'f4', 'dimensions': ('obs',),
                                          'description': 'times saved for analysis'}
                                 }
        else:
            variable_key_dict = {'dom_paths': {'datatype': 'i4', 'dimensions': ('n_lon_bins', 'n_lat_bins'),
                                               'description': 'unique particle visits'},
                                 'retain_SG_north': {'datatype': 'i4', 'dimensions': ('trajectory', 'transit_info'),
                                                      'description': 'transit hours to polygon'},
                                 'CG_lon': {'datatype': 'f4', 'dimensions': ('obs',),
                                            'description': 'average longitude over time'},
                                 'CG_lat': {'datatype': 'f4', 'dimensions': ('obs',),
                                            'description': 'average latitude over time'},
                                 'o2_exp': {'datatype': 'f4', 'dimensions': ('obs',),
                                            'description': 'experienced oxygen concentration along trajectory'},
                                 'chl_exp': {'datatype': 'f4', 'dimensions': ('obs',),
                                             'description': 'experienced chlorophyll concentration along trajectory'},
                                 'temp_exp': {'datatype': 'f4', 'dimensions': ('obs',),
                                              'description': 'experienced temperature along trajectory'},
                                 'lat_init': {'datatype': 'f4', 'dimensions': ('trajectory',),
                                              'description': 'initial latitude'},
                                 'lon_init': {'datatype': 'f4', 'dimensions': ('trajectory',),
                                              'description': 'initial longitude'},
                                 'time': {'datatype': 'f4', 'dimensions': ('obs',),
                                          'description': 'times saved for analysis'}
                                 }


        for variable in variable_key_dict:
            if variable not in self.analysis_df.variables.keys():
                self.analysis_df.createVariable(variable, variable_key_dict[variable]['datatype'],
                                                variable_key_dict[variable]['dimensions'])
                self.analysis_df[variable].description = variable_key_dict[variable]['description']
        return














