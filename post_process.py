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
        tr_file = fpath.trajectory_path + self.analysis_df.trajectory_file
        self.trajectory_df = nc.Dataset(tr_file, mode='r')


        times = self.trajectory_df.variables['time']
        self.dates = num2date(times, times.units)

        # extract data from trajectory file
        self.lat = self.trajectory_df['lat']
        self.lon = self.trajectory_df['lon']
        self.shp_p = np.shape(self.lat)[0]
        self.shp_t = np.shape(self.lat)[1]
        # self.z = self.nc_file['z']

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
            if p_i == 0:
                self.init_polygon()
            self.recruit(self.lon[p_i, :], self.lat[p_i, :])

        self.analysis_df.variables['dom_paths'][:] = dom_paths
        self.analysis_df.variables['recruit_SG_north'][:, 0] = self.transit_hours
        self.analysis_df.variables['recruit_SG_north'][:, 1] = self.visit_index

        print('Closing: ' + self.analysis_file)
        os.system("echo closing analysis file")
        self.analysis_df.close()
        return

    def recruit(self, lon_p, lat_p):
        gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries.from_xy(lon_p, lat_p))  # Create Geo dataframe from positional data
        polya = gpd.GeoDataFrame(geometry=self.polygon)  # Geo dataframe from polygon data
        in_pol = gpd.tools.sjoin(gdf, polya, predicate="within", how='left')  # Check points are in polygon , predicate="within", how='left'
        in_idx = np.where(in_pol.index_right > -1)  # Find all values in the polygon
        p_id = np.zeros(np.shape(lon_p)[0])  # Initialize vector for storage
        p_id[in_idx[0]] = 1  # Store as ones
        if np.sum(p_id) > 0:
            visits = np.where(p_id == 1)
            self.visit_index[self.p_i] = visits[0][0]
            t_to_poly = self.dates[visits[0][0]] - self.dates[0]
            self.transit_hours[self.p_i] = (t_to_poly.days * 24) + np.floor(t_to_poly.seconds * 1 / (60 * 60)).astype(int)
        return

    def init_polygon(self):
        import matplotlib.pyplot as plt
        # lon_1, lat_1 are bottom left coordinates
        lon_lims = [-39.5, -35]
        lat_lims = [-54, -53]
        self.lon_1 = lon_lims[0]
        self.lon_2 = lon_lims[1]
        self.lat_1 = lat_lims[0]
        self.lat_2 = lat_lims[1]
        coords = ((self.lon_1, self.lat_1), (self.lon_1, self.lat_2), (self.lon_2, self.lat_2),
                  (self.lon_2, self.lat_1), (self.lon_1, self.lat_1))  # tuple item;
        polygon1 = Polygon(coords)
        self.polygon = gpd.GeoSeries(polygon1)
        return

    def init_ncfile(self):
        # As I am appending to a file, I should check if dimensions or variables already exist
        dimension_key_dict = {'transit_info': 2, 'n_parts': self.shp_p}
        for dimension in dimension_key_dict:
            if dimension not in self.analysis_df.dimensions.keys():
                self.analysis_df.createDimension(dimension, dimension_key_dict[dimension])

        variable_key_dict = {'dom_paths': {'datatype':'i4','dimensions':('n_lon_bins', 'n_lat_bins'), 'description': 'unique particle visits'},
                             'recruit_SG_north': {'datatype':'i4','dimensions':('n_parts', 'transit_info'), 'description': 'transit hours to polygon'}}
        for variable in variable_key_dict:
            if variable not in self.analysis_df.variables.keys():
                self.analysis_df.createVariable(variable, variable_key_dict[variable]['datatype'], variable_key_dict[variable]['dimensions'])
                self.analysis_df[variable].description = variable_key_dict[variable]['description']

        return












