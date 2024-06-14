import os
import numpy as np
import netCDF4 as nc
from netCDF4 import num2date
import geopandas as gpd
from shapely.geometry import Polygon
import xarray as xr
import matplotlib.pyplot as plt


class PostProcess:

    def __init__(self, analysis_file, fpath):
        # files for analysis
        self.fpath = fpath
        self.analysis_file = analysis_file
        self.analysis_df = nc.Dataset(analysis_file, mode='r+')
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

        return

    def trajectory_analysis(self, test):
        # First, check if we are running a test analysis or full analysis:
        if test:
            self.p_stride = 10
            self.t_stride = 10
        else:
            self.p_stride = 1
            self.t_stride = 1

        # Set up a netcdf file for storing output from the analysis:
        self.init_ncfile()
        self.analysis_df.variables['dom_paths'][:] = 0  # Initialize dom_path matrix
        self.analysis_df.variables['recruit'][:] = 0
        self.analysis_df.variables['CG'][:, :] = 0

        id_lat = np.digitize(self.lat, self.lat_bin_vals)
        id_lon = np.digitize(self.lon, self.lon_bin_vals)
        dom_paths = np.zeros([np.shape(self.lon_bin_vals)[0], np.shape(self.lat_bin_vals)[0]])
        self.transit_hours = np.zeros(self.shp_p)
        self.visit_index = np.zeros(self.shp_p)

        for p_i in range(0, self.shp_p, self.p_stride):
            self.p_i = p_i
            print("Particle analysis : " + str(self.p_i + 1) + " of " + str(self.shp_p))  # print particle id number in analysis
            idx_p = id_lat[p_i, :] < np.shape(self.lat_bin_vals)[0]
            dom_paths[id_lon[p_i, idx_p], id_lat[p_i, idx_p]] = dom_paths[id_lon[p_i, idx_p], id_lat[p_i, idx_p]] + 1
            if p_i == 0:
                self.init_square_polygon()
            self.recruit(self.lon[p_i, :], self.lat[p_i, :])

        self.analysis_df.variables['dom_paths'][:] = dom_paths
        self.analysis_df.variables['recruit'][:, 0, 0] = self.transit_hours
        self.analysis_df.variables['recruit'][:, 1, 0] = self.visit_index
        #todo: analysis- generic retention, z dom_paths, transit_times

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

    def init_square_polygon(self):
        import matplotlib.pyplot as plt
        # lon_1, lat_1 are bottom left coordinates
        self.lon_1 = self.analysis_df.variables['square_polygons'][0, 0]
        self.lon_2 = self.analysis_df.variables['square_polygons'][0, 1]
        self.lat_1 = self.analysis_df.variables['square_polygons'][0, 2]
        self.lat_2 = self.analysis_df.variables['square_polygons'][0, 3]
        coords = ((self.lon_1, self.lat_1), (self.lon_1, self.lat_2), (self.lon_2, self.lat_2),
                  (self.lon_2, self.lat_1), (self.lon_1, self.lat_1))  # tuple item;
        polygon1 = Polygon(coords)
        self.polygon = gpd.GeoSeries(polygon1)
        return

    def unique_visits(self, dom_points):
        unique_rows = np.unique(dom_points, axis=0).astype(int)
        unique_rows = unique_rows[unique_rows[:, 0] > -1]
        id1 = unique_rows[:, 0]
        id2 = unique_rows[:, 1]

        return


    def init_ncfile(self):
        # As I am appending to a file, I should check if dimensions or variables already exist
        try:
            self.analysis_df.dimensions['time']
        except:
            self.analysis_df.createDimension('time', self.shp_t)
        try:
            self.analysis_df.dimensions['transit_info']
        except:
            self.analysis_df.createDimension('transit_info', 2) # transit time and index of visit
        try:
            self.analysis_df.dimensions['particles']
        except:
            self.analysis_df.createDimension('particles', self.shp_p)
        try:
            self.analysis_df.dimensions['lon_lat']
        except:
            self.analysis_df.createDimension('lon_lat', 2)
        try:
            self.analysis_df.variables['dom_paths']
        except:
            self.analysis_df.createVariable('dom_paths', 'i4', ('n_lon_bins', 'n_lat_bins'))
        try:
            self.analysis_df.variables['recruit']
        except:
            self.analysis_df.createVariable('recruit', 'i4', ('particles','transit_info', 'n_polygons'))
            self.analysis_df['recruit'].description = 'Transit hours to polygon'
        try:
            self.analysis_df.variables['CG']
        except:
            self.analysis_df.createVariable('CG', 'f4', ('time', 'lon_lat'))
        return


    def print_pstep(self, descriptor):
        if self.p_i == 0:
            print("First step in calculation")
        else:
            print(descriptor + " : " + str(self.p_i + 1) + " of " + str(self.shp_p))
        return







