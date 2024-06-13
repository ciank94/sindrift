import os

import numpy as np
import netCDF4 as nc
from netCDF4 import num2date
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
            self.t_stride = 6
        else:
            self.p_stride = 1
            self.t_stride = 6

        # Set up a netcdf file for storing output from the analysis:
        self.init_ncfile()
        self.analysis_df.variables['dom_paths'][:] = 0  # Initialize dom_path matrix
        self.analysis_df.variables['recruit'][:] = 0
        self.analysis_df.variables['CG'][:, :] = 0

        #for t_i in range(0, self.shp_t, 1):
         #   print('Time analysis : ' + str(t_i + 1) + " of " + str(self.shp_t))
          #  self.analysis_df.variables['CG'][t_i, 0] = np.nanmedian(self.lon[:, t_i])
           # self.analysis_df.variables['CG'][t_i, 1] = np.nanmedian(self.lat[:, t_i])

        # master loop for all analysis calculations
        for p_i in range(0, self.shp_p, self.p_stride):
            # general steps for each analysis
            self.p_i = p_i
            print("Particle analysis : " + str(self.p_i + 1) + " of " + str(self.shp_p))  # print particle id number in analysis
            if p_i == 0:
                if self.fpath.node == 'server':
                    cmd = "echo first particle in analysis"
                    os.system(cmd)

            self.lon_p = self.lon[p_i, :]
            self.lat_p = self.lat[p_i, :]

            # specific variables stored
            dom_points = np.ones([self.shp_t, 2]) * -1  # dominant paths, use -1 for inactive individuals

            # loop over each time step t_i
            for t_i in range(0, self.shp_t, self.t_stride):
                self.lon_t = self.lon_p[t_i]
                self.lat_t = self.lat_p[t_i]
                self.get_closest_point()
                # assign to the closest bin
                dom_points[t_i, 0] = self.lon_id
                dom_points[t_i, 1] = self.lat_id

            # Then, look at the unique lat-long ids for particle p_i
            self.unique_visits(dom_points)

            if p_i == 0:
                self.init_square_polygon()
            self.recruit()

        #todo: analysis- generic retention, z dom_paths, transit_times

        print('Closing: ' + self.analysis_file)
        os.system("echo closing analysis file")
        self.summarise_file()
        self.analysis_df.close()
        return

    def recruit(self):
        import geopandas as gpd
        #crs = "EPSG:4326"
        gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries.from_xy(self.lon_p, self.lat_p))  # Create Geo dataframe from positional data
        polya = gpd.GeoDataFrame(geometry=self.polygon)  # Geo dataframe from polygon data
        #gdf.crs = polya.crs  # Make sure we use the same projections for both
        in_pol = gpd.tools.sjoin(gdf, polya, predicate="within", how='left')  # Check points are in polygon , predicate="within", how='left'
        in_idx = np.where(in_pol.index_right > -1)  # Find all values in the polygon
        p_id = np.zeros(np.shape(self.lon_p)[0])  # Initialize vector for storage
        p_id[in_idx[0]] = 1  # Store as ones
        if np.sum(p_id) > 0:
            visits = np.where(p_id == 1)
            visit_index = visits[0][0]
            t_to_poly = self.dates[visits[0][0]] - self.dates[0]
            transit_hours = (t_to_poly.days * 24) + np.floor(t_to_poly.seconds * 1 / (60 * 60)).astype(int)
            self.analysis_df.variables['recruit'][self.p_i, 0, 0] = transit_hours
            self.analysis_df.variables['recruit'][self.p_i, 1, 0] = visit_index
        return

    def init_square_polygon(self):
        from shapely.geometry import Polygon
        import geopandas as gpd
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
        self.analysis_df.variables['dom_paths'][id1, id2] = self.analysis_df.variables['dom_paths'][
                                                                id1, id2] + 1  # write to file
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

    def get_closest_point(self):
        self.lon_id = np.argmin(np.sqrt((self.lon_t - self.lon_bin_vals[:]) ** 2))
        self.lat_id = np.argmin(np.sqrt((self.lat_t - self.lat_bin_vals[:]) ** 2))
        return

    def print_pstep(self, descriptor):
        if self.p_i == 0:
            print("First step in calculation")
        else:
            print(descriptor + " : " + str(self.p_i + 1) + " of " + str(self.shp_p))
        return

    def summarise_file(self):
        savefile = self.fpath.trajectory_path + self.analysis_df.trajectory_file_prefix + 'summary' + '.txt'
        with open(savefile, "w") as text_file:
            print('Writing data to file: ' + savefile)
            text_file.writelines(str(self.analysis_df))







