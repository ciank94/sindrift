import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


class PostProcess:

    def __init__(self, analysis_file):
        # files for analysis
        self.analysis_file = analysis_file
        self.analysis_df = nc.Dataset(analysis_file, mode='r+')
        self.trajectory_df = nc.Dataset(self.analysis_df.trajectory_file, mode='r')

        # extract data from trajectory file
        self.lat = self.trajectory_df['lat']
        self.lon = self.trajectory_df['lon']
        self.shp_p = np.shape(self.lat)[0]
        self.shp_t = np.shape(self.lat)[1]
        # self.z = self.nc_file['z']

        #add from analysis file
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

        # master loop for all analysis calculations
        for p_i in range(0, self.shp_p, self.p_stride):
            # general steps for each analysis
            self.print_pstep(descriptor="Particle analysis", p=p_i)  # print particle id number in analysis
            lon_p = self.lon[p_i, :]
            lat_p = self.lat[p_i, :]

            # specific variables stored
            dom_points = np.ones([self.shp_t, 2]) * -1  # dominant paths, use -1 for inactive individuals

            # loop over each time step t_i
            for t_i in range(0, self.shp_t, self.t_stride):
                self.lon_t = lon_p[t_i]
                self.lat_t = lat_p[t_i]
                self.get_closest_point()
                # assign to the closest bin
                dom_points[t_i, 0] = self.lon_id
                dom_points[t_i, 1] = self.lat_id

            # Then, look at the unique lat-long ids for particle p_i
            self.unique_visits(dom_points)


            #todo: analysis- generic recruitment, generic retention, z dom_paths, transit_times

        print('Closing: ' + self.analysis_file)
        self.analysis_df.close()
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
            self.analysis_df.dimensions['particles']
        except:
            self.analysis_df.createDimension('particles', self.shp_p)
        try:
            self.analysis_df.variables['dom_paths']
        except:
            self.analysis_df.createVariable('dom_paths', 'i4', ('n_lon_bins', 'n_lat_bins'))
        return

    def get_closest_point(self):
        self.lon_id = np.argmin(np.sqrt((self.lon_t - self.lon_bin_vals[:]) ** 2))
        self.lat_id = np.argmin(np.sqrt((self.lat_t - self.lat_bin_vals[:]) ** 2))
        return

    def print_pstep(self, descriptor, p):
        if p == 0:
            print("First step in calculation")
        else:
            print(descriptor + " : " + str(p + 1) + " of " + str(self.shp_p))
        if p == self.shp_p:
            print("Final particle calculation")
        return





