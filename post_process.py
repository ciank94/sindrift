import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


class PostProcess:

    def __init__(self, analysis_file):
        self.analysis_df = nc.Dataset(analysis_file, mode='r+')
        self.trajectory_df = nc.Dataset(self.analysis_df.trajectory_file, mode='r')
        breakpoint()
        # Organise file data
        #self.outfile_path = outfile_path
        #self.key = key
        #self.year = y_i
        #self.release_n = r_i + 1
        #self.nc_file = nc.Dataset(trajectory_file)

        # extract data from trajectory file
        self.lat = self.nc_file['lat']
        self.lon = self.nc_file['lon']
        self.shp_p = np.shape(self.lat)[0]
        self.shp_t = np.shape(self.lat)[1]
        #self.z = self.nc_file['z']

        # Initialize settings for dominant pathway calculation:
        self.init_region(key)

        # output file with data from simulation
        self.analysis_file = (self.outfile_path + self.key + '_' + str(self.year) + '_R'
                         + str(self.release_n) + '_trajectory_analysis.nc')
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

        # master loop for all analysis calculations
        # loop over each particle p_i
        for p_i in range(0, self.shp_p, self.p_stride):
            # general steps for each analysis
            self.print_pstep(descriptor="Particle analysis", p=p_i)  # print particle id number in analysis
            lon_p = self.lon[p_i, :]
            lat_p = self.lat[p_i, :]

            # specific variables stored
            temp_points = np.ones([self.shp_t, 2]) * -1 # dominant paths, use -1 for inactive individuals

            # loop over each time step t_i
            for t_i in range(0, self.shp_t, self.t_stride):
                self.lon_t = lon_p[t_i]
                self.lat_t = lat_p[t_i]
                self.get_closest_point()
                temp_points[t_i, 0] = self.lon_id
                temp_points[t_i, 1] = self.lat_id

            # Then, look at the unique lat-long ids for particle p_i
            unique_rows = np.unique(temp_points, axis=0).astype(int)
            unique_rows = unique_rows[unique_rows[:, 0] > -1]
            id1 = unique_rows[:, 0]
            id2 = unique_rows[:, 1]
            self.write_dom_paths(id1, id2)  # write to file

            #todo: analysis- generic recruitment, generic retention, z dom_paths, transit_times

        print('Closing: ' + self.analysis_file)
        self.outfile.close()
        return


    def init_ncfile(self):
        self.outfile = nc.Dataset(self.analysis_file, 'r+')
        self.outfile.createDimension('time', self.shp_t)
        self.outfile.createDimension('particles', self.shp_p)
        self.outfile.createDimension('lon_bins', self.shp_lon_bins)
        self.outfile.createDimension('lat_bins', self.shp_lat_bins)
        self.outfile.createVariable('dom_paths', 'i4', ('lon_bins', 'lat_bins'))
        self.outfile.variables['dom_paths'][:] = 0
        return

    def write_dom_paths(self, id1, id2):
        self.outfile.variables['dom_paths'][id1, id2] = self.outfile.variables['dom_paths'][id1, id2] + 1
        return

    def diagnostics(self):
        self.init_plot()
        import cartopy.feature as cfeature
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor='face',
                                                facecolor='lightgrey')
        self.ax.add_feature(land_10m)
        self.ax.coastlines(resolution='10m', linewidth=0.7)
        step_p = np.floor(np.shape(self.lon)[0] / 1000).astype(int)
        lon_1 = self.lon[0:-1:step_p, :]
        lat_1 = self.lat[0:-1:step_p, :]
        c_vals = np.arange(0, np.shape(lat_1)[1])
        c_vals = c_vals * np.ones([np.shape(lat_1)[0], np.shape(lat_1)[1]])
        plt.scatter(lon_1, lat_1, c=c_vals, s=1.3, cmap='YlOrRd', alpha=1, linewidth=1.3, linestyle='-', marker='_')
        plt.scatter(lon_1[:, 0], lat_1[:, 0], s=18, facecolor='yellow', edgecolors='k', alpha=0.9, linewidth=0.6)
        plt.scatter(lon_1[:, -1], lat_1[:, -1], s=18, facecolor='red', edgecolors='k', alpha=0.9, linewidth=0.6)
        self.init_region("SG")
        self.ax.set_extent([self.min_lon, self.max_lon, self.min_lat, self.max_lat])
        plt_name = 'SG_diagnose_test'
        self.save_plot(plt_name)

        #breakpoint()

        #self.ax.set_extent([min_lon, max_lon, min_lat, max_lat])

    def init_plot(self):
        import cartopy.crs as ccrs
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(projection=ccrs.PlateCarree())
        return

    def save_plot(self, plt_name):
        savefile = self.outfile_path + plt_name + '.png'
        print('Saving file: ' + savefile)
        plt.savefig(savefile, dpi=400)
        plt.close()
        return


    def dominant_paths(self):
        #self.init_ncfile()
        for p_i in range(0, self.shp_p, self.p_stride):
            lon_p = self.lon[p_i, :]
            lat_p = self.lat[p_i, :]
            #self.write_ncfile(lon_p, lat_p, p_i)
            temp_points = np.ones([self.shp_t, 2])*-1
            self.print_pstep(descriptor="Particle number: ", p=p_i)

            for t in range(0, self.shp_t, self.t_stride):
                self.lon_t = lon_p[t]
                self.lat_t = lat_p[t]
                self.get_closest_point()
                temp_points[t, 0] = self.lon_id
                temp_points[t, 1] = self.lat_id

            # Then, look at the unique lat-long ids for particle p_i
            unique_rows = np.unique(temp_points, axis=0).astype(int)
            unique_rows = unique_rows[unique_rows[:, 0] > -1]
            if np.shape(unique_rows)[0] > 0:
                self.dom_matrix[unique_rows[:, 0], unique_rows[:, 1]] = self.dom_matrix[
                                                                              unique_rows[:, 0], unique_rows[:, 1]] + 1
        #print('Closing transposed file')
        #self.outfile.close()

        print('Saving: ' + self.dom_file)
        np.save(self.dom_file, self.dom_matrix)
        return



    def print_pstep(self, descriptor, p):
        if p == 0:
            print("First step in calculation")
        else:
            print(descriptor + " : " + str(p + 1) + " of " + str(self.shp_p))
        if p == self.shp_p:
            print("Final particle calculation")
        return

    def get_closest_point(self):
        self.lon_id = np.argmin(np.sqrt((self.lon_t - self.lon_range[:]) ** 2))
        self.lat_id = np.argmin(np.sqrt((self.lat_t - self.lat_range[:]) ** 2))
        return





