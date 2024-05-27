import os.path
import matplotlib.pyplot as plt
import matplotlib.collections as mcol
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import netCDF4 as nc


class Analyse:

    def __init__(self, path, key, year):
        # Organise file data
        self.lon_id = None
        self.lat_id = None
        self.lat_t = None
        self.lon_t = None
        self.bin_res = 0.5
        self.filepath = path
        self.results = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/results/'
        self.key = key
        self.tr_file = self.filepath + self.key + '_' + str(year) + '_trajectory.nc'
        self.nc_file = nc.Dataset(self.tr_file)
        self.bath_file = self.results + 'bath.npy'
        self.bath_file_lon = self.results + 'bath_lon.npy'
        self.bath_file_lat = self.results + 'bath_lat.npy'
        self.bath = np.load(self.bath_file)
        self.bath_lon = np.load(self.bath_file_lon)
        self.bath_lat = np.load(self.bath_file_lat)

        # extract data from trajectory file
        self.lat = self.nc_file['lat']
        self.lon = self.nc_file['lon']
        self.shp_p = np.shape(self.lat)[0]
        self.shp_t = np.shape(self.lat)[1]
        #self.z = self.nc_file['z']

        #Initialize matrix:
        self.init_region(key)
        self.domin_matrix = np.zeros([self.shp_lon_bins, self.shp_lat_bins])
        self.d_scale = 40
        self.max_scale = 0.5

        # plotting parameters
        self.bath_contours = np.linspace(0, 3000, 10)
        self.bath_cmap = plt.get_cmap('Blues')
        self.dom_cmap = plt.get_cmap('Reds')
        self.depth_colors = np.arange(0, 4500, 200)
        return

    def dominant_paths(self):

        for p_i in range(0, self.shp_p, 1):
            lon_p = self.lon[p_i, :]
            lat_p = self.lat[p_i, :]
            temp_points = np.ones([self.shp_t, 2])*-1
            self.print_pstep(descriptor="Dominant path calculation", p = p_i)
            if lat_p.mask.shape == ():
                shp_mask = 0
            else:
                shp_mask = 1

            for t in range(0, self.shp_t, 1):
                self.lon_t = lon_p[t]
                self.lat_t = lat_p[t]
                self.get_closest_point()
                if shp_mask == 1:
                    if lat_p.mask[t]:  # check if individual is inactive
                        break
                    else:
                        temp_points[t, 0] = self.lon_id
                        temp_points[t, 1] = self.lat_id

            # Then, look at the unique latlong ids for particle p_i
            unique_rows = np.unique(temp_points, axis=0).astype(int)
            unique_rows = unique_rows[unique_rows[:, 0] > -1]
            if np.shape(unique_rows)[0] > 0:
                self.domin_matrix[unique_rows[:, 0], unique_rows[:, 1]] = self.domin_matrix[
                                                                              unique_rows[:, 0], unique_rows[:, 1]] + 1
        return

    def plot_dom_paths(self):
        self.init_plot()
        self.plot_background()
        n_levels = np.arange(np.min(self.domin_matrix), self.max_scale*np.max(self.domin_matrix),
                             np.max(self.domin_matrix) / self.d_scale)
        self.plot1 = plt.contourf(self.lon_range, self.lat_range, self.domin_matrix.T, levels=n_levels, cmap=self.dom_cmap,
                     transform=ccrs.PlateCarree(), extend='both')
        self.c_max = self.max_scale*np.max(self.domin_matrix)
        self.caxis_title = 'Probability (%)'
        self.add_cbar()
        plt_name = self.key + "_dom_paths"
        self.save_plot(plt_name)
        return

    def get_closest_point(self):
        self.lon_id = np.argmin(np.sqrt((self.lon_t - self.lon_range[:]) ** 2))
        self.lat_id = np.argmin(np.sqrt((self.lat_t - self.lat_range[:]) ** 2))
        return

    def print_pstep(self, descriptor, p):
        if p == 0:
            print("First step in calculation")
        else:
            print(descriptor + " : " + str(p + 1) + " of " + str(self.shp_p))
        if p == self.shp_p:
            print("Final particle calculation")
        return


    def plot_trajectory(self, region):
        self.init_region(region)
        self.init_plot()
        self.plot_background()

        self.plot_depth()
        self.c_max = 4500
        self.caxis_title = 'Depth (m)'
        self.add_cbar()

        step_v = np.floor(np.shape(self.lon)[0]/1000).astype(int)
        step_v2 = np.floor(np.shape(self.lon)[1]/(np.shape(self.lon)[1]*0.99)).astype(int)

        lon_1 = self.lon[0:-1:step_v, 0:-1:step_v2]
        lat_1 = self.lat[0:-1:step_v, 0:-1:step_v2]
        c_vals = np.arange(0,np.shape(lat_1)[1])
        c_vals = c_vals*np.ones([np.shape(lat_1)[0], np.shape(lat_1)[1]])
        #

        plt.scatter(lon_1, lat_1, c=c_vals,s=1.3, cmap='YlOrRd', alpha=1, linewidth=1.3, linestyle='-', marker='_')
        plt.scatter(lon_1[:, 0], lat_1[:, 0], s=18, facecolor='yellow', edgecolors='k', alpha=0.9, linewidth=0.6)
        plt.scatter(lon_1[:, -1], lat_1[:, -1], s=18, facecolor='red', edgecolors='k', alpha=0.9, linewidth=0.6)

        plt_name = self.key + "_worms"
        self.save_plot(plt_name)
        return



    def plot_background(self):
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor='face',
                                                facecolor='lightgrey')
        self.ax.add_feature(land_10m)
        self.ax.coastlines(resolution='10m', linewidth=0.7)
        plt.contour(self.bath_lon, self.bath_lat, self.bath, self.bath_contours, colors='k', alpha=0.2, linewidths=0.7,
                    transform=ccrs.PlateCarree())

        # set extent and grid lines;
        gl = self.ax.gridlines(draw_labels=True, alpha=0.4)
        gl.top_labels = False
        gl.right_labels = False
        self.ax.set_extent([self.min_lon, self.max_lon, self.min_lat, self.max_lat])
        return

    def plot_depth(self):
        #levels= self.depth_colors
        self.plot1 = plt.contourf(self.bath_lon, self.bath_lat, self.bath, levels = self.depth_colors, cmap = self.bath_cmap,
                  transform=ccrs.PlateCarree(), extend='both')
        return


    def init_plot(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(projection=ccrs.PlateCarree())
        return

    def add_cbar(self):
        cbar = plt.colorbar(self.plot1, extend = 'both')
        cbar.ax.set_ylabel(self.caxis_title, loc='center', size=9, weight='bold')
        cbar.ax.tick_params(labelsize=10, rotation=0)
        plt.clim(0, self.c_max)
        return

    def save_plot(self, plt_name):
        savefile = self.results + plt_name + '.png'
        print('Saving file: ' + savefile)
        plt.savefig(savefile, dpi=400)
        plt.close()
        return

    def init_region(self, region):


        if region == "SG":
            self.min_lon = -40.5
            self.max_lon = -33.8
            self.min_lat = -57.5
            self.max_lat = -51.5


        if region == "AP":
            self.min_lon = -65.3
            self.max_lon = -51
            self.min_lat = -69
            self.max_lat = -56

        if region == "SO":
            self.min_lon = -50
            self.max_lon = -41
            self.min_lat = -65
            self.max_lat = -57

        if region == "full":
            self.min_lon = -65
            self.max_lon = -31
            self.min_lat = -70
            self.max_lat = -50

        if region == "APSO":
            self.min_lon = -75
            self.max_lon = -30
            self.min_lat = -70
            self.max_lat = -50

        self.lat_range = np.arange(self.min_lat - 20, self.max_lat + 15, self.bin_res)
        self.lon_range = np.arange(self.min_lon - 20, self.max_lon + 15, self.bin_res)
        self.shp_lon_bins = np.shape(self.lon_range)[0]
        self.shp_lat_bins = np.shape(self.lat_range)[0]

        return

class sinRead:

    def __init__(self, var):
        self.sinmod_path = 'E:/fromIngrid/'
        f_name_start = 'samplesNSEW_2020'
        f_ext = '.nc'
        self.f_name = self.sinmod_path + f_name_start + var + f_ext
        self.nc_file = nc.Dataset(self.f_name)
        #file_save = 'C:/Users/ciank/PycharmProjects/Krill_data/SINdrift/animation.mp4'
        return


