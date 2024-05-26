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
        self.bin_res = 0.6
        self.filepath = path
        self.results = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/results/'
        self.key = key
        self.tr_file = self.filepath + self.key + '_' + str(year) +  '_trajectory.nc'
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
        #self.z = self.nc_file['z']

        # plotting parameters
        self.bath_contours = np.linspace(0, 3000, 10)
        self.bath_cmap = plt.get_cmap('Blues')
        self.depth_colors = np.arange(0, 4500, 200)
        return

    def dominant_paths(self):
        breakpoint()
        lat = self.lat[0, :]
        lon = self.lon[0, :]

        arr = [lon[:], lat[:]]

        unique_rows = np.unique(arr, axis=0)

        lon_mesh, lat_mesh = np.meshgrid(self.lon_range, self.lat_range)

        shp_lat = np.shape(self.lat)[0]
        shp_lon_range = np.shape(self.lon_range)[0]
        shp_lat_range = np.shape(self.lat_range)[0]
        dens_m = np.zeros([shp_lon_range, shp_lat_range])
        n_catches = np.zeros([shp_lon_range, shp_lat_range])
        dens_f = np.zeros([shp_lon_range, shp_lat_range])
        # best option: loop over first individual and save latlong bin values the individual has visited, then I can
        # use unique rows to find the unique bin values in each case, start with low bin resolution for clarity
        for ij in range(0, shp_lat):
            lat_id = np.argmin(np.sqrt((lat.iloc[ij] - self.lat_range[:]) ** 2))
            lon_id = np.argmin(np.sqrt((lon.iloc[ij] - self.lon_range[:]) ** 2))
            dens_m[lon_id, lat_id] = dens_m[lon_id, lat_id] + 1
            n_catches[lon_id, lat_id] = n_catches[lon_id, lat_id] + 1

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

        self.lat_range = np.arange(self.min_lat - 10, self.max_lat + 6, self.bin_res)
        self.lon_range = np.arange(self.min_lon - 10, self.max_lon + 6, self.bin_res)

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


