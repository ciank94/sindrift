import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import netCDF4 as nc

class PlotData:

    def __init__(self, fuse):
        # Settings for getting bathymetry files
        # todo: document how to get the bathymetry files as below
        # todo: look into saving information from simulation in post_process as class/ dictionary
        self.bath_file = fuse.outfile_path + 'bath.npy'
        self.bath_file_lon = fuse.outfile_path + 'bath_lon.npy'
        self.bath_file_lat = fuse.outfile_path + 'bath_lat.npy'
        self.bath = np.load(self.bath_file)
        self.bath_lon = np.load(self.bath_file_lon)
        self.bath_lat = np.load(self.bath_file_lat)
        self.min_lon = -75
        self.max_lon = -30
        self.min_lat = -70
        self.max_lat = -50

        #SG
        self.min_lon = -40.5
        self.max_lon = -33.8
        self.min_lat = -57.5
        self.max_lat = -51.5

        self.bin_res = 0.02
        self.lat_range = np.arange(self.min_lat - 20, self.max_lat + 15, self.bin_res)
        self.lon_range = np.arange(self.min_lon - 20, self.max_lon + 15, self.bin_res)
        self.results = fuse.outfile_path




        self.save_prefix = fuse.key_id + '_' + str(fuse.year_id)

        #self.min_lon = -50
        #self.max_lon = -41
        #self.min_lat = -65
        #self.max_lat = -57

        #self.min_lon = -65.3
        #self.max_lon = -51
        #self.min_lat = -69
        #self.max_lat = -56

        # Dominant pathways color scaling:
        self.d_scale = 40
        self.max_scale = 0.8

        # plotting parameters
        self.bath_contours = np.linspace(0, 3000, 10)
        self.bath_cmap = plt.get_cmap('Blues')
        #self.dom_cmap = plt.get_cmap('Reds')
        self.dom_cmap = plt.get_cmap('OrRd')
        self.depth_colors = np.arange(0, 4500, 200)

    def plot_dom_paths(self, fuse):
        self.init_plot()
        self.plot_background()
        n_levels = np.arange(np.min(fuse.dom_path_i), self.max_scale * np.max(fuse.dom_path_i),
                            np.max(fuse.dom_path_i)*self.max_scale / self.d_scale)
        #breakpoint()
        self.plot1 = plt.contourf(self.lon_range, self.lat_range, fuse.dom_path_i.T, levels=n_levels,
                                cmap=self.dom_cmap,
                                 transform=ccrs.PlateCarree(), extend='both')
        #self.plot1 = plt.pcolormesh(self.lon_range, self.lat_range, fuse.dom_path_i.T,
         #                         cmap=self.dom_cmap,
          #                        transform=ccrs.PlateCarree())
        self.c_max = self.max_scale * np.max(fuse.dom_path_i)
        self.caxis_title = 'Probability (%)'
        self.add_cbar()
        plt_name = self.save_prefix + "_dom_paths"
        self.save_plot(plt_name)
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

    def plot_trajectory(self, region):
        self.init_region(region)
        self.init_plot()
        self.plot_background()

        self.plot_depth()
        self.c_max = 4500
        self.caxis_title = 'Depth (m)'
        self.add_cbar()

        step_v = np.floor(np.shape(self.lon)[0] / 1000).astype(int)
        step_v2 = np.floor(np.shape(self.lon)[1] / (np.shape(self.lon)[1] * 0.99)).astype(int)

        lon_1 = self.lon[0:-1:step_v, 0:-1:step_v2]
        lat_1 = self.lat[0:-1:step_v, 0:-1:step_v2]
        c_vals = np.arange(0, np.shape(lat_1)[1])
        c_vals = c_vals * np.ones([np.shape(lat_1)[0], np.shape(lat_1)[1]])
        #

        plt.scatter(lon_1, lat_1, c=c_vals, s=1.3, cmap='YlOrRd', alpha=1, linewidth=1.3, linestyle='-', marker='_')
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
        # levels= self.depth_colors
        self.plot1 = plt.contourf(self.bath_lon, self.bath_lat, self.bath, levels=self.depth_colors,
                                  cmap=self.bath_cmap,
                                  transform=ccrs.PlateCarree(), extend='both')
        return

    def init_plot(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(projection=ccrs.PlateCarree())
        return

    def add_cbar(self):
        cbar = plt.colorbar(self.plot1, extend='both')
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




class FuseData:
    def __init__(self, infile_path, outfile_path, year_ids, release_ids, key_ids):
        self.dom_path_i = None
        self.outfile_path = outfile_path
        self.infile_path = infile_path
        self.year_id = year_ids
        self.release_ids = release_ids
        self.key_id = key_ids
        self.n_parts = 10000

        # Fusion and statistics of datasets for plotting:
        self.fuse_dominant_paths()

    def fuse_dominant_paths(self):
        c_i = 0
        year = self.year_id
        key = self.key_id
        for release in self.release_ids:
            c_i = c_i + 1
            file = self.infile_path + key + '_' + str(year) + '_R' + str(release) + '_trajectory_analysis.nc'
            f1 = nc.Dataset(file)
            if c_i == 1:
                self.dom_path_i = f1['dom_paths'][:]
                f1.close()
            else:
                self.dom_path_i = self.dom_path_i + f1['dom_paths'][:]
                f1.close()

        #self.dom_path_i = self.dom_path_i/(c_i*self.n_parts)
        #breakpoint()
        self.dom_path_i = self.dom_path_i / np.nanmax(self.dom_path_i)
        if self.dom_path_i[0,0] > 3*np.mean(np.unique(self.dom_path_i)):
            self.dom_path_i[0, 0] = 0

        self.dom_path_i = self.dom_path_i / np.nanmax(self.dom_path_i)


