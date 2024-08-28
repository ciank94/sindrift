import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
#from matplotlib.patches import Polygon
from shapely.geometry.polygon import LinearRing, Polygon
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import netCDF4 as nc
import os
from datetime import datetime
from netCDF4 import num2date
class PlotData:

    def __init__(self, key, year, compile_folder, analysis_folder):
        self.figures_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/figures/'
        self.key = key
        self.year_n = year
        self.year = str(year)
        self.compile_folder = compile_folder
        self.analysis_folder = analysis_folder
        self.file_prefix = key + '_' + self.year + '_'
        self.save_prefix = self.figures_path + self.file_prefix
        self.analysis_file = self.analysis_folder + self.file_prefix + 'R1_trajectory_analysis.nc'
        self.df = nc.Dataset(self.analysis_file)

        # polygon for analysis
        self.n_poly = 0
        self.poly = self.get_polygon(self.n_poly)

        # For creating the bathymetry file
        self.load_bathymetry()

        # parameters for plotting:
        self.load_plot_params()
        #self.lon_lat_extent()
        return

    def read_catch(self):
        catch_file = self.figures_path + 'C1_680.csv'
        self.time_file = self.figures_path + 'time.npy'

        csv_file = pd.read_csv(catch_file, sep=',')
        self.area_rule = csv_file.asd_code == 483  # key code for area
        self.get_time(csv_file)
        self.df = csv_file[self.area_rule]

        id_c = self.year_c == self.year_n + 1
        avg_catch = np.nanmean(self.df.krill_greenweight_kg[id_c]/1000)
        return avg_catch



    def get_time(self, csv_file):
        if not os.path.exists(self.time_file):  # if file doesn't exist
            time_array = csv_file.datetime_haul_start
            shp_t = np.shape(time_array)[0]
            time_store = np.zeros([shp_t, 4])

            for i in range(0, shp_t):
                d1 = datetime.strptime(time_array[i], "%Y-%m-%d %H:%M:%S")
                time_store[i, 0] = d1.hour
                time_store[i, 1] = d1.day
                time_store[i, 2] = d1.month
                time_store[i, 3] = d1.year

            np.save(self.time_file, time_store)
            time_array = np.load(self.time_file)
            print('Saving file: ' + self.time_file)
        else:
            time_array = np.load(self.time_file)
            print('File exists: ' + self.time_file)
        time_sub = time_array[self.area_rule]
        self.hour = time_sub[:, 0].astype(int)
        self.day = time_sub[:, 1].astype(int)
        self.month = time_sub[:, 2].astype(int)
        self.year_c = time_sub[:, 3].astype(int)
        return

    def get_recruits(self):
        filename = self.compile_folder + self.file_prefix + 'recruit_SG.csv'
        r_table = pd.read_csv(filename)
        recruit_mean = np.mean(r_table.recruit_number * 100 / 10000)
        return recruit_mean

    def get_recruit_times(self):
        filename = self.compile_folder + self.file_prefix + 'recruit_SG.csv'
        r_table = pd.read_csv(filename)
        recruit_time = np.mean(r_table.recruit_time)/24
        return recruit_time

    def get_temp_exp(self):
        filename = self.compile_folder + self.file_prefix + 'temp_exp.npy'
        temp_exp = np.load(filename)
        temp_mean = np.nanmean(temp_exp)
        return temp_mean

    def get_chl_exp(self):
        filename = self.compile_folder + self.file_prefix + 'chl_exp.npy'
        chl_exp = np.load(filename)
        chl_mean = np.nanmean(chl_exp)
        return chl_mean

    def get_o2_exp(self):
        filename = self.compile_folder + self.file_prefix + 'o2_exp.npy'
        o2_exp = np.load(filename)
        o2_mean = np.nanmean(o2_exp)
        return o2_mean


    def load_plot_params(self):

        # Dominant pathways color scaling:
        self.d_scale = 40
        if self.key == 'SOIN':
            self.max_val = 5

        # plotting parameters
        self.bath_contours = np.linspace(0, 3000, 10)
        self.bath_cmap = plt.get_cmap('Blues')
        #self.dom_cmap = plt.get_cmap('Reds')
        self.dom_cmap = plt.get_cmap('OrRd')
        self.depth_colors = np.arange(0, 4500, 200)
        self.site_recruit_cmap = plt.get_cmap('hot')

        # offsets for figure boundaries:
        self.lon_offset = 1
        self.lat_offset = 3
        return

    def plot_dom_paths(self, release_n):
        self.init_plot()
        self.plot_background(background=self.key)
        filename = self.compile_folder + self.file_prefix + 'dom_paths.npy'
        dom_paths = np.load(filename)
        dom_paths = dom_paths.astype(float)
        #dom_paths[dom_paths == 0] = np.nan
        dom_paths = (dom_paths / ((release_n) * 10000)) * 100
        n_levels = np.arange(np.nanmin(dom_paths), self.max_val, self.max_val / 25)
        #self.plot1 = plt.contourf(self.df['lon_bin_vals'][:], self.df['lat_bin_vals'][:], dom_paths.T, levels=n_levels, cmap=self.dom_cmap,
                                  #transform=ccrs.PlateCarree(), extend='both')
        self.plot1 = plt.pcolormesh(self.df['lon_bin_vals'][:], self.df['lat_bin_vals'][:], dom_paths.T,
                                  cmap=self.dom_cmap,
                                  transform=ccrs.PlateCarree())
        self.add_cbar(c_max=self.max_val, caxis_title='probability (%)')
        plt_name = self.file_prefix + "dom_paths"
        self.save_plot(plt_name)
        return

    def plot_site_recruit_t(self, release_n):
        self.init_plot()
        self.plot_background(background='SOI')
        filename = self.compile_folder + self.file_prefix + 'site_recruits.npy'
        site_recruits = np.load(filename)
        filename2 = self.compile_folder + self.file_prefix + 'recruit_SG.csv'
        r_table = pd.read_csv(filename2)
        print('recruit number mean ' + self.year + ': ' + str(np.mean(r_table.recruit_number*100/10000)))
        print('recruit number std ' + self.year + ': ' + str(np.std(r_table.recruit_number*100/10000)))
        print('recruit time mean ' + self.year + ': ' + str(np.mean(r_table.recruit_time/24)))
        print('recruit time std ' + self.year + ': ' + str(np.std(r_table.recruit_time / 24)))
        self.ax.scatter(self.df.variables['lon_init'][:], self.df.variables['lat_init'][:], s=10, facecolor='none', edgecolors='gray',
                        alpha=0.3, linewidth=0.2)
        c_vals = site_recruits[2, :] / (release_n) * 100
        #print('recruit mean values for ' + self.year + ': ' + str(np.mean(c_vals)))
        #print('recruit std values for ' + self.year + ': ' + str(np.std(c_vals)))
        self.plot1 = plt.scatter(self.df.variables['lon_init'][:], self.df.variables['lat_init'][:], c=c_vals, s=10, edgecolors='gray',
                                 vmin=np.nanmean(c_vals)/2,
                                 vmax=100, linewidth=0.2, cmap=self.site_recruit_cmap)
        self.add_cbar(c_max=60, caxis_title='Percentage (%)')
        self.save_plot(plt_name=self.file_prefix + 'site_recruits')
        return

    def plot_hist_environment(self):
        filename = self.compile_folder + self.file_prefix + 'chl_exp.npy'
        chl_exp = np.load(filename)
        print('chlorophyll mean: ' + self.year + ': ' + str(np.mean(chl_exp)))
        print('chlorophyll std: ' + self.year + ': ' + str(np.std(chl_exp)))
        plt.plot(chl_exp, c='k')
        plt.ylabel('mg m-3')
        plt.ylim([0.03, 0.26])
        self.save_plot(plt_name=self.file_prefix + 'chl_exp')

        filename = self.compile_folder + self.file_prefix + 'o2_exp.npy'
        o2_exp = np.load(filename)
        print('o2 mean: ' + self.year + ': ' + str(np.mean(o2_exp)))
        print('o2 std: ' + self.year + ': ' + str(np.std(o2_exp)))
        plt.plot(o2_exp, c='k')
        plt.ylabel('mmol m-3')
        plt.ylim([320, 365])
        self.save_plot(plt_name=self.file_prefix + 'o2_exp')

        filename = self.compile_folder + self.file_prefix + 'temp_exp.npy'
        temp_exp = np.load(filename)
        print('temp mean: ' + self.year + ': ' + str(np.nanmean(temp_exp)))
        print('temp std: ' + self.year + ': ' + str(np.nanstd(temp_exp)))
        plt.plot(temp_exp, c='k')
        plt.ylabel('C')
        plt.ylim([-2, 2])
        self.save_plot(plt_name=self.file_prefix + 'temp_exp')

        filename = self.compile_folder + self.file_prefix + 'CG_lat.npy'
        CG_lat = np.load(filename)
        plt.plot(CG_lat, c='k')
        plt.ylabel('latitude')
        plt.ylim([-60.4, -57])
        self.save_plot(plt_name=self.file_prefix + 'CG_lat')
        return



    def plot_CG_paths(self, df):
        self.init_plot()
        self.plot_background()
        [self.ax.plot(df.CG[:, 0, i], df.CG[:, 1, i], alpha=0.9, linewidth=0.5, c='k') for i in range(0, np.shape(df.CG)[2])]
        plt_name = self.save_prefix + "CG_plot"
        self.save_plot(plt_name)
        return

    def plot_site_recruits(self, df):
        self.init_plot()
        self.plot_background(background="SOI")
        c_vals = df.site_recruits[:]/(df.counter_r + 1) * 100
        self.ax.scatter(df.lon_init, df.lat_init, s=10, facecolor='none', edgecolors='gray', alpha=0.3, linewidth=0.2)
        self.plot1 = plt.scatter(df.lon_init, df.lat_init, c=c_vals,s=10, edgecolors='gray', vmin=np.nanmean(c_vals)/2,
                               vmax=np.nanmean(c_vals)*2, linewidth=0.2, cmap=self.site_recruit_cmap)
        self.add_cbar(c_max=np.nanmax(c_vals), caxis_title='Percentage (%)')
        self.save_plot(plt_name=self.save_prefix + 'site_recruits')
        return






    def plot_init(self):
        self.init_plot()
        self.plot_background()
        plt.scatter(self.lon[:,0], self.lat[:,0])
        plt_name = self.save_prefix + "init_plot"
        self.save_plot(plt_name)
        return


    def load_bathymetry(self):
        self.bath_res = 0.04  # bathymetry resolution
        self.bath_file = self.figures_path + 'bath.npy'
        self.bath_file_lon = self.figures_path + 'bath_lon.npy'
        self.bath_file_lat = self.figures_path + 'bath_lat.npy'
        if not os.path.exists(self.bath_file):
            self.get_bathfile()
            print('Creating ' + self.bath_file)
        self.bath = np.load(self.bath_file)
        self.bath_lon = np.load(self.bath_file_lon)
        self.bath_lat = np.load(self.bath_file_lat)
        return

    def get_bathfile(self):
        bathfile_gebco = self.fpath.figures_path + 'gebco_2023.nc'

        # take in the gebco bathymetry file
        bath_f = nc.Dataset(bathfile_gebco)
        e = np.array(bath_f['elevation'][:]).astype(float)
        lat_e = np.array(bath_f['lat'][:])
        lon_e = np.array(bath_f['lon'][:])

        e[e > 0] = np.nan
        e *= -1

        new_lat = np.arange(np.min(lat_e), np.max(lat_e), self.bath_res)
        new_lon = np.arange(np.min(lon_e), np.max(lon_e), self.bath_res)

        shp_lat = np.shape(new_lat)[0]
        shp_lon = np.shape(new_lon)[0]
        e_new = np.zeros([shp_lat, shp_lon])
        # todo: use digitize instead of argmin;
        for i in range(0, shp_lat):
            for j in range(0, shp_lon):
                lat_id = np.argmin(np.sqrt((lat_e[:] - new_lat[i]) ** 2))
                lon_id = np.argmin(np.sqrt((lon_e[:] - new_lon[j]) ** 2))
                e_new[i, j] = e[lat_id, lon_id]

        np.save(self.bath_file, e_new)
        np.save(self.bath_file_lat, new_lat)
        np.save(self.bath_file_lon, new_lon)
        bath_f.close()
        return



    def plot_time_recruits(self, df):
        dates = df.recruit_table.date
        recruit_n = df.recruit_table.recruit_number
        recruit_t = df.recruit_table.recruit_time/24

        fig, ax1 = plt.subplots()
        plt.xticks(rotation=45, fontsize=10)
        ax2 = ax1.twinx()
        ax1.set_ylabel('recruit time (days)', color='b', fontsize=11)
        ax2.set_ylabel('recruit number', color='r', fontsize=11)
        ax1.plot(dates, recruit_t, c = 'b')
        ax2.plot(dates, recruit_n, c = 'r')
        #ax1.set_xlabel('date', fontsize=8)
        plt.grid(alpha=0.45)  # nice and clean grid
        plt_name = self.save_prefix + "recruit_times"
        self.save_plot(plt_name)



    def plot_recruits(self):
        recruit = self.recruit[:,:]
        recruit_t = recruit[:,0]
        recruit_i = recruit[:,1]
        id1 = np.where(recruit_t>0)
        self.init_plot()
        self.plot_background()
        i_vals = recruit_i[id1[0][:]].astype(int)
        particle_n = id1[0][:].astype(int)
        x = np.array(self.lon[particle_n, :])
        y = np.array(self.lat[particle_n, :])
        [self.ax.plot(x[i, 0:i_vals[i]], y[i, 0:i_vals[i]], alpha=0.1, linewidth=0.5, c='r') for i in range(0, np.shape(id1)[1])]
        self.ax.scatter(x[:, 0], y[:, 0], alpha=0.7, c='r', s=1)
        x,y = self.poly.exterior.xy
        self.ax.plot(x, y, color='y', alpha=0.9,
                linewidth=2, solid_capstyle='round', zorder=2)
        plt_name = self.save_prefix + "recruit"
        self.save_plot(plt_name)

        n_recruits = np.shape(id1)[1]


        #self.ax.plot(polygon1)
        #self.ax.add_patch(polygon1, facecolor='r', alpha=0.4)
        #
        #breakpoint()
        return n_recruits



    def get_polygon(self, poly_n):
        lon_lims = [-39.5, -35]
        lat_lims = [-54, -53]
        self.lon_1 = lon_lims[0]
        self.lon_2 = lon_lims[1]
        self.lat_1 = lat_lims[0]
        self.lat_2 = lat_lims[1]
        coords = ((self.lon_1, self.lat_1), (self.lon_1, self.lat_2), (self.lon_2, self.lat_2),
                  (self.lon_2, self.lat_1), (self.lon_1, self.lat_1))  # tuple item;
        poly = Polygon(coords)
        return poly

    def plot_trajectory(self):
        self.init_plot()
        self.plot_background()
        self.plot_depth()
        self.add_cbar(caxis_title='depth (m)', c_max=4500)

        step_v = np.floor(np.shape(self.lon)[0] / 1000).astype(int)
        step_v2 = np.floor(np.shape(self.lon)[1] / (np.shape(self.lon)[1] * 0.99)).astype(int)

        lon_1 = self.lon[0:-1:step_v, 0:-1:step_v2]
        lat_1 = self.lat[0:-1:step_v, 0:-1:step_v2]
        c_vals = np.arange(0, np.shape(lat_1)[1])
        c_vals = c_vals * np.ones([np.shape(lat_1)[0], np.shape(lat_1)[1]])
        #

        plt.scatter(lon_1, lat_1, c=c_vals, s=1.3, cmap='YlOrRd', alpha=1, linewidth=1.3, linestyle='-', marker='_')
        plt.scatter(lon_1[:, -1], lat_1[:, -1], s=18, facecolor='red', edgecolors='k', alpha=0.9, linewidth=0.6)
        plt.scatter(lon_1[:, 0], lat_1[:, 0], s=5, facecolor='yellow', edgecolors='g', alpha=0.3, linewidth=0.5)

        plt_name = self.save_prefix + "worms"
        self.save_plot(plt_name)
        return

    def plot_background(self, background, ax_name):
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor='face',
                                                facecolor='lightgrey')
        ax_name.add_feature(land_10m)
        ax_name.coastlines(resolution='10m', linewidth=0.7)
        ax_name.contour(self.bath_lon, self.bath_lat, self.bath, self.bath_contours, colors='k', alpha=0.2, linewidths=0.7,
                    transform=ccrs.PlateCarree())

        # set extent and grid lines;
        gl = ax_name.gridlines(draw_labels=True, alpha=0.4)
        gl.top_labels = False
        gl.right_labels = False
        #self.ax.set_extent(
            #[self.min_lon, self.max_lon, self.min_lat, self.max_lat])
        if background == "BSSI":
            #self.AP_lon_lat_extent()
            self.gen_lon_lat_extent()
            ax_name.set_extent(
                [self.min_lon, self.max_lon, self.min_lat,
                 self.max_lat])
            #self.ax.set_extent(
             #   [self.min_lon - self.lon_offset, self.max_lon + self.lon_offset, self.min_lat - self.lat_offset,
              #   self.max_lat + self.lat_offset])
        elif background == "SOI":
            self.SOI_lon_lat_extent()
            ax_name.set_extent(
                [self.min_lon, self.max_lon, self.min_lat,
                 self.max_lat])
        elif background == "SOIN":
            self.SOIN_lon_lat_extent()
            ax_name.set_extent(
                [self.min_lon, self.max_lon, self.min_lat,
                 self.max_lat])
        elif background == 'SG_ex':
            lon_lim = [-39.4, -35]
            lat_lim = [-55, -53.2]
            ax_name.set_extent(
                [lon_lim[0], lon_lim[1], lat_lim[0],
                 lat_lim[1]])
        else:
            self.gen_lon_lat_extent()
            ax_name.set_extent([self.min_lon, self.max_lon, self.min_lat,
                            self.max_lat])
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

    def add_cbar(self, c_max, caxis_title, axis_name):
        cbar = axis_name.colorbar(self.plot1, extend='both', pad=0.01)
        cbar.ax.set_ylabel(caxis_title, loc='center', size=9, weight='bold')
        cbar.ax.tick_params(labelsize=10, rotation=0)
        axis_name.clim(0, c_max)
        return

    def save_plot(self, plt_name):
        savefile = self.figures_path + plt_name + '.png'
        print('Saving file: ' + savefile)
        plt.savefig(savefile, dpi=400)
        plt.close()
        return

    def AP_lon_lat_extent(self):
        self.min_lon = -64
        self.max_lon = -34
        self.min_lat = -70
        self.max_lat = -50

        return

    def SOI_lon_lat_extent(self):
        #todo: put all site_lims in same function and add ifelse statements
        self.min_lon = -49
        self.max_lon = -41
        self.min_lat = -64.5
        self.max_lat = -58
        return

    def gen_lon_lat_extent(self):
        # SG extent;
        self.min_lon = -45
        self.max_lon = -32
        self.min_lat = -61
        self.max_lat = -50
        return

    def SOIN_lon_lat_extent(self):
        # SG extent;
        self.min_lon = -50
        self.max_lon = -31
        self.min_lat = -65
        self.max_lat = -50
        return


class CatchData:
    def __init__(self):
        self.figures_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/figures/'
        catch_file = self.figures_path + 'C1_680.csv'
        self.time_file = self.figures_path + 'time.npy'
        self.csv_file = pd.read_csv(catch_file, sep=',')
        return

    def catch_facts(self):
        self.get_area(482)
        np.nanmean(self.df.krill_greenweight_kg)/1000

        np.sum((self.month == 12) | (self.month == 1)| (self.month == 2)| (self.month == 3)) / np.sum(self.month>0)

        self.get_area(481)
        np.sum((self.month == 4) | (self.month == 5) | (self.month == 6) | (self.month == 3)) / np.sum(self.month > 0)

        self.get_area(483)
        np.sum((self.year_c == 2019))
        breakpoint()

        np.sum(self.csv_file.asd_code == 481)  # ap
        np.sum(self.csv_file.asd_code == 481) / self.csv_file.shape[0]
        np.sum(self.csv_file.asd_code == 482) / self.csv_file.shape[0]
        np.sum(self.csv_file.asd_code == 483) / self.csv_file.shape[0]


    def plot_lat_lon(self):
        self.get_area(481)
        lat = self.df.latitude_haul_start
        lon = self.df.longitude_haul_start
        b_depth = self.df.depth_bottom_haul_start_m
        g_depth = self.df.depth_gear_haul_start_m

        fig, axs = plt.subplots(3, 4, figsize=(24, 14))

        axs[0, 0].hist(lon, facecolor='red', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)
        axs[0, 1].hist(lat, facecolor='red', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)
        axs[0, 2].hist(g_depth, facecolor='red', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)
        axs[0, 3].hist(b_depth, facecolor='red', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)

        for i in range(0, 4):
            axs[0, i].xaxis.set_tick_params(labelsize=14)
            axs[0, i].yaxis.set_tick_params(labelsize=14)

        axs[0, 2].set_xlim([0, 350])
        axs[0, 3].set_xlim([0, 4000])


        self.get_area(482)
        lat = self.df.latitude_haul_start
        lon = self.df.longitude_haul_start
        b_depth = self.df.depth_bottom_haul_start_m
        g_depth = self.df.depth_gear_haul_start_m

        axs[1, 0].hist(lon, facecolor='blue', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)
        axs[1, 1].hist(lat, facecolor='blue', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)
        axs[1, 2].hist(g_depth, facecolor='blue', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)
        axs[1, 3].hist(b_depth, facecolor='blue', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)

        for i in range(0, 4):
            axs[1, i].xaxis.set_tick_params(labelsize=14)
            axs[1, i].yaxis.set_tick_params(labelsize=14)

        axs[1, 2].set_xlim([0, 350])
        axs[1, 3].set_xlim([0, 4000])


        self.get_area(483)
        lat = self.df.latitude_haul_start
        lon = self.df.longitude_haul_start
        b_depth = self.df.depth_bottom_haul_start_m
        g_depth = self.df.depth_gear_haul_start_m

        axs[2, 0].hist(lon, facecolor='olive', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)
        axs[2, 1].hist(lat, facecolor='olive', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)
        axs[2, 2].hist(g_depth, facecolor='olive', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)
        axs[2, 3].hist(b_depth, facecolor='olive', edgecolor='gray', linewidth=0.5, alpha=0.75, bins=30)

        for i in range(0, 4):
            axs[2, i].xaxis.set_tick_params(labelsize=14)
            axs[2, i].yaxis.set_tick_params(labelsize=14)

        axs[2, 1].set_xlim([-55, -53])
        axs[2, 2].set_xlim([0, 350])
        axs[2, 3].set_xlim([0, 4000])

        axs[0, 0].set_ylabel('Catch events', fontsize=15)
        axs[1, 0].set_ylabel('Catch events', fontsize=15)
        axs[2, 0].set_ylabel('Catch events', fontsize=15)

        axs[2, 0].set_xlabel('Longitude', fontsize=15)
        axs[2, 1].set_xlabel('Latitude', fontsize=15)
        axs[2, 2].set_xlabel('Gear depth (m)', fontsize=15)
        axs[2, 3].set_xlabel('Bottom depth (m)', fontsize=15)

        axs[0, 0].grid(alpha=0.45)  # nice and clean grid
        axs[0, 1].grid(alpha=0.45)  # nice and clean grid
        axs[1, 0].grid(alpha=0.45)  # nice and clean grid
        axs[1, 1].grid(alpha=0.45)  # nice and clean grid
        axs[2, 0].grid(alpha=0.45)  # nice and clean grid
        axs[2, 1].grid(alpha=0.45)  # nice and clean grid
        axs[2, 2].grid(alpha=0.45)  # nice and clean grid
        axs[2, 3].grid(alpha=0.45)  # nice and clean grid
        axs[1, 2].grid(alpha=0.45)  # nice and clean grid
        axs[1, 3].grid(alpha=0.45)  # nice and clean grid
        self.save_plot(plt_name='fishing_lon_lat')
        return



    def plot_fishing_season(self):
        self.get_area(481)
        ap_c = np.zeros(12)
        ap_av = np.zeros(12)
        std_ap = np.zeros(12)
        c = -1
        for i in range(1, 13):
            c = c + 1
            ap_c[c] = np.sum(self.month == i)
            ap_av[c] = np.nanmean(self.df.krill_greenweight_kg[self.month==i]) / 1000
            std_ap[c] = np.nanstd(self.df.krill_greenweight_kg[self.month==i]) / 1000

        self.get_area(482)
        so_c = np.zeros(12)
        so_av = np.zeros(12)
        std_so = np.zeros(12)
        c = -1
        for i in range(1, 13):
            c = c + 1
            so_c[c] = np.sum(self.month == i)
            so_av[c] = np.nanmean(self.df.krill_greenweight_kg[self.month == i]) / 1000
            std_so[c] = np.nanstd(self.df.krill_greenweight_kg[self.month == i]) / 1000

        self.get_area(483)
        sg_c = np.zeros(12)
        sg_av = np.zeros(12)
        std_sg = np.zeros(12)
        c = -1
        for i in range(1, 13):
            c = c + 1
            sg_c[c] = np.sum(self.month == i)
            sg_av[c] = np.nanmean(self.df.krill_greenweight_kg[self.month == i]) / 1000
            std_sg[c] = np.nanstd(self.df.krill_greenweight_kg[self.month == i]) / 1000

        fig, axs = plt.subplots(2, 2, figsize=(20, 8))
        x_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        axs[0,0].bar(x_name, ap_c, color='red', alpha=0.75)
        axs[0,0].bar(x_name, so_c, bottom=ap_c, color='blue', alpha=0.75)
        axs[0,0].bar(x_name, sg_c, bottom=ap_c + so_c, color='olive', alpha=0.75)
        axs[0,0].set_ylabel("frequency", fontsize=13)
        axs[0,0].set_xlabel("month", fontsize=13)
        axs[0,0].legend(["AP", "SO", "SG"], fontsize=13)
        axs[0,0].grid(alpha=0.45)
        axs[0, 0].set_ylim([0, 18000])

        axs[1, 0].plot(x_name, ap_av, 'r-o', linewidth=4, alpha=0.75)
        #axs[1, 0].scatter(x_name, ap_av, c='r', linewidth=5, s=3)
        #axs[1, 0].errorbar(x_name, ap_av, yerr=std_ap, color="r", alpha=0.5)
        axs[1, 0].plot(x_name, so_av, 'b-o', linewidth=4, alpha=0.75)
        axs[1, 0].plot(x_name, sg_av,'-o', c='olive', linewidth=4, alpha=0.75)
        axs[1, 0].set_ylabel("weight (tonnes)", fontsize=13)
        axs[1, 0].set_xlabel("month", fontsize=13)
        axs[1, 0].legend(["AP", "SO", "SG"], fontsize=13)
        axs[1, 0].grid(alpha=0.45)
        axs[1, 0].set_ylim([0, 50])


        self.get_area(481)
        ap_c = np.zeros(18)
        ap_av = np.zeros(18)
        std_ap = np.zeros(18)
        c = -1
        for i in range(2006, 2024):
            c = c + 1
            ap_c[c] = np.sum(self.year_c == i)
            ap_av[c] = np.nanmean(self.df.krill_greenweight_kg[self.year_c == i]) / 1000
            std_ap[c] = np.nanstd(self.df.krill_greenweight_kg[self.year_c == i]) / 1000

        self.get_area(482)
        so_c = np.zeros(18)
        so_av = np.zeros(18)
        std_so = np.zeros(18)
        c = -1
        for i in range(2006, 2024):
            c = c + 1
            so_c[c] = np.sum(self.year_c == i)
            so_av[c] = np.nanmean(self.df.krill_greenweight_kg[self.year_c == i]) / 1000
            std_so[c] = np.nanstd(self.df.krill_greenweight_kg[self.year_c == i]) / 1000

        self.get_area(483)
        sg_c = np.zeros(18)
        sg_av = np.zeros(18)
        std_sg = np.zeros(18)
        c = -1
        for i in range(2006, 2024):
            c = c + 1
            sg_c[c] = np.sum(self.year_c == i)
            sg_av[c] = np.nanmean(self.df.krill_greenweight_kg[self.year_c == i]) / 1000
            std_sg[c] = np.nanstd(self.df.krill_greenweight_kg[self.year_c == i]) / 1000

        x_name = ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017',
                  '2018', '2019', '2020', '2021', '2022', '2023']

        axs[0, 1].bar(x_name, ap_c, color='red', alpha=0.75)
        axs[0, 1].bar(x_name, so_c, bottom=ap_c, color='blue', alpha=0.75)
        axs[0, 1].bar(x_name, sg_c, bottom=ap_c + so_c, color='olive', alpha=0.75)
        axs[0, 1].set_ylabel("frequency", fontsize=13)
        axs[0, 1].set_xlabel("year", fontsize=13)
        axs[0, 1].legend(["AP", "SO", "SG"], fontsize=13)
        axs[0, 1].grid(alpha=0.45)
        axs[0, 1].set_ylim([0, 18000])

        axs[1, 1].plot(x_name, ap_av, 'r-o', linewidth=4, alpha=0.75)
        axs[1, 1].plot(x_name, so_av, 'b-o', linewidth=4, alpha=0.75)
        axs[1, 1].plot(x_name, sg_av, '-o', c='olive', linewidth=4, alpha=0.75)
        axs[1, 1].set_ylabel("weight (tonnes)", fontsize=13)
        axs[1, 1].set_xlabel("year", fontsize=13)
        axs[1, 1].legend(["AP", "SO", "SG"], fontsize=13)
        axs[1, 1].grid(alpha=0.45)
        axs[1, 1].set_ylim([0, 50])

        plt.tight_layout()
        self.save_plot(plt_name='fishing_season')
        return

    def get_area(self, asd_code):
        self.area_rule = self.csv_file.asd_code == asd_code  # key code for area
        self.get_time(self.csv_file)
        self.df = self.csv_file[self.area_rule]
        return


    def get_time(self, csv_file):
        if not os.path.exists(self.time_file):  # if file doesn't exist
            time_array = csv_file.datetime_haul_start
            shp_t = np.shape(time_array)[0]
            time_store = np.zeros([shp_t, 4])

            for i in range(0, shp_t):
                d1 = datetime.strptime(time_array[i], "%Y-%m-%d %H:%M:%S")
                time_store[i, 0] = d1.hour
                time_store[i, 1] = d1.day
                time_store[i, 2] = d1.month
                time_store[i, 3] = d1.year

            np.save(self.time_file, time_store)
            time_array = np.load(self.time_file)
            print('Saving file: ' + self.time_file)
        else:
            time_array = np.load(self.time_file)
            print('File exists: ' + self.time_file)
        time_sub = time_array[self.area_rule]
        self.hour = time_sub[:, 0].astype(int)
        self.day = time_sub[:, 1].astype(int)
        self.month = time_sub[:, 2].astype(int)
        self.year_c = time_sub[:, 3].astype(int)
        return

    def save_plot(self, plt_name):
        savefile = self.figures_path + plt_name + '.png'
        print('Saving file: ' + savefile)
        plt.savefig(savefile, dpi=400)
        plt.close()
        return

def plot_SOIN_recruit_dom_paths(compile_folder, analysis_folder):
    fig, axs = plt.subplots(figsize=(24, 12), nrows=3, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()}, sharex=True,
                            layout='constrained', sharey=True)

    years = np.arange(2005, 2009+1,1)
    release_n = 20
    max_v = 7.84
    l_max = max_v
    skip_v = 0.02
    offset = 0.15

    idx = -1
    key_name = 'SOIN'
    idy = 0
    for y in years:
        idx = idx + 1
        p_plot = PlotData(key=key_name, year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        p_plot.plot_background(background='BSSI', ax_name=axs[idy, idx])
        filename = compile_folder + p_plot.file_prefix + 'dom_paths.npy'
        dom_paths = np.load(filename)
        dom_paths = dom_paths.astype(float)
        # dom_paths[dom_paths == 0] = np.nan
        dom_paths = (dom_paths / ((release_n) * 10000)) * 100
        dom_paths[dom_paths>max_v] = max_v- offset
        #max_val = np.nanmax(dom_paths) / 50
        #n_levels = np.arange(np.nanmin(dom_paths), np.nanmax(dom_paths), max_val)
        n_levels = np.arange(0, l_max, skip_v)
    # self.plot1 = plt.contourf(self.df['lon_bin_vals'][:], self.df['lat_bin_vals'][:], dom_paths.T, levels=n_levels, cmap=self.dom_cmap,
    # transform=ccrs.PlateCarree(), extend='both')
        d_map = axs[idy, idx].contourf(p_plot.df['lon_bin_vals'][:], p_plot.df['lat_bin_vals'][:], dom_paths.T,
                                cmap=plt.get_cmap('jet'), levels=n_levels,
                                transform=ccrs.PlateCarree(), vmin =0, vmax = max_v)
        axs[idy, idx].set_title(str(y+1))


    years = np.arange(2010, 2014 + 1, 1)
    idx = -1
    idy = 1
    for y in years:
        idx = idx + 1
        p_plot = PlotData(key=key_name, year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        p_plot.plot_background(background='BSSI', ax_name=axs[idy, idx])
        filename = compile_folder + p_plot.file_prefix + 'dom_paths.npy'
        dom_paths = np.load(filename)
        dom_paths = dom_paths.astype(float)
        # dom_paths[dom_paths == 0] = np.nan
        dom_paths = (dom_paths / ((release_n) * 10000)) * 100
        dom_paths[dom_paths > max_v] = max_v - offset
        # max_val = np.nanmax(dom_paths) / 50
        # n_levels = np.arange(np.nanmin(dom_paths), np.nanmax(dom_paths), max_val)
        n_levels = np.arange(0, l_max, skip_v)
        # self.plot1 = plt.contourf(self.df['lon_bin_vals'][:], self.df['lat_bin_vals'][:], dom_paths.T, levels=n_levels, cmap=self.dom_cmap,
        # transform=ccrs.PlateCarree(), extend='both')
        # im = axs[idy, idx].imshow([p_plot.df['lon_bin_vals'][:], p_plot.df['lat_bin_vals'][:], dom_paths.T],
        #                           cmap=plt.get_cmap('jet'), vmin=0, vmax=max_v)
        #levels = n_levels
        d_map = axs[idy, idx].contourf(p_plot.df['lon_bin_vals'][:], p_plot.df['lat_bin_vals'][:], dom_paths.T,
                                        cmap=plt.get_cmap('jet'), levels = n_levels,
                                        transform=ccrs.PlateCarree(), vmin=0, vmax=max_v)
        axs[idy, idx].set_title(str(y + 1))

    cbar = fig.colorbar(d_map, ax=axs[idy, idx], shrink=0.99, extend='both')
    cbar.ax.tick_params(labelsize=15)
    #ticks = range(0, max_v + 1, 1)
    #im.set_clim(0, 5)




    years = np.arange(2015, 2019 + 1, 1)
    idx = -1
    idy = 2
    for y in years:
        idx = idx + 1
        p_plot = PlotData(key=key_name, year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        p_plot.plot_background(background='BSSI', ax_name=axs[idy, idx])
        filename = compile_folder + p_plot.file_prefix + 'dom_paths.npy'
        dom_paths = np.load(filename)
        dom_paths = dom_paths.astype(float)
        # dom_paths[dom_paths == 0] = np.nan
        dom_paths = (dom_paths / ((release_n) * 10000)) * 100
        dom_paths[dom_paths > max_v] = max_v - offset
        # max_val = np.nanmax(dom_paths) / 50
        # n_levels = np.arange(np.nanmin(dom_paths), np.nanmax(dom_paths), max_val)
        n_levels = np.arange(0, l_max, skip_v)
        # self.plot1 = plt.contourf(self.df['lon_bin_vals'][:], self.df['lat_bin_vals'][:], dom_paths.T, levels=n_levels, cmap=self.dom_cmap,
        # transform=ccrs.PlateCarree(), extend='both')
        d_map = axs[idy, idx].contourf(p_plot.df['lon_bin_vals'][:], p_plot.df['lat_bin_vals'][:], dom_paths.T,
                                       cmap=plt.get_cmap('jet'), levels=n_levels,
                                       transform=ccrs.PlateCarree(), vmin=0, vmax=max_v)
        axs[idy, idx].set_title(str(y+1))


    name_p = 'dom_paths_' + key_name
    p_plot.save_plot(plt_name=name_p)
    return

def plot_BSSI_recruit_dom_paths(compile_folder, analysis_folder):
    fig, axs = plt.subplots(figsize=(24, 12), nrows=3, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()},
                            sharex=True,
                            layout='constrained', sharey=True)

    years = np.arange(2005, 2009 + 1, 1)
    release_n = 20
    max_v = 4.84
    l_max = max_v
    skip_v = 0.02
    offset = 0.15

    idx = -1
    key_name = 'BSSI'
    idy = 0
    for y in years:
        idx = idx + 1
        p_plot = PlotData(key=key_name, year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        p_plot.plot_background(background='BSSI', ax_name=axs[idy, idx])
        filename = compile_folder + p_plot.file_prefix + 'dom_paths.npy'
        dom_paths = np.load(filename)
        dom_paths = dom_paths.astype(float)
        # dom_paths[dom_paths == 0] = np.nan
        dom_paths = (dom_paths / ((release_n) * 10000)) * 100
        dom_paths[dom_paths > max_v] = max_v - offset
        # max_val = np.nanmax(dom_paths) / 50
        # n_levels = np.arange(np.nanmin(dom_paths), np.nanmax(dom_paths), max_val)
        n_levels = np.arange(0, l_max, skip_v)
        # self.plot1 = plt.contourf(self.df['lon_bin_vals'][:], self.df['lat_bin_vals'][:], dom_paths.T, levels=n_levels, cmap=self.dom_cmap,
        # transform=ccrs.PlateCarree(), extend='both')
        d_map = axs[idy, idx].contourf(p_plot.df['lon_bin_vals'][:], p_plot.df['lat_bin_vals'][:], dom_paths.T,
                                       cmap=plt.get_cmap('jet'), levels=n_levels,
                                       transform=ccrs.PlateCarree(), vmin=0, vmax=max_v)
        axs[idy, idx].set_title(str(y + 1))

    years = np.arange(2010, 2014 + 1, 1)
    idx = -1
    idy = 1
    for y in years:
        idx = idx + 1
        p_plot = PlotData(key=key_name, year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        p_plot.plot_background(background='BSSI', ax_name=axs[idy, idx])
        filename = compile_folder + p_plot.file_prefix + 'dom_paths.npy'
        dom_paths = np.load(filename)
        dom_paths = dom_paths.astype(float)
        # dom_paths[dom_paths == 0] = np.nan
        dom_paths = (dom_paths / ((release_n) * 10000)) * 100
        dom_paths[dom_paths > max_v] = max_v - offset
        # max_val = np.nanmax(dom_paths) / 50
        # n_levels = np.arange(np.nanmin(dom_paths), np.nanmax(dom_paths), max_val)
        n_levels = np.arange(0, l_max, skip_v)
        # self.plot1 = plt.contourf(self.df['lon_bin_vals'][:], self.df['lat_bin_vals'][:], dom_paths.T, levels=n_levels, cmap=self.dom_cmap,
        # transform=ccrs.PlateCarree(), extend='both')
        d_map = axs[idy, idx].contourf(p_plot.df['lon_bin_vals'][:], p_plot.df['lat_bin_vals'][:], dom_paths.T,
                                       cmap=plt.get_cmap('jet'), levels=n_levels,
                                       transform=ccrs.PlateCarree(), vmin=0, vmax=max_v)
        axs[idy, idx].set_title(str(y + 1))

    cbar = fig.colorbar(d_map, ax=axs[idy, idx], shrink=0.99, extend='both')
    cbar.ax.tick_params(labelsize=15)

    years = np.arange(2015, 2019 + 1, 1)
    idx = -1
    idy = 2
    for y in years:
        idx = idx + 1
        p_plot = PlotData(key=key_name, year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        p_plot.plot_background(background='BSSI', ax_name=axs[idy, idx])
        filename = compile_folder + p_plot.file_prefix + 'dom_paths.npy'
        dom_paths = np.load(filename)
        dom_paths = dom_paths.astype(float)
        # dom_paths[dom_paths == 0] = np.nan
        dom_paths = (dom_paths / ((release_n) * 10000)) * 100
        dom_paths[dom_paths > max_v] = max_v - offset
        # max_val = np.nanmax(dom_paths) / 50
        # n_levels = np.arange(np.nanmin(dom_paths), np.nanmax(dom_paths), max_val)
        n_levels = np.arange(0, l_max, skip_v)
        # self.plot1 = plt.contourf(self.df['lon_bin_vals'][:], self.df['lat_bin_vals'][:], dom_paths.T, levels=n_levels, cmap=self.dom_cmap,
        # transform=ccrs.PlateCarree(), extend='both')
        d_map = axs[idy, idx].contourf(p_plot.df['lon_bin_vals'][:], p_plot.df['lat_bin_vals'][:], dom_paths.T,
                                       cmap=plt.get_cmap('jet'), levels=n_levels,
                                       transform=ccrs.PlateCarree(), vmin=0, vmax=max_v)
        axs[idy, idx].set_title(str(y + 1))

    name_p = 'dom_paths_' + key_name
    p_plot.save_plot(plt_name=name_p)
    return




def plot_retain(compile_folder, analysis_folder):
    years = np.arange(2006, 2020 + 1, 1)
    counter=-1
    flush_time = np.zeros(np.shape(years))
    for y in years:
        counter = counter + 1
        p_plot = PlotData(key='SOIN', year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        filename = p_plot.compile_folder + p_plot.file_prefix + 'retain_SG.csv'
        r_table = pd.read_csv(filename)
        flush_time[counter] = np.nanmean((r_table.retain_time/24))

    fig, ax1 = plt.subplots(1, 1, figsize=(26, 14))
    ax1.bar(np.arange(0, np.shape(flush_time)[0]), flush_time, color='r', alpha=.75)
    ax1.set_ylabel('time (days)', color='r', fontsize=15)
    ax1.xaxis.set_tick_params(labelsize=14)
    ax1.yaxis.set_tick_params(labelsize=14)
    uniq_years = np.unique(years)
    plt.xticks(np.arange(0, np.shape(uniq_years)[0]),
               ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
                '2019', '2020'])
    plt.grid(alpha=0.45)  # nice and clean grid
    p_plot.save_plot('flushing_time')

    # filename = p_plot.compile_folder + p_plot.file_prefix + 'site_retention.npy'
    # ret_sites = np.load(filename)
    # p_plot.init_plot()
    # p_plot.plot_background(background='SG')
    # plt.scatter(ret_sites[0,:],ret_sites[1,:])
    return

def plot_arrivals(compile_folder, analysis_folder):
    import datetime
    years = np.arange(2005, 2019 + 1, 1)
    # the catch dataset is now in figures folder; use the data to compare against simulations;
    recruit_v = np.zeros(np.shape(years))
    catch_v = np.zeros(np.shape(years))

    counter = -1
    for y in years:
        counter = counter + 1
        p_plot = PlotData(key='BSSI', year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        filename = p_plot.compile_folder + p_plot.file_prefix + 'recruit_SG.csv'
        r_table = pd.read_csv(filename)
        r_table.date = pd.to_datetime(r_table.date, format="%d/%m/%Y")
        df = r_table.sort_values(by='date')
        dates = [datetime.datetime.strptime(ts, "%d/%m/%Y") for ts in r_table.date]
        idx = dates.sort()
        breakpoint()
        dates.sort()
        sorteddates = [datetime.datetime.strftime(ts, "%d/%m/%Y") for ts in dates]
        sorteddates

        # catch_v[counter] = p_plot.read_catch()
        # recruit_v[counter] = p_plot.get_recruits()


def plot_linreg(compile_folder, analysis_folder):
    from scipy import stats
    years = np.arange(2005, 2019 + 1, 1)
    # the catch dataset is now in figures folder; use the data to compare against simulations;
    recruit_v = np.zeros(np.shape(years))
    catch_v = np.zeros(np.shape(years))

    counter = -1
    for y in years:
        counter = counter + 1
        p_plot = PlotData(key='BSSI', year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        catch_v[counter] = p_plot.read_catch()
        recruit_v[counter] = p_plot.get_recruits()

    fig, ax1 = plt.subplots(2, 1, figsize=(12, 8))

    varx = recruit_v
    vary = catch_v
    mask = ~np.isnan(varx) & ~np.isnan(vary)
    res = stats.linregress(varx[mask], vary[mask])
    print(stats.pearsonr(varx[mask], vary[mask]))

    ax1[0].plot(varx, vary, 'r.', label='original data', markersize=12)
    ax1[0].plot(varx, res.intercept + res.slope * varx, 'k', label='fitted line', linewidth=4)
    ax1[0].grid(alpha=0.45)
    ax1[0].set_ylabel('weight (tonnes)', fontsize=15)
    ax1[0].set_xlabel('AP fraction recruited (%)', fontsize=15)


    counter = -1
    for y in years:
        counter = counter + 1
        p_plot = PlotData(key='SOIN', year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        catch_v[counter] = p_plot.read_catch()
        recruit_v[counter] = p_plot.get_recruits()

    varx = recruit_v
    vary = catch_v
    mask = ~np.isnan(varx) & ~np.isnan(vary)
    res = stats.linregress(varx[mask], vary[mask])
    print(stats.pearsonr(varx[mask], vary[mask]))
    ax1[1].grid(alpha=0.45)
    ax1[1].plot(varx, vary, 'r.', label='original data', markersize=12)
    ax1[1].plot(varx, res.intercept + res.slope * varx, 'k', label='fitted line', linewidth=4)
    ax1[1].set_ylabel('weight (tonnes)', fontsize=15)
    ax1[1].set_xlabel('SO fraction recruited (%)', fontsize=15)


    ax1[0].xaxis.set_tick_params(labelsize=14)
    ax1[0].yaxis.set_tick_params(labelsize=14)
    ax1[1].xaxis.set_tick_params(labelsize=14)
    ax1[1].yaxis.set_tick_params(labelsize=14)

    p_plot.save_plot('recruit_correlation')
    return


def plot_recruit_stat(compile_folder, analysis_folder):
    years = np.arange(2005, 2019 + 1, 1)
    release_number = 10

    # the catch dataset is now in figures folder; use the data to compare against simulations;
    recruit_v = np.zeros(np.shape(years))
    recruit_t = np.zeros(np.shape(years))
    catch_v = np.zeros(np.shape(years))
    temp_v = np.zeros(np.shape(years))
    chl_v = np.zeros(np.shape(years))
    o2_v = np.zeros(np.shape(years))

    counter = -1
    for y in years:
        counter = counter + 1
        p_plot = PlotData(key='BSSI', year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        catch_v[counter] = p_plot.read_catch()
        recruit_v[counter] = p_plot.get_recruits()
        recruit_t[counter] = p_plot.get_recruit_times()
        temp_v[counter] = p_plot.get_temp_exp()
        chl_v[counter] = p_plot.get_chl_exp()
        o2_v[counter] = p_plot.get_o2_exp()


    fig, ax1 = plt.subplots(2, 2, figsize=(26, 14))

    axis1_title = 'weight (tonnes)'
    axis2_title = 'time (days)'
    ax2 = ax1[0, 0].twinx()
    ax2.set_ylabel(axis1_title, color='k', fontsize=15)
    ax2.plot(catch_v, 'k--', linewidth=4, alpha=0.75)
    ax1[0, 0].set_xlabel('year', fontsize=15)
    ax1[0, 0].set_ylabel(axis2_title, color='r', fontsize=15)
    ax1[0, 0].plot(np.arange(0, np.shape(recruit_t)[0]), recruit_t, color='r', linewidth=4)
    ax1[0,0].xaxis.set_tick_params(labelsize=14)
    ax1[0,0].yaxis.set_tick_params(labelsize=14)
    ax2.yaxis.set_tick_params(labelsize=14)



    uniq_years = np.unique(years)
    plt.xticks(np.arange(0, np.shape(uniq_years)[0]),
               ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
                '2019', '2020'])
    plt.grid(alpha=0.45)  # nice and clean grid

    axis1_title = 'weight (tonnes)'
    axis2_title = 'percentage recruited (%)'
    ax2 = ax1[0, 1].twinx()
    ax2.set_ylabel(axis1_title, color='k', fontsize=15)
    ax2.plot(catch_v, 'k--', linewidth=4, alpha=0.75)
    ax1[0, 1].set_ylabel(axis2_title, color='r', fontsize=15)
    ax1[0, 1].set_xlabel('year', fontsize=15)
    ax1[0, 1].plot(np.arange(0, np.shape(recruit_v)[0]), recruit_v, color='r', linewidth=4)
    ax1[0, 1].xaxis.set_tick_params(labelsize=14)
    ax1[0, 1].yaxis.set_tick_params(labelsize=14)
    ax2.yaxis.set_tick_params(labelsize=14)
    uniq_years = np.unique(years)
    plt.xticks(np.arange(0, np.shape(uniq_years)[0]),
               ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
                '2019', '2020'])
    plt.grid(alpha=0.45)  # nice and clean grid

    years = np.arange(2005, 2019 + 1, 1)
    release_number = 10

    # the catch dataset is now in figures folder; use the data to compare against simulations;
    recruit_v = np.zeros(np.shape(years))
    recruit_t = np.zeros(np.shape(years))
    catch_v = np.zeros(np.shape(years))
    temp_v = np.zeros(np.shape(years))
    chl_v = np.zeros(np.shape(years))
    o2_v = np.zeros(np.shape(years))

    counter = -1
    for y in years:
        counter = counter + 1
        p_plot = PlotData(key='SOIN', year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
        catch_v[counter] = p_plot.read_catch()
        recruit_v[counter] = p_plot.get_recruits()
        recruit_t[counter] = p_plot.get_recruit_times()
        temp_v[counter] = p_plot.get_temp_exp()
        chl_v[counter] = p_plot.get_chl_exp()
        o2_v[counter] = p_plot.get_o2_exp()

    axis1_title = 'weight (tonnes)'
    axis2_title = 'time (days)'
    ax2 = ax1[1, 0].twinx()
    ax2.set_ylabel(axis1_title, color='k', fontsize=15)
    ax2.plot(catch_v, 'k--', linewidth=4, alpha=0.75)
    ax1[1, 0].set_ylabel(axis2_title, color='blue', fontsize=15)
    ax1[1, 0].set_xlabel('year', fontsize=15)
    ax1[1, 0].plot(np.arange(0, np.shape(recruit_t)[0]), recruit_t, color='b' , linewidth=4)
    ax1[1, 0].xaxis.set_tick_params(labelsize=14)
    ax1[1, 0].yaxis.set_tick_params(labelsize=14)
    ax2.yaxis.set_tick_params(labelsize=14)
    uniq_years = np.unique(years)
    plt.xticks(np.arange(0, np.shape(uniq_years)[0]),
               ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
                '2019', '2020'])
    plt.grid(alpha=0.45)  # nice and clean grid

    axis1_title = 'weight (tonnes)'
    axis2_title = 'percentage recruited (%)'
    ax2 = ax1[1, 1].twinx()
    ax2.set_ylabel(axis1_title, color='k', fontsize=15)
    ax2.plot(catch_v, 'k--', linewidth=4, alpha=0.75)
    ax1[1, 1].set_ylabel(axis2_title, color='b', fontsize=15)
    ax1[1, 1].set_xlabel('year', fontsize=15)
    ax1[1, 1].plot(np.arange(0, np.shape(recruit_v)[0]), recruit_v, color='b', linewidth=4)
    ax1[1, 1].xaxis.set_tick_params(labelsize=14)
    ax1[1, 1].yaxis.set_tick_params(labelsize=14)
    ax1[0, 0].set_ylim([140, 240])
    ax1[1, 0].set_ylim([140, 240])
    ax1[0, 1].set_ylim([0, 12])
    ax1[1, 1].set_ylim([0, 12])
    ax2.yaxis.set_tick_params(labelsize=14)
    uniq_years = np.unique(years)
    plt.xticks(np.arange(0, np.shape(uniq_years)[0]),
               ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
                '2019', '2020'])
    plt.grid(alpha=0.45)  # nice and clean grid


    p_plot.save_plot('recruit_stat')
    return

def plot_temp_SG(compile_folder, analysis_folder):
    from scipy import stats
    filename = compile_folder + 'CMEMS_TEMP_SGfull_2006.nc'
    nc_file = nc.Dataset(filename)
    dp = nc_file.variables['depth']
    theta = nc_file['thetao']
    lat = nc_file.variables['latitude']
    lon = nc_file.variables['longitude']
    times = num2date(nc_file['time'], nc_file['time'].units)
    cdata = CatchData()
    cdata.get_area(483)
    counter = -1

    # box 1 = WCB (north-western shelf of south georgia)
    lon1=-55
    lon2=-30
    lat1=-56
    lat2=-50
    id_lon = (lon[:] > lon1) & (lon[:] < lon2)
    id_lat = (lat[:] > lat1) & (lat[:] < lat2)
    id_dp = 0
    month =np.zeros(times.shape[0])
    year = np.zeros(times.shape[0])
    for i in range(0, times.shape[0]):
        year[i] = times[i].year
        month[i] = times[i].month

    theta_sub = theta[:, id_dp, id_lat, id_lon].data
    theta_sub[theta_sub < -20000] = np.nan
    theta_box = np.nanmean(theta_sub, axis=(1, 2))

    lat_c = cdata.df.latitude_haul_start
    lon_c = cdata.df.longitude_haul_start
    id_lat_c = (lat_c[:] > lat1) & (lat_c[:] < lat2)
    id_lon_c = (lon_c[:] > lon1) & (lon_c[:] < lon2)
    catches = np.zeros(15)
    temp = np.zeros(15)
    data_list = []
    for y in range(2006, 2020+1):
        counter = counter + 1
        id_y = cdata.year_c == y
        catches[counter] = np.nanmean(cdata.df.krill_greenweight_kg[id_lat_c & id_lon_c & id_y]) / 1000
        id_time = year == y
        vals = theta_box[id_time]
        temp[counter] = np.nanmean(vals)
        mask = ~np.isnan(vals)
        data_list.append(vals[mask])

    fig, ax1 = plt.subplots(2, 1, figsize=(26, 14))

    axis1_title = 'weight (tonnes)'
    axis2_title = 'temperature ($^\circ$C)'
    idx = 0
    ax2 = ax1[idx].twinx()
    ax1[idx].plot(np.arange(1, 16), temp, 'b', linewidth=4)
    ax2.set_ylabel(axis1_title, color='k', fontsize=15)
    ax2.plot(np.arange(1, 16), catches, 'k--', linewidth=4, alpha=0.75)
    ax1[idx].set_ylabel(axis2_title, color='b', fontsize=15)

    ax1[idx].xaxis.set_tick_params(labelsize=14)
    ax1[idx].yaxis.set_tick_params(labelsize=14)
    # ax1.set_ylim([0, 250])
    # ax1.set_ylim([0, 250])
    # ax1.set_ylim([0, 14])
    # ax1.set_ylim([0, 14])
    ax2.yaxis.set_tick_params(labelsize=14)
    plt.xticks(np.arange(1, 16),
               ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
                '2019', '2020'])
    plt.grid(alpha=0.45)  # nice and clean grid

    varx = temp
    vary = catches
    mask = ~np.isnan(varx) & ~np.isnan(vary)
    res = stats.linregress(varx[mask], vary[mask])
    print(stats.pearsonr(varx[mask], vary[mask]))


    counter = -1
    catches = np.zeros(4)
    temp = np.zeros(4)
    data_list = []
    for m in range(6, 9+1):
        counter = counter + 1
        id_m = cdata.month == m
        catches[counter] = np.nanmean(cdata.df.krill_greenweight_kg[id_lat_c & id_lon_c & id_m]) / 1000
        id_time = month == m
        vals = theta_box[id_time]
        temp[counter] = np.nanmean(vals)
        mask = ~np.isnan(vals)
        data_list.append(vals[mask])

    axis1_title = 'weight (tonnes)'
    axis2_title = 'temperature ($^\circ$C)'
    idx = 1
    ax2 = ax1[idx].twinx()
    ax1[idx].plot(np.arange(1, 5), temp, 'b', linewidth=4)
    ax2.set_ylabel(axis1_title, color='k', fontsize=15)
    ax2.plot(np.arange(1, 5), catches, 'k--', linewidth=4, alpha=0.75)
    ax1[idx].set_ylabel(axis2_title, color='b', fontsize=15)
    #ax1[idx].bar(np.arange(0, np.shape(temp)[0]), temp, color='b', alpha=0.75)

    ax1[idx].xaxis.set_tick_params(labelsize=14)
    ax1[idx].yaxis.set_tick_params(labelsize=14)
    # ax1.set_ylim([0, 250])
    # ax1.set_ylim([0, 250])
    # ax1.set_ylim([0, 14])
    # ax1.set_ylim([0, 14])
    ax2.yaxis.set_tick_params(labelsize=14)
    plt.xticks(np.arange(1, 5),
               ['Jun', 'Jul', 'Aug', 'Sep'])
    plt.grid(alpha=0.45)  # nice and clean grid
    varx = temp
    vary = catches
    mask = ~np.isnan(varx) & ~np.isnan(vary)
    res = stats.linregress(varx[mask], vary[mask])
    print(stats.pearsonr(varx[mask], vary[mask]))


    plt_name = 'temp_by_year'
    figures_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/figures/'
    savefile = figures_path + plt_name + '.png'
    print('Saving file: ' + savefile)
    plt.savefig(savefile, dpi=400)
    plt.close()
    return

def plot_catch_points(compile_folder, analysis_folder):
    cdata = CatchData()
    lon = cdata.csv_file.longitude_haul_start
    lat = cdata.csv_file.latitude_haul_start
    figures_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/figures/'
    bath_file = figures_path + 'bath.npy'
    bath_file_lon = figures_path + 'bath_lon.npy'
    bath_file_lat = figures_path + 'bath_lat.npy'
    bath_contours = np.arange(0, 5750, 300)
    bath = np.load(bath_file)
    bath_lon = np.load(bath_file_lon)
    bath_lat = np.load(bath_file_lat)
    fig = plt.figure(figsize=(12, 8))
    ax_name = fig.add_subplot(projection=ccrs.PlateCarree())
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                            edgecolor='face',
                                            facecolor='lightgrey')
    ax_name.add_feature(land_10m)
    ax_name.coastlines(resolution='10m', linewidth=0.7)
    ax_name.contour(bath_lon, bath_lat, bath, bath_contours, colors='k', alpha=0.2, linewidths=0.7,
                    transform=ccrs.PlateCarree())
    d_map = ax_name.contourf(bath_lon, bath_lat, bath, levels=bath_contours,
                  transform=ccrs.PlateCarree(), cmap= plt.get_cmap('Blues'), vmin=0, vmax=4000)

    #d_map = ax_name.imshow(bath_lon, bath_lat, bath, levels=bath_contours,
                            # transform=ccrs.PlateCarree(), cmap= plt.get_cmap('Blues'))

    cbar = plt.colorbar(d_map, extend='both', pad=0.01, ax = ax_name)
    cbar.ax.set_ylabel('depth (m)', loc='center', size=9, weight='bold')
    cbar.ax.tick_params(labelsize=10, rotation=0)

    # set extent and grid lines;
    gl = ax_name.gridlines(draw_labels=True, alpha=0.4)
    gl.top_labels = False
    gl.right_labels = False
    ax_name.scatter(lon, lat, c='r', s=1)

    ax_name.set_extent(
        [-64, -34, -70, -50])
    cdata.save_plot(plt_name='fishing_points')

    return

def plot_SG_rec_area(compile_folder, analysis_folder):
    figures_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/figures/'
    bath_file = figures_path + 'bath.npy'
    bath_file_lon = figures_path + 'bath_lon.npy'
    bath_file_lat = figures_path + 'bath_lat.npy'
    bath_contours = np.arange(0, 5750, 300)
    bath = np.load(bath_file)
    bath_lon = np.load(bath_file_lon)
    bath_lat = np.load(bath_file_lat)
    p_plot = PlotData(key='BSSI', year=2006, compile_folder=compile_folder, analysis_folder=analysis_folder)
    fig = plt.figure(figsize=(12, 8))
    ax_name = fig.add_subplot(projection=ccrs.PlateCarree())
    p_plot.plot_background(background='SG_ex', ax_name=ax_name)
    d_map = ax_name.contourf(bath_lon, bath_lat, bath, levels=bath_contours,
                             transform=ccrs.PlateCarree(), cmap=plt.get_cmap('Blues'), vmin=0, vmax=4000)

    # d_map = ax_name.imshow(bath_lon, bath_lat, bath, levels=bath_contours,
    # transform=ccrs.PlateCarree(), cmap= plt.get_cmap('Blues'))

    cbar = plt.colorbar(d_map, extend='both', pad=0.01, ax=ax_name, shrink=0.5)
    cbar.ax.set_ylabel('depth (m)', loc='center', size=9, weight='bold')
    cbar.ax.tick_params(labelsize=10, rotation=0)

    gl = ax_name.gridlines(draw_labels=True, alpha=0.4)
    gl.top_labels = False
    gl.right_labels = False
    p_plot.save_plot(plt_name='SG_area_rec')
    return


def plot_worms(compile_folder, analysis_folder, trajectory_folder):
    idx = -1
    key_name = 'BSSI'
    idy = 0
    y = 2016
    idx = idx + 1
    p_plot = PlotData(key=key_name, year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
    trajectory_file = trajectory_folder + p_plot.file_prefix + 'R1_trajectory.nc'
    nc_file = nc.Dataset(trajectory_file)
    filename = p_plot.compile_folder + p_plot.file_prefix + 'site_recruits.npy'
    r_table = np.load(filename)
    idx = r_table[2,:]>0
    lon = nc_file['lon'][idx, :]
    lat = nc_file['lat'][idx, :]
    lon = lon[::6,::2]
    lat = lat[::6,::2]
    figures_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/figures/'
    bath_file = figures_path + 'bath.npy'
    bath_file_lon = figures_path + 'bath_lon.npy'
    bath_file_lat = figures_path + 'bath_lat.npy'
    bath_contours = np.arange(0, 5750, 300)
    bath = np.load(bath_file)
    bath_lon = np.load(bath_file_lon)
    bath_lat = np.load(bath_file_lat)
    fig, ax_name = plt.subplots(figsize=(10, 12), nrows=2, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()},
                            layout='constrained')

    p_plot.plot_background(background='AP', ax_name=ax_name[0])
    ax_name[0].set_extent(
        [-64, -34, -70, -50])
    d_map = ax_name[0].contourf(bath_lon, bath_lat, bath, levels=bath_contours,
                             transform=ccrs.PlateCarree(), cmap=plt.get_cmap('Blues'), vmin=0, vmax=4000)
    site_recruit_cmap=plt.get_cmap('OrRd')
    c_vals = np.arange(0, lon.shape[1], 1)*np.ones(lon.shape)
    [ax_name[0].plot(lon[i, :], lat[i,:], color='r', linewidth=2, alpha=0.3, markersize=1) for i in range(0, lon.shape[0])]

    #ax_name.scatter(lon, lat, c=c_vals,s=10, edgecolors='gray', vmin=np.nanmean(c_vals)/2,
                              # vmax=np.nanmean(c_vals)*2, linewidth=0.2, cmap=site_recruit_cmap)



    key_name = 'SOIN'
    p_plot = PlotData(key=key_name, year=y, compile_folder=compile_folder, analysis_folder=analysis_folder)
    trajectory_file = trajectory_folder + p_plot.file_prefix + 'R1_trajectory.nc'
    nc_file = nc.Dataset(trajectory_file)
    filename = p_plot.compile_folder + p_plot.file_prefix + 'site_recruits.npy'
    r_table = np.load(filename)
    idx = r_table[2, :] > 0
    lon = nc_file['lon'][idx, :]
    lat = nc_file['lat'][idx, :]
    lon = lon[::6, ::2]
    lat = lat[::6, ::2]
    p_plot.plot_background(background='AP', ax_name=ax_name[1])
    ax_name[1].set_extent(
        [-64, -34, -70, -50])
    d_map = ax_name[1].contourf(bath_lon, bath_lat, bath, levels=bath_contours,
                                transform=ccrs.PlateCarree(), cmap=plt.get_cmap('Blues'), vmin=0, vmax=4000)
    site_recruit_cmap = plt.get_cmap('OrRd')
    c_vals = np.arange(0, lon.shape[1], 1) * np.ones(lon.shape)
    [ax_name[1].plot(lon[i, :], lat[i, :], color='r', linewidth=2, alpha=0.3, markersize=1) for i in
     range(0, lon.shape[0])]

    p_plot.save_plot(plt_name='AP_SO_worms')
    return
    #p_plot.plot_background(background='BSSI', ax_name=axs[idy, idx])





