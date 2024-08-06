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

    def plot_background(self, background):
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
        #self.ax.set_extent(
            #[self.min_lon, self.max_lon, self.min_lat, self.max_lat])
        if background == "AP":
            self.AP_lon_lat_extent()
            self.ax.set_extent(
                [self.min_lon, self.max_lon, self.min_lat,
                 self.max_lat])
            #self.ax.set_extent(
             #   [self.min_lon - self.lon_offset, self.max_lon + self.lon_offset, self.min_lat - self.lat_offset,
              #   self.max_lat + self.lat_offset])
        elif background == "SOI":
            self.SOI_lon_lat_extent()
            self.ax.set_extent(
                [self.min_lon, self.max_lon, self.min_lat,
                 self.max_lat])
        elif background == "SOIN":
            self.SOIN_lon_lat_extent()
            self.ax.set_extent(
                [self.min_lon, self.max_lon, self.min_lat,
                 self.max_lat])
        else:
            self.gen_lon_lat_extent()
            self.ax.set_extent([self.min_lon, self.max_lon, self.min_lat,
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

    def add_cbar(self, c_max, caxis_title):
        cbar = plt.colorbar(self.plot1, extend='both', pad=0.01)
        cbar.ax.set_ylabel(caxis_title, loc='center', size=9, weight='bold')
        cbar.ax.tick_params(labelsize=10, rotation=0)
        plt.clim(0, c_max)
        return

    def save_plot(self, plt_name):
        savefile = self.figures_path + plt_name + '.png'
        print('Saving file: ' + savefile)
        plt.savefig(savefile, dpi=400)
        plt.close()
        return

    def AP_lon_lat_extent(self):
        self.min_lon = -65.3
        self.max_lon = -51
        self.min_lat = -69
        self.max_lat = -56
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
        self.min_lon = -43
        self.max_lon = -33
        self.min_lat = -58
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

        fig, axs = plt.subplots(3, 4, figsize=(20, 14))

        axs[0, 0].hist(lon, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)

        axs[0, 1].hist(lat, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)
        axs[0, 2].hist(g_depth, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)
        axs[0, 3].hist(b_depth, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)

        axs[0, 2].set_xlim([0, 350])
        axs[0, 3].set_xlim([0, 4000])

        self.get_area(482)
        lat = self.df.latitude_haul_start
        lon = self.df.longitude_haul_start
        b_depth = self.df.depth_bottom_haul_start_m
        g_depth = self.df.depth_gear_haul_start_m

        axs[1, 0].hist(lon, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)
        axs[1, 1].hist(lat, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)
        axs[1, 2].hist(g_depth, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)
        axs[1, 3].hist(b_depth, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)

        axs[1, 2].set_xlim([0, 350])
        axs[1, 3].set_xlim([0, 4000])


        self.get_area(483)
        lat = self.df.latitude_haul_start
        lon = self.df.longitude_haul_start
        b_depth = self.df.depth_bottom_haul_start_m
        g_depth = self.df.depth_gear_haul_start_m

        axs[2, 0].hist(lon, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)
        axs[2, 1].hist(lat, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)
        axs[2, 2].hist(g_depth, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)
        axs[2, 3].hist(b_depth, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.75, bins=30)

        axs[2, 2].set_xlim([0, 350])
        axs[2, 3].set_xlim([0, 4000])

        axs[0, 0].set_ylabel('AP', fontsize=13)
        axs[1, 0].set_ylabel('SO', fontsize=13)
        axs[2, 0].set_ylabel('SG', fontsize=13)

        axs[2, 0].set_xlabel('Longitude', fontsize=13)
        axs[2, 1].set_xlabel('Latitude', fontsize=13)
        axs[2, 2].set_xlabel('Gear depth (m)', fontsize=13)
        axs[2, 3].set_xlabel('Bottom depth (m)', fontsize=13)

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
        breakpoint()



    def plot_fishing_season(self):
        self.get_area(481)
        ap_c = np.zeros(12)
        c = -1
        for i in range(1, 13):
            c = c + 1
            ap_c[c] = np.sum(self.month == i)

        self.get_area(482)
        so_c = np.zeros(12)
        c = -1
        for i in range(1, 13):
            c = c + 1
            so_c[c] = np.sum(self.month == i)

        self.get_area(483)
        sg_c = np.zeros(12)
        c = -1
        for i in range(1, 13):
            c = c + 1
            sg_c[c] = np.sum(self.month == i)

        fig, axs = plt.subplots(2, 1, figsize=(12, 8))
        x_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        axs[0].bar(x_name, ap_c, color='r')
        axs[0].bar(x_name, so_c, bottom=ap_c, color='b')
        axs[0].bar(x_name, sg_c, bottom=ap_c + so_c, color='y')
        axs[0].set_ylabel("Catch events", fontsize=13)
        axs[0].set_xlabel("Month", fontsize=13)
        axs[0].legend(["AP", "SO", "SG"], fontsize=13)

        self.get_area(481)
        ap_c = np.zeros(18)
        c = -1
        for i in range(2006, 2024):
            c = c + 1
            ap_c[c] = np.sum(self.year_c == i)

        self.get_area(482)
        so_c = np.zeros(18)
        c = -1
        for i in range(2006, 2024):
            c = c + 1
            so_c[c] = np.sum(self.year_c == i)

        self.get_area(483)
        sg_c = np.zeros(18)
        c = -1
        for i in range(2006, 2024):
            c = c + 1
            sg_c[c] = np.sum(self.year_c == i)

        x_name = ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017',
                  '2018', '2019', '2020', '2021', '2022', '2023']
        axs[1].bar(x_name, ap_c, color='r')
        axs[1].bar(x_name, so_c, bottom=ap_c, color='b')
        axs[1].bar(x_name, sg_c, bottom=ap_c + so_c, color='y')
        axs[1].set_ylabel("Catch events", fontsize=13)
        axs[1].set_xlabel("Year", fontsize=13)
        axs[1].legend(["AP", "SO", "SG"], fontsize=13, loc='upper left')
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


def check_retain(compile_folder, analysis_folder):
    p_plot = PlotData(key='SGCM', year=2005, compile_folder=compile_folder, analysis_folder=analysis_folder)
    filename = p_plot.compile_folder + p_plot.file_prefix + 'retain_SG.csv'
    r_table = pd.read_csv(filename)

    filename = p_plot.compile_folder + p_plot.file_prefix + 'site_retention.npy'
    ret_sites = np.load(filename)
    p_plot.init_plot()
    p_plot.plot_background(background='SG')
    plt.scatter(ret_sites[0,:],ret_sites[1,:])
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

    fig, ax1 = plt.subplots(2, 2, figsize=(20, 14))

    axis1_title = 'catch'
    axis2_title = 'time (days)'
    ax2 = ax1[0, 0].twinx()
    ax2.set_ylabel(axis1_title, color='r', fontsize=13)
    ax2.plot(catch_v, 'r')
    ax1[0, 0].set_ylabel(axis2_title, color='b', fontsize=13)
    ax1[0, 0].bar(np.arange(0, np.shape(recruit_t)[0]), recruit_t, color='b')
    uniq_years = np.unique(years)
    plt.xticks(np.arange(0, np.shape(uniq_years)[0]),
               ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
                '2019', '2020'])
    plt.grid(alpha=0.45)  # nice and clean grid

    axis1_title = 'catch'
    axis2_title = 'fraction recruited (%)'
    ax2 = ax1[0, 1].twinx()
    ax2.set_ylabel(axis1_title, color='r', fontsize=13)
    ax2.plot(catch_v, 'r')
    ax1[0, 1].set_ylabel(axis2_title, color='b', fontsize=13)
    ax1[0, 1].bar(np.arange(0, np.shape(recruit_v)[0]), recruit_v, color='b')
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

    axis1_title = 'catch'
    axis2_title = 'time (days)'
    ax2 = ax1[1, 0].twinx()
    ax2.set_ylabel(axis1_title, color='r', fontsize=13)
    ax2.plot(catch_v, 'r')
    ax1[1, 0].set_ylabel(axis2_title, color='b', fontsize=13)
    ax1[1, 0].bar(np.arange(0, np.shape(recruit_t)[0]), recruit_t, color='b')
    uniq_years = np.unique(years)
    plt.xticks(np.arange(0, np.shape(uniq_years)[0]),
               ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
                '2019', '2020'])
    plt.grid(alpha=0.45)  # nice and clean grid

    axis1_title = 'catch'
    axis2_title = 'fraction recruited (%)'
    ax2 = ax1[1, 1].twinx()
    ax2.set_ylabel(axis1_title, color='r', fontsize=13)
    ax2.plot(catch_v, 'r')
    ax1[1, 1].set_ylabel(axis2_title, color='b', fontsize=13)
    ax1[1, 1].bar(np.arange(0, np.shape(recruit_v)[0]), recruit_v, color='b')
    ax1[0, 0].set_ylim([0, 250])
    ax1[1, 0].set_ylim([0, 250])
    ax1[0, 1].set_ylim([0, 14])
    ax1[1, 1].set_ylim([0, 14])
    uniq_years = np.unique(years)
    plt.xticks(np.arange(0, np.shape(uniq_years)[0]),
               ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
                '2019', '2020'])
    plt.grid(alpha=0.45)  # nice and clean grid


    p_plot.save_plot('recruit_stat')
    return

