# todo: add functionality for plotting diagnostics from simulation- add options so that the worms plot can be done alone as a check
from configure import FileExplorer
from compile_releases import StoreReleases
from plot_data import (
    PlotData,
    CatchData,
    plot_recruit_stat,
    plot_retain,
    plot_SOIN_recruit_dom_paths,
    plot_BSSI_recruit_dom_paths,
    plot_temp_SG,
    plot_catch_points,
    plot_linreg,
    plot_arrivals,
    plot_SG_rec_area,
    plot_worms,
    plot_poster_dom_paths,
    plot_recruit_dom_paths,
    plot_ant_sub,
    plot_transit_distributions,
    plot_particles,
    plot_seasonal_rec,
)

#folder = "A:/Cian_sinmod/opendrift/"
# folder = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/'
folder = "S:/opendrift/"
#folder = "S:/nird/projects/NS9828K/opendrift/"
compile_folder = folder + "compile/"
analysis_folder = folder + "analysis/"
trajectory_folder = folder + "trajectory/"
phys_folder = folder + "phys_states/"


# figure 1:first figure in paper that shows the position of fronts relative to catch data ##
# plot_ant_sub()
#plot_catch_points(compile_folder, analysis_folder)

# figure 2:
#cdata = CatchData()
# cdata.catch_facts()
#cdata.plot_fishing_season()

# figure 3: show particles backtracked from SG for one year;
# plot_particles(compile_folder, analysis_folder, trajectory_folder, phys_folder)

# figure 4:
# plot_worms(compile_folder, analysis_folder, trajectory_folder)

# figure 5:
plot_seasonal_rec(compile_folder)

# figure 6:
# plot_transit_distributions(compile_folder, analysis_folder)

# figure 7:
# plot_linreg(compile_folder, analysis_folder)
