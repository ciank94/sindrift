from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic, reader_global_landmask
from configure import Case

file_path = '/cluster/projects/nn9828k/Cian_sinmod/opendrift'
file_name = '/CMEMS_GLPHYS_D_norm_2000_2000.nc'
phys_states = file_path + file_name

# Simulation settings
time_step_hours = 6  # negative time is backwards stepping of model
duration_days = 120  # look into (time=start_time + i*time_step) for setting the start and end time of simulations;

key_list = ["SG_NW"]
for key in key_list:

    o = OceanDrift(loglevel=20)  # log_level= 0 for full diagnostics, 50 for none
    reader_samples = reader_netCDF_CF_generic.Reader(phys_states)  # read forcing variables
    reader_landmask = reader_global_landmask.Reader()  # high resolution coast for particle beaching etc.
    o.add_reader([reader_landmask, reader_samples])  # add readers to model instance

    case = Case(reader_samples, file_path, time_step_hours, duration_days)
    case.get_scenarios(key=key)

    # Simulation seeding
    o.disable_vertical_motion()
    o.seed_elements(lon=case.lon_init,
                    lat=case.lat_init,
                    time=case.t_init,
                    number=case.n_part,
                    radius=case.radius)

    # Simulation running
    o.run(duration=case.duration,
          time_step=case.time_step,
          outfile=case.trajectory_file,
          export_variables=case.export_variables)






#print(reader_samples)






