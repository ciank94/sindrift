#file_path = '/cluster/projects/nn9828k/Cian_sinmod/copernicus_client/results/'
from post_process import Analyse
file_path = 'A:/Cian_sinmod/copernicus_client/results/'
year = 2000
key = "APSO"
x = Analyse(file_path, key, year)
x.plot_trajectory("APSO")
x.dominant_paths()
x.plot_dom_paths()
breakpoint()
#y_end = 2011