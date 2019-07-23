## Haven't figured out how to install pyroms
# looking at https://github.com/ESMG/pyroms

import pyroms
import pyroms.pyroms_toolbox as pyroms_toolbox
#import pyroms_toolbox

#import pyroms.pyroms_toolbox.pyroms_toolbox.quiver as pyquiver
#import pyroms.pyroms.pyroms as pyroms


# bring in the data
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
inputfile = dir+'/upper_ches_his.nc'
gridfile = dir + '/ROMS_grid.nc'
gridid = pyroms.pyroms.pyroms.grid.get_ROMS_grid(gridfile)

pyroms_toolbox.pyroms_toolbox.quiver('u', 'v', 10, 3, gridid, filename=inputfile)
