from stochpy_radiometer.do_pproc import do_pproc
#import do_pproc
from stochpy_radiometer.utils import read_pipeline_ini,fix_params
import stochpy_radiometer.detgeom as dg
import numpy as np
import optparse

def parse_command_line():
    """
    parse command line
    """
    parser = optparse.OptionParser()
    parser.add_option(
            "--ini-file", "-i", help="params file",
            dest="ini", type=str)
    params, args = parser.parse_args()
    return params

# parse cmd line
cmd_params = parse_command_line()
# read params
params = read_pipeline_ini(cmd_params.ini)['postpostproc']
binary_params = read_pipeline_ini(cmd_params.ini)['binary_params']
# fix param variable types
params = fix_params(params)
binary_params = fix_params(binary_params)
# get direction and unit vector
Direction = dg.SkyPosition(float(params['decDeg']), float(params['raHr']))
params['k_hat'] = Direction.k
#run it!
#print('Doing full post processing...')
do_pproc(params, binary_params)
