
#
# To use this script, make sure you set python path enviroment variable to include the following paths
#   ${SE2WAVE_ROOT}/utils/python:${PETSC_DIR}/lib/petsc/bin
# e.g.
#   export PYTHONPATH=${SE2WAVE_ROOT}/utils/python:${PETSC_DIR}/lib/petsc/bin
#

import PetscBinaryIO as pio
import json as j


def se2wave_load_coordinates(filename,debug=False):
  io = pio.PetscBinaryIO()
  data = dict()
  with open(filename) as fp:
    v = io.readInteger(fp)
    data['mx'] = v
    
    v = io.readInteger(fp)
    data['my'] = v
    
    v = io.readInteger(fp)
    data['nx'] = v
    
    v = io.readInteger(fp)
    data['ny'] = v
    
    objecttype = io.readObjectType(fp)
    v = io.readVec(fp)
    data['coor'] = v
  
  if debug:
    print('<debug> se2wave_load_coordinates()')
    print('<debug> #elements-x',data['mx'])
    print('<debug> #elements-y',data['my'])
    print('<debug> #basis-x',data['nx'])
    print('<debug> #basis-y',data['ny'])
  
  return data


def se2wave_load_wavefield(filename,has_displacement,has_velocity,debug=False):
  io = pio.PetscBinaryIO()
  data = dict()
  with open(filename) as fp:
    v = io.readInteger(fp)
    data['mx'] = v
    
    v = io.readInteger(fp)
    data['my'] = v
    
    v = io.readInteger(fp)
    data['nx'] = v
    
    v = io.readInteger(fp)
    data['ny'] = v
    
    v = io.readInteger(fp)
    data['step'] = v
    
    v = io.readReal(fp)
    data['time'] = v
    
    if has_displacement:
      objecttype = io.readObjectType(fp)
      v = io.readVec(fp)
      data['displ'] = v
    
    if has_velocity:
      objecttype = io.readObjectType(fp)
      v = io.readVec(fp)
      data['vel'] = v

  if debug:
    print('<debug> se2wave_load_wavefield()')
    print('<debug> #elements-x',data['mx'])
    print('<debug> #elements-y',data['my'])
    print('<debug> #basis-x',data['nx'])
    print('<debug> #basis-y',data['ny'])
    print('<debug> step',data['step'])
    print('<debug> time',data['time'])
  
  return data


def se2wave_load_json(filename,debug=False):
  jdata = dict()
  with open(filename, "r") as fp:
    jdata = j.load(fp)
    
  if debug:
    wdata = jdata['se2wave']
    print('<debug> se2wave_load_json()')
    print('<debug> spatial_dimension:',wdata['spatial_dimension'])
    print('<debug> nx:',wdata['nx'])
    print('<debug> time:',wdata['time'])
    print('<debug> fields:',wdata['data']['fields'])
    print('<debug> datafilename:',wdata['data']['filename'])
  return jdata

