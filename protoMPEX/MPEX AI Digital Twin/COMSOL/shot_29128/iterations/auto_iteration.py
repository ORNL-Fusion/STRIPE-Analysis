import subprocess
import os
from Sheath_impedance_MPEXAI import Sheath_impedance_MPEXAI
from itertools import product

N_ITERATIONS = 5

def runComsol(comsolFile):
 commandList=['comsol','batch','-inputfile',comsolFile,'-methodcall','main'];
 #print(commandList)
 subprocess.run(commandList)


def doCase(caseDir):
 if not os.path.exists(caseDir):
  os.makedirs(caseDir)
 wkdir=os.path.join(caseDir,'initial')
 if not os.path.exists(wkdir):
  os.makedirs(wkdir)
 nextDirs=[]
 for i in range(1,N_ITERATIONS):
  nextDirs.append(os.path.join(caseDir,'it'+str(i)))
  if not os.path.exists(nextDirs[-1]):
   os.makedirs(nextDirs[-1])
 # make sure comsol file for initial case exists
 initialComsolFile=os.path.join(wkdir,'initial.mph')
 if not os.path.isfile(initialComsolFile):
  os.system('cp initial.mph '+initialComsolFile)
 # run comsol for initial case
 runComsol(initialComsolFile)
 # check if output files exist
 eFieldFile=os.path.join(wkdir,'Efield.csv')
 wpFile=os.path.join(wkdir,'window_parameters.csv')
 loadingFile=os.path.join(wkdir,'loading.csv')
 if os.path.isfile(eFieldFile) and os.path.isfile(wpFile) and os.path.isfile(loadingFile):
  print('Comsol output exists');
 else:
  print('Comsol output NOT FOUND');
 # run sheath model
 Sheath_impedance_MPEXAI(eFieldFile,wpFile)
 for i,nextDir in enumerate(nextDirs):
  if not os.path.isdir(nextDir):
   os.mkdir(nextDir);
  comsolfile=os.path.join(nextDir,'loop_it'+str(i+1)+'.mph');
  # move epsilon and sigma files
  os.replace('epsilon_py.csv',os.path.join(nextDir,'epsilon.csv'))
  os.replace('sigma_py.csv',os.path.join(nextDir,'sigma.csv'))
  if nextDir==nextDirs[-1]:
   break;
  # make sure comsolfile exists
  if not os.path.isfile(comsolfile):
   os.system('cp iteration.mph '+comsolfile)
  # run comsol again
  runComsol(comsolfile)
  # check if output files exist
  eFieldFile=os.path.join(nextDir,'Efield.csv')
  wpFile=os.path.join(nextDir,'window_parameters.csv')
  loadingFile=os.path.join(nextDir,'loading.csv')
  if os.path.isfile(eFieldFile) and os.path.isfile(wpFile) and os.path.isfile(loadingFile):
   print('Comsol output exists');
  else:
   print('Comsol output NOT FOUND');
  # run sheath model
  Sheath_impedance_MPEXAI(eFieldFile,wpFile)

if os.path.isfile('initial.mph') and os.path.isfile('iteration.mph'):
    pass
else:
    print('This script needs two comsol files in the working directory:');
    print('    initial.mph (for vacuum sheaths)')
    print('    iteration.mph (which reads sheath parameters from epsilon.csv and from sigma.csv)')
    exit()
doCase('MPEXAI')

