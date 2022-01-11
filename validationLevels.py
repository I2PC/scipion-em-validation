# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import sys

import pyworkflow.plugin as pwplugin
from pyworkflow.project import Manager
from pyworkflow.utils.path import makePath, copyFile, cleanPath

def usage(message=''):
    print("\nMake a Structure Based Drug Design for a collection of PDBs"
          "\n\n   -> scipion3 python validationLevels.py [opt1=val1 opt2=val2 ...] "
          "\n         project:  myProject"
          "\n         LEVEL 0 ====="
          "\n            map:  mymap.mrc"
          "\n            sampling:  1 [A]"
          "\n            threshold:  0.03"
          "\n            resolution:  2 [A]"
          "\n         LEVEL 1 ====="
          "\n            map1:  mymap1.mrc"
          "\n            map2:  mymap2.mrc"
          "\n         LEVEL 2 ====="
          "\n            avgs:  myaverages.mrcs"
          "\n            avgSampling:  2 [A]"
          "\n            symmetry: c1"
          "\n         LEVEL 3 ====="
          "\n            particles:  myparticles.mrcs"
          "\n            ptclSampling:  1 [A]"
          "\n            kV:  300 [kV]"
          "\n            Cs:  2.7 [mm]"
          "\n            Q0:  0.1"
          "\n         LEVEL 4 ====="
          "\n            hasAngles:  yes"
          )
    message = "\n\n  >>  %s\n" % message if message != '' else ''
    print(message)
    sys.exit(1)
    # ~/scipion3/scipion3 python ~/data/Dropbox/H/scipion-em-validation/validationLevels.py project=TestValidation map=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010948_XmippProtLocSharp/extra/sharpenedMap_1.mrc sampling=0.74 threshold=0.0025 resolution=2.6 map1=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010450_XmippProtReconstructHighRes/extra/Iter001/volume01.vol map2=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010450_XmippProtReconstructHighRes/extra/Iter001/volume02.vol avgs=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/011849_XmippProtCropResizeParticles/extra/output_images.stk avgSampling=1.24 symmetry=o particles=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010450_XmippProtReconstructHighRes/particles.sqlite ptclSampling=0.74 kV=300 Cs=2.7 Q0=0.1 hasAngles=yes

if any(i in sys.argv for i in ['-h', '-help', '--help', 'help']):
    usage()

# Manager will halp as to find paths, create project...
manager = Manager()

# Fixing some parameters depending on the arguments or taking the default ones
FNMAP1 = None
FNMAP2 = None
for arg in sys.argv:
    if arg.startswith('project='):
        PROJECT_NAME = arg.split('project=')[1]
    if arg.startswith('map='):
        FNMAP = arg.split("map=")[1]
    if arg.startswith('sampling='):
        TS = float(arg.split("sampling=")[1])
    if arg.startswith('threshold='):
        MAPTHRESHOLD = float(arg.split("threshold=")[1])
    if arg.startswith('resolution='):
        MAPRESOLUTION = float(arg.split("resolution=")[1])
    if arg.startswith('map1='):
        FNMAP1 = arg.split("map1=")[1]
    if arg.startswith('map2='):
        FNMAP2 = arg.split("map2=")[1]
    if arg.startswith('avgs='):
        FNAVGS = arg.split("avgs=")[1]
    if arg.startswith('avgSampling='):
        TSAVG = float(arg.split("avgSampling=")[1])
    if arg.startswith('symmetry='):
        SYM = arg.split("symmetry=")[1]
    if arg.startswith('particles='):
        FNPARTICLES = arg.split("particles=")[1]
    if arg.startswith('ptclSampling='):
        TSPARTICLES = float(arg.split("ptclSampling=")[1])
    if arg.startswith('kV='):
        KV = float(arg.split("kV=")[1])
    if arg.startswith('Cs='):
        CS = float(arg.split("Cs=")[1])
    if arg.startswith('Q0='):
        Q0 = float(arg.split("Q0=")[1])
    if arg.startswith('hasAngles='):
        HASANGLES = arg.split("hasAngles=")[1]

# Detect level
argsPresent = [x.split('=')[0] for x in sys.argv]
LEVEL0 = ["map", "sampling", "threshold", "resolution"]
LEVEL1 = ["map1", "map2"]
LEVEL2 = ["avgs", "avgSampling", "symmetry"]
LEVEL3 = ["particles", "ptclSampling", "kV", "Cs", "Q0"]
LEVEL4 = ["hasAngles"]

def detectLevel(labels, args):
    retval = True
    for label in labels:
        if not label in args:
            retval = False
            break
    return retval
levels = []
if detectLevel(LEVEL0, argsPresent):
    levels.append(0)
if detectLevel(LEVEL1, argsPresent):
    levels.append(1)
if detectLevel(LEVEL2, argsPresent):
    levels.append(2)
if detectLevel(LEVEL3, argsPresent):
    levels.append(3)
if detectLevel(LEVEL4, argsPresent):
    levels.append(4)

if len(levels)==0 or not 0 in levels:
    usage()

# Creating the project
projectDir = manager.getProjectPath(PROJECT_NAME)
if os.path.exists(projectDir):
    cleanPath(projectDir)
project = manager.createProject(PROJECT_NAME)
fnProjectDir = project.getPath()
os.chdir(fnProjectDir)

# Create report
from validationReport import ValidationReport
report = ValidationReport(fnProjectDir, levels)

# Level 0
from validationLevel0 import level0
protImportMap, protCreateMask, bfactor = level0(project, report, FNMAP, FNMAP1, FNMAP2, TS, MAPTHRESHOLD, MAPRESOLUTION,
                                       skipAnalysis = False)

# Level 1
if 1 in levels:
    from validationLevel1 import level1
    protImportMap1, protImportMap2 = level1(project, report, FNMAP1, FNMAP2, TS, MAPRESOLUTION,
                                            protImportMap, protCreateMask, skipAnalysis = False)

# Level 2
if 2 in levels:
    from validationLevel2 import level2
    protImportAvgs, protAvgsResizeMap = level2(project, report, protImportMap, FNAVGS, TSAVG, SYM, skipAnalysis = False)

# Level 3
if 3 in levels:
    from validationLevel3 import level3
    protImportParticles, protResizeMap, protResizeAvgs = level3(project, report, protImportMap, protImportAvgs,
                                                                FNPARTICLES, TSPARTICLES, KV, CS, Q0,
                                                                skipAnalysis = False)

# Level 4
if 4 in levels:
    from validationLevel4 import level4
    level4(project, report, protImportMap, protCreateMask, protResizeMap, SYM, MAPRESOLUTION, bfactor,
           skipAnalysis = False)

# Close report
report.closeReport()