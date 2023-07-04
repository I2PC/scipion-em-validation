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
import math
import numpy as np

import pyworkflow.plugin as pwplugin
from pyworkflow.project import Manager
from pyworkflow.utils.path import makePath, copyFile, cleanPath
from resourceManager import sendToSlurm, waitOutput, waitUntilFinishes
from pwem.convert.atom_struct import AtomicStructHandler
from validationReport import readMap
import json

import configparser

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'config.yaml'))
useSlurm = config['QUEUE'].getboolean('USE_SLURM')

class OutOfChainsError(Exception): #TODO: remove it when updating pwem repo
    pass

class OutOfAtomsError(Exception): #TODO: remove it when updating pwem repo
    pass
class UpdatedAtomicStructHandler(AtomicStructHandler): #TODO: remove it when updating pwem repo
    """
    Class that contain utilities to handle pdb/cif files.
    Updates: get the number of atoms in the structure and raise an error 
    when trying to write a structure with more tha 99999 atoms as a PDB file
    """

    def numberAtomsInStructure(self, structure, writeAsPdb=False):
        """
        Get number of atoms in the given structure.

        When using this function to write a PDB file (writeAsPdb=True) 
        and the number of atoms is greater than 99999 raises an OutOfAtomsError
        """
        atom_records = structure.get_atoms()
        n_atoms = len(list(atom_records))

        if writeAsPdb and n_atoms > 99999:
            raise OutOfAtomsError
        else:
            return n_atoms
        
    def writeAsPdb(self, pdbFile):
        """ 
        Save structure as PDB. Be aware that this is not a lossless conversion
        Returns False if conversion is not possible. 
        True otherwise.
        Updates: check that the number of atoms present in the structure
        is not greater than 99999.
        """
        # check input is not PDB
        if self.type == self.PDB:
            pass
        else:
            # rename long chains
            try:
                chainmap = self.renameChains(self.structure)
            except OutOfChainsError:
                print("Too many chains to represent in PDB format")
                return False

            for new, old in chainmap.items():
                # for new, old in list(chainmap.items()):
                if new != old:
                    print("Renaming chain {0} to {1}".format(old, new))
            
            # Get number of atoms
            try:
                atoms = self.numberAtomsInStructure(self.getStructure(), writeAsPdb=True)
            except OutOfAtomsError:
                print("Too many atoms to represent in PDB format")
                return False

        self._write(pdbFile)

        return True

def usage(message=''):
    print("\nMake a Map Validation Report"
          "\n\n   -> scipion3 python validationLevels.py [opt1=val1 opt2=val2 ...] "
          "\n         inputParams: [uploads_path]/inputParams.json"
          "\n         or"
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
          "\n         LEVEL 5 ====="
          '\n            micrographs:  dir/*.mrc'
          '\n            micSampling:  1 [A]'
          "\n         LEVEL A ====="
          '\n            atomicModel:  mymodel.pdb [or .cif]'
          '\n            doMultimodel:  yes # This method is computationally very expensive'
          "\n         LEVEL W ====="
          '\n            workflow:  workflow.json'
          "\n         LEVEL O ====="
          '\n            xlm:  xlm.txt'
          '\n'
          '\n            saxs: saxsCurve.dat'
          '\n'
          '\n            untiltedMic: untilted.mrc'
          '\n            tiltedMic: tilted.mrc'
          '\n            tiltkV: 300 [kV]'
          '\n            tiltCs: 2.7 [mm]'
          '\n            tiltQ0: 0.1'
          '\n            tiltSampling: 1 [A]'
          '\n            tiltAngle: 15 [degrees]'
          '\n            untiltedCoords: coords [.xmd from Xmipp or .json from Eman]'
          '\n            tiltedCoords: coords [.xmd from Xmipp or .json from Eman]'
          )
    message = "\n\n  >>  %s\n" % message if message != '' else ''
    print(message)
    sys.exit(1)
    # ~/scipion3/scipion3 python ~/data/Dropbox/H/scipion-em-validation/validationLevels.py project=TestValidation map=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010948_XmippProtLocSharp/extra/sharpenedMap_1.mrc sampling=0.74 threshold=0.0025 resolution=2.6 map1=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010450_XmippProtReconstructHighRes/extra/Iter001/volume01.vol map2=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010450_XmippProtReconstructHighRes/extra/Iter001/volume02.vol avgs=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/012458_XmippProtCropResizeParticles/extra/output_images.stk avgSampling=1.24 symmetry=o particles=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010450_XmippProtReconstructHighRes/particles.sqlite ptclSampling=0.74 kV=300 Cs=2.7 Q0=0.1 hasAngles=yes micrographs="/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/006458_XmippProtMovieCorr/extra/*mic.mrc" micSampling=0.49 atomicModel=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/centered4V1W.pdb doMultimodel=yes workflow=http://nolan.cnb.csic.es/cryoemworkflowviewer/workflow/637ca2bbcd57e45e88f6fabb7f6b1095a3ca0de6 xlm=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/xlm.txt saxs=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/SASDE55.dat untiltedMic=/home/coss/scipion3/data/tests/eman/mics/ip3r10252011-0005_0-2.hdf tiltedMic=/home/coss/scipion3/data/tests/eman/mics/ip3r10252011-0005_10.hdf tiltkV=200 tiltCs=2 tiltQ0=0.1 tiltSampling=1.88 tiltAngle=60 untiltedCoords=/home/coss/scipion3/data/tests/eman/coords/ip3r10252011-0005_0-2_info.json tiltedCoords=/home/coss/scipion3/data/tests/eman/coords/ip3r10252011-0005_10_info.json
    # ~/scipion3/scipion3 python ~/data/Dropbox/H/scipion-em-validation/validationLevels.py project=TestValidation map=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010948_XmippProtLocSharp/extra/sharpenedMap_1.mrc sampling=0.74 threshold=0.0025 resolution=2.6 map1=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010450_XmippProtReconstructHighRes/extra/Iter001/volume01.vol map2=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010450_XmippProtReconstructHighRes/extra/Iter001/volume02.vol avgs=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/012458_XmippProtCropResizeParticles/extra/output_images.stk avgSampling=1.24 symmetry=o particles=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/010450_XmippProtReconstructHighRes/particles.sqlite ptclSampling=0.74 kV=300 Cs=2.7 Q0=0.1 hasAngles=yes micrographs="/home/coss/ScipionUserData/projects/Example_10248_Scipion3/Runs/006458_XmippProtMovieCorr/extra/*mic.mrc" micSampling=0.49 atomicModel=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/centered4V1W.pdb doMultimodel=yes workflow=http://nolan.cnb.csic.es/cryoemworkflowviewer/workflow/637ca2bbcd57e45e88f6fabb7f6b1095a3ca0de6 saxs=/home/coss/ScipionUserData/projects/Example_10248_Scipion3/SASDE55.dat
    # ~/scipion3/scipion3 python ~/data/Dropbox/H/scipion-em-validation/validationLevels.py project=Validation11337 map=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB11337/emd_11337.map sampling=1.047 threshold=0.165 resolution=3.3 symmetry=c1  atomicModel=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB11337/6zp7_updated.pdb doMultimodel=no map1=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB11337/emd_11337_half_map_1.map map2=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB11337/emd_11337_half_map_2.map
    # ~/scipion3/scipion3 python ~/data/Dropbox/H/scipion-em-validation/validationLevels.py project=Validation22301 map=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB22301/emd_22301.map sampling=0.52 threshold=0.1 resolution=3.7 symmetry=c1  atomicModel=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB22301/6xs6_updated.cif doMultimodel=no
    # ~/scipion3/scipion3 python ~/data/Dropbox/H/scipion-em-validation/validationLevels.py project=Validation22838 map=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB22838/emd_22838.map sampling=1.058 threshold=0.2 resolution=3.84 symmetry=c1  atomicModel=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB22838/7kec_updated.pdb doMultimodel=no
    # ~/scipion3/scipion3 python ~/data/Dropbox/H/scipion-em-validation/validationLevels.py project=Validation11668 map=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB11668/emd_11668.map sampling=0.492 threshold=0.15 resolution=1.15 symmetry=o  atomicModel=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB11668/7a6a_updated.cif doMultimodel=no
    # ~/scipion3/scipion3 python ~/data/Dropbox/H/scipion-em-validation/validationLevels.py project=Validation11668_05 map=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB11668/emd_11668.map sampling=0.492 threshold=0.05 resolution=1.15 symmetry=o  atomicModel=/home/coss/data/Dropbox/Aplicaciones/ShareLaTeX/MapValidation/EMDB11668/7a6a_updated.cif doMultimodel=no

if any(i in sys.argv for i in ['-h', '-help', '--help', 'help']):
    usage()

# Manager will help as to find paths, create project...
manager = Manager()

# Fixing some parameters depending on the arguments or taking the default ones
PROJECT_NAME = None
PATH_UPLOADS = None
PATH_LOGS = None
FNMAP = None
FN_METADATA = None
TS = None
MAPTHRESHOLD = None
MAPRESOLUTION = None
HASANGLES = False
doMultimodel = False
FNMODEL = None
JOB_NAME = None
JOB_DESCRIPTION= None
EMDB_ID = None
PDB_ID = None
IS_EMDB_ENTRY = False

MAPCOORDX = 0
MAPCOORDY = 0
MAPCOORDZ = 0
FNMAP1 = None
FNMAP2 = None
XLM = None
SAXS = None

UNTILTEDMIC = None
TILTEDMIC = None
TILTKV = None
TILTCS = None
TILTQ0 = None
TILTTS = None
TILTANGLE = None
UNTILTEDCOORDS = None
TILTEDCOORDS = None

def readParamsFromFile(filename):
    try:
        with open(filename,mode="r") as f:
            jdata = json.load(f)

        global PROJECT_NAME, PATH_UPLOADS, PATH_LOGS, FNMAP, FN_METADATA, TS, MAPTHRESHOLD, MAPRESOLUTION, MAPCOORDX, MAPCOORDY, MAPCOORDZ, FNMAP1, FNMAP2, HASANGLES, doMultimodel, FNMODEL

        IS_EMDB_ENTRY = jdata["isEMDBentry"]
        EMDB_ID = jdata["emdbId"]
        PDB_ID = jdata["pdbId"]
        JOB_NAME = jdata["jobName"]
        JOB_DESCRIPTION = jdata["jobDescription"]
        PROJECT_NAME = jdata["emdbId"]
        PATH_UPLOADS = jdata["jobUploadsDir"]
        PATH_LOGS = jdata["jobLogsDir"]
        FNMAP = os.path.join(jdata["jobUploadsDir"], jdata["mapFileName"]) 
        FN_METADATA = os.path.join(jdata["jobUploadsDir"], jdata["metaFileName"])
        TS =  jdata["sampling"]
        MAPTHRESHOLD =  jdata["threshold"]
        MAPRESOLUTION =  jdata["resolution"]
        MAPCOORDX =  jdata["mapCoordX"]
        MAPCOORDY =  jdata["mapCoordY"]
        MAPCOORDZ =  jdata["mapCoordZ"]
        FNMAP1 =  os.path.join(jdata["jobUploadsDir"], jdata["map1"])
        FNMAP2 =  os.path.join(jdata["jobUploadsDir"], jdata["map2"])
        HASANGLES =  jdata["hasAngles"]
        doMultimodel =  jdata["doMultimodel"]        
        FNMODEL =  os.path.join(jdata["jobUploadsDir"], jdata["modelFileName"])
    except Exception as ex:
        print("ERROR:", ex)

levels = []

for arg in sys.argv:
    # if this param is present, the rest of params will be readed from the metadata file
    # there should not be more params in the command line
    if arg.startswith('inputParams'):
        FN_PARAMS = arg.split('inputParams=')[1]
        readParamsFromFile(FN_PARAMS)
        if FNMAP:
            levels.append("0")
        if FNMAP and FNMODEL:
            levels.append("A")

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
    if arg.startswith('mapCoordX='):
        MAPCOORDX = float(arg.split("mapCoordX=")[1])
    if arg.startswith('mapCoordY='):
        MAPCOORDY = float(arg.split("mapCoordY=")[1])
    if arg.startswith('mapCoordZ='):
        MAPCOORDZ = float(arg.split("mapCoordZ=")[1])
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
    if arg.startswith('micrographs='):
        MICPATTERN = arg.split("micrographs=")[1]
    if arg.startswith('micSampling='):
        TSMIC = float(arg.split("micSampling=")[1])
    if arg.startswith('atomicModel='):
        FNMODEL = arg.split("atomicModel=")[1]
    if arg.startswith('doMultimodel='):
        doMultimodel = arg.split("doMultimodel=")[1]=="yes"
    if arg.startswith('workflow='):
        WORKFLOW = arg.split("workflow=")[1]
    if arg.startswith('xlm='):
        XLM = arg.split("xlm=")[1]
    if arg.startswith('saxs='):
        SAXS = arg.split("saxs=")[1]
    if arg.startswith('untiltedMic='):
        UNTILTEDMIC = arg.split("untiltedMic=")[1]
    if arg.startswith('tiltedMic='):
        TILTEDMIC = arg.split("tiltedMic=")[1]
    if arg.startswith('tiltkV='):
        TILTKV = float(arg.split("tiltkV=")[1])
    if arg.startswith('tiltCs='):
        TILTCS = float(arg.split("tiltCs=")[1])
    if arg.startswith('tiltQ0='):
        TILTQ0 = float(arg.split("tiltQ0=")[1])
    if arg.startswith('tiltSampling='):
        TILTTS = float(arg.split("tiltSampling=")[1])
    if arg.startswith('tiltAngle='):
        TILTANGLE = float(arg.split("tiltAngle=")[1])
    if arg.startswith('untiltedCoords='):
        UNTILTEDCOORDS = arg.split("untiltedCoords=")[1]
    if arg.startswith('tiltedCoords='):
        TILTEDCOORDS = arg.split("tiltedCoords=")[1]

# Detect level
argsPresent = [x.split('=')[0] for x in sys.argv]
LEVEL0 = ["map", "sampling", "threshold", "resolution"]
LEVEL1 = ["map1", "map2"]
LEVEL2 = ["avgs", "avgSampling", "symmetry"]
LEVEL3 = ["particles", "ptclSampling", "kV", "Cs", "Q0"]
LEVEL4 = ["hasAngles"]
LEVEL5 = ["micrographs", "micSampling"]
LEVELA = ["atomicModel", "doMultimodel"]
LEVELW = ["workflow"]
LEVELOa = ["xlm"]
LEVELOb = ["saxs"]
LEVELOc = ["untiltedMic","tiltedMic","tiltkV","tiltCs","tiltQ0","tiltSampling","tiltAngle","untiltedCoords",
           "tiltedCoords"]

def detectLevel(labels, args):
    retval = True
    for label in labels:
        if not label in args:
            retval = False
            break
    return retval
if detectLevel(LEVEL0, argsPresent):
    levels.append("0")
if detectLevel(LEVEL1, argsPresent) and detectLevel(LEVEL0, argsPresent):
    levels.append("1")
if detectLevel(LEVEL2, argsPresent) and detectLevel(LEVEL0, argsPresent):
    levels.append("2")
if detectLevel(LEVEL3, argsPresent) and detectLevel(LEVEL2, argsPresent):
    levels.append("3")
if detectLevel(LEVEL4, argsPresent) and detectLevel(LEVEL3, argsPresent):
    levels.append("4")
if detectLevel(LEVEL5, argsPresent) and detectLevel(LEVEL3, argsPresent):
    levels.append("5")
if detectLevel(LEVELA, argsPresent) and detectLevel(LEVEL0, argsPresent):
    levels.append("A")
if detectLevel(LEVELW, argsPresent):
    levels.append("W")
if detectLevel(LEVELOa, argsPresent) and detectLevel(LEVELA, argsPresent):
    levels.append("O")
if detectLevel(LEVELOb, argsPresent) and detectLevel(LEVELA, argsPresent): #TODO tiene sentido que chequee el nivel A cuando este nivel no depende del atomic model?? Mejor detectLevel(LEVEL0, ...)?
    if not "O" in levels:
        levels.append("O")
if detectLevel(LEVELOc, argsPresent) and detectLevel(LEVELA, argsPresent): #TODO tiene sentido que chequee el nivel A cuando este nivel no depende del atomic model?? Mejor detectLevel(LEVEL0, ...)?
    if not "O" in levels:
        levels.append("O")

if len(levels)==0 or not "0" in levels:
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

# Validate inputs formats
wrongInputs = {'errors':[], 'warnings':[]}
# check 'map' arg
fnDir, fnBase = os.path.split(FNMAP)
protImportMapChecker = project.newProtocol(pwplugin.Domain.importFromPlugin('pwem.protocols', 'ProtImportVolumes', doRaise=True),
                                           objLabel='check format - import map',
                                           filesPath=os.path.join(fnDir,FNMAP),
                                           samplingRate=TS,
                                           setOrigCoord=True,
                                           x=MAPCOORDX,
                                           y=MAPCOORDY,
                                           z=MAPCOORDZ)
if useSlurm:
    sendToSlurm(protImportMapChecker)
project.launchProtocol(protImportMapChecker)
#waitOutput(project, protImportMapChecker, 'outputVolume')
waitUntilFinishes(project, protImportMapChecker)
if protImportMapChecker.isFailed():
    wrongInputs['errors'].append({'param': 'map', 'value': FNMAP, 'cause': 'There is a problem reading the volume map file'})

else:
    # check if we can have a proper mask with the thresold specified
    protCreateMaskChecker = project.newProtocol(pwplugin.Domain.importFromPlugin('xmipp3.protocols.protocol_preprocess', 'XmippProtCreateMask3D', doRaise=True),
                                                objLabel='check proper mask',
                                                inputVolume=protImportMapChecker.outputVolume,
                                                threshold=MAPTHRESHOLD,
                                                doBig=True,
                                                doMorphological=True,
                                                elementSize=math.ceil(2/TS)) # Dilation by 2A
    if useSlurm:
        sendToSlurm(protCreateMaskChecker)
    project.launchProtocol(protCreateMaskChecker)
    waitUntilFinishes(project, protCreateMaskChecker)

    M = readMap(protCreateMaskChecker.outputMask.getFileName()).getData()
    totalMass = np.sum(M)
    if not totalMass > 0:
        wrongInputs['errors'].append({'param': 'threshold', 'value': MAPTHRESHOLD, 'cause': 'The mask obtained from the volume map is empty, try to lower the threshold value'})

if "1" in levels:
    # check 'map1' and 'map2' arg
    # 'map1'
    fnDir, fnBase = os.path.split(FNMAP1)
    protImportMap1Checker = project.newProtocol(pwplugin.Domain.importFromPlugin('pwem.protocols', 'ProtImportVolumes', doRaise=True),
                                                objLabel='check format - import half1',
                                                filesPath=fnDir,
                                                filesPattern=FNMAP1,
                                                samplingRate=TS,
                                                setOrigCoord=True,
                                                x=MAPCOORDX,
                                                y=MAPCOORDY,
                                                z=MAPCOORDZ)

    if useSlurm:
        sendToSlurm(protImportMap1Checker)
    project.launchProtocol(protImportMap1Checker)
    #waitOutput(project, protImportMap1Checker, 'outputVolume')
    waitUntilFinishes(project, protImportMap1Checker)
    if protImportMap1Checker.isFailed():
        wrongInputs['errors'].append({'param': 'map1', 'value': FNMAP1, 'cause': 'There is a problem reading the half-map1 file'})

    # 'map2'
    fnDir, fnBase = os.path.split(FNMAP2)
    protImportMap2Checker = project.newProtocol(pwplugin.Domain.importFromPlugin('pwem.protocols', 'ProtImportVolumes', doRaise=True),
                                                objLabel='check format - import half2',
                                                filesPath=fnDir,
                                                filesPattern=FNMAP2,
                                                samplingRate=TS,
                                                setOrigCoord=True,
                                                x=MAPCOORDX,
                                                y=MAPCOORDY,
                                                z=MAPCOORDZ)

    if useSlurm:
        sendToSlurm(protImportMap2Checker)
    project.launchProtocol(protImportMap2Checker)
    #waitOutput(project, protImportMap2Checker, 'outputVolume')
    waitUntilFinishes(project, protImportMap2Checker)
    if protImportMap2Checker.isFailed():
        wrongInputs['errors'].append({'param': 'map2', 'value': FNMAP2, 'cause': 'There is a problem reading the half-map2 file'})

if "2" in levels:
    # Check 'avgs' arg
    protImportAvgsChecker = project.newProtocol(pwplugin.Domain.importFromPlugin('pwem.protocols', 'ProtImportAverages', doRaise=True),
                                                objLabel='check format - import averages',
                                                filesPath=FNAVGS,
                                                samplingRate=TSAVG)
    if useSlurm:
        sendToSlurm(protImportAvgsChecker)
    project.launchProtocol(protImportAvgsChecker)
    #waitOutput(project, protImportAvgsChecker, 'outputAverages')
    waitUntilFinishes(project, protImportAvgsChecker)
    if protImportAvgsChecker.isFailed():
        wrongInputs['errors'].append({'param': 'avgs', 'value': FNAVGS, 'cause': 'There is a problem reading the 2D Classes file'})

if "3" in levels:
    # Check 'particles' arg
    protImportParticlesChecker = project.newProtocol(pwplugin.Domain.importFromPlugin('pwem.protocols', 'ProtImportParticles', doRaise=True),
                                                     objLabel='check format - import particles',
                                                     filesPath=FNPARTICLES,
                                                     samplingRate=TSPARTICLES,
                                                     voltage=KV,
                                                     sphericalAberration=CS,
                                                     amplitudeContrast=Q0)
    if FNPARTICLES.endswith(".sqlite"):
        protImportParticlesChecker.importFrom.set(protImportParticlesChecker.IMPORT_FROM_SCIPION)
        protImportParticlesChecker.sqliteFile.set(FNPARTICLES)
    elif FNPARTICLES.endswith(".xmd"):
        protImportParticlesChecker.importFrom.set(protImportParticlesChecker.IMPORT_FROM_XMIPP)
        protImportParticlesChecker.mdFile.set(FNPARTICLES)
    elif FNPARTICLES.endswith(".star"):
        protImportParticlesChecker.importFrom.set(protImportParticlesChecker.IMPORT_FROM_RELION)
        protImportParticlesChecker.starFile.set(FNPARTICLES)
    if useSlurm:
        sendToSlurm(protImportParticlesChecker)
    project.launchProtocol(protImportParticlesChecker)
    #waitOutput(project, protImportParticlesChecker, 'outputParticles')
    waitUntilFinishes(project, protImportParticlesChecker)
    # check if particles has alignment
    if protImportParticlesChecker.isFailed() or not protImportParticlesChecker.outputParticles.hasAlignment():
        wrongInputs['errors'].append({'param': 'particles', 'value': FNPARTICLES, 'cause': 'There is a problem reading the particles file'})

if "5" in levels:
    # Check 'micrographs' arg
    protImportMicrographsChecker = project.newProtocol(pwplugin.Domain.importFromPlugin('pwem.protocols', 'ProtImportMicrographs', doRaise=True),
                                     objLabel='check format - import mics',
                                     samplingRate=TSMIC,
                                     voltage=KV,
                                     sphericalAberration=CS,
                                     amplitudeContrast=Q0)
    if MICPATTERN.endswith(".sqlite"):
        protImportMicrographsChecker.importFrom.set(protImportMicrographsChecker.IMPORT_FROM_SCIPION)
        protImportMicrographsChecker.sqliteFile.set(MICPATTERN)
    else:
        protImportMicrographsChecker.filesPattern.set(MICPATTERN)
    if useSlurm:
        sendToSlurm(protImportMicrographsChecker)
    project.launchProtocol(protImportMicrographsChecker)
    #waitOutput(project, protImportMicrographsChecker, 'outputMicrographs')
    waitUntilFinishes(project, protImportMicrographsChecker)
    if protImportMicrographsChecker.isFailed():
        wrongInputs['errors'].append({'param': 'micrographs', 'value': MICPATTERN, 'cause': 'There is a problem reading the micrographs file'})

if "A" in levels and not protImportMapChecker.isFailed():
    # Check 'atomicModel' arg

    writeAtomicModelFailed = False
    try: # Check if biopython can read atomic model file
        h = UpdatedAtomicStructHandler()
        h.read(FNMODEL)
        try:  # Check if biopython can convert atomic file to PDB
            # Get structure ID
            structure_id = os.path.basename(FNMODEL)
            structure_id = structure_id[:4] if len(structure_id) > 4 else "1xxx"
            pdbFile = '%s.pdb' % (structure_id)
            
            # Get tmp pdb  from imput atomic model to work on
            fnPdb = os.path.join(project.getTmpPath(), pdbFile) #TODO: poner en otra carpet??
            h.writeAsPdb(fnPdb)
        except OutOfChainsError:
            wrongInputs['warnings'].append({'param': 'atomicModel', 'value': FNMODEL, 'cause': 'Atomic model file not valid. Some programs cannot handle it due to size: Too many chains to represent in PDB format'})
            writeAtomicModelFailed = True
        except OutOfAtomsError:
            wrongInputs['warnings'].append({'param': 'atomicModel', 'value': FNMODEL, 'cause': 'Atomic model file not valid. Some programs cannot handle it due to size: Too many atoms to represent in PDB format'})
            writeAtomicModelFailed = True
        except:
            wrongInputs['warnings'].append({'param': 'atomicModel', 'value': FNMODEL, 'cause': 'Atomic model file not valid. Some programs cannot handle it because it cannot be safely written in PDB format'})
            writeAtomicModelFailed = True
    except:
        wrongInputs['error'].append({'param': 'atomicModel', 'value': FNMODEL, 'cause': 'Atomic model file not valid. It cannot be safely read'})

    if not writeAtomicModelFailed:
        protImportAtomicModelChecker = project.newProtocol(pwplugin.Domain.importFromPlugin('pwem.protocols', 'ProtImportPdb', doRaise=True),
                                                        objLabel='check format - import atomic',
                                                        inputPdbData=1,
                                                        pdbFile=fnPdb)
        protImportAtomicModelChecker.inputVolume.set(protImportMapChecker.outputVolume)
        if useSlurm:
            sendToSlurm(protImportAtomicModelChecker)
        project.launchProtocol(protImportAtomicModelChecker)
        #waitOutput(project, protImportAtomicModelChecker, 'outputPdb')
        waitUntilFinishes(project, protImportAtomicModelChecker)
        if protImportAtomicModelChecker.isFailed():
            wrongInputs['errors'].append({'param': 'atomicModel', 'value': fnPdb, 'cause': 'There is a problem reading the atomic model file'})


if "O" in levels and not protImportMapChecker.isFailed():
    # Check 'xlm', 'saxs', 'untiltedMic', 'tiltedMic', 'untiltedCoords', 'tiltedCoords' args
    # 'xlm'
    if "A" in levels and not writeAtomicModelFailed and not protImportAtomicModelChecker.isFailed():  #TODO: add  ... and XLM is not None:
        protImportXLMChecker = project.newProtocol(pwplugin.Domain.importFromPlugin('xlmtools.protocols', 'ProtWLM', doRaise=True),
                                                   objLabel="check format - XLM",
                                                   xlList=XLM)
        protImportXLMChecker.pdbs.set([protImportAtomicModelChecker.outputPdb])
        if useSlurm:
            sendToSlurm(protImportXLMChecker)
        project.launchProtocol(protImportXLMChecker)
        #waitOutput(project, protImportXLMChecker, 'crosslinkStruct_1')
        waitUntilFinishes(project, protImportXLMChecker)

        if protImportXLMChecker.isFailed():
            wrongInputs['errors'].append({'param': 'xlm', 'value': XLM, 'cause': 'There is a problem reading the XML file'})
    #TODO: por que se ejecutan todos los checkers si yo solo he especificado los argumentos del nivel O.a?? No deberia puesto que esta fallando el job porque ejecuta los chequers de postseudo sin tener suficientes parametros (no se especifican en el comando de lanzamiento). REVISAR!!!!!!!!!! 
    # 'sax'
    #TODO: Add 'if SAXS is not None:'
    protCreateMask = project.newProtocol(pwplugin.Domain.importFromPlugin('xmipp3.protocols.protocol_preprocess', 'XmippProtCreateMask3D', doRaise=True),
                                         objLabel='check format - create mask',
                                         inputVolume=protImportMapChecker.outputVolume,
                                         threshold=MAPTHRESHOLD,
                                         doBig=True,
                                         doMorphological=True,
                                         elementSize=math.ceil(2/TS)) # Dilation by 2A
    if useSlurm:
        sendToSlurm(protCreateMask)
    project.launchProtocol(protCreateMask)
    #waitOutput(project, protCreateMask, 'outputMask')
    waitUntilFinishes(project, protCreateMask)


    protPseudo = project.newProtocol(pwplugin.Domain.importFromPlugin('continuousflex.protocols', 'FlexProtConvertToPseudoAtoms', doRaise=True),
                                     objLabel="check format - convert Map to Pseudo",
                                     maskMode=2,
                                     pseudoAtomRadius=1.5)
    protPseudo.inputStructure.set(protImportMapChecker.outputVolume)
    protPseudo.volumeMask.set(protCreateMask.outputMask)
    if useSlurm:
        sendToSlurm(protPseudo)
    project.launchProtocol(protPseudo)
    #waitOutput(project, protPseudo, 'outputVolume')
    #waitOutput(project, protPseudo, 'outputPdb')
    waitUntilFinishes(project, protPseudo)


    protImportSaxsChecker = project.newProtocol(pwplugin.Domain.importFromPlugin('atsas.protocols',
                                                                                 'AtsasProtConvertPdbToSAXS', doRaise=True),
                                                objLabel="check format - SAXS",
                                                experimentalSAXS=SAXS)
    protImportSaxsChecker.inputStructure.set(protPseudo.outputPdb)
    if useSlurm:
        sendToSlurm(protImportSaxsChecker)
    project.launchProtocol(protImportSaxsChecker)
    if protImportSaxsChecker.isFailed():
        wrongInputs['errors'].append({'param': 'saxs', 'value': SAXS, 'cause': 'There is a problem reading the SAXS file'})

    # 'untiltedMic' and 'tiltedMic'
    #TODO: Add 'if not [x for x in (UNTILTEDMIC, TILTEDMIC, TILTKV, TILTCS, TILTQ0, TILTTS, TILTANGLE, UNTILTEDCOORDS, TILTEDCOORDS) if x is None]: # Checks that none of the variables are None'
    protImportTiltPairsChecker = project.newProtocol(pwplugin.Domain.importFromPlugin('pwem.protocols', 'ProtImportMicrographsTiltPairs', doRaise=True),
                                                     objLabel="check format - import tilt pairs",
                                                     patternUntilted=UNTILTEDMIC,
                                                     patternTilted=TILTEDMIC,
                                                     voltage=TILTKV,
                                                     ampContrast=TILTQ0,
                                                     sphericalAberration=TILTCS,
                                                     samplingRate=TILTTS)
    if useSlurm:
        sendToSlurm(protImportTiltPairsChecker)
    project.launchProtocol(protImportTiltPairsChecker)
    #waitOutput(project, protImportTiltPairsChecker, 'outputMicrographsTiltPair')
    waitUntilFinishes(project, protImportTiltPairsChecker)

    if protImportTiltPairsChecker.isFailed():
        wrongInputs['errors'].append( {'param': 'untiltedMic', 'value': UNTILTEDMIC, 'cause': 'There is a problem reading the untilted mic file'})
        wrongInputs['errors'].append({'param': 'tiltedMic', 'value': TILTEDMIC, 'cause': 'There is a problem reading the tilted mic file'})

    # 'untiltedCoords' and 'tiltedCoords'
    x, y, z = protImportMapChecker.outputVolume.getDimensions()
    Ts = protImportMapChecker.outputVolume.getSamplingRate()
    dMap = x * Ts
    boxSize = int(dMap / TILTTS)
    protImportCoordsChecker = project.newProtocol(pwplugin.Domain.importFromPlugin('pwem.protocols', 'ProtImportCoordinatesPairs', doRaise=True),
                                                  objLabel="check format - import paired coordinates",
                                                  patternUntilted=UNTILTEDCOORDS,
                                                  patternTilted=TILTEDCOORDS,
                                                  boxSize=boxSize)
    if UNTILTEDCOORDS.endswith('.json'):
        protImportCoordsChecker.importFrom.set(1)
    protImportCoordsChecker.inputMicrographsTiltedPair.set(protImportTiltPairsChecker.outputMicrographsTiltPair)
    if useSlurm:
        sendToSlurm(protImportCoordsChecker)
    project.launchProtocol(protImportCoordsChecker)
    #waitOutput(project, protImportCoordsChecker, 'outputCoordinatesTiltPair')
    waitUntilFinishes(project, protImportCoordsChecker)
    if protImportCoordsChecker.isFailed():
        wrongInputs['errors'].append({'param': 'untiltedCoords', 'value': UNTILTEDCOORDS, 'cause': 'There is a problem reading the untilted coords file'})
        wrongInputs['errors'].append({'param': 'tiltedCoords', 'value': TILTEDCOORDS, 'cause': 'There is a problem reading the tilted coords file'})

with open (os.path.join(report.fnReportDir, 'wrongInputs.json'), 'w') as f:
        json.dump(wrongInputs, f)
# if some input data was wrong do whatever we want: inform the user, write error msg in report, etc.
#if protImportMapChecker.isFailed() or protImportMap1Checker.isFailed() or protImportMap2Checker.isFailed() or \
if protImportMapChecker.isFailed() or \
        (protImportMap1Checker.isFailed() if "1" in levels and 'protImportMap1Checker' in locals() else None) or \
        (protImportMap2Checker.isFailed() if "1" in levels and 'protImportMap2Checker' in locals() else None) or \
        (protImportAvgsChecker.isFailed() if "2" in levels and 'protImportAvgsChecker' in locals() else None) or \
        (protImportParticlesChecker.isFailed() if "3" in levels and 'protImportParticlesChecker' in locals() else None) or \
        (protImportMicrographsChecker.isFailed() if "5" in levels and 'protImportMicrographsChecker' in locals() else None) or \
        (protImportAtomicModelChecker.isFailed() if "A" in levels and 'protImportAtomicModelChecker' in locals() else None) or \
        (protImportXLMChecker.isFailed() if "O" in levels and 'protImportXLMChecker' in locals() else None) or \
        (protImportSaxsChecker.isFailed() if "O" in levels and 'protImportSaxsChecker' in locals() else None) or \
        (protImportTiltPairsChecker.isFailed() if "O" in levels and 'protImportTiltPairsChecker' in locals() else None) or \
        (protImportCoordsChecker.isFailed() if "O" in levels and 'protImportCoordsChecker' in locals() else None) or \
        len(wrongInputs['errors'])>0:
    print("Some input data was not correct")

else: # go ahead
    print("All inputs were correct, let's process them!")
    # Level 0
    from validationLevel0 import level0
    protImportMap, protCreateMask, bfactor, protResizeMap, protResizeMask = level0(
        project, report, FNMAP, FNMAP1, FNMAP2, TS, MAPTHRESHOLD, MAPRESOLUTION, MAPCOORDX, MAPCOORDY, MAPCOORDZ, skipAnalysis = False)

    # Level 1
    if "1" in levels:
        from validationLevel1 import level1
        protImportMap1, protImportMap2 = level1(project, report, FNMAP1, FNMAP2, TS, MAPRESOLUTION, MAPCOORDX, MAPCOORDY, MAPCOORDZ,
                                                protImportMap, protResizeMap, protCreateMask, protResizeMask,
                                                skipAnalysis = False)

    # Level 2
    if "2" in levels:
        from validationLevel2 import level2
        protImportAvgs, protAvgsResizeMap = level2(project, report, protImportMap, FNAVGS, TSAVG, SYM, skipAnalysis = False)

    # Level 3
    if "3" in levels:
        from validationLevel3 import level3
        protImportParticles, protResizeParticlesMap, protResizeAvgs = level3(project, report, protImportMap, protImportAvgs,
                                                                    FNPARTICLES, TSPARTICLES, KV, CS, Q0,
                                                                    skipAnalysis = False)

    # Level 4
    if "4" in levels:
        from validationLevel4 import level4
        protResizeParticles = level4(project, report, protImportMap, protCreateMask, protResizeParticlesMap, SYM,
                                     MAPRESOLUTION, bfactor, protResizeMap, protResizeMask, skipAnalysis = False)

    # Level 5
    if "5" in levels:
        from validationLevel5 import level5
        level5(project, report, protImportParticles, KV, CS, Q0, MICPATTERN, TSMIC, skipAnalysis = False)

    # Level A
    #TODO: pass writeAtomicModelFailed to levelA() to write the warning in the report
    if "A" in levels:
        from validationLevelA import levelA
        protAtom = levelA(project, report, protImportMap, FNMODEL, fnPdb, writeAtomicModelFailed, MAPRESOLUTION, doMultimodel, MAPCOORDX, MAPCOORDY, MAPCOORDZ, skipAnalysis = False)
    else:
        protAtom = None

    # Level W
    if "W" in levels:
        from validationLevelW import levelW
        levelW(project, report, WORKFLOW, skipAnalysis = False)

    # Level O
    #TODO: pass writeAtomicModelFailed to levelO() to write the warning in the report (section O.a.)
    if "O" in levels:
        from validationLevelO import levelO
        levelO(project, report, protImportMap, protCreateMask, protAtom, XLM, SAXS,
               UNTILTEDMIC, TILTEDMIC, TILTKV, TILTCS, TILTQ0, TILTTS, TILTANGLE, UNTILTEDCOORDS, TILTEDCOORDS, SYM,
               skipAnalysis = False)

    # Close report
    report.abstractResolution(MAPRESOLUTION)
    # Check if there are warnings
    with open(os.path.join(fnProjectDir, "validationReport", "summaryWarnings.tex")) as summaryWarnings:
        content = summaryWarnings.read()
        if "STATUS" not in content and "WARNINGS" not in content:
            report.fhSummaryWarnings.write("No warnings.")
    report.closeReport()