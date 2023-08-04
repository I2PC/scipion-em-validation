#!/usr/bin/env python
import os
import sys
import re
import argparse
import logging
from csv import reader
from utils import *


def downloadEmdbFile(url, fileName, localPath, gunzip=False):
    """
    Download halfmaps & additional map files from EMDB
    """
    localFileName = fileName.replace('.gz', '') if gunzip else fileName
    fullFileName = os.path.join(localPath, localFileName)
    if not os.path.exists(fullFileName) or os.stat(fullFileName).st_size == 0:
        logger.info('Getting %s' % localFileName)
        fileName = downloadFile(url, localPath, fileName, raw=True)
        if gunzip:
            logger.info('Uncompress gziped file %s' % fileName)
            fileName = ungzipFile(os.path.join(
                localPath, fileName), fullFileName)
    else:
        logger.warning('File %s already exists' % localFileName)
    assert fileName, "ERROR: %s file not found" % fileName
    return fileName


def downloadPdbModel(pdbId, workDir, force=False):
    """
    Download model file from PDB
    """
    # pdbFileName = pdbId + '.pdb'
    # destFile = os.path.join(workDir, pdbFileName)
    cifFileName = pdbId + '.cif'
    destFile = os.path.join(workDir, cifFileName)
    try:
        if force or not os.path.exists(destFile) or os.stat(destFile).st_size == 0:
            # download PDB model file
            # url = PDB_REPOSITORY + fileName
            # https://files.rcsb.org/download/7BB7.pdb
            # url = URL_PDB_RCSB_REPOSITORY + pdbId + '.pdb'
            # pdbFileName = downloadFile(url, workDir, pdbFileName)

            # download mmCIF model file
            # https://files.rcsb.org/download/7BB7.cif
            url = URL_PDB_RCSB_REPOSITORY + cifFileName
            logger.info('Getting %s' % cifFileName)
            destFile = downloadFile(url, workDir, cifFileName)

        else:
            logger.warning('File %s already exists' % destFile)
    except Exception as ex:
        logger.warning(ex)
        return None
    # assert destFile, "ERROR: %s atomic model not found" % destFile
    return destFile


def getJobRunCommand(inputParams):
    """
        compose a command to run a job with 'inputDataObj'
    """
    # script -c '
    cmmd = "script -c \'"
    # export DISPLAY=:1.0
    cmmd = cmmd + SERVER_DISPLAY_SET
    # /home/vrs/scipion3/scipion3 python
    cmmd = cmmd + ' && ' + PATH_PYTHON + ' ' + CMD_PYTHON
    # /home/vrs/scipion3/scipion-em-validation/validationLevels.py
    cmmd = cmmd + ' ' + PATH_SCRIPT + CMD_SCRIPT

    # EMDBentryId=EMD-35930
    params = "EMDBid=" + inputParams["emdbId"] + "'"

    logging = os.path.join(
        inputParams["jobLogsDir"], "runJob_%s.log" % inputParams["emdbId"])

    cmmdLine = cmmd + " " + params + " -f " + logging
    return cmmdLine


def getJobRunCommand_old(inputParams):
    """
        compose a command to run a job with 'inputDataObj'
    """

    uploadsPath = inputParams["jobUploadsDir"]
    # script -c '
    cmmd = "script -c \'"
    # export DISPLAY=:1.0
    cmmd = cmmd + SERVER_DISPLAY_SET
    # /home/vrs/scipion3/scipion3 python
    cmmd = cmmd + ' && ' + PATH_PYTHON + ' ' + CMD_PYTHON
    # /home/vrs/scipion3/scipion-em-validation/validationLevels.py
    cmmd = cmmd + ' ' + PATH_SCRIPT + CMD_SCRIPT
    # project=' + prjName
    params = 'project=' + inputParams["emdb_id"]
    # map=/home/vrs/ScipionUserData/projects/679f4c1e128e/map.mrc
    params = params + ' ' + 'map=' + \
        os.path.join(uploadsPath, inputParams["mapFileName"])
    # sampling=0.74
    params = params + ' ' + 'sampling=' + str(inputParams["sampling"])
    # threshold=0.0025
    params = params + ' ' + 'threshold=' + str(inputParams["threshold"])
    # resolution=2.6
    params = params + ' ' + 'resolution=' + str(inputParams["resolution"])
    # mapCoordX=0 mapCoordY=0 mapCoordZ=0
    params = params + ' ' + 'mapCoordX=' + str(inputParams["mapCoordX"]) + ' ' + \
        'mapCoordY=' + str(inputParams["mapCoordY"]) + \
        ' ' + 'mapCoordZ=' + str(inputParams["mapCoordZ"])
    # map1=/home/vrs/ScipionUserData/projects/679f4c1e128e/volume01.vol
    # map2=/home/vrs/ScipionUserData/projects/679f4c1e128e/volume02.vol
    if "map1" in inputParams.keys():
        params = params + ' ' + 'map1=' + \
            os.path.join(uploadsPath, inputParams["map1"]) + ' ' + \
            'map2=' + os.path.join(uploadsPath, inputParams["map2"])
    # "hasAngles": False,
    params = params + ' ' + 'hasAngles=' + str(False)
    # atomicModel=/home/vrs/ScipionUserData/projects/Example_10248_Scipion3/pdb8dm8.ent
    if "modelFileName" in inputParams.keys():
        params = params + ' ' + 'atomicModel=' + \
            os.path.join(uploadsPath,  inputParams["modelFileName"])
    # doMultimodel=False
    params = params + ' ' + 'doMultimodel=' + str(False)

    params = params + "\'"

    logging = inputParams["jobLogsDir"] + \
        "job_%s.log" % inputParams["emdb_id"]

    cmmdLine = cmmd + ' ' + params + ' -f ' + logging
    return cmmdLine


def main(argv):

    parser = argparse.ArgumentParser(description="Prepare Job to Run")
    parser.add_argument(
        "-m", "--map", help="Volume map EMDB ID (EMD-1234 or EMD-12345)", required=True)
    parser.add_argument(
        "-f", "--force-download", help="download files without checking if they already exist", required=False, action='store_true')
    parser.add_argument(
        "-l", "--logFile", help="log file. By default 'prepareJob.log' in a dedicated 'logs' folder.", required=False)
    parser.add_argument(
        "-t", "--test", help="perform a trial run with no changes made", required=False, action='store_true')
    args = parser.parse_args()
    if args.logFile:
        logFile = args.logFile
        logsDir = os.path.dirname(logFile)
        logFilename = os.path.basename(logFile)
    else:
        logsDir = os.path.join(PATH_TOOLS_DIR, DIR_TOOLS_LOGS)
        logFilename = os.path.join(
            logsDir, os.path.splitext(Path(__file__).name)[0] + '.log')
    logger = logSetup(__name__, logsDir, logFilename)

    force_download = args.force_download
    if force_download:
        logger.warning(
            'Files will be downloaded regardless of whether they already exist')
    test_only = args.test
    if test_only:
        logger.warning(
            'Performing a trial run. No permanent changes will be made')

    emdbId = args.map
    logger.info('Prepare to run Job %s' % emdbId)
    rgx = re.compile(r'^EMD-\d{4,5}$')
    if not rgx.match(emdbId):
        logger.error(
            "Not a valid EMDB ID %s. Should start with EMD- followed by 4-5 digits)" % emdbId)
        exit(1)

    logger.info('Processing %s' % emdbId)
    mapNum = emdbId.replace('EMD-', '')
    jobName = emdbId

    # Create (if not exists) a directory for the Job
    logger.info('Creating Job dir %s' % jobName)
    jobDir = createDir(PATH_WORK_DIR, jobName)

    # Get EMDB metadata
    jData = getEmdbMetadata(mapNum, jobDir)
    if not jData:
        removeDir(jobDir)
        logger.error('EMDB data not found for %s' % emdbId)
        exit(2)

    jobUploadsDir = os.path.join(PATH_APP, DIR_PROJ_UPLOADS, emdbId)
    jobLogsDir = os.path.join(PATH_APP, DIR_PROJ_LOGS)

    if "map" in jData:
        emdbId = jData["emdb_id"]
        gzMapFileName = jData["map"]["file"]
        mapFileName = gzMapFileName.replace('.gz', '')
        logger.info('Map file name: %s' % gzMapFileName)
        metaFileName = emdbId + ".json"
        logger.info('Metadata file name: %s' % metaFileName)
        size_x = int(jData["map"]["dimensions"]["col"])
        size_y = int(jData["map"]["dimensions"]["col"])
        size_z = int(jData["map"]["dimensions"]["col"])
        logger.info("Map size (number of grid points): %g x %g x %g" %
                    (size_x, size_y, size_z))

        org_x = float(jData["map"]["origin"]["col"])
        org_y = float(jData["map"]["origin"]["row"])
        org_z = float(jData["map"]["origin"]["sec"])
        logger.info('Map origin (x,y,z) : (%g, %g, %g)' %
                    (org_x, org_y, org_z))

        sampling_x = float(jData["map"]["pixel_spacing"]["x"]["valueOf_"])
        sampling_y = float(jData["map"]["pixel_spacing"]["y"]["valueOf_"])
        sampling_z = float(jData["map"]["pixel_spacing"]["z"]["valueOf_"])
        if sampling_x != sampling_y != sampling_z:
            logger.warning('pixel spacing values are different')
        logger.info('Map sampling (pixel_spacing / voxel size): %g x %g x %g Å' %
                    (sampling_x, sampling_y, sampling_z))
        threshold = float(jData["map"]["contour_list"]["contour"][0]["level"])
        logger.info('Map threshold (recommended contour level): %f' %
                    threshold)

    if "structure_determination_list" in jData:
        jData["structure_determination_list"]
        resolution = float(jData["structure_determination_list"]["structure_determination"][0]
                           ["image_processing"][0]["final_reconstruction"]["resolution"]["valueOf_"])
        logger.info('Map resolution: %f' % resolution)

    halfmapFiles = []
    additionalMapFiles = []
    maskFiles = []
    if "interpretation" in jData:
        if "half_map_list" in jData["interpretation"]:
            halfmaps = jData["interpretation"]["half_map_list"]["half_map"]
            for halfmapFile in halfmaps:
                halfmapFiles.append(halfmapFile["file"])
            logger.info('Half-maps: %s' % halfmapFiles)

        if "additional_map_list" in jData["interpretation"]:
            additionalMaps = jData["interpretation"]["additional_map_list"]["additional_map"]
            for additionalMap in additionalMaps:
                additionalMapFiles.append(additionalMap["file"])
            logger.info('Additional maps: $s' % additionalMapFiles)

        if "segmentation_list" in jData["interpretation"]:
            masks = jData["interpretation"]["segmentation_list"]["segmentation"]
            for mask in masks:
                maskFiles.append(mask["file"])
            logger.info('masks: %s' % maskFiles)

    atomModelIDs = []
    atomModelFiles = []
    if "pdb_list" in jData["crossreferences"]:
        atomModels = jData["crossreferences"]["pdb_list"]["pdb_reference"]
        for model in atomModels:
            atomModelIDs.append(model["pdb_id"])
        if len(atomModels) > 1:
            logger.warning("There are %d atomic models: %s" %
                           (len(atomModels), atomModelIDs))
        logger.info("Atomic model %s" % atomModelIDs[0])
        aModelId = atomModelIDs[0]

    # Save job params to a json file
    inputParams = {
        "isEMDBentry": True,
        "jobName": "",
        "jobDescription": "",
        "emdbId": emdbId,
        "jobUploadsDir": jobUploadsDir,
        "jobLogsDir": jobLogsDir,
        "mapFileName": mapFileName,
        "metaFileName": metaFileName,
        "sampling": sampling_x,
        "threshold": threshold,
        "resolution": resolution,
        "mapCoordX": round(abs(org_x) * sampling_x, 4),
        "mapCoordY": round(abs(org_y) * sampling_y, 4),
        "mapCoordZ": round(abs(org_z) * sampling_z, 4),
        "hasAngles": False,
        "doMultimodel": False
    }
    if halfmapFiles:
        inputParams["map1"] = halfmapFiles[0].replace('.gz', '')
        inputParams["map2"] = halfmapFiles[1].replace('.gz', '')

    inputParams["pdbId"] = aModelId

    save_json(data=inputParams, path=jobDir, filename=FN_LOCAL_JSON_PARAMS)

    # Save command line to a bash file (*.sh)
    cmmd = getJobRunCommand(inputParams)
    save2file(data=cmmd, path=jobDir, filename="runJob_%s.sh" %
              emdbId, append=False)


if __name__ == '__main__':
    main(sys.argv[1:])
