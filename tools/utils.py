#
# utility methods & functions
#
import json
import os
import os.path
import sys
import requests
import ftplib
import shutil
import gzip
import re
import logging
from io import BytesIO
from pathlib import Path

PARAM_FIELDS = ['map', 'sampling', 'threshold', 'resolution', 'mapCoordX', 'mapCoordY', 'mapCoordZ',
                'map1', 'map2',
                'avgs', 'avgSampling', 'symmetry',
                'particles', 'ptclSampling', 'kV', 'Cs', 'Q0',
                'hasAngles',
                'micrographs', 'micSampling',
                'atomicModel', 'doMultimodel',
                'workflow',
                'xlm',
                'saxs',
                'untiltedMic', 'tiltedMic', 'tiltkV', 'tiltCs', 'tiltQ0', 'tiltSampling', 'tiltAngle',
                'untiltedCoords', 'tiltedCoords']

PARAM_FILE_FIELDS = ['map', 'map1', 'map2',
                     'avgs',
                     'particles',
                     'micrographs',
                     'atomicModel',
                     'workflow',
                     'xlm',
                     'saxs',
                     'untiltedMic', 'tiltedMic',
                     'untiltedCoords', 'tiltedCoords']

#
# utility paths & urls
#
URL_EMDB_EBI_REPOSITORY = 'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/'
URL_EMDB_WWPDB_REPOSITORY = 'https://ftp.wwpdb.org/pub/emdb/structures/'
URL_EMDB_RCSB_REPOSITORY = 'https://ftp.rcsb.org/pub/emdb/structures/'
URL_EMDB_FTP_SERVER = "ftp.ebi.ac.uk"
URL_EMDB_FTP_DIR = "pub/databases/emdb/structures/%s/other"
URL_EMDB_REST_API = 'https://www.ebi.ac.uk/emdb/api/entry/'
URL_PDB_EBI_REPOSITORY = 'https://www.ebi.ac.uk/pdbe/entry-files/download/'
URL_PDB_RCSB_REPOSITORY = 'https://files.rcsb.org/download/'
PATH_HOME = Path.home()
PATH_TOOLS_DIR = Path.cwd()
PATH_APP = os.path.join(PATH_HOME, 'ScipionUserData')
DIR_PROJ_UPLOADS = 'uploads'
DIR_PROJ_LOGS = 'logs'
DIR_TOOLS_LOGS = 'logs'
PATH_WORK_DIR = os.path.join(PATH_APP, DIR_PROJ_UPLOADS)
FN_LOCAL_JSON_PARAMS = 'inputParams.json'


SERVER_DISPLAY_SET = 'export DISPLAY=:1.0'
PATH_PYTHON = os.path.join(PATH_HOME, 'scipion3/scipion3')
CMD_PYTHON = 'python'
PATH_SCRIPT = os.path.join(PATH_HOME, 'scipion3/scipion-em-validation/')
CMD_SCRIPT = 'validationLevels.py'


logger = logging.getLogger(__name__)


def logSetup(name, path, filename):
    """
    Configure the logging files
    """
    # create the dir and file if not exists
    save2file("", path, filename)
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fileHandler = logging.FileHandler(
        filename=os.path.join(path, filename), mode='a')
    fileHandler.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setLevel(logging.ERROR)
    formatter = logging.Formatter(
        fmt='%(asctime)s %(levelname)s %(threadName)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fileHandler.setFormatter(formatter)
    consoleHandler.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fileHandler)
    logger.addHandler(consoleHandler)
    return logger


def save_json(data,  path, filename, createIfNotExist=True):
    """
    Save the data as json file
    """
    logger.info("Save json %s %s" % (path, filename))
    if path and createIfNotExist:
        os.makedirs(path, exist_ok=True)
    full_path = os.path.join(path, filename)
    with open(full_path, 'w') as f:
        json.dump(data, f)
    return f.name


def read_json(fileName):
    try:
        with open(fileName) as json_file:
            jData = json.load(json_file)
        return jData
    except Exception as ex:
        logger.error('Could not read file %s: %s' % (fileName, ex))


def save2file(data, path, filename, createIfNotExist=True, append=True):
    """
    Save data into a text file
    if data is a list of items, will be concatenated in a single string
    file will be created if not exists, by default
    data will be appended, by default
    """
    logger.info("Save data: %s %s" % (path, filename))
    if path and createIfNotExist:
        os.makedirs(path, exist_ok=True)
    full_path = os.path.join(path, filename)
    strData = ''
    if isinstance(data, list):
        strData = ' '.join([str(item) for item in data if item])
        strData += '\n'
    else:
        strData = data
    if append:
        with open(full_path, 'a') as f:
            f.write(strData)
    else:
        with open(full_path, 'w') as f:
            f.write(strData)
    return f.name


def deleteFile(fileName, filePath):
    """
    Delete file
    """
    try:
        logger.info("Delete file %s" % fileName)
        os.remove(os.path.join(filePath, fileName))
    except FileNotFoundError as ex:
        logger.error("Could not delete file %s %s" % (fileName, ex.strerror))
        logger.warning(
            "Seems like %s wasn't there to be deleted. Never mind!", fileName)


def createDir(path, name):
    """
    Create (if not exists) a working directory
    """
    fullPath = os.path.join(path, name)
    pathObj = Path(fullPath)
    if not pathObj.exists():
        pathObj.mkdir(parents=True)
        logger.info('Created directory %s' % fullPath)
    else:
        logger.warning('Directory already exists %s' % fullPath)
    return fullPath


def removeDir(fullPath, force=False):
    """
    Remove a directory, must be empty.
    """
    try:
        Path(fullPath).rmdir()
        logger.info('Deleted directory %s' % fullPath)
    except Exception as ex:
        raise


def ungzipFile(srcFile, destFile, remove=True):
    """
    Uncompress a gz file
    by default (remove=True) the gzipped file will be
    removed after
    """
    logger.info('Unzipping %s -> %s' % (srcFile, destFile))
    if srcFile.endswith('.gz'):
        with open(destFile, 'wb') as fout, gzip.GzipFile(srcFile) as fgzip:
            shutil.copyfileobj(fgzip, fout)
    if remove and os.path.exists(srcFile):
        os.remove(srcFile)
        logger.info('Removing gz file %s' % srcFile)
    return destFile


def downloadFile(url, dirPath, fileName, raw=False):
    """
    Download a file

    If raw=True, use Response.raw and shutil.copyfileobj()
    with binary files in order to stream to disk directly 
    from the raw socket response from the server
    without using excessive memory
    """
    fullfname = os.path.join(dirPath, fileName)
    if not os.path.exists(fullfname):
        logger.info('Downloading %s to %s' % (url, fullfname))
        try:
            with requests.get(url, stream=True) as req:
                if req.status_code == 404:
                    logger.error('Could not download file %s' % url)
                    return None
                if raw:
                    with open(fullfname, 'wb') as fd:
                        shutil.copyfileobj(req.raw, fd, length=8192)
                else:
                    with open(fullfname, "w") as fd:
                        fd.write(req.text)

        except Exception as ex:
            logger.error('Could not download file %s: %s' % (url, ex))
    else:
        logger.warning('File %s already exists.' % fullfname)
    return fileName


def getEmdbMetadata(mapNum, dirPath):
    """
    Get map metadata from EMDB API
    """
    entry = 'EMD-' + mapNum
    jsonMapfile = 'EMD-' + mapNum + '.json'
    # check parameters json
    url = URL_EMDB_REST_API + entry
    logger.info('Getting required parameters from %s' % url)
    fileName = downloadFile(url, dirPath, jsonMapfile)
    if fileName:
        return read_json(os.path.join(dirPath, fileName))


def downloadEMDB_halfMaps(emdb_id, path):
    """
    Download half-maps files from EMDB
    """
    outName = os.path.join(path, emdb_id)
    if not os.path.exists(outName):
        try:
            ftp = ftplib.FTP(URL_EMDB_FTP_SERVER)
            ftp.login()
            ftp.cwd(URL_EMDB_FTP_DIR % emdb_id)
            fnames = ftp.nlst()
            half_names_inServer = [None, None]

            for fname in fnames:
                match_objs = re.match(".*half[-_\.]*(map_?)*([12])", fname)
                if match_objs:
                    half_names_inServer[int(match_objs.group(2)) - 1] = fname

            logger.info('Half maps' % half_names_inServer)

            if None in half_names_inServer:
                raise ValueError("Half map not available")

            tmp_halfs = []
            for halfNum in range(1, 3):
                try:
                    tmp_half_name = outName+"_half_%d.mrc" % halfNum
                    tmp_halfs.append(tmp_half_name)

                    with BytesIO() as flo:
                        ftp.retrbinary(
                            'RETR ' + half_names_inServer[halfNum-1], flo.write)
                        flo.seek(0)

                        with open(tmp_half_name, 'wb') as fout, gzip.GzipFile(fileobj=flo) as fgzip:
                            shutil.copyfileobj(fgzip, fout)

                except ftplib.all_errors as e:
                    deleteFile(tmp_half_name)
                    msg = "Error downloading half-maps for emdb_id: %s to: %s" % (
                        emdb_id, outName)
                    logger.error(msg)
                    raise ValueError(msg)

        except ftplib.all_errors as e:
            msg = "Error downloading half-maps for emdb_id: %s to: %s" % (
                emdb_id, outName)
            logger.error(msg)
            raise ValueError(msg)
    else:
        logger.warning('File already exists %s' % outName)
    return outName


def saveIntermediateData(fnProjectDir, protocol, isFile, fileOrParamKey, value, details=None):
    """
    Save intermediate interesting data (files and params) from each method to be stored when the validation finishes
    'fnProjectDir': report path
    'protocol': protocol name
    'isFile': True if is a file, False if is a param
    'fileOrParamKey': measure name
    'value': measure value (a number, a string, a list, etc.)
    'details': extra info
    """
    file = os.path.join(fnProjectDir, 'intermediateData.json')
    if os.path.exists(file):
        with open(file, 'r') as json_file:
            intermediateData = json.load(json_file)
            if protocol in intermediateData:
                if 'files' in intermediateData[protocol] if isFile else 'params' in intermediateData[protocol]:
                    intermediateData[protocol]['files' if isFile else 'params'].append({'name': fileOrParamKey, 'value': value, 'details': details})
                else:
                    intermediateData[protocol]['files' if isFile else 'params'] = [{'name': fileOrParamKey, 'value': value, 'details': details}]
            else:
                intermediateData[protocol] = {'files' if isFile else 'params': [{'name': fileOrParamKey, 'value': value, 'details': details}]}

        with open(file, 'w') as json_file:
            json.dump(intermediateData, json_file)
            json_file.close()
    else:
        with open(file, 'w') as json_file:
            intermediateData = {protocol: {'files' if isFile else 'params': [{'name': fileOrParamKey, 'value': value, 'details': details}]}}
            json.dump(intermediateData, json_file)
            json_file.close()

    return intermediateData