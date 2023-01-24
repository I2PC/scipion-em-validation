import os
import time
import json
from pyworkflow.object import Boolean
from pyworkflow.protocol import getProtocolFromDb

#TODO: add function to set whether the user wants to use slurm or not

def waitOutput(project, prot, outputAttributeName, sleepTime=10, timeOut=432000):
  """ Wait until the output is being generated by the protocol. """

  def _loadProt():
    # Load the last version of the protocol from its own database
    loadedProt = getProtocolFromDb(prot.getProject().path,
                                   prot.getDbPath(),
                                   prot.getObjId())
    # Close DB connections
    loadedProt.getProject().closeMapper()
    loadedProt.closeMappers()
    return loadedProt

  counter = 1
  prot2 = _loadProt()

  numberOfSleeps = timeOut / sleepTime

  while (not prot2.hasAttribute(outputAttributeName)) and prot2.isActive():
    time.sleep(sleepTime)
    prot2 = _loadProt()
    if counter > numberOfSleeps:
      print("Timeout (%s) reached waiting for %s at %s" % (timeOut, outputAttributeName, prot))
      break
    counter += 1

  # Update the protocol instance to get latest changes
  project._updateProtocol(prot)

def waitOutputFile(project, prot, outputFileName, sleepTime=10, timeOut=432000):
  """ Wait until the output file is being generated by the protocol. """
  def _loadProt():
    # Load the last version of the protocol from its own database
    loadedProt = getProtocolFromDb(prot.getProject().path,
                                   prot.getDbPath(),
                                   prot.getObjId())
    # Close DB connections
    loadedProt.getProject().closeMapper()
    loadedProt.closeMappers()
    return loadedProt

  counter = 1
  prot2 = _loadProt()

  numberOfSleeps = timeOut / sleepTime
  while (not os.path.exists(os.path.join(project.getPath(), prot2._getExtraPath(outputFileName)))) and prot2.isActive():
    time.sleep(sleepTime)
    prot2 = _loadProt()
    if counter > numberOfSleeps:
      print("Timeout (%s) reached waiting for %s at %s" % (timeOut, outputFileName, prot))
      break
    counter += 1

def waitUntilFinishes(project, prot, sleepTime=10, timeOut=432000):
  """ Wait until the protocol finishes. """

  def _loadProt():
    # Load the last version of the protocol from its own database
    loadedProt = getProtocolFromDb(prot.getProject().path,
                                   prot.getDbPath(),
                                   prot.getObjId())
    # Close DB connections
    loadedProt.getProject().closeMapper()
    loadedProt.closeMappers()
    return loadedProt

  counter = 1
  prot2 = _loadProt()

  numberOfSleeps = timeOut / sleepTime

  while not prot2.isFinished():
    time.sleep(sleepTime)
    prot2 = _loadProt()
    if counter > numberOfSleeps:
      print("Timeout (%s) reached waiting for %s to finish" % (timeOut, prot))
      break
    counter += 1

  # Update the protocol instance to get latest changes
  project._updateProtocol(prot)

#TODO: allow the user to add parameters (e.g. prot.gpuList.set(7))

def sendToSlurm(prot, GPU=False, nGPUs=1):
    prot._useQueue.set(Boolean(True))
    QUEUE_PARAMS = (u'myslurmqueue', {u'JOB_TIME': u'120', u'JOB_MEMORY': u'8192', u'GPU_COUNT': u'%s' % (nGPUs if GPU else 0)})
    prot._queueParams.set(json.dumps(QUEUE_PARAMS))

def skipSlurm(prot):
    prot._useQueue.set(Boolean(False))
    prot.gpuList.set(7)
