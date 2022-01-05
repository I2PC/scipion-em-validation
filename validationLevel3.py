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

import math
import numpy as np
import os

import pyworkflow.plugin as pwplugin
from pyworkflow.utils.path import cleanPath
from xmipp3.convert import writeSetOfParticles
import xmipp3

from validationReport import reportHistogram

def importParticles(project, label, protImportMap, fnParticles, TsParticles, kV, Cs, Q0):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportParticles', doRaise=True)
    protImport = project.newProtocol(Prot,
                                     objLabel=label,
                                     filesPath=fnParticles,
                                     samplingRate=TsParticles,
                                     voltage=kV,
                                     sphericalAberration=Cs,
                                     amplitudeContrast=Q0)
    if fnParticles.endswith(".sqlite"):
        protImport.importFrom.set(protImport.IMPORT_FROM_SCIPION)
        protImport.sqliteFile.set(fnParticles)
    elif fnParticles.endswith(".xmd"):
        protImport.importFrom.set(protImport.IMPORT_FROM_XMIPP)
        protImport.mdFile.set(fnParticles)
    elif fnParticles.endswith(".star"):
        protImport.importFrom.set(protImport.IMPORT_FROM_RELION)
        protImport.starFile.set(fnParticles)
    project.launchProtocol(protImport, wait=True)
    if protImport.isFailed():
        raise Exception("Import averages did not work")
    return protImport

    XdimMap = protImportMap.outputVolume.getDim()[0]
    TsMap = protImportMap.outputVolume.getSamplingRate()
    AMap = XdimMap * TsMap

    XdimPtcls = protImport.outputParticles.getDim()[0]

    if XdimPtcls==XdimMap and TsAvg==TsMap:
        return protImport

    XdimPtclsp = int(AMap/TsParticles)
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCropResizeParticles', doRaise=True)
    protResize1 = project.newProtocol(Prot,
                                      objLabel="Resize Ptcls",
                                      doWindow=True,
                                      windowOperation=1 if XdimPtclsp>XdimPtcls else 0,
                                      windowSize=XdimPtclsp)
    protResize1.inputParticles.set(protImport.outputParticles)
    project.launchProtocol(protResize1, wait=True)

    protResize2 = project.newProtocol(Prot,
                                      objLabel="Resize and resample Ptcls",
                                      doResize=True,
                                      resizeSamplingRate=TsMap,
                                      doWindow=True,
                                      windowOperation=1 if XdimMap>XdimPtclsp else 0,
                                      windowSize=XdimMap)
    protResize2.inputParticles.set(protResize1.outputParticles)
    project.launchProtocol(protResize2, wait=True)
    return protResize2

def reportInput(project, report, fnParticles, protParticles):
    particlesStack = os.path.join(report.getReportDir(),"particles.xmd")
    writeSetOfParticles(protParticles.outputParticles, particlesStack)
    toWrite = \
"""
\\section{Particles}
Set of Particles classes: %s \\\\
\\\\
The first 32 can be seen in Fig. \\ref{fig:particles}.\\\\
""" % (fnParticles.replace('_','\_').replace('/','/\-'))
    report.write(toWrite)

    report.setOfImages(particlesStack, xmipp3.MDL_IMAGE, "First particles of the set of particles provided by the user",
                       "fig:particles", os.path.join(report.getReportDir(),"particles_"), "1.5cm", 8, imgMax=31)
    cleanPath(particlesStack)

def level3(project, report, protImportMap, fnParticles, TsParticles, kV, Cs, Q0, skipAnalysis = False):
    # Import particles
    protParticles = importParticles(project, "import particles", protImportMap, fnParticles, TsParticles, kV, Cs, Q0)
    reportInput(project, report, fnParticles, protParticles)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 3 analysis')
    return protParticles