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

def importAvgs(project, label, protImportMap, fnAvgs, TsAvg):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportAverages', doRaise=True)
    protImport = project.newProtocol(Prot,
                               objLabel=label,
                               filesPath=fnAvgs,
                               samplingRate=TsAvg)
    project.launchProtocol(protImport, wait=True)
    if protImport.isFailed():
        raise Exception("Import averages did not work")

    XdimMap = protImportMap.outputVolume.getDim()[0]
    TsMap = protImportMap.outputVolume.getSamplingRate()
    AMap = XdimMap * TsMap

    XdimAvgs = protImport.outputAverages.getDim()[0]

    if XdimAvgs==XdimMap and TsAvg==TsMap:
        return prot

    XdimAvgsp = int(AMap/TsAvg)
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCropResizeParticles', doRaise=True)
    protResize1 = project.newProtocol(Prot,
                                      objLabel="Resize Avgs",
                                      doWindow=True,
                                      windowOperation=1 if XdimAvgsp>XdimAvgs else 0,
                                      windowSize=XdimAvgsp)
    protResize1.inputParticles.set(protImport.outputAverages)
    project.launchProtocol(protResize1, wait=True)

    protResize2 = project.newProtocol(Prot,
                                      objLabel="Resize and resample Avgs",
                                      doResize=True,
                                      resizeSamplingRate=TsMap,
                                      doWindow=True,
                                      windowOperation=1 if XdimMap>XdimAvgsp else 0,
                                      windowSize=XdimMap)
    protResize2.inputParticles.set(protResize1.outputAverages)
    project.launchProtocol(protResize2, wait=True)
    return protResize2

def reportInput(project, report, fnAvgs, protAvgs):
    avgStack = os.path.join(report.getReportDir(),"avgs.xmd")
    writeSetOfParticles(protAvgs.outputAverages, avgStack)
    toWrite = \
"""
\\section{2D Classes}
Set of 2D classes: %s \\\\
\\\\
The classes can be seen in Fig. \\ref{fig:classes2D}.\\\\
""" % (fnAvgs.replace('_','\_').replace('/','/\-'))
    report.write(toWrite)

    report.setOfImages(avgStack, xmipp3.MDL_IMAGE, "Set of 2D classes provided by the user",
                       "fig:classes2D", os.path.join(report.getReportDir(),"avg2D_"), "1.5cm", 8)
    cleanPath(avgStack)

def level2(project, report, importMap, fnAvgs, TsAvg, skipAnalysis = False):
    # Import averages
    protAvgs = importAvgs(project, "import averages", importMap, fnAvgs, TsAvg)
    reportInput(project, report, fnAvgs, protAvgs)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 2 analysis')
    return protAvgs