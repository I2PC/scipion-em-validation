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

import numpy as np
import os

import pyworkflow.plugin as pwplugin
from pyworkflow.utils.path import cleanPath
from xmipp3.convert import writeSetOfParticles
import xmipp3

from validationReport import reportHistogram

from resourceManager import waitOutput, sendToSlurm, waitUntilFinishes

import configparser

from resources.constants import *

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'config.yaml'))
useSlurm = config['QUEUE'].getboolean('USE_SLURM')

def importAvgs(project, label, protImportMap, fnAvgs, TsAvg):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportAverages', doRaise=True)
    protImport = project.newProtocol(Prot,
                               objLabel=label,
                               filesPath=fnAvgs,
                               samplingRate=TsAvg)
    if useSlurm:
        sendToSlurm(protImport)
    project.launchProtocol(protImport)
    #waitOutput(project, protImport, 'outputAverages')
    waitUntilFinishes(project, protImport)
    if protImport.isFailed():
        raise Exception("Import averages did not work")
    if protImport.isAborted():
        raise Exception("Import averages was MANUALLY ABORTED")

    XdimMap = protImportMap.outputVolume.getDim()[0]
    TsMap = protImportMap.outputVolume.getSamplingRate()
    AMap = XdimMap * TsMap

    XdimAvgs = protImport.outputAverages.getDim()[0]

    if XdimAvgs==XdimMap and TsAvg==TsMap:
        return protImport

    XdimAvgsp = int(AMap/TsAvg)
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCropResizeParticles', doRaise=True)
    protResize1 = project.newProtocol(Prot,
                                      objLabel="Resize Avgs",
                                      doWindow=True,
                                      windowOperation=1,
                                      windowSize=XdimAvgsp)
    protResize1.inputParticles.set(protImport.outputAverages)
    if useSlurm:
        sendToSlurm(protResize1)
    project.launchProtocol(protResize1)
    #waitOutput(project, protResize1, 'outputAverages')
    waitUntilFinishes(project, protResize1)

    protResize2 = project.newProtocol(Prot,
                                      objLabel="Resize and resample Avgs",
                                      doResize=True,
                                      resizeSamplingRate=TsMap,
                                      doWindow=True,
                                      windowOperation=1,
                                      windowSize=XdimMap)
    protResize2.inputParticles.set(protResize1.outputAverages)
    if useSlurm:
        sendToSlurm(protResize2)
    project.launchProtocol(protResize2)
    waitUntilFinishes(project, protResize2)

    return protImport, protResize2


def compareReprojections(project, report, protImportMap, protAvgs, symmetry):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCompareReprojections', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="2.a Compare reprojections",
                               optimizeGray=True,
                               doEvaluateResiduals=True,
                               symmetryGroup=symmetry)
    prot.inputSet.set(protAvgs.outputAverages)
    prot.inputVolume.set(protImportMap.outputVolume)
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'reprojections')
    waitUntilFinishes(project, prot)

    secLabel = "sec:fsc3d"
    msg = \
"""
\\subsection{Level 2.a Reprojection consistency}
\\label{%s}
\\textbf{Explanation}:\\\\ 
The 2D classes can be aligned against the reconstructed map, then the correlation between reprojections of the map 
and the 2D classes can be analyzed. Also, analyzing the residuals (2D class minus the corresponding reprojection) can 
reveal systematic differences between them.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("2.a Compare reprojections", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_COMPARE_PROJECTIONS)
        report.writeSummary("2.a Compare reprojections", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    md = xmipp3.MetaData(prot._getExtraPath("anglesCont.xmd"))
    md.sort(xmipp3.MDL_MAXCC)

    fnHist = os.path.join(report.getReportDir(),"reprojectionCCHist.png")
    cc = md.getColumnValues(xmipp3.MDL_MAXCC)
    reportHistogram(cc,'Cross-correlation',fnHist)

    toWrite =\
"""Fig. \\ref{fig:reprojectionCChist} shows the histogram of the cross-correlation between the 2D classes and 
the map reprojections. The average correlation is %f, and its range is [%f,%f]. Now we show
the 2D classes, the corresponding reprojection, the difference between both (residual), the covariance matrix of the
residual image, and the correlation between the 2D class and the reprojection. For a perfect match, the residual
would be just noise, and its covariance matrix should be a diagonal. Rows have been sorted by correlation so that
the worse correlating images are displayed at the beginning.

Also, 2D classes of a high-resolution map should also be of high resolution. This cannot, for the moment, be
automatically assessed. But a visual inspection should confirm that the resolution of the 2D classes match the
reported resolution of the map.

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Histogram of the correlation coefficient between the 2D classes provided by the user and the
    corresponding reprojections.}
    \\label{fig:reprojectionCChist}
\\end{figure}
"""%(np.mean(cc), np.min(cc), np.max(cc), fnHist)
    report.write(toWrite)

    report.showj(md,
                 [xmipp3.MDL_IMAGE, xmipp3.MDL_IMAGE_REF, xmipp3.MDL_IMAGE_RESIDUAL, xmipp3.MDL_IMAGE_COVARIANCE, \
                  xmipp3.MDL_MAXCC],
                 [True, True, True, True, False],
                 ['', '', '', '', '%5.2f '],
                 ["2D Class","Reprojection","Residual","Covariance",'Correlation'],
                 os.path.join(report.getReportDir(),"reproj_"), "2cm")

    # Warnings
    warnings=[]
    testWarnings = False
    if np.sum(np.array(cc)<0.7)/len(cc)>0.2 or testWarnings:
        warnings.append("{\\color{red} \\textbf{A large fraction of the 2D classes, %4.1f\\%%, correlate less "\
                        "than 0.7 with reprojections of the input map}}"%(np.sum(np.array(cc)<0.7)/len(cc)*100))
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the proportion of classes for which the correlation is below
0.7 is smaller than 20\\%.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "2.a Reprojection consistency", secLabel)

    if len(warnings)==0:
        report.writeAbstract("The 2D classes provided by the user do not seem to correlate well with the "\
                             "reprojections of the map (see Sec. \\ref{%s}). "%secLabel)


def reportInput(project, report, fnAvgs, protAvgs):
    avgStack = os.path.join(report.getReportDir(),"avgs.xmd")
    writeSetOfParticles(protAvgs.outputAverages, avgStack)

    # Get file basename to write it in the report
    basenameFnAvgs = os.path.basename(fnAvgs)

    toWrite = \
"""
\\section{2D Classes}
Set of 2D classes: %s \\\\
\\\\
The classes can be seen in Fig. \\ref{fig:classes2D}.\\\\
""" % (basenameFnAvgs.replace('_','\_').replace('/','/\-'))
    report.write(toWrite)

    report.setOfImages(avgStack, xmipp3.MDL_IMAGE, "Set of 2D classes provided by the user",
                       "fig:classes2D", os.path.join(report.getReportDir(),"avg2D_"), "1.5cm", 8)
    cleanPath(avgStack)

def level2(project, report, protImportMap, fnAvgs, TsAvg, symmetry, skipAnalysis = False):
    # Import averages
    protImportAvgs, protAvgs = importAvgs(project, "import averages", protImportMap, fnAvgs, TsAvg)
    reportInput(project, report, fnAvgs, protAvgs)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 2 analysis')
        compareReprojections(project, report, protImportMap, protAvgs, symmetry)
    return protImportAvgs, protAvgs
