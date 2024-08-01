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
import scipy

import pyworkflow.plugin as pwplugin
import xmipp3

from validationReport import reportHistogram, reportPlot, reportMultiplePlots
from resourceManager import waitOutput, sendToSlurm, skipSlurm, waitUntilFinishes, waitOutputFile

import configparser

from resources.constants import *

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'config.yaml'))
useSlurm = config['QUEUE'].getboolean('USE_SLURM')
gpuIdSkipSlurm = config['QUEUE'].getint('GPU_ID_SKIP_SLURM')

def resizeProject(project, protMap, protParticles, resolution):
    Xdim = protMap.outputVolume.getDim()[0]
    Ts = protMap.outputVolume.getSamplingRate()
    AMap = Xdim * Ts

    TsTarget = resolution/2
    Xdimp = AMap/TsTarget
    Xdimp = int(2*math.floor(Xdimp/2))

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCropResizeParticles', doRaise=True)
    protResizeParticles = project.newProtocol(Prot,
                                              objLabel="Resize Ptcls Ts=%2.1f"%TsTarget,
                                              doResize=True,
                                              resizeSamplingRate=TsTarget,
                                              doWindow=True,
                                              windowOperation=1,
                                              windowSize=Xdimp)
    protResizeParticles.inputParticles.set(protParticles.outputParticles)
    if useSlurm:
        sendToSlurm(protResizeParticles)
    project.launchProtocol(protResizeParticles)
    #waitOutput(project, protResizeParticles, 'outputParticles')
    waitUntilFinishes(project, protResizeParticles)

    return protResizeParticles


def similarityMeasures(project, report, protMap, protMask, protParticles, symmetry, resolution):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtReconstructHighRes', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="4.a Highres local",
                               symmetryGroup=symmetry,
                               alignmentMethod=1,
                               numberOfIterations=1,
                               maximumTargetResolution=max(10,resolution),
                               contShift=False,
                               contAngles=False,
                               useGpu=False)
    prot.inputVolumes.set(protMap.outputVolume)
    prot.inputParticles.set(protParticles.outputParticles)
    prot.nextMask.set(protMask.outputMask)
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputVolume')
    #waitOutput(project, prot, 'outputParticles')
    waitUntilFinishes(project, prot)

    secLabel = "sec:outlierDetection"
    msg = \
"""
\\subsection{Level 4.a Similarity criteria}
\\label{%s}
\\textbf{Explanation}:\\\\ 
We measured the similarity between the experimental images, with the angles and shifts provided by the user, and
reprojections of the input map along the same direction. We measured the correlation and IMED (IMage Euclidean Distance,
which is a generalized measure of the standard Euclidean Distance in which nearby pixels also contribute to the
calculation of the final distance between the image at a given point)
(see this \\href{%s}{link} for more details) between both sets of images. If the set of particles is properly assigned there should be a 
single peak in the
1D histograms of these two similarity measures, and a single cloud in their joint scatter plot. The presence of 
multiple peaks could reveal a mixture of different conformations, the presence of misaligned particles or 
contaminations, or the difference between isolated particles and particles with other objects around.\\\\ 
\\\\
It must be noted that there is a dependence between similarity metrics and defocus. Typically this dependence
is such that lower defoci imply lower similarity due to the smallest contrast. You have to be sure that the groups
seen in the similarity measures are not caused by defocus groups.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, SIMILARITY_MEASURES_DOI)
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.a Similarity criteria", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_SIMILARITY_CRITERIA)
        report.writeSummary("4.a Similarity criteria", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    md=xmipp3.MetaData(prot._getExtraPath(os.path.join('Iter001','angles.xmd')))
    imed = md.getColumnValues(xmipp3.MDL_IMED)
    cc = md.getColumnValues(xmipp3.MDL_CORRELATION_IDX)
    defocus = 0.5*(np.array(md.getColumnValues(xmipp3.MDL_CTF_DEFOCUSU))+\
                   np.array(md.getColumnValues(xmipp3.MDL_CTF_DEFOCUSV)))

    fnImedHist = os.path.join(report.getReportDir(),"reprojectionImedHist.png")
    reportHistogram(imed,'IMED',fnImedHist)
    fnCCHist = os.path.join(report.getReportDir(),"reprojectionCCHist.png")
    fnImedCC = os.path.join(report.getReportDir(),"reprojectionImedCC.png")
    fnDefocusCC = os.path.join(report.getReportDir(),"reprojectionDefocusCC.png")
    reportHistogram(cc,'Cross-correlation',fnCCHist)
    reportPlot(imed, cc, "IMED", "Cross-correlation", fnImedCC, plotType='scatter')
    reportPlot(defocus, cc, "Defocus (A)", "Cross-correlation", fnDefocusCC, plotType='scatter')

    msg=\
"""Fig. \\ref{fig:reprojectionSimilarity} shows the histograms of the cross-correlation and IMED, a joint scatter
plot and the dependence of the cross-correlation with the defocus.\\\\

\\begin{figure}[H]
    \centering
    \includegraphics[width=6.5cm]{%s}
    \includegraphics[width=6.5cm]{%s} \\\\
    \includegraphics[width=6.5cm]{%s}
    \includegraphics[width=6.5cm]{%s}
    \\caption{Top: Histogram of the cross-correlation (CC) and IMED between the experimental images and their corresponding
             reprojections. Bottom: Scatter plots of CC vs IMED, and CC vs defocus.}
    \\label{fig:reprojectionSimilarity}
\\end{figure}

"""%(fnCCHist, fnImedHist, fnImedCC, fnDefocusCC)

    # def bic(X,K):
    #     from sklearn.mixture import GaussianMixture
    #     gmm = GaussianMixture(n_components=K).fit(X)
    #     return gmm.bic(X)
    # X = np.zeros((len(cc),2))
    # X[:,0]=imed
    # X[:,1]=cc
    # msg+="We calculated the BIC of a clustering of the IMED-CC scatter plot with a Gaussian Mixture" \
    #      "with K=1, 2, and 3 components. The clustering was repeated 20 times and we report here the mean "\
    #      "Bayesian Information Criterion (BIC) and its standard error. \\\\"
    # bestS = -1
    # bestK = -1
    # for K in range(1,4):
    #     S = []
    #     for n in range(20):
    #         S.append(bic(X,K))
    #     Savg = np.mean(S)
    #     Sstd = np.std(S)
    #     msg+="For K=%d, the average BIC was %f and its standard deviation %f.\\\\"%(K, Savg, Sstd)
    #     if Savg>bestS:
    #         bestS=Savg
    #         bestK=K
    # msg+="\\\\From these values, we conclude the best number of classes is %d.\n\n"%bestK

    report.write(msg)

    # Warnings
    report.writeWarningsAndSummary(None, "4.a Similarity criteria", secLabel)

def alignabilitySmoothness(project, report, protMap, protMask, protParticles, symmetry, resolution):

    secLabel = "sec:alignabilitySmoothness"
    msg = \
"""
\\subsection{Level 4.b Alignability smoothness}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This algorithm (see this \\href{%s}{link} for more details) analyzes the smoothness of the correlation function over the projection sphere and
the stability of its maximum. Ideally, the angular assignment given by the user should coincide with the maximum
of the smoothed cross-correlation landscape.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, ALIGNABILITY_SMOOTHNESS_DOI)
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtAngularGraphConsistency', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="4.b Alignability smoothness",
                               inputVolume=protMap.outputVolume,
                               inputParticles=protParticles.outputParticles,
                               symmetryGroup=symmetry,
                               maximumTargetResolution=resolution,
                               numberOfMpi=8)
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputParticles')
    #waitOutput(project, prot, 'outputParticlesAux')
    waitUntilFinishes(project, prot)

    if prot.isFailed():
        report.writeSummary("4.b Alignability smoothness", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_ALIGN_SMOOTHNESS)
        report.writeSummary("4.b Alignability smoothness", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    md = xmipp3.MetaData(prot._getExtraPath("Iter1/anglesDisc.xmd"))
    dist2Max = np.array(md.getColumnValues(xmipp3.MDL_GRAPH_DISTANCE2MAX_PREVIOUS))

    fnDistHist = os.path.join(report.getReportDir(), "graphDistHist.png")
    reportHistogram(dist2Max, 'Angular distance to maximum', fnDistHist)

    avgDist = np.mean(dist2Max)

    fDist = np.sum(dist2Max > 10) / dist2Max.size * 100

    msg = \
"""Fig. \\ref{fig:graphValidation} shows the histogram of the angular distance 
between the angular assignment given by the user and the maximum of the smoothed landscape of cross-correlations. 
plot. The average angular distance %4.3f. The percentage of images whose distance is larger than 10 is %4.1f\\%%.\\\\

\\begin{figure}[H]
    \centering
    \includegraphics[width=6.5cm]{%s}
    \\caption{Histogram of the angular distance between the angular assignment given
              by the user and the maximum of the smoothed landscape of cross-correlation.}
    \\label{fig:graphValidation}
\\end{figure}
""" % (avgDist, fDist, fnDistHist)
    report.write(msg)

    warnings = []
    testWarnings = False
    if fDist > 30 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images whose angular assignment is significantly "\
                          "away from the smoothed maximum is too high, %4.1f\\%%}}" % fDist)
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if less than 30\\% of the images have their angular assignment 
is further than 10 degrees from the smoothed cross-correlation maximum. 
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "4.b Alignability smoothness", secLabel)

    if len(warnings) > 0:
        report.writeAbstract("It seems that the input particles cannot be easily aligned (see Sec. \\ref{%s}). " % \
                             secLabel)


def multirefAlignability(project, report, protMap, protMask, protParticles, symmetry):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtMultiRefAlignability', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="4.c Alignability",
                               useGpu=False,
                               symmetryGroup=symmetry,
                               numberOfMpi=8)
    prot.inputVolumes.set(protMap.outputVolume)
    prot.inputParticles.set(protParticles.outputParticles)
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputParticles')
    #waitOutput(project, prot, 'outputVolumes')
    waitUntilFinishes(project, prot)

    secLabel = "sec:multirefAlignability"
    msg = \
"""
\\subsection{Level 4.c Alignability precision and accuracy}
\\label{%s}
\\textbf{Explanation}:\\\\ 
The precision (see this \\href{%s}{link} for more details) analyzes the orientation distribution of the
best matching reprojections from the reference volume. If the high values are clustered around the same orientation,
then the precision is close to 1. Otherwise, it is closer to -1. Below 0.5 the best directions tend to be scattered.
The alignability accuracy (see this \\href{%s}{link} for more details) compares the final angular assignment with the result of a new angular 
assignment. The similarity between both is again encoded between -1 and 1.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, MULTIREF_ALIGNABILITY_1_DOI, MULTIREF_ALIGNABILITY_2_DOI)
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.c Alignability", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot
    
    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_ALIGNABILITY)
        report.writeSummary("4.c Alignability", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    md=xmipp3.MetaData(prot._getExtraPath("vol001_pruned_particles_alignability.xmd"))
    prec = np.array(md.getColumnValues(xmipp3.MDL_SCORE_BY_ALIGNABILITY_PRECISION))
    acc = np.array(md.getColumnValues(xmipp3.MDL_SCORE_BY_ALIGNABILITY_ACCURACY))
    idxPrec = np.logical_and(prec>=-1, prec<=1)
    idxAcc = np.logical_and(acc>=-1, acc<=1)
    idx = np.logical_and(idxPrec, idxAcc)

    prec = prec[idx]
    acc = acc[idx]

    fnPrecHist = os.path.join(report.getReportDir(),"alignabilityPrecisionHist.png")
    reportHistogram(prec,'Precision',fnPrecHist)

    fnAccHist = os.path.join(report.getReportDir(),"alignabilityAccuracyHist.png")
    reportHistogram(acc,'Accuracy',fnAccHist)

    fnAccPrec = os.path.join(report.getReportDir(),"alignabilityAccuracyPrecision.png")

    reportPlot(prec, acc, "Precision", "Accuracy", fnAccPrec, plotType='scatter')

    avgAcc = np.mean(acc)
    avgPrec = np.mean(prec)

    fAcc = np.sum(acc<0.5)/acc.size*100
    fPrec = np.sum(prec<0.5)/prec.size*100

    msg=\
"""Fig. \\ref{fig:multirefAlignability} shows the histograms of the accuracy and precision, and a joint scatter
plot. The average accuracy was %4.3f and the average precision %4.3f. The percentage of images whose accuracy
is below 0.5 is %4.1f\\%%, and the percentage of images whose precision is below 0.5 is %4.1f\\%%.\\\\

\\begin{figure}[H]
    \centering
    \includegraphics[width=6.5cm]{%s}
    \includegraphics[width=6.5cm]{%s} \\\\
    \includegraphics[width=9cm]{%s}
    \\caption{Top: Histogram of the accuracy and precision. Bottom: Scatter plot of both measures.}
    \\label{fig:multirefAlignability}
\\end{figure}
"""%(avgAcc, avgPrec, fAcc, fPrec, fnAccHist, fnPrecHist, fnAccPrec)
    report.write(msg)

    warnings=[]
    testWarnings = False
    if fAcc>30 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with low alignability accuracy is too high, "\
                        "%4.1f\\%%}}"%fAcc)
    if fPrec>30 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with low alignability precision is too high, "\
                        "%4.1f\\%%}}"%fPrec)
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if less than 30\\% of the images have an 1) accuracy and 2)
precision smaller than 0.5. 
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "4.c Alignability precision and accuracy", secLabel)

    if len(warnings)>0:
        report.writeAbstract("It seems that the input particles cannot be easily aligned (see Sec. \\ref{%s}). "%\
                             secLabel)

def compareAlignment(project, report, refmap, protRefParticles, protReconstruction, symmetry, label,
                     agent1, agent2, fnRoot, generateSlices=True):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtAlignVolumeParticles', doRaise=True)
    protAlign = project.newProtocol(Prot,
                                    objLabel="4.d"+label+" align",
                                    symmetryGroup=symmetry)
    protAlign.inputReference.set(refmap)
    protAlign.inputVolume.set(protReconstruction.outputVolume)
    protAlign.inputParticles.set(protReconstruction.outputParticles)
    if useSlurm:
        sendToSlurm(protAlign)
    project.launchProtocol(protAlign)
    #waitOutput(project, protAlign, 'outputVolume')
    #waitOutput(project, protAlign, 'outputParticles')
    waitUntilFinishes(project, protAlign)

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCompareAngles', doRaise=True)
    protCompare = project.newProtocol(Prot,
                                    objLabel="4.d"+label+" compare",
                                    symmetryGroup=symmetry)
    protCompare.inputParticles1.set(protRefParticles.outputParticles)
    protCompare.inputParticles2.set(protReconstruction.outputParticles)
    if useSlurm:
        sendToSlurm(protCompare)
    project.launchProtocol(protCompare)
    #waitOutput(project, protCompare, 'outputParticles')
    waitUntilFinishes(project, protCompare)

    allShiftDiffs = []
    allAngleDiffs = []
    for particle in protCompare.outputParticles:
        allShiftDiffs.append(particle.getAttributeValue("_xmipp_shiftDiff"))
        allAngleDiffs.append(particle.getAttributeValue("_xmipp_angleDiff"))
    allShiftDiffs.sort()
    allAngleDiffs.sort()

    medShift = np.median(allShiftDiffs)
    medAngle = np.median(allAngleDiffs)
    madShift = scipy.stats.median_absolute_deviation(allShiftDiffs)
    madAngle = scipy.stats.median_absolute_deviation(allAngleDiffs)

    thresholdShift = 5 # medShift+3*madShift
    thresholdAngle = 5 # medAngle+3*madAngle

    outliersShift = np.sum(np.array(allShiftDiffs)>thresholdShift)/len(allShiftDiffs)*100
    outliersAngles = np.sum(np.array(allAngleDiffs)>thresholdAngle)/len(allAngleDiffs)*100

    fnShift = os.path.join(report.getReportDir(),fnRoot+"_shiftDiff.png")
    fnAngles = os.path.join(report.getReportDir(),fnRoot+"_anglesDiff.png")

    nparticles = [x+1 for x in range(len(allShiftDiffs))]
    reportPlot(nparticles, allShiftDiffs, "Particle No.", "Shift difference (A)", fnShift)
    reportPlot(nparticles, allAngleDiffs, "Particle No.", "Angular difference difference (degrees)", fnAngles)

    simplifiedLabel = label.split()[-1]

    if generateSlices:
        msg=\
    """Fig. \\ref{fig:check%s} shows some representative slices of the reconstruction performed by %s for checking its
    correctness. """%(simplifiedLabel, simplifiedLabel)
        fnMap = protReconstruction.outputVolume.getFileName()
        if fnMap.endswith(".mrc"):
            fnMap+=":mrc"
        report.orthogonalSlices("check%s"%simplifiedLabel, msg,
                                "Slices of maximum variation in the three dimensions of the "\
                                "map reconstructed by %s"%simplifiedLabel, fnMap,
                                "fig:check%s"%simplifiedLabel, maxVar=True)

    msg=\
"""Fig. \\ref{fig:comparison%s} shows the shift and angular difference between the alignment given by %s
and the one calculated by %s. The median shift difference was %4.1f\\AA, and the median angular difference
%5.1f degrees. Their corresponding median absolute deviations were %4.1f and %4.1f, respectively. Particles
with a shift difference larger than 5\\AA~or an angular difference larger than 5 degrees would be considered as 
incorrectly assigned in one of the two assignments (the user's or the new one). %4.1f\\%% of particles were considered 
to have an uncertain shift, and %4.1f\\%% of particles were considered to have an uncertain alignment. 
\\\\

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \includegraphics[width=9cm]{%s}
    \\caption{Top: Shift difference between the alignment given by %s and the one calculated by %s. Bottom:
    Angular difference. The X-axis represents all particles sorted by their difference.}
    \\label{fig:comparison%s}
\\end{figure}
"""%(simplifiedLabel, agent1, agent2, medShift, medAngle, madShift, madAngle, outliersShift, outliersAngles,
     fnShift, fnAngles, agent1, agent2, simplifiedLabel)
    report.write(msg)

    return outliersShift, outliersAngles

def relionAlignment(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry, resolution):
    Prot = pwplugin.Domain.importFromPlugin('relion.protocols',
                                            'ProtRelionRefine3D', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="4.d Relion Refine",
                               initialLowPassFilterA=10,
                               symmetryGroup=symmetry
                               )
    prot.referenceVolume.set(protResizeMap.outputVol)
    prot.inputParticles.set(protResizeParticles.outputParticles)
    prot.referenceMask.set(protResizeMask.outputVol)
    if useSlurm:
        sendToSlurm(prot, GPU=True)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputVolume')
    #waitOutput(project, prot, 'outputParticles')
    #waitOutput(project, prot, 'outputFSC')
    waitUntilFinishes(project, prot)

    secLabel = "sec:relionAlignment"
    msg = \
"""
\\subsection{Level 4.d1 Relion alignment}
\\label{%s}
\\textbf{Explanation}:\\\\ 
We have performed an independent angular assignment using Relion autorefine (see this \\href{%s}{link} for more details). Images were
downsampled to a pixel size of 3\\AA. Then, we measured the difference between the angular assignment of the 
particles given by the user and the one done by Relion.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, RELION_ALIGNMENT_AND_CLASSIFICATION_DOI)
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.d1 Relion alignment", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_RELION_ALIGN)
        report.writeSummary("4.d1 Relion alignment", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    shiftOutliers, angleOutliers = compareAlignment(project, report, protResizeMap.outputVol,
                                                    protResizeParticles, prot, symmetry,
                                                    "1. Relion", "the user", "Relion", "alignmentRelion")
    warnings = []
    testWarnings = False
    if shiftOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain shift is larger than 20\\%}}")
    if angleOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain angles is larger than 20\\%}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if less than 20\\% of the images have 1) a shift difference
larger than 5\\AA, and 2) an angular difference larger than 5 degrees.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "4.d1 Relion alignment", secLabel)

    if len(warnings)>0:
        report.writeAbstract("It seems that the angular assignment given by the user does not match with the one "\
                             "produced by Relion (see Sec. \\ref{%s}). "%secLabel)
    return prot

def cryosparcAlignment(project, report, protMap, protMask, protParticles, symmetry):
    Prot = pwplugin.Domain.importFromPlugin('cryosparc2.protocols',
                                            'ProtCryoSparcNonUniformRefine3D', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="4.d Cryosparc Refine")
    prot.referenceVolume.set(protMap.outputVol)
    prot.inputParticles.set(protParticles.outputParticles)
    prot.refMask.set(protMask.outputVol)
    if symmetry[0]=="c":
        prot.symmetryGroup.set(0)
        prot.symmetryOrder.set(int(symmetry[1:]))
    elif symmetry[0]=="d":
        prot.symmetryGroup.set(1)
        prot.symmetryOrder.set(int(symmetry[1:]))
    elif symmetry[0]=="o":
        prot.symmetryGroup.set(3)
    elif symmetry=="i1":
        prot.symmetryGroup.set(4)
    elif symmetry=="i2":
        prot.symmetryGroup.set(5)
    if useSlurm:
        skipSlurm(prot, gpuIdSkipSlurm)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputVolume')
    #waitOutput(project, prot, 'outputParticles')
    #waitOutput(project, prot, 'outputFSC')
    waitUntilFinishes(project, prot)

    secLabel = "sec:cryosparcAlignment"
    msg = \
"""
\\subsection{Level 4.d2 CryoSparc alignment}
\\label{%s}
\\textbf{Explanation}:\\\\ 
We have performed an independent angular assignment using CryoSparc non-homogeneous refinement (see this \\href{%s}{link} for more details).
Images were downsampled to a pixel size of 3\\AA.  Then, we measured the difference between the angular assignment
of the particles given by the user and the one done by CryoSparc.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, CRYOSPARC_ALIGNMENT_DOI)
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.d2 CryoSparc alignment", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_CRYOSPARC_ALIGN)
        report.writeSummary("4.d2 CryoSparc alignment", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    shiftOutliers, angleOutliers = compareAlignment(project, report, protMap.outputVol, protParticles, prot, symmetry,
                                                    "2. Cryosparc", "the user", "CryoSparc", "alignmentCryosparc")

    warnings = []
    testWarnings = False
    if shiftOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain shift is larger than 20\\%}}")
    if angleOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain angles is larger than 20\\%}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if less than 20\\% of the images have 1) a shift difference
larger than 5\\AA, and 2) an angular difference larger than 5 degrees.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "4.d2 CryoSparc alignment", secLabel)
    if len(warnings)>0:
        report.writeAbstract("It seems that the angular assignment given by the user does not match with the one "\
                             "produced by CryoSparc (see Sec. \\ref{%s}). "%secLabel)

    return prot

def compareRelionAndCryosparc(project, report, protRelion, protCryoSparc, symmetry):
    secLabel = "sec:relionCryosparc"
    msg = \
"""
\\subsection{Level 4.d3 Relion/CryoSparc alignments}
\\label{%s}
\\textbf{Explanation}:\\\\ 
In Secs. \\ref{sec:relionAlignment} and \\ref{sec:cryosparcAlignment} we compared the angular assignment given by
the user to the angular assignment of Relion and CryoSparc, respectively. We now compare these two alignments as a way
to measure the ``intrinsic'' uncertainty in the angular assignment. This comparison gives an estimate of the 
alignability of the input images.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if protRelion.isFailed() or protCryoSparc.isFailed():
        report.writeSummary("4.d3 Relion/CryoSparc alignments", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return 

    if protRelion.isAborted() or protCryoSparc.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_RELION_CRYOSPARC_ALIGN)
        report.writeSummary("4.d3 Relion/CryoSparc alignments", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return

    shiftOutliers, angleOutliers = compareAlignment(project, report, protRelion.outputVolume, protRelion, protCryoSparc,
                                                    symmetry, "3. Relion/Cryosparc", "Relion", "CryoSparc",
                                                    "alignmentRelionCryosparc", False)

    warnings = []
    testWarnings = False
    if shiftOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain shift is larger than 20\\%}}")
    if angleOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain angles is larger than 20\\%}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if less than 20\\% of the images have 1) a shift difference
larger than 5\\AA, and 2) an angular difference larger than 5 degrees.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "4.d3 Relion/CryoSparc alignments", secLabel)
    if len(warnings)>0:
        report.writeAbstract("It seems that the angular assignment produced by Relion does not match with the one "\
                             "produced by Cryosparc (see Sec. \\ref{%s}). This is probably a sign of the difficulty "\
                             "to align these particles. " %secLabel)

def relionClassification(project, report, protMap, protMask, protParticles, symmetry):
    Prot = pwplugin.Domain.importFromPlugin('relion.protocols',
                                            'ProtRelionClassify3D', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="4.e Relion classify",
                               copyAlignment=True,
                               initialLowPassFilterA=10,
                               symmetryGroup=symmetry,
                               numberOfClasses=2,
                               doImageAlignment=False,
                               doGpu=False,
                               numberOfMpi=8)
    prot.referenceVolume.set(protMap.outputVol)
    prot.inputParticles.set(protParticles.outputParticles)
    prot.referenceMask.set(protMask.outputVol)
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputClasses')
    #waitOutput(project, prot, 'outputVolumes')
    waitUntilFinishes(project, prot)

    secLabel = "sec:relionClassification"
    msg = \
"""
\\subsection{Level 4.e Classification without alignment}
\\label{%s}
\\textbf{Explanation}:\\\\ 
We have performed a 3D classification of the input particles in two classes without aligning them using Relion
(see this \\href{%s}{link} for more details) to confirm they belong to a single state. Images were downsampled to a pixel size of 3\\AA. 
A valid result would be: 1) a class attracting most particles and an almost empty class, 2) two classes with an 
arbitrary number of images in each one, but without any significant structural difference between the two.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, RELION_ALIGNMENT_AND_CLASSIFICATION_DOI)
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.e Classification without alignment", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_CLASSIFICATION_WITHOUT_ALIGN)
        report.writeSummary("4.e Classification without alignment", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    totalImgs = protParticles.outputParticles.getSize()
    Nimgs = []
    representative = []
    for cls in prot.outputClasses.iterItems():
        Nimgs.append(len(cls.getIdSet()))
        representative.append(cls.getRepresentative().getFileName())
    warnings = None
    if len(Nimgs)==0:
        report.write("{\\color{red} The classification did not produce any output}")
        report.writeSummary("4.e Classification without alignment", secLabel, "{\\color{red} Could not be measured}\\\\")
        return prot
    elif len(Nimgs)==1:
        report.write("The classification converged to a single class with %d out of %d images in it.\\\\"%\
                     (Nimgs[0], totalImgs))
        warnings = []
        msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the classification converged to a single class.
\\\\

"""
        report.write(msg)
    elif len(Nimgs)==2:
        #waitOutputFile(project, prot, representative[0].split("/")[-1])
        #waitOutputFile(project, prot, representative[1].split("/")[-1])
        V1 = xmipp3.Image(representative[0]+":mrc")
        V2 = xmipp3.Image(representative[1]+":mrc")
        Vdiff = V1-V2
        msg=\
"""The classification converged to two classses with %d and %d images, out of a total of %d input images. To allow for 
visual inspection we show Class 1 in Fig. \\ref{fig:class1}, Class2 in Fig. \\ref{fig:class2}, and their difference in
Fig. \\ref{fig:classDiff}. There should not be any difference between both classes, the slices of the difference volume
should be just noise without any visible structure related to the macromolecule.\\\\"""%(Nimgs[0], Nimgs[1], totalImgs)
        report.write(msg)

        report.orthogonalSlices("class1", "",
                                "Slices of maximum variation in Class 1.", V1, "fig:class1", maxVar=True)
        report.orthogonalSlices("class2", "",
                                "Slices of maximum variation in Class 2.", V2, "fig:class2", maxVar=True)
        report.orthogonalSlices("classDiff", "",
                                "Slices of maximum variation in Class 1-Class 2", Vdiff, "fig:classDiff", maxVar=True)

    report.writeWarningsAndSummary(warnings, "4.e Classification without alignment", secLabel)
    return prot

def validateOverfitting(project, report, protMap, protMask, protParticles, symmetry, resolution):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtValidateOverfitting', doRaise=True)

    Nptcls = protParticles.outputParticles.getSize()

    samplePtcls=""
    for x in [0.03, 0.1, 0.3, 0.75, 1.0]:
        samplePtcls+="%d "%int(0.5*x*Nptcls)
    prot = project.newProtocol(Prot,
                               objLabel="4.f Overfitting detection",
                               numberOfIterations=5,
                               numberOfParticles=samplePtcls,
                               doNoise=True,
                               numberOfMpi=8,
                               symmetry=symmetry)

    prot.input3DReference.set(protMap.outputVol)
    prot.inputParticles.set(protParticles.outputParticles)
    if useSlurm:
        sendToSlurm(prot, GPU=True)
    project.launchProtocol(prot)
    waitUntilFinishes(project, prot)

    secLabel = "sec:relionClassification"
    msg = \
"""
\\subsection{Level 4.f Overfitting detection}
\\label{%s}
\\textbf{Explanation}:\\\\ 
The detection of overfitting can be performed through a series of 5 reconstructions with an increasing number
of experimental particles and the same number of pure noise particles (see this \\href{%s}{link} for more details). The resolution of
the reconstructions with experimental particles should always be better than those from noise. \\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, VALIDATE_OVERFITTING_DOI)
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.f Overfitting detection", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot
    
    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_OVERFITTING_DETECT)
        report.writeSummary("4.f Overfitting detection", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    def readResults(fn):
        xVal = []
        yVal = []
        yErr = []
        fileValues = open(fn,'r')
        for line in fileValues:
            values = line.split()
            xVal.append(float(values[0]))
            yVal.append(float(values[1]))
            yErr.append(float(values[2]))
        return xVal, yVal, yErr

    x, y, ye = readResults(prot._defineResultsTxt())
    xn, yn, yen = readResults(prot._defineResultsNoiseTxt())

    fnPlot = os.path.join(report.getReportDir(),"overfitting.png")
    reportMultiplePlots(x,[np.reciprocal(y), np.reciprocal(yn)], "No. Particles", "Resolution Freq. (1/A)", fnPlot,
                        ['Experimental particles', 'Noise particles'])

    overfitting = (np.sum(np.array(yn)<np.array(y))>0)
    msg=\
"""We tested with subsets of %s particles. Fig. \\ref{fig:overfitting} shows the inverse of the resolution as a function
of the number of particles.
"""%(', '.join(samplePtcls.split()[:-1]) + ' and ' + samplePtcls.split()[-1])

    if overfitting:
        msg+="We have detected that the resolution of pure noise particles is sometimes better than the one of "\
             "true particles."
    msg+=\
"""\\\\

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Inverse of the resolution as a function of the number of particles.}
    \\label{fig:overfitting}
\\end{figure}
"""%fnPlot
    report.write(msg)

    warnings=[]
    testWarnings = False
    if overfitting or testWarnings:
        warnings.append("{\\color{red} \\textbf{the resolution of pure noise particles is sometimes better than the "\
                        "one of true particles.}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the resolution of noise particles is never better than the
resolution of true particles.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "4.f Overfitting detection", secLabel)
    if len(warnings)>0:
        report.writeAbstract("It seems that there might be some overfitting (see Sec. \\ref{%s}). "%secLabel)

def angularDistributionEfficiency(project, report, protResizeParticles, symmetry, resolution, bfactor):
    Xdim = protResizeParticles.outputParticles.getDim()[0]
    Ts = protResizeParticles.outputParticles.getSamplingRate()
    APtcls = Xdim * Ts

    Prot = pwplugin.Domain.importFromPlugin('cryoef.protocols',
                                            'ProtCryoEF', doRaise=True)

    prot = project.newProtocol(Prot,
                               objLabel="4.g CryoEF",
                               symmetryGroup=symmetry,
                               diam=0.9*APtcls,
                               FSCres=resolution,
                               angAcc=10,
                               Bfact=-bfactor)

    prot.inputParticles.set(protResizeParticles.outputParticles)
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputVolume1')
    #waitOutput(project, prot, 'outputVolume2')
    waitUntilFinishes(project, prot)

    secLabel = "sec:AngularDistributionEfficiency"
    msg = \
"""
\\subsection{Level 4.g Angular distribution efficiency}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method (see this \\href{%s}{link} for more details) evaluates the ability of the angular distribution to fill the Fourier space. 
It determines a resolution
per direction based on the number of particles in each direction and reports the distribution efficiency, a 
number between 0 (inefficient) and 1 (total efficiency).\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, ANGULAR_DISTR_EFFICIENCY_DOI)
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.g Angular distribution efficiency", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot
    
    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_ANGULAR_DISTR_EFFICIENCY)
        report.writeSummary("4.g Angular distribution efficiency", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    fh=open(prot._getExtraPath("input_angles_PSFres.dat"))
    R=[]
    for line in fh.readlines():
        R.append(float(line))
    fh.close()
    fnRHist = os.path.join(report.getReportDir(),"angEffHist.png")
    reportHistogram(R,'Resolution (A)', fnRHist)

    avgDirResolution=np.mean(R)
    msg=\
"""Fig. \\ref{fig:efficiency} shows the histogram of the measured resolutions per direction. The average resolution
was %4.1f \\AA, and its range [%4.1f,%4.1f].

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Histogram of the directional resolution according to the angular distribution efficiency.}
    \\label{fig:efficiency}
\\end{figure}
"""%(avgDirResolution, np.min(R), np.max(R), fnRHist)
    report.write(msg)

    # Warnings
    warnings=[]
    testWarnings = False
    if resolution<0.8*avgDirResolution or testWarnings:
        warnings.append("{\\color{red} \\textbf{The resolution reported by the user, %5.2f\\AA, is at least 80\\%% "\
                        "smaller than the average directional resolution, %5.2f \\AA.}}" %\
                        (resolution, avgDirResolution))
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the resolution reported by the user is larger than 0.8 times
the average directional resolution.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "4.g Angular distribution efficiency", secLabel)

def samplingCompensationFactor(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry):

    secLabel = "sec:SCF"
    msg = \
"""
\\subsection{Level 4.h Sampling compensation factor}
\\label{%s}
\\textbf{Explanation}:\\\\ 
The SCF (see these links [\\href{%s}{link 1}, \\href{%s}{link 2}] for more details) is a measure of how effectively projections have been arranged on the
projection sphere in order to maximize
the resulting SNR. It is evaluated over the particle assignments and estimates how much the global SSNR/FSC has been
decremented due to deviations from uniformity in the projection distribution. The SCF calculation is based solely on
the concept of Fourier sampling, which decouples from other parameters in simple estimates. The SCF does not take
into account the CTF or any of the microscope effects. A low value of the SCF indicates an ineffective set of Euler
angles. A poor distribution of sampling can additionally result in other issues, such as mis-assignment of Eulerian
angles to the particles in the data, but that is not measured directly by the SCF. 

The SCF is calculated by first evaluating the effective number of measurements that take place on a spiral grid on a 
Fourier space sphere at some intermediate value of Fourier radius. Over these measurements, it is the ratio of the 
harmonic mean of the sampling divided by the mean sampling; the resultant is a quantity equal to or less than 1, 
where a value of 1 corresponds to a set of uniform views. For pure side views, the SCF value is theoretically 
evaluated at 0.81. Pure side views give rise to a fully sampled and largely isotropic reconstruction, and the ~20%% 
deviation from unity means that the spectral SNR is attenuated proportionally by that amount. Reconstructions 
characterized by an SCF value between 0.81 and 1 are thus healthy reconstructions. Between ~0.5 and 0.81, one
typically begins to observe small problems in the map and slightly elongated features due to non-uniform sampling and 
resolution anisotropy. When the SCF drops below 0.5, there tend to be more serious issues in the reconstruction and 
artifacts from resolution anisotropy. Values of SCF lower than 0.5 should encourage the experimentalist to reconsider 
the preparation of the sample or the design of the experiment, for example by tilting the stage. In all cases, the 
SCF provides a direct measure of the deviation from unity for the global spectral SNR. An SCF value of 0.5 and 0.25 
means that two and four times as much data would need to be collected, respectively, to reach an equivalent 
resolution as a fully sampled isotopic map characterized by an SCF of 1.0. However, as noted above, low SCF values 
typically go hand in hand with issues that are more severe and not directly measured by sampling, such as Euler 
angle mis-assignment.

\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, SAMPLING_COMPENS_FACTOR_1_DOI, SAMPLING_COMPENS_FACTOR_2_DOI)
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('scf.protocols',
                                            'ScfProtAnalysis', doRaise=True)

    symStr='C1'
    if symmetry.startswith('c') or symmetry.startswith('d'):
        symStr=symmetry.upper()
    elif symmetry.startswith("i"):
        symStr = 'Icos'
    elif symmetry=="o":
        symStr='Oct'

    prot = project.newProtocol(Prot,
                               objLabel="4.h SCF",
                               inParticles=protResizeParticles.outputParticles,
                               numberToUse=-1,
                               sym=symStr)

    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    waitUntilFinishes(project, prot)

    msg='The results of the SCF analysis was:\\\\ \n'
    fh=open(prot._getPath('extra.txt'))
    results={}
    for line in fh.readlines():
        tokens = line.replace('INFO:root:','').split('=')
        key = tokens[0].strip()
        try:
            value = float(tokens[1].strip())
            msg+="%s=%8.4f\\\\ \n"%(key, value)
            results[key] = value
        except:
            pass
    fh.close()
    msg+="\n\n"
    report.write(msg)

    msg=\
"""Fig. \\ref{fig:scf} shows the SCF plot for this angular distribution.

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{SCF plot.}
    \\label{fig:scf}
\\end{figure}
"""%os.path.join(project.getPath(),prot._getExtraPath("particleAnglesTilt0.jpg"))
    report.write(msg)

    # Warnings
    warnings=[]
    testWarnings = False
    if results['SCFStar']<0.5 or testWarnings:
        warnings.append("{\\color{red} \\textbf{SCF*} is smaller than 0.5}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the SCF*$>$0.5.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "4.h SCF", secLabel)


def ctfStability(project, report, protRefinement, protResizeParticles, protResizeMask):
    Prot = pwplugin.Domain.importFromPlugin('relion.protocols',
                                            'ProtRelionPostprocess', doRaise=True)

    protPostprocess = project.newProtocol(Prot,
                                          objLabel="4.i PostProcess")
    protPostprocess.protRefine.set(protRefinement)
    protPostprocess.solventMask.set(protResizeMask.outputVol)
    if useSlurm:
        sendToSlurm(protPostprocess)
    project.launchProtocol(protPostprocess)
    #waitOutput(project, protPostprocess, 'outputVolume')
    waitUntilFinishes(project, protPostprocess)

    secLabel = "sec:ctfStability"
    msg = \
"""
\\subsection{Level 4.i CTF stability}
\\label{%s}
\\textbf{Explanation}:\\\\ 
We estimated the per-particle defocus, B-factor, astigmatism, and phase shift using Relion's ctf refine (see this \\href{%s}{link} for more details). Ideally,
the differences in defoci cannot be larger than the ice thickness. We also estimated the local magnification offsets 
(which should be around 0) and the B-factor.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, CTF_STABILITY_DOI)
    report.write(msg)
    if protPostprocess.isFailed():
        report.writeSummary("4.i CTF stability", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return protPostprocess
    
    if protPostprocess.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_CTF_STABILITY)
        report.writeSummary("4.i CTF stability", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return protPostprocess

    Prot = pwplugin.Domain.importFromPlugin('relion.protocols',
                                            'ProtRelionCtfRefinement', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="4.i CTF Refinement",
                               doCtfFitting=True,
                               fitDefocus=2,
                               fitAstig=2,
                               fitBfactor=2,
                               fitPhaseShift=2,
                               numberOfMpi=8)
    prot.inputParticles.set(protResizeParticles.outputParticles)
    prot.inputPostprocess.set(protPostprocess)
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputParticles')
    waitUntilFinishes(project, prot)

    md1 = xmipp3.MetaData("particles@"+prot._getPath("input_particles.star"))
    md2 = xmipp3.MetaData("particles@"+prot._getExtraPath("particles_ctf_refine.star"))
    md1.sort(xmipp3.RLN_PARTICLE_ID)
    md2.sort(xmipp3.RLN_PARTICLE_ID)

    du1=md1.getColumnValues(xmipp3.RLN_CTF_DEFOCUSU)
    dv1=md1.getColumnValues(xmipp3.RLN_CTF_DEFOCUSV)
    d1=0.5*(np.array(du1)+np.array(dv1))
    a1=np.abs(np.array(du1)-np.array(dv1))

    du2=md2.getColumnValues(xmipp3.RLN_CTF_DEFOCUSU)
    dv2=md2.getColumnValues(xmipp3.RLN_CTF_DEFOCUSV)
    d2 = 0.5 * (np.array(du2) + np.array(dv2))
    a2 = np.abs(np.array(du2) - np.array(dv2))

    ddiff = d1-d2
    adiff = a1-a2

    phase=md2.getColumnValues(xmipp3.RLN_CTF_PHASESHIFT)
    B=md2.getColumnValues(xmipp3.RLN_CTF_BFACTOR)
    s=md2.getColumnValues(xmipp3.RLN_CTF_SCALEFACTOR)

    pD = np.percentile(d2,[2.5, 50, 97.5])
    pA = np.percentile(a2,[2.5, 50, 97.5])
    pDdiff = np.percentile(ddiff,[2.5, 50, 97.5])
    pAdiff = np.percentile(adiff,[2.5, 50, 97.5])
    pPhase = np.percentile(phase,[2.5, 50, 97.5])
    pB = np.percentile(B,[2.5, 50, 97.5])
    pS = np.percentile(s,[2.5, 50, 97.5])

    msg=\
"""The following list shows the median, confidence intervals and links to the histograms for the refined parameters.
Ideally these should all concentrate around 0, except for the defocus and the phase shift that 
must be centered around their true values.

\\begin{center}
    \\begin{tabular}{cccc}
        \\textbf{Parameter} & \\textbf{Median} & \\textbf{95\\%% Confidence interval} & \\textbf{Histogram} \\\\
        \\hline
        Defocus (\\AA) & %5.2f & [%5.1f,%5.1f] & Fig. \\ref{fig:defocus} \\\\
        Astigmatism (\\AA) & %5.2f & [%5.1f,%5.1f] & Fig. \\ref{fig:astigmatism} \\\\
        Defocus difference (\\AA) & %5.2f & [%5.1f,%5.1f] & Fig. \\ref{fig:ddiffHist} \\\\
        Astigmatism difference (\\AA) & %5.2f & [%5.1f,%5.1f] & Fig. \\ref{fig:adiffHist} \\\\
        Phase shift (degrees) & %5.2f & [%5.1f,%5.1f] & Fig. \\ref{fig:phaseHist} \\\\
        B-factor (\\AA$^2$) & %5.2f & [%5.1f,%5.1f] & Fig. \\ref{fig:BHist} \\\\
        Scale factor & %5.3f & [%5.3f,%5.3f] & Fig. \\ref{fig:SHist} \\\\
    \\end{tabular}
\\end{center}

"""%(pD[1], pD[0], pD[2],
     pA[1], pA[0], pA[2],
     pDdiff[1], pDdiff[0], pDdiff[2],
     pAdiff[1], pAdiff[0], pAdiff[2],
     pPhase[1], pPhase[0], pPhase[2],
     pB[1], pB[0], pB[2],
     pS[1], pS[0], pS[2])

    def addHistogram(report, y, ylabel, fnHist, caption, figLabel):
        fnHist = os.path.join(report.getReportDir(),fnHist)
        reportHistogram(y, ylabel, fnHist)
        msg = \
"""
\\begin{figure}[H]
    \centering
    \includegraphics[width=7cm]{%s}
    \\caption{%s}
    \\label{%s}
\\end{figure}

"""%(fnHist, caption, figLabel)
        return msg

    msg += addHistogram(report, d2, "Defocus", "localDHist.png",
                        "Histogram of the defocus after local refinement (\\AA).", "fig:defocus")
    msg += addHistogram(report, a2, "Astigmatism", "localAHist.png",
                        "Histogram of the astigmatism after local refinement (\\AA).", "fig:astigmatism")
    msg += addHistogram(report, ddiff, "Defocus difference", "localDdiffHist.png",
                        "Histogram of the difference in defocus after local refinement (\\AA).", "fig:ddiffHist")
    msg += addHistogram(report, adiff, "Astigmatism difference", "localAdiffHist.png",
                        "Histogram of the difference in astigmatism after local refinement (\\AA).", "fig:adiffHist")
    msg += addHistogram(report, phase, "CTF Phase shift (degrees)", "localPhaseHist.png",
                        "Histogram of the CTF phase shift (degrees).", "fig:phaseHist")
    msg += addHistogram(report, B, "B-factor (\\AA$^2$)", "localBHist.png",
                        "Histogram of the B-factor (\\AA$^2$).", "fig:BHist")
    msg += addHistogram(report, s, "Scale factor", "localSHist.png",
                        "Histogram of the Scale factor.", "fig:SHist")
    report.write(msg)

    warnings = []
    def addWarning(p, ylabel, lowerLimit, upperLimit, check0=True):
        if p[0]<lowerLimit:
            warnings.append("{\\color{red} \\textbf{The lower limit of %s is smaller than accepted lower bound "\
                            "(%5.1f), it is %5.1f.}}"%(ylabel, p[0],lowerLimit))
        if p[2]>upperLimit:
            warnings.append("{\\color{red} \\textbf{The upper limit of %s is smaller than accepted upper bound "\
                            "(%5.1f), it is %5.1f.}}"%(ylabel, p[2],upperLimit))
        if check0 and (p[0]>0 or p[2]<0):
            warnings.append("{\\color{red} \\textbf{The 95\\%% confidence interval of %s is not centered.}}"%ylabel)
    addWarning(pA, "astigmatism", -5000, 5000, check0=False)
    addWarning(pDdiff, "defocus difference", -5000, 5000)
    addWarning(pAdiff, "astigmatism difference", -5000, 5000)
    addWarning(pB, "B-factor", -5, 5)
    addWarning(pS, "scale factor", -0.1, 0.1)
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if 1) astigmatism is between -5000 and 5000; 2) the defocus
difference is between -5000 and 5000; 3) the astigmatism difference is between -5000 and 5000; 4) the B-factor is
between -5 and 5; and 5) the scale factor is between -0.1 and 0.1. 
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "4.i CTF stability", secLabel)
    if len(warnings)>0:
        report.writeAbstract("It seems that there is some problem with the CTF (see Sec. \\ref{%s}). "%secLabel)
    return prot

def level4(project, report, protMap, protMask, protParticles, symmetry, resolution, bfactor,
           protResizeMap, protResizeMask, skipAnalysis = False):

    protResizeParticles = resizeProject(project, protMap, protParticles, resolution)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 4 analysis')
        msg = \
"""This analysis compares the experimental images provided along with their angular assignment to the reconstructed
map."""
        report.write(msg)
        similarityMeasures(project, report, protMap, protMask, protParticles, symmetry, resolution)
        alignabilitySmoothness(project, report, protMap, protMask, protParticles, symmetry, resolution)
        multirefAlignability(project, report, protMap, protMask, protParticles, symmetry)
        protRelion = relionAlignment(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry, resolution)
        protCryoSparc = cryosparcAlignment(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry)
        compareRelionAndCryosparc(project, report, protRelion, protCryoSparc, symmetry)
        relionClassification(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry)
        validateOverfitting(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry, resolution)
        angularDistributionEfficiency(project, report, protResizeParticles, symmetry, resolution, bfactor)
        samplingCompensationFactor(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry)
        ctfStability(project, report, protCryoSparc, protResizeParticles, protResizeMask)
    return protResizeParticles
