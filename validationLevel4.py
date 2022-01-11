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

import glob
import math
import numpy as np
import os
import scipy

from pwem.emlib.metadata import iterRows
import pyworkflow.plugin as pwplugin
from pyworkflow.utils.path import cleanPath
from xmipp3.convert import writeSetOfParticles
import xmipp3

from validationReport import reportHistogram, reportPlot, reportMultiplePlots

def resizeProject(project, protMap, protMask, protParticles, resolution):
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
    project.launchProtocol(protResizeParticles, wait=True)

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCropResizeVolumes', doRaise=True)
    protResizeMap = project.newProtocol(Prot,
                                        objLabel="Resize Volume Ts=%2.1f"%TsTarget,
                                        doResize=True,
                                        resizeSamplingRate=TsTarget,
                                        doWindow=True,
                                        windowOperation=1,
                                        windowSize=Xdimp)
    protResizeMap.inputVolumes.set(protMap.outputVolume)
    project.launchProtocol(protResizeMap, wait=True)

    protResizeMask = project.newProtocol(Prot,
                                         objLabel="Resize Mask Ts=%2.1f"%TsTarget,
                                         doResize=True,
                                         resizeSamplingRate=TsTarget,
                                         doWindow=True,
                                         windowOperation=1,
                                         windowSize=Xdimp)
    protResizeMask.inputVolumes.set(protMask.outputMask)
    project.launchProtocol(protResizeMask, wait=True)

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtPreprocessVolumes', doRaise=True)
    protPreprocessMask = project.newProtocol(Prot,
                                             objLabel="Binarize",
                                             doThreshold=True,
                                             thresholdType=1,
                                             threshold=0.5,
                                             fillType=1)
    protPreprocessMask.inputVolumes.set(protResizeMask.outputVol)
    project.launchProtocol(protPreprocessMask, wait=True)

    return protResizeParticles, protResizeMap, protPreprocessMask


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
                               contAngles=False)
    prot.inputVolumes.set(protMap.outputVolume)
    prot.inputParticles.set(protParticles.outputParticles)
    prot.nextMask.set(protMask.outputMask)
    project.launchProtocol(prot, wait=True)

    bblCitation = \
"""\\bibitem[Sorzano et~al., 2015]{Sorzano2015b}
Sorzano, C. O.~S., Vargas, J., {de la Rosa-Trev\\'{\i}n}, J.~M., Ot{\\'o}n, J.,
  {\\'A}lvarez-Cabrera, A.~L., Abrishami, V., Sesmero, E., Marabini, R., and
  Carazo, J.~M. (2015).
\\newblock A statistical approach to the initial volume problem in single
  particle analysis by electron microscopy.
\\newblock {\em J. Structural Biology}, 189:213--219."""
    report.addCitation("Sorzano2015b", bblCitation)

    secLabel = "sec:outlierDetection"
    msg = \
"""
\\subsection{Level 4.a Similarity criteria}
\\label{%s}
\\textbf{Explanation}:\\\\ 
We measured the similarity between the experimental images, with the angles and shifts provided by the user, and
reprojections of the input map along the same direction. We measured the correlation and IMED \\cite{Sorzano2015b}
between both sets of images. If the set of particles is properly assigned there should be a single peak in the
1D histograms of these two similarity measures, and a single cloud in their joint scatter plot. The presence of 
multiple peaks could reveal a mixture of different conformations, the presence of misaligned particles or 
contaminations, or the difference between isolated particles and particles with other objects around.\\\\ 
\\\\
It must be noted that there is a dependence between similarity metrics and defocus. Typically this dependence
is such that lower defoci imply lower similarity due to the smallest contrast. You have to be sure that the groups
seen in the similarity measures are not caused by defocus groups.
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.a Similarity criteria", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
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

def alignabilitySmoothness(project, report, protMap, protMask, protParticles, symmetry):
    bblCitation = \
"""\\bibitem[M{\\'e}ndez et~al., 2021]{Mendez2021b}
M{\\'e}ndez, J., Gardu{\\~n}o, E., Carazo, J.~M., and Sorzano, C. O.~S. (2021).
\\newblock Identification of incorrectly oriented particles in {Cryo-EM} single
  particle analysis.
\\newblock {\em J. Structural Biology}, 213:107771."""
    report.addCitation("Mendez2021b", bblCitation)

    secLabel = "sec:alignabilitySmoothness"
    msg = \
"""
\\subsection{Level 4.b Alignability smoothness}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This algorithm \\cite{Mendez2021b} analyzes the smoothness of the correlation function over the projection sphere and
the stability of its maximum.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    report.writeSummary("4.b Alignability smoothness", secLabel, "{\\color{red} Not in Scipion}")
    report.write("{\\color{red} \\textbf{ERROR: Not in Scipion.}}\\\\ \n")

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
    project.launchProtocol(prot, wait=True)

    bblCitation = \
"""\\bibitem[Vargas et~al., 2017]{Vargas2017}
Vargas, J., Melero, R., G{\\'o}mez-Blanco, J., Carazo, J.~M., and Sorzano, C.
  O.~S. (2017).
\\newblock Quantitative analysis of {3D} alignment quality: its impact on
  soft-validation, particle pruning and homogeneity analysis.
\\newblock {\em Scientific Reports}, 7:6307."""
    report.addCitation("Vargas2017", bblCitation)

    bblCitation = \
"""\\bibitem[Vargas et~al., 2016]{Vargas2016}
Vargas, J., Ot{\\'o}n, J., Marabini, R., Carazo, J.~M., and Sorzano, C. O.~S.
  (2016).
\\newblock Particle alignment reliability in single particle electron
  cryomicroscopy: a general approach.
\\newblock {\em Scientific Reports}, 6:21626."""
    report.addCitation("Vargas2016", bblCitation)

    secLabel = "sec:multirefAlignability"
    msg = \
"""
\\subsection{Level 4.b Alignability}
\\label{%s}
\\textbf{Explanation}:\\\\ 
Alignability precision and accuracy. The precision \\cite{Vargas2016} analyzes the orientation distribution of the
best matching reprojections from the reference volume. If the high values are clustered around the same orientation,
then the precision is close to 1. Otherwise, it is closer to -1. Below 0.5 the best directions tend to be scattered.
The alignability accuracy \\cite{Vargas2017} compares the final angular assignment with the result of a new angular 
assignment. The similarity between both is again encoded between -1 and 1.
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.c Alignability", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
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
    report.writeWarningsAndSummary(warnings, "4.c Alignability", secLabel)

def compareAlignment(project, report, refmap, protRefParticles, protReconstruction, symmetry, label, fnRoot,
                     generateSlices=True):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtAlignVolumeParticles', doRaise=True)
    protAlign = project.newProtocol(Prot,
                                    objLabel="4.d"+label+" align",
                                    symmetryGroup=symmetry)
    protAlign.inputReference.set(refmap)
    protAlign.inputVolume.set(protReconstruction.outputVolume)
    protAlign.inputParticles.set(protReconstruction.outputParticles)
    project.launchProtocol(protAlign, wait=True)

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCompareAngles', doRaise=True)
    protCompare = project.newProtocol(Prot,
                                    objLabel="4.d"+label+" compare",
                                    symmetryGroup=symmetry)
    protCompare.inputParticles1.set(protRefParticles.outputParticles)
    protCompare.inputParticles2.set(protReconstruction.outputParticles)
    project.launchProtocol(protCompare, wait=True)

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
"""Fig. \\ref{fig:comparison%s} shows the shift and angular difference between the alignment given by the user
and the one calculated by %s. The median shift difference was %4.1f\\AA, and the median angular difference
%5.1f degrees. Their corresponding median absolute deviations were %4.1f and %4.1f, respectively. Particles
with a shift difference larger than 5\\AA~or an angular difference larger than 5 degrees would be considered as 
incorrectly assigned in one of the two assignments (the user's or the new one). %4.1f\\%% particles were considered 
to have an uncertain shift, and %4.1f\\%% particles were considered to have an uncertain alignment. 
\\\\

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \includegraphics[width=9cm]{%s}
    \\caption{Top: Shift difference between the alignment given by the user and the one calculated by %s. Bottom:
    Angular difference. The X-axis represents all particles sorted by their difference.}
    \\label{fig:comparison%s}
\\end{figure}
"""%(simplifiedLabel, simplifiedLabel, medShift, medAngle, madShift, madAngle, outliersShift, outliersAngles,
     fnShift, fnAngles, simplifiedLabel, simplifiedLabel)
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
    project.launchProtocol(prot, wait=True)

    bblCitation = \
        """\\bibitem[Scheres, 2012]{Scheres2012}
Scheres, S. H.~W. (2012).
\\newblock A {B}ayesian view on cryo-{EM} structure determination.
\\newblock {\em J. {M}olecular {B}iology}, 415:406--418."""
    report.addCitation("Scheres2012", bblCitation)

    secLabel = "sec:relionAlignment"
    msg = \
"""
\\subsection{Level 4.d1 Relion alignment}
\\label{%s}
\\textbf{Explanation}:\\\\ 
We have performed an independent angular assignment using Relion autorefine \\cite{Scheres2012}. Images were
downsampled to a pixel size of 3\\AA. Then, we measured the difference between the angular assignment of the 
particles given by the user and the one done by Relion.
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.d1 Relion alignment", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

    shiftOutliers, angleOutliers = compareAlignment(project, report, protResizeMap.outputVol,
                                                    protResizeParticles, prot, symmetry,
                                                    "1. Relion", "alignmentRelion")
    warnings = []
    testWarnings = False
    if shiftOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain shift is larger than 20\\%}}")
    if angleOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain angles is larger than 20\\%}}")
    report.writeWarningsAndSummary(warnings, "4.d1 Relion alignment", secLabel)
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
    project.launchProtocol(prot, wait=True)

    bblCitation = \
        """\\bibitem[Punjani et~al., 2017]{Punjani2017b}
Punjani, A., Brubaker, M.~A., and Fleet, D.~J. (2017).
\\newblock Building proteins in a day: Efficient {3D} molecular structure
  estimation with electron cryomicroscopy.
\\newblock {\em {IEEE} Trans. Pattern Analysis \& Machine Intelligence},
  39:706--718."""
    report.addCitation("Punjani2017b", bblCitation)

    secLabel = "sec:cryosparcAlignment"
    msg = \
"""
\\subsection{Level 4.d2 CryoSparc alignment}
\\label{%s}
\\textbf{Explanation}:\\\\ 
We have performed an independent angular assignment using CryoSparc non-homogeneous refinement \\cite{Punjani2017b}.
Images were downsampled to a pixel size of 3\\AA.  Then, we measured the difference between the angular assignment
of the particles given by the user and the one done by CryoSparc.
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.d2 CryoSparc alignment", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

    shiftOutliers, angleOutliers = compareAlignment(project, report, protMap.outputVol, protParticles, prot, symmetry,
                                                    "2. Cryosparc", "alignmentCryosparc")

    warnings = []
    testWarnings = False
    if shiftOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain shift is larger than 20\\%}}")
    if angleOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain angles is larger than 20\\%}}")
    report.writeWarningsAndSummary(warnings, "4.d2 CryoSparc alignment", secLabel)

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
alignability of the input images.
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    shiftOutliers, angleOutliers = compareAlignment(project, report, protRelion.outputVolume, protRelion, protCryoSparc,
                                                    symmetry, "3. Relion/Cryosparc", "alignmentRelionCryosparc",
                                                    False)

    warnings = []
    testWarnings = False
    if shiftOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain shift is larger than 20\\%}}")
    if angleOutliers>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of images with uncertain angles is larger than 20\\%}}")
    report.writeWarningsAndSummary(warnings, "4.d3 Relion/CryoSparc alignments", secLabel)

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
    project.launchProtocol(prot, wait=True)

    bblCitation = \
        """\\bibitem[Scheres, 2012]{Scheres2012}
Scheres, S. H.~W. (2012).
\\newblock A {B}ayesian view on cryo-{EM} structure determination.
\\newblock {\em J. {M}olecular {B}iology}, 415:406--418."""
    report.addCitation("Scheres2012", bblCitation)

    secLabel = "sec:relionClassification"
    msg = \
"""
\\subsection{Level 4.e Classification without alignment}
\\label{%s}
\\textbf{Explanation}:\\\\ 
We have performed a 3D classification of the input particles in two classes without aligning them using Relion
\\cite{Scheres2012}. Images were downsampled to a pixel size of 3\\AA. A valid result would be: 1) a class attracting
most particles and an almost empty class, 2) two classes with an arbitrary number of images in each one, but without
any significant structural difference between the two.
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.e Classification without alignment", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
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
    elif len(Nimgs)==2:
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
    project.launchProtocol(prot, wait=True)

    bblCitation = \
"""\\bibitem[Heymann, 2015]{Heymann2015}
Heymann, B. (2015).
\\newblock Validation of {3DEM} reconstructions: The phantom in the noise.
\\newblock {\em AIMS Biophysics}, 2:21--35."""
    report.addCitation("Heymann2015", bblCitation)

    secLabel = "sec:relionClassification"
    msg = \
"""
\\subsection{Level 4.f Overfitting detection}
\\label{%s}
\\textbf{Explanation}:\\\\ 
The detection of overfitting can be performed through a series of 5 reconstructions with an increasing number
of experimental particles and the same number of pure noise particles \\cite{Heymann2015}. The resolution of
the reconstructions with experimental particles should always be better than those from noise. \\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.f Overfitting detection", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
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
"""%(', '.join(samplePtcls.split()))

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
    if overfitting>30 or testWarnings:
        warnings.append("{\\color{red} \\textbf{the resolution of pure noise particles is sometimes better than the "\
                        "one of true particles.}}")
    report.writeWarningsAndSummary(warnings, "4.f Overfitting detection", secLabel)

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
    project.launchProtocol(prot, wait=True)

    bblCitation = \
"""\\bibitem[Naydenova and Russo, 2017]{Naydenova2017}
Naydenova, K. and Russo, C.~J. (2017).
\\newblock Measuring the effects of particle orientation to improve the
  efficiency of electron cryomicroscopy.
\\newblock {\em Nature communications}, 8:629."""
    report.addCitation("Naydenova2017", bblCitation)

    secLabel = "sec:AngularDistributionEfficiency"
    msg = \
"""
\\subsection{Level 4.g Angular distribution efficiency}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method \\cite{Naydenova2017} evaluates the ability of the angular distribution to fill the Fourier space. 
It determines a resolution
per direction based on the number of particles in each direction and reports the distribution efficiency, a 
number between 0 (inefficient) and 1 (total efficiency).\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("4.g Angular distribution efficiency", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
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
    report.writeWarningsAndSummary(warnings, "4.g Angular distribution efficiency", secLabel)

def samplingCompensationFactor(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry):
    bblCitation = \
"""\\bibitem[Baldwin and Lyumkis, 2020]{Baldwin2020}
Baldwin, P.~R. and Lyumkis, D. (2020).
\\newblock Non-uniformity of projection distributions attenuates resolution in
  {Cryo-EM}.
\\newblock {\em Progress in Biophysics and Molecular Biology}, 150:160--183."""
    report.addCitation("Baldwin2020", bblCitation)

    secLabel = "sec:SCF"
    msg = \
        """
        \\subsection{Level 4.h Sampling compensation factor}
        \\label{%s}
        \\textbf{Explanation}:\\\\ 
        SCF \\cite{Baldwin2020} measures the ability of the angular distribution to fill the Fourier space.\\\\
        \\\\
        \\textbf{Results:}\\\\
        \\\\
        """ % secLabel
    report.write(msg)
    report.writeSummary("4.h SCF", secLabel, "{\\color{red} Not in Scipion}")
    report.write("{\\color{red} \\textbf{ERROR: Not in Scipion.}}\\\\ \n")


def level4(project, report, protMap, protMask, protParticles, symmetry, resolution, bfactor, skipAnalysis = False):
    # Resize to the given resolution
    protResizeParticles, protResizeMap, protResizeMask = resizeProject(project, protMap, protMask, protParticles,
                                                                       resolution)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 4 analysis')
        similarityMeasures(project, report, protMap, protMask, protParticles, symmetry, resolution)
        alignabilitySmoothness(project, report, protMap, protMask, protParticles, symmetry)
        multirefAlignability(project, report, protMap, protMask, protParticles, symmetry)
        protRelion = relionAlignment(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry, resolution)
        protCryoSparc = cryosparcAlignment(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry)
        compareRelionAndCryosparc(project, report, protRelion, protCryoSparc, symmetry)
        relionClassification(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry)
        validateOverfitting(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry, resolution)
        angularDistributionEfficiency(project, report, protResizeParticles, symmetry, resolution, bfactor)
        samplingCompensationFactor(project, report, protResizeMap, protResizeMask, protResizeParticles, symmetry)