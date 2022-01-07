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

from pwem.emlib.metadata import iterRows
import pyworkflow.plugin as pwplugin
from pyworkflow.utils.path import cleanPath
from xmipp3.convert import writeSetOfParticles
import xmipp3

from validationReport import reportHistogram, reportPlot, reportMultiplePlots, readStack

def importParticles(project, label, protImportMap, protImportClasses, fnParticles, TsParticles, kV, Cs, Q0):
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

    XdimMap = protImportMap.outputVolume.getDim()[0]
    TsMap = protImportMap.outputVolume.getSamplingRate()
    AMap = XdimMap * TsMap

    XdimPtcls = protImport.outputParticles.getDim()[0]

    if XdimPtcls==XdimMap and TsParticles==TsMap:
        protResizeMap = protImport
    else:
        XdimPtclsp = int(AMap/TsParticles)
        Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                                'XmippProtCropResizeParticles', doRaise=True)
        protResize1 = project.newProtocol(Prot,
                                          objLabel="Resize Ptcls 1 ",
                                          doWindow=True,
                                          windowOperation=1,
                                          windowSize=XdimPtclsp)
        protResize1.inputParticles.set(protImport.outputParticles)
        project.launchProtocol(protResize1, wait=True)

        protResize2 = project.newProtocol(Prot,
                                          objLabel="Resize and resample Ptcls 2",
                                          doResize=True,
                                          resizeSamplingRate=TsMap,
                                          doWindow=True,
                                          windowOperation=1,
                                          windowSize=XdimMap)
        protResize2.inputParticles.set(protResize1.outputParticles)
        project.launchProtocol(protResize2, wait=True)
        protResizeMap = protResize2

    XdimClasses = protImportClasses.outputAverages.getDim()[0]
    TsClasses = protImportClasses.outputAverages.getSamplingRate()

    if XdimPtcls==XdimClasses and TsParticles==TsClasses:
        protResizeAvgs = protImport
    else:
        AClasses = XdimClasses * TsClasses

        XdimPtclsp = int(AClasses/TsParticles)
        Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                                'XmippProtCropResizeParticles', doRaise=True)
        protResize1 = project.newProtocol(Prot,
                                          objLabel="Resize Ptcls 3",
                                          doWindow=True,
                                          windowOperation=1,
                                          windowSize=XdimPtclsp)
        protResize1.inputParticles.set(protImport.outputParticles)
        project.launchProtocol(protResize1, wait=True)

        protResize2 = project.newProtocol(Prot,
                                          objLabel="Resize and resample Ptcls 4",
                                          doResize=True,
                                          resizeSamplingRate=TsClasses,
                                          doWindow=True,
                                          windowOperation=1,
                                          windowSize=XdimClasses)
        protResize2.inputParticles.set(protResize1.outputParticles)
        project.launchProtocol(protResize2, wait=True)
        protResizeAvgs = protResize2

    return protImport, protResizeMap, protResizeAvgs

def classAnalysis(project, report, protParticles, protClasses):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtStrGpuCrrSimple', doRaise=True)
    protGL2D = project.newProtocol(Prot,
                                   objLabel="3.ab GL2D")
    protGL2D.inputRefs.set(protClasses.outputAverages)
    protGL2D.inputParticles.set(protParticles.outputParticles)
    project.launchProtocol(protGL2D, wait=True)

    bblCitation = \
"""\\bibitem[Sorzano et~al., 2014]{Sorzano2014}
Sorzano, C. O.~S., Vargas, J., de~la Rosa-Trev\\'{\i}n, J.~M.,
  Zald\\'{\i}var-Peraza, A., Ot{\\'o}n, J., Abrishami, V., Foche, I., Marabini,
  R., Caffarena, G., and Carazo, J.~M. (2014).
\\newblock Outlier detection for single particle analysis in electron
  microscopy.
\\newblock In {\em Proc. Intl. Work-Conference on Bioinformatics and Biomedical
  Engineering, IWBBIO}, page 950."""
    report.addCitation("Sorzano2014", bblCitation)

    secLabel = "sec:outlierDetection"
    msg = \
"""
\\subsection{Level 3.a Outlier detection}
\\label{%s}
\\textbf{Explanation}:\\\\ 
The set of particles is classified into the input set of 2D classes of Level 2. The number of particles that are
considered to be outliers in those classes is reported. A particle is an outlier if its Mahalanobis distance to the
centroid of the class is larger than 3 \\cite{Sorzano2014}\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if protGL2D.isFailed():
        report.writeSummary("3.a Outlier detection", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return protGL2D

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCoreAnalysis', doRaise=True)
    protCore = project.newProtocol(Prot,
                                   objLabel="3.ab Core analysis")
    protCore.inputClasses.set(protGL2D.outputClasses)
    project.launchProtocol(protCore, wait=True)

    md = xmipp3.MetaData(protCore._getExtraPath(os.path.join("level_00","level_classes_core.xmd")))
    md.sort(xmipp3.MDL_MODELFRAC)

    fnFracHist = os.path.join(report.getReportDir(),"coreFracHist.png")
    frac = md.getColumnValues(xmipp3.MDL_MODELFRAC)
    reportHistogram(frac,'Core fraction',fnFracHist)

    fnCountHist = os.path.join(report.getReportDir(),"class2DCountHist.png")
    count = md.getColumnValues(xmipp3.MDL_CLASS_COUNT)
    reportHistogram(count,'Image count',fnCountHist)

    toWrite=\
"""The following table shows the input classes, the number of particles assigned to them, and the fraction of these
particles that are considered to be part of the core (the closer to 1, the better). Fig. \\ref{fig:coreFracHist}
shows the histogram of the core fraction of the classes. Fig. \\ref{fig:classCountHist}
shows the histogram of the size of the classes.
 
\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Histogram of the core fraction of the 2D classes.}
    \\label{fig:coreFracHist}
\\end{figure}
 
\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Histogram of the number of particles assigned to the 2D classes.}
    \\label{fig:classCountHist}
\\end{figure}

"""%(fnFracHist, fnCountHist)
    report.write(toWrite)

    report.showj(md,
                 [xmipp3.MDL_IMAGE, xmipp3.MDL_CLASS_COUNT, xmipp3.MDL_MODELFRAC],
                 [True, False, False],
                 ["", "%d ", "%4.3f "],
                 ["2D Class", "No. Particles", "Core fraction"],
                 os.path.join(report.getReportDir(), "core_"), "2cm")

    # Warnings
    warnings=[]
    testWarnings = False
    if np.sum(np.array(frac)<0.7)/len(frac)>0.2 or testWarnings:
        warnings.append("{\\color{red} \\textbf{A large fraction of the 2D classes are rather unstable. In particular, "\
                        "%d classes have a core that is smaller than 70\\%% of the images assigned}}"%\
                        (np.sum(np.array(frac)<0.7)))
    report.writeWarningsAndSummary(warnings, "3.a Outlier detection", secLabel)

    # 3.b internal consistency
    Ts = protParticles.outputParticles.getSamplingRate()
    f05 = np.array(md.getColumnValues(xmipp3.MDL_CLASSIFICATION_FRC_05))
    R05 = Ts*np.reciprocal(f05[f05>0])
    fnR05Hist = os.path.join(report.getReportDir(),"class2DFRC05Hist.png")
    reportHistogram(R05,'Resolution at FRC=0.5 (A)',fnR05Hist)

    fnR05 = os.path.join(report.getReportDir(),"class2DFRC05.png")
    reportPlot(count, f05/Ts, "No. images in the class", "Freq. at which FRC=0.5 (A$^{-1}$)", fnR05,
               plotType='scatter')
    secLabel = "sec:internalConsistency"
    msg = \
"""
\\subsection{Level 3.b Classification internal consistency}
\\label{%s}
\\textbf{Explanation}:\\\\ 
The input particles are classified in 2D clusters. The quality of the 2D clusters is assessed through Fourier Ring Correlation.\\\\
\\textbf{Results:}\\\\
Fig. \\ref{fig:classFRC05Hist} shows the histogram of the resolution of each one of the classes. This resolution
strongly depends on the number of particles assigned to the class, and this server only sees a small fraction of the
particles. Fig. \\ref{fig:classFRC05} shows a scatter plot of the resolution (in \\AA$^{-1}$ in the classes versus the 
number of particles.
\\\\
The following table shows each class, the number of particles assigned to it, and its resolution as measured by
the FRC=0.5.
\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Histogram of the resolution at FRC=0.5 of the different classes.}
    \\label{fig:classFRC05Hist}
\\end{figure}

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Scatter plot of the frequency at which FRC=0.5 (\\AA$^{-1}$) vs the number of particles assigned to
              each class.}
    \\label{fig:classFRC05}
\\end{figure}

""" % (secLabel,fnR05Hist, fnR05)
    report.write(msg)

    for objId in md:
        frc05 = md.getValue(xmipp3.MDL_CLASSIFICATION_FRC_05, objId)
        if frc05>0:
            md.setValue(xmipp3.MDL_CLASSIFICATION_FRC_05, Ts/frc05, objId)
        else:
            md.setValue(xmipp3.MDL_CLASSIFICATION_FRC_05, 999.0, objId)

    md.sort(xmipp3.MDL_CLASSIFICATION_FRC_05)
    report.showj(md,
                 [xmipp3.MDL_IMAGE, xmipp3.MDL_CLASS_COUNT, xmipp3.MDL_CLASSIFICATION_FRC_05],
                 [True, False, False],
                 ["", "%d ", "%5.1f "],
                 ["2D Class", "No. Particles", "Resolution (\\AA)"],
                 os.path.join(report.getReportDir(), "classfrc_"), "2cm")

    # Warnings
    report.writeWarningsAndSummary(None, "3.b 2D Classification internal consistency", secLabel)

def newClassification(project, report, protParticles, protClasses):
    Prot = pwplugin.Domain.importFromPlugin('cryosparc2.protocols',
                                            'ProtCryo2D', doRaise=True)
    protClassif2D = project.newProtocol(Prot,
                                        objLabel="3.c CryoSparc 2D",
                                        numberOfClasses=protClasses.outputAverages.getSize())
    protClassif2D.inputParticles.set(protParticles.outputParticles)
    project.launchProtocol(protClassif2D, wait=True)

    bblCitation = \
"""\\bibitem[Punjani et~al., 2017]{Punjani2017b}
Punjani, A., Brubaker, M.~A., and Fleet, D.~J. (2017).
\\newblock Building proteins in a day: Efficient {3D} molecular structure
  estimation with electron cryomicroscopy.
\\newblock {\em {IEEE} Trans. Pattern Analysis \& Machine Intelligence},
  39:706--718."""
    report.addCitation("Punjani2017b", bblCitation)

    secLabel = "sec:externalConsistency"
    msg = \
"""
\\subsection{Level 3.c Classification external consistency}
\\label{%s}
\\textbf{Explanation}:\\\\ 
The input particles were classified with CryoSparc \\cite{Punjani2017b} using the same number of classes
as the ones provided by the user. Except for the difference in number of particles between the original classification
and the number of particles available to the server, the new classes should resemble the old ones\\\\
\\textbf{Results:}\\\\
""" % (secLabel)
    report.write(msg)

    if protClassif2D.isFailed():
        report.writeSummary("3.c Classification external consistency", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return protClassif2D

    fileList = glob.glob(protClassif2D._getExtraPath("cryosparc_*_class_averages_scaled.mrcs"))
    if len(fileList)==0:
        report.writeSummary("3.c Classification external consistency", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: Cannot find the output classes.}}\\\\ \n")
        return protClassif2D

    fnDensity = os.path.join(report.getReportDir(),"corrDensity.png")
    msg = \
"""
Fig. \\ref{fig:newClassification} shows the new classification. The classification provided by the user is in Fig. 
\\ref{fig:classes2D}.

"""
    report.write(msg)

    report.setOfImages(fileList[0], xmipp3.MDL_IMAGE, "Set of 2D classes calculated by CryoSparc. "\
                       "These should be compared to those in Fig. \\ref{fig:classes2D}.",
                       "fig:newClassification", os.path.join(report.getReportDir(),"newAvg2D_"), "1.5cm", 8)

    avgStack = os.path.join(report.getReportDir(),"avgs.xmd")
    writeSetOfParticles(protClasses.outputAverages, avgStack)
    classesOld = readStack(os.path.join(report.getReportDir(),"avgs.xmd"))
    cleanPath(avgStack)

    classesNew = xmipp3.Image(fileList[0]+":mrc").getData()

    def arrayToList(classes):
        C = []
        Z,_,_ = classes.shape
        for i in range(Z):
            I = xmipp3.Image()
            I.setData(classes[i,:,:])
            C.append(I)
        return C

    classesOld = arrayToList(classesOld)
    classesNew = arrayToList(classesNew)

    def locateMatch(I,classes):
        bestCorr = -1e38
        besti = None
        for i in range(len(classes)):
            corr = I.correlationAfterAlignment(classes[i])
            if corr>bestCorr:
                bestCorr = corr
                besti = i
        return bestCorr, besti

    def compareC1C2(report, C1, C2, C1label, C2label, show):
        msg = None
        if show:
            msg=\
    """The following table shows for each class in %s which is the best match in %s and its correlation coefficient.
    \\begin{longtable}{ccc}
      \\centering
      \\textbf{%s class} & \\textbf{%s class} & \\textbf{Correlation} \\\\ 
    """%(C1label, C2label, C1label, C2label)

        corr = []
        newIdx = []
        for i in range(len(C1)):
            bestCorr, besti = locateMatch(C1[i], C2)

            if show:
                img1 = os.path.join(report.getReportDir(),"classes_%s_%s_%d.jpg"%(C1label,C2label,i))
                img2 = os.path.join(report.getReportDir(),"classes_%s_%s_%d_match.jpg"%(C1label,C2label,i))
                C1[i].write(img1)
                C2[besti].write(img2)
                bestCorrText = "%5.3f"%bestCorr
                if bestCorr<0.8:
                    bestCorrText="{\\color{red} %s}"%bestCorrText
                msg+="  \includegraphics[width=2cm]{%s} & \includegraphics[width=2cm]{%s} & %s \\\\ \n"%\
                     (img1, img2, bestCorrText)

            corr.append(bestCorr)
            newIdx.append(besti)
        if show:
            msg += \
"""\\end{longtable}

"""
        return corr, msg

    corrOld, table = compareC1C2(report, classesOld, classesNew, "User", "New", True)
    corrNew, _     = compareC1C2(report, classesNew, classesOld, "New", "User", False)

    from scipy.stats import gaussian_kde, ks_2samp
    densityOld = gaussian_kde(corrOld)
    densityNew = gaussian_kde(corrNew)
    densityOld.covariance_factor = lambda: .25
    densityOld._compute_covariance()
    densityNew.covariance_factor = lambda: .25
    densityNew._compute_covariance()
    allCorrs = corrOld+corrNew

    fnDensity = os.path.join(report.getReportDir(), "corrDensity.png")
    corraxis = np.linspace(np.min(allCorrs),np.max(allCorrs),200)
    reportMultiplePlots(corraxis, [densityOld(corraxis), densityNew(corraxis)], 'Correlation', 'Prob. Density Function',
                        fnDensity, ['User vs New', 'New vs User'])

    D, pvalue = ks_2samp(corrOld, corrNew)
    msg = \
"""
Fig. \\ref{fig:correlationDensity} shows the probability density function of the correlation of the user classes
compared to the newly computed and viceversa. Ideally, these two distributions should be similar. We compared these
two distributions with a Kolmogorov-Smirnov (KS) two-sample test. The KS statistic was %f and the p-value %f.

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Probability density function of the correlation of the user classes compared to the newly computed and
    viceversa.}
    \\label{fig:correlationDensity}
\\end{figure}

%s
"""%(D, pvalue, fnDensity,table)
    report.write(msg)

    warnings=[]
    testWarnings = False
    if np.min(corrOld)<0.8 or testWarnings:
        warnings.append("{\\color{red} \\textbf{Some user classes correlate less than 0.8 with the newly calculated "\
                        "classes}}")
    if pvalue<0.05 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The equality of the two correlation " \
                        "distributions (user vs new, new vs user) was rejected with a p-value of %f}}"%pvalue)
    report.writeWarningsAndSummary(warnings, "3.c 2D Classification external consistency", secLabel)

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

def level3(project, report, protImportMap, protClasses, fnParticles, TsParticles, kV, Cs, Q0,
           skipAnalysis = False):
    # Import particles
    protParticles, protResizeMap, protResizeAvgs = importParticles(project, "import particles", protImportMap,
                                                                   protClasses, fnParticles, TsParticles, kV, Cs, Q0)
    reportInput(project, report, fnParticles, protParticles)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 3 analysis')
        classAnalysis(project, report, protResizeAvgs, protClasses)
        newClassification(project, report, protResizeAvgs, protClasses)
    return protParticles, protResizeMap, protResizeAvgs