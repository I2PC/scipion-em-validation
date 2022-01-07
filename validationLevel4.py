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

from validationReport import reportHistogram, reportPlot, reportMultiplePlots

def similarityMeasures(project, report, protMap, protMask, protParticles, symmetry, resolution):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtReconstructHighRes', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="4.a Highres local",
                               symmetryGroup=symmetry,
                               alignmentMethod=1,
                               numberOfIterations=1,
                               maximumTargetResolution=max(3,resolution),
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


def level4(project, report, protMap, protMask, protParticles, symmetry, resolution, skipAnalysis = False):
    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 4 analysis')
        similarityMeasures(project, report, protMap, protMask, protParticles, symmetry, resolution)
