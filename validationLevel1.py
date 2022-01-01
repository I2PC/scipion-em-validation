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
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy

import pyworkflow.plugin as pwplugin
from validationReport import readMap, latexEnumerate, calculateSha256, CDFFromHistogram, CDFpercentile, reportPlot
import xmipp3

def importMap(project, label, fnMap, Ts):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportVolumes', doRaise=True)
    fnDir, fnBase = os.path.split(fnMap)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               filesPath=fnDir,
                               filesPattern=fnMap,
                               samplingRate=Ts)
    project.launchProtocol(prot, wait=True)
    return prot

def blocres(project, report, label, map, mask):
    bblCitation = \
"""\\bibitem[Cardone et~al., 2013]{Cardone2013}
Cardone, G., Heymann, J.~B., and Steven, A.~C. (2013).
\\newblock One number does not fit all: Mapping local variations in resolution
  in cryo-em reconstructions.
\\newblock {\em J. Structural Biology}, 184:226--236."""
    report.addCitation("Cardone2013", bblCitation)

    secLabel = "sec:blocres"
    msg = \
"""
\\subsection{Level 1.a Local resolution with Blocres}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method \\cite{Cardone2013} computes a local Fourier Shell Correlation (FSC) between the two half maps.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    report.writeSummary("1.a Blocres", secLabel, "{\\color{red} Binary installation fails}")
    report.write("{\\color{red} \\textbf{ERROR: Binary installation fails.}}\\\\ \n")

def resmap(project, report, label, map, mask):
    bblCitation = \
"""\\bibitem[Kucukelbir et~al., 2014]{Kucukelbir2014}
Kucukelbir, A., Sigworth, F.~J., and Tagare, H.~D. (2014).
\\newblock Quantifying the local resolution of cryo-{EM} density maps.
\\newblock {\em Nature Methods}, 11:63--65."""
    report.addCitation("Kucukelbir2014", bblCitation)

    secLabel = "sec:resmap"
    msg = \
"""
\\subsection{Level 1.b Local resolution with Resmap}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method \\cite{Kucukelbir2014} is based on a test hypothesis testing of the superiority of signal over noise at different frequencies.\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    report.writeSummary("1.b Resmap", secLabel, "{\\color{red} Not fully automatic}")
    report.write("{\\color{red} \\textbf{ERROR: Not fully automatic.}}\\\\ \n")

def monores(project, report, label, protImportMap, protCreateMask, resolution):
    Ts = protImportMap.outputVolume.getSamplingRate()

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtMonoRes', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               useHalfVolumes=True,
                               minRes=2*Ts,
                               maxRes=max(10,5*resolution))
    prot.associatedHalves.set(protImportMap.outputVolume)
    prot.maskExcl.set(protCreateMask.outputMask)
    project.launchProtocol(prot, wait=True)

    bblCitation = \
"""\\bibitem[Vilas et~al., 2018]{Vilas2018}
Vilas, J.~L., G{\\'o}mez-Blanco, J., Conesa, P., Melero, R., de~la
  Rosa~Trev\\'{\i}n, J.~M., Ot{\\'o}n, J., Cuenca, J., Marabini, R., Carazo,
  J.~M., Vargas, J., and Sorzano, C. O.~S. (2018).
\\newblock {MonoRes}: automatic and unbiased estimation of local resolution for
  electron microscopy maps.
\\newblock {\em Structure}, 26:337--344."""
    report.addCitation("Vilas2018", bblCitation)

    secLabel = "sec:monores"
    msg = \
"""
\\subsection{Level 1.c Local resolution with MonoRes}
\\label{%s}
\\textbf{Explanation}:\\\\ 
MonoRes \\cite{Vilas2018} This methods evaluates the local energy of a point with respect to the distribution of 
energy in the noise. This comparison is performed at multiple frequencies and for each one, the monogenic 
transformation separates the amplitude and phase of the input map.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("1.c MonoRes", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

    md = xmipp3.MetaData()
    md.read(prot._getExtraPath("hist.xmd"))
    x_axis = md.getColumnValues(xmipp3.MDL_X)
    y_axis = md.getColumnValues(xmipp3.MDL_COUNT)

    fnHistMonoRes = os.path.join(report.getReportDir(), "histMonoRes.png")
    reportPlot(x_axis[:-2], y_axis[:-2], 'Resolution (A)', '# of voxels', fnHistMonoRes, plotType="bar",
               barWidth=(x_axis[-1] - x_axis[0]) / len(x_axis))

    R, RCDF=CDFFromHistogram(x_axis[:-2], y_axis[:-2])
    Rpercentiles = CDFpercentile(R, RCDF, Fp=[0.025, 0.25, 0.5, 0.75, 0.975])
    resolutionP = CDFpercentile(R, RCDF, xp=resolution)

    toWrite=\
"""
Fig. \\ref{fig:histMonores} shows the histogram of the local resolution according to MonoRes. Some representative
percentiles are:

\\begin{center}
    \\begin{tabular}{|c|c|}
        \\hline
        \\textbf{Percentile} & Resolution(\AA) \\\\
        \\hline
        2.5\\%% & %5.2f \\\\
        \\hline
        25\\%% & %5.2f \\\\
        \\hline
        50\\%% & %5.2f \\\\
        \\hline
        75\\%% & %5.2f \\\\
        \\hline
        97.5\\%% & %5.2f \\\\
        \\hline
    \\end{tabular}
\\end{center}

The reported resolution, %5.2f \AA, is at the percentile %5.2f. 
Fig. \\ref{fig:monoresColor} shows some representative views of the local resolution

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Histogram of the local resolution according to MonoRes.}
    \\label{fig:histMonores}
\\end{figure}

"""%(Rpercentiles[0], Rpercentiles[1], Rpercentiles[2], Rpercentiles[3], Rpercentiles[4], resolution, resolutionP,
     fnHistMonoRes)
    report.write(toWrite)

    report.colorIsoSurfaces("", "Local resolution according to Monores.", "fig:monoresColor",
                            project, "monoresViewer", protImportMap.outputVolume.getFileName(),
                            Ts, prot._getExtraPath("monoresResolutionChimera.mrc"), Rpercentiles[0], Rpercentiles[-1])

    # Warnings
    warnings=[]
    testWarnings = False
    if resolutionP<0.05 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \AA, is particularly with respect "\
                        "to the local resolution distribution. It occupies the %5.2f percentile}}"%\
                        (resolution,resolutionP))
    report.writeWarningsAndSummary(warnings, "1.c Monores", secLabel)
    return prot

def reportInput(project, report, fnMap1, fnMap2, protImportMap1, protImportMap2):
    toWrite=\
"""
\\section{Half maps}
Half map 1: %s \\\\
SHA256 hash: %s \\\\ 
\\\\
Half map 2: %s \\\\
SHA256 hash: %s \\\\ 
\\\\
Slices of the first half map can be seen in Fig. \\ref{fig:maxVarHalf1}.\\\\
Slices of the second half map can be seen in Fig. \\ref{fig:maxVarHalf2}.\\\\
Slices of the difference between both maps can be seen in Fig. \\ref{fig:maxVarHalfDiff}. There should not be 
any structure in this difference. Sometimes some patterns are seen if the map is symmetric.\\\\
"""%(fnMap1.replace('_','\_').replace('/','/\-'), calculateSha256(fnMap1),\
     fnMap2.replace('_','\_').replace('/','/\-'), calculateSha256(fnMap2))
    report.write(toWrite)

    fnMap1 = os.path.join(project.getPath(), protImportMap1.outputVolume.getFileName())
    fnMap2 = os.path.join(project.getPath(), protImportMap2.outputVolume.getFileName())
    report.orthogonalSlices("half1", "", "Slices of maximum variation in the three dimensions of Half 1", fnMap1,
                            "fig:maxVarHalf1", maxVar=True)
    report.orthogonalSlices("half2", "", "Slices of maximum variation in the three dimensions of Half 2", fnMap2,
                            "fig:maxVarHalf2", maxVar=True)

    V1 = readMap(fnMap1).getData()
    V2 = readMap(fnMap2).getData()
    Vdiff = V1-V2
    report.orthogonalSlices("halfDiff", "",
                            "Slices of maximum variation in the three dimensions of the difference Half1-Half2.", Vdiff,
                            "fig:maxVarHalfDiff", maxVar=True)


def level1(project, report, fnMap1, fnMap2, Ts, resolution, protImportMap, protCreateMask, skipAnalysis = False):
    # Import maps
    protImportMap1 = importMap(project, "import half1", fnMap1, Ts)
    if protImportMap1.isFailed():
        raise Exception("Import map did not work")
    protImportMap2 = importMap(project, "import half2", fnMap2, Ts)
    if protImportMap2.isFailed():
        raise Exception("Import map did not work")
    reportInput(project, report, fnMap1, fnMap2, protImportMap1, protImportMap2)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 1 analysis')
        blocres(project, report, "1.a Blocres", protImportMap, protCreateMask)
        resmap(project, report, "1.b Resmap", protImportMap, protCreateMask)
        monores(project, report, "1.c Monores", protImportMap, protCreateMask, resolution)

    return protImportMap1, protImportMap2