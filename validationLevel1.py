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
from validationReport import readMap, latexEnumerate, calculateSha256

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


def level1(project, report, fnMap1, fnMap2, Ts, protImportMap, protCreateMask, skipAnalysis = False):
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

    return protImportMap1, protImportMap2