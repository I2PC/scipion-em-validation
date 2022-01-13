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

from validationReport import plotMicrograph

def importMicrographs(project, label, fnMics, TsMics, kV, Cs, Q0):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportMicrographs', doRaise=True)
    protImport = project.newProtocol(Prot,
                                     objLabel=label,
                                     samplingRate=TsMics,
                                     voltage=kV,
                                     sphericalAberration=Cs,
                                     amplitudeContrast=Q0)
    if fnMics.endswith(".sqlite"):
        protImport.importFrom.set(protImport.IMPORT_FROM_SCIPION)
        protImport.sqliteFile.set(fnMics)
    else:
        protImport.filesPattern.set(fnMics)
    project.launchProtocol(protImport, wait=True)
    if protImport.isFailed():
        raise Exception("Import micrographs did not work")

    return protImport

def extractCoords(project, label, protImportParticles, protMics):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtExtractCoords', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label)
    prot.inputParticles.set(protImportParticles.outputParticles)
    prot.inputMicrographs.set(protMics.outputMicrographs)
    project.launchProtocol(prot, wait=True)
    if prot.isFailed():
        raise Exception("Extract coordinates did not work")

    return prot

def reportInput(project, report, MICPATTERN, protCoords, protMics):
    setOfCoords = protCoords.outputCoordinates
    fnMics = [(x.getObjId(), x.getFileName()) for x in setOfCoords.iterMicrographs()]

    boxSize = setOfCoords.getBoxSize()
    Ts = protMics.outputMicrographs.getSamplingRate()

    msg = \
"""
\\section{Micrographs}
Set of Micrographs: %s \\\\
\\\\
%d micrographs were provided by the user. The first 2 can be seen in Fig. \\ref{fig:micrographs}.

\\begin{figure}[H]
    \centering
"""%(MICPATTERN.replace('_','\_').replace('/','/\-'), len(fnMics))

    for i in range(min(2,len(fnMics))):
        micId, fnMic = fnMics[i]
        coordsmic = [(c.getX(), c.getY()) for c in setOfCoords.iterCoordinates(micId)]
        fnOut = os.path.join(report.getReportDir(),"mic%d.png"%i)
        plotMicrograph(fnMic, fnOut, coordsmic, boxSize, Ts)
        msg+="\includegraphics[width=8cm]{%s}\\\\ \n"%fnOut

    msg+=\
"""    \\caption{Probability density function of the correlation of the user classes compared to the newly computed and
    viceversa.}
    \\label{fig:micrographs}
\\end{figure}
"""
    report.write(msg)

def micCleaner(project, report, label, protCoords):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtDeepMicrographScreen', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label)
    prot.inputCoordinates.set(protCoords.outputCoordinates)
    project.launchProtocol(prot, wait=True)

    bblCitation = \
"""\\bibitem[Sanchez-Garcia et~al., 2020]{Sanchez2020}
Sanchez-Garcia, R., Segura, J., Maluenda, D., Sorzano, C. O.~S., and Carazo,
  J.~M. (2020).
\\newblock {MicrographCleaner}: A python package for cryo-{EM} micrograph
  cleaning using deep learning.
\\newblock {\em J. Structural Biology}, 210:107498."""
    report.addCitation("Sanchez2020", bblCitation)

    secLabel = "sec:micCleaner"
    msg = \
"""
\\subsection{Level 5.a Micrograph cleaner}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method assigns a score between 0 and 1 reflecting the probability that the coordinate is outside a region with
aggregations, ice crystals, carbon edges, etc. \\cite{Sanchez2020} \\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("5.a Micrograph cleaner", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

    return prot

def level5(project, report, protImportParticles, KV, CS, Q0, MICPATTERN, TSMIC, skipAnalysis=False):
    protMics = importMicrographs(project, "Import mics", MICPATTERN, TSMIC, KV, CS, Q0)
    protCoords = extractCoords(project, "Extract coords", protImportParticles, protMics)

    reportInput(project, report, MICPATTERN, protCoords, protMics)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 5 analysis')
        micCleaner(project, report, "5.a MicCleaner", protCoords)
    return protMics, protCoords