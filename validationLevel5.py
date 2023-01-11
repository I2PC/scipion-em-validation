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

import os

import pyworkflow.plugin as pwplugin

from validationReport import plotMicrograph
from resourceManager import waitOutput, sendToSlurm, waitUntilFinishes

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
    sendToSlurm(protImport)
    project.launchProtocol(protImport)
    #waitOutput(project, protImport, 'outputMicrographs')
    waitUntilFinishes(project, protImport)
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
    sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputCoordinates')
    waitUntilFinishes(project, prot)
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
        msg+="\includegraphics[width=7cm]{%s}\\\\ \n"%fnOut

    msg+=\
"""    \\caption{Two example micrographs with their coordinates.}
    \\label{fig:micrographs}
\\end{figure}
"""
    report.write(msg)

def micCleaner(project, report, label, protCoords):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtDeepMicrographScreen', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               threshold=0.9)
    prot.inputCoordinates.set(protCoords.outputCoordinates)

    sendToSlurm(prot, GPU=True)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputCoordinates_Auto_090')
    waitUntilFinishes(project, prot)

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
This method assigns a score between 0 (bad coordinate) and 1 (good coordinate) reflecting the probability that the
coordinate is outside a region with aggregations, ice crystals, carbon edges, etc. \\cite{Sanchez2020} \\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("5.a Micrograph cleaner", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

    Ninput = protCoords.outputCoordinates.getSize()
    Noutput = prot.outputCoordinates_Auto_090.getSize()

    msg=\
"%d coordinates out of %d (%4.1f \\%%) were scored below 0.9 by MicrographCleaner.\n\n"%(Ninput-Noutput, Ninput,
                                                                                         (Ninput-Noutput)/Ninput*100)
    report.write(msg)

    # Warnings
    warnings = []
    testWarnings = False
    if Noutput < 0.8 * Ninput or testWarnings:
        warnings.append("{\\color{red} \\textbf{More than 20\\% of the particles were flagged as being suspicious.}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if less than 20\\% of the coordinates are suspected to
lie in aggregations, contaminations, ice crystals, etc.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "5.a Micrograph cleaner", secLabel)

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
