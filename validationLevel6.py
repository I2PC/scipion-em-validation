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

def importModel(project, label, protImportMap, FNMODEL):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportPdb', doRaise=True)
    protImport = project.newProtocol(Prot,
                                     objLabel=label,
                                     inputPdbData=1,
                                     pdbFile=FNMODEL)
    protImport.inputVolume.set(protImportMap.outputVolume)
    project.launchProtocol(protImport, wait=True)
    if protImport.isFailed():
        raise Exception("Import micrographs did not work")

    return protImport

def mapq(project, report, protImportMap, protAtom):
    bblCitation = \
"""\\bibitem[Pintilie and Chiu, 2021]{Pintilie2021}
Pintilie, G. and Chiu, W. (2021).
\\newblock Validation, analysis and annotation of cryo-{EM} structures.
\\newblock {\em Acta Crystallographica D, Struct. Biol.}, 77:1142--1152.
"""
    report.addCitation("Pintilie2021", bblCitation)

    secLabel = "sec:mapq"
    msg = \
        """
        \\subsection{Level 6.a MAP-Q}
        \\label{%s}
        \\textbf{Explanation}:\\\\ 
        MAP-Q \\cite{Pintilie2021} computes the local correlation between the map and each one of its atoms assumed to
        have a Gaussian shape.\\\\
        \\\\
        \\textbf{Results:}\\\\
        \\\\
        """ % secLabel
    report.write(msg)
    report.writeSummary("6.a MAP-Q", secLabel, "{\\color{red} Not in Scipion}")
    report.write("{\\color{red} \\textbf{ERROR: Not in Scipion.}}\\\\ \n")

def reportInput(project, report, FNMODEL):
    msg = \
"""
\\section{Atomic model}
Atomic model: %s \\\\
\\\\"""%FNMODEL.replace('_','\_').replace('/','/\-')
    report.write(msg)

    msg = "See Fig. \\ref{fig:modelInput}.\\\\"
    report.atomicModel("modelInput", msg, "Input atomic model", FNMODEL, "fig:modelInput")

def level6(project, report, protImportMap, FNMODEL, skipAnalysis=False):
    protAtom = importModel(project, "Import atomic", protImportMap, FNMODEL)

    reportInput(project, report, FNMODEL)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 6 analysis')
        mapq(project, report, protImportMap, protAtom)
