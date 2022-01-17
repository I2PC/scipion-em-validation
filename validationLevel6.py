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

from validationReport import reportHistogram

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

def fscq(project, report, protImportMap, protAtom):
    bblCitation = \
"""\\bibitem[Ram{\\'i}rez-Aportela et~al., 2021]{Ramirez2021}
Ram{\\'i}rez-Aportela, E., Maluenda, D., Fonseca, Y.~C., Conesa, P., Marabini,
  R., Heymann, J.~B., Carazo, J.~M., and Sorzano, C. O.~S. (2021).
\\newblock Fsc-q: A cryoem map-to-atomic model quality validation based on the
  local fourier shell correlation.
\\newblock {\em Nature Communications}, 12(1):1--7.
"""
    report.addCitation("Ramirez2021", bblCitation)

    secLabel = "sec:fscq"
    msg = \
"""
\\subsection{Level 6.b FSC-Q}
\\label{%s}
\\textbf{Explanation}:\\\\ 
FSC-Q \\cite{Ramirez2021} compares the local FSC between the map and the atomic model to the local FSC of the two 
half maps.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtConvertPdb', doRaise=True)
    protConvert = project.newProtocol(Prot,
                                      objLabel="6.b Convert",
                                      inputPdbData=1,
                                      sampling=protImportMap.outputVolume.getSamplingRate(),
                                      vol=True)
    protConvert.pdbObj.set(protAtom.outputPdb)
    protConvert.volObj.set(protImportMap.outputVolume)
    project.launchProtocol(protConvert, wait=True)

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtValFit', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="6.b FSC-Q",
                               inputPDB=protAtom.outputPdb.getFileName())
    prot.inputVolume.set(protImportMap.outputVolume)
    prot.pdbMap.set(protConvert.outputVolume)
    project.launchProtocol(prot, wait=True)

    V = xmipp3.Image(prot._getExtraPath("pdb_volume.map"))
    FSCQr = xmipp3.Image(prot._getExtraPath("diferencia_norm.map"))
    fscqr = FSCQr.getData()[V.getData()>0.97]
    avgFSCQr = np.mean(fscqr)
    ci = np.percentile(fscqr, [2.5, 97.5])
    f15=float(np.sum(abs(fscqr)>1.5))/fscqr.size*100

    fnHist = os.path.join(report.getReportDir(),"fscqrHist.png")
    reportHistogram(np.clip(fscqr, -1.5, 1.5), "FSC-Qr", fnHist)

    msg =\
"""Fig. \\ref{fig:fscqHist} shows the histogram of the normalized FSC-Qr (between -1.5 and 1.5) and Fig. 
\\ref{fig:fscq} the colored isosurface of the atomic model converted to map. The
average FSC-Qr is %5.2f, its 95\\%% confidence interval would be [%5.2f,%5.2f]. The percentage of values
whose FSC-Qr absolute value is beyond 1.5 is %5.1f \\%%.

\\begin{figure}[H]
  \\centering
  \\includegraphics[width=8cm]{%s}
  \\caption{Histogram of the FSC-Qr limited to -1.5 and 1.5.}
  \\label{fig:fscqHist}
\\end{figure}

"""%(avgFSCQr, ci[0], ci[1], f15, fnHist)
    report.colorIsoSurfaces(msg, "Isosurface of the atomic model colored by FSC-Qr between -1.5 and 1.5",
                            "fig:fscq", project, "fscq", prot._getExtraPath("pdb_volume.map"),
                            protImportMap.outputVolume.getSamplingRate(),
                            prot._getExtraPath("diferencia_norm.map"), -1.5, 1.5)
    warnings = []
    testWarnings = False
    if f15>10 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of voxels that have a FSC-Qr larger than 1.5 in "\
                        "absolute value is %5.1f, that is larger than 10\\%%}}"%f15)
    report.writeWarningsAndSummary(warnings, "6.b FSC-Q", secLabel)

def multimodel(project, report, protImportMap, protAtom):
    bblCitation = \
"""\\bibitem[Herzik et~al., 2019]{Herzik2019}
Herzik, M.~A., Fraser, J.~S., and Lander, G.~C. (2019).
\\newblock A multi-model approach to assessing local and global cryo-{EM} map
  quality.
\\newblock {\em Structure}, 27(2):344--358.e3.
"""
    report.addCitation("Herzik2019", bblCitation)

    secLabel = "sec:multimodel"
    msg = \
"""
\\subsection{Level 6.c Multimodel stability}
\\label{%s}
\\textbf{Explanation}:\\\\ 
The method of \\cite{Herzik2019} estimates the ambiguity of the atomic model in each region of the CryoEM map due to 
the different local resolutions or local heterogeneity.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    report.writeSummary("6.c Multimodel", secLabel, "{\\color{red} Not in Scipion}")
    report.write("{\\color{red} \\textbf{ERROR: Not in Scipion.}}\\\\ \n")

def guinierModel(project, report, protImportMap, protAtom, resolution):
    map = protImportMap.outputVolume
    fnMap = map.getFileNamme()
    Ts = map.getSamplingRate()

    fnAtom = protAtom.outputVolume.get()
    fnOut = os.path.join(report.getReportDir(), "sharpenedModel.mrc")
    args = "-i %s -o %s --sampling %f --maxres %s --auto"%(fnAtom, fnOut, Ts, resolution)

    scipionHome = getScipionHome()
    scipion3 = os.path.join(scipionHome, 'scipion3')
    output = subprocess.check_output([scipion3, 'run xmipp_volume_correct_bfactor %s' % args])

    p = subprocess.Popen('%s run xmipp_volume_correct_bfactor %s' % (scipion3, args), shell=True,
                         stderr=subprocess.PIPE)

    dinv2, lnFMap, _ = readGuinier(os.path.join(report.getReportDir(), "sharpenedModel.mrc") + ".guinier")
    _, lnFAtom, _ = readGuinier(fnOut + ".guinier")
    p = np.polyfit(lnFmap, lnfAtom, 1)
    lnFMapp = p[0]*lnFMap+p[1]

    fnPlot = os.path.join(report.getReportDir(), 'BfactorAtom.png')
    reportMultiplePlots(dinv2, [lnFAtom, lnFMapp], '1/Resolution^2 (1/A^2)', 'log Structure factor', fnPlot,
                        ['Atomic model', 'Experimental map'])

    secLabel = "sec:bfactorModel"
    msg = \
"""\\subsection{Level 6.d Map-Model Guinier analysis}
\\label{%s}
\\textbf{Explanation:}\\\\
We compared the Guinier plot \\cite{Rosenthal2003} of the atomic and the experimental map model. We made a 
polynomial fitting between both to make sure that they had comparable scales. Ideally, this polynomial should be the 
identity function. \\\\
\\\\
\\textbf{Results:}\\\\
Fig. \\ref{fig:BfactorModel} shows the logarithm (in natural units) of the structure factor (the module squared of the
Fourier transform) of the atom model and the experimental map. The polynomial that fits the experimental map into the
atomic model was $lnF_{fittedMap}=%f lnF_{map} +(%f)$

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Guinier plot of the atom model and experimental map. 
              The X-axis is the square of the inverse of the resolution in \\AA.}
    \\label{fig:BfactorModel}
\\end{figure}

""" % (secLabel, p[0], p[1], fnPlot)
    report.write(msg)

def phenix(project, report, protImportMap, protAtom, resolution):
    bblCitation = \
"""\\bibitem[Barad et~al., 2015]{Barad2015}
Barad, B.~A., Echols, N., Wang, R. Y.-R., Cheng, Y., DiMaio, F., Adams, P.~D.,
  and Fraser, J.~S. (2015).
\\newblock {EMR}inger: side chain-directed model and map validation for 3{D}
  cryo-electron microscopy.
\\newblock {\em Nature Methods}, 12(10):943--946.
"""
    report.addCitation("Barad2015", bblCitation)

    bblCitation = \
"""\\bibitem[Afonine et~al., 2018]{Afonine2018}
Afonine, P.~V., Klaholz, B.~P., Moriarty, N.~W., Poon, B.~K., Sobolev, O.~V.,
  Terwilliger, T.~C., Adams, P.~D., and Urzhumtsev, A. (2018).
\newblock New tools for the analysis and validation of cryo-{EM} maps and
  atomic models.
\newblock {\em Acta Crystallographica D, Struct. Biol.}, 74:814--840.
"""
    report.addCitation("Afonine2018", bblCitation)

    secLabel = "sec:phenix"
    msg = \
"""
\\subsection{Level 6.e Phenix validation}
\\label{%s}
\\textbf{Explanation}:\\\\ 
Phenix provides a number of tools to assess the agreement between the experimental map and its atomic model
\\cite{Afonine2018}. Among them, EMringer \\cite{Barad2015} compares the side chains of the atomic model to the 
CryoEM map.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)


def reportInput(project, report, FNMODEL):
    msg = \
"""
\\section{Atomic model}
Atomic model: %s \\\\
\\\\"""%FNMODEL.replace('_','\_').replace('/','/\-')
    report.write(msg)

    msg = "See Fig. \\ref{fig:modelInput}.\\\\"
    report.atomicModel("modelInput", msg, "Input atomic model", FNMODEL, "fig:modelInput")

def level6(project, report, protImportMap, FNMODEL, resolution, doMultimodel, skipAnalysis=False):
    protAtom = importModel(project, "Import atomic", protImportMap, FNMODEL)

    reportInput(project, report, FNMODEL)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 6 analysis')
        mapq(project, report, protImportMap, protAtom)
        fscq(project, report, protImportMap, protAtom)
        if doMultimodel:
            multimodel(project, report, protImportMap, protAtom)
        guinierModel(project, report, protImportMap, protAtom, resolution)
        phenix(project, report, protImportMap, protAtom, resolution)
