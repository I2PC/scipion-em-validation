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
import pickle
import scipy
import subprocess

from scipion.utils import getScipionHome
from pwem.emlib.metadata import iterRows
import pyworkflow.plugin as pwplugin
from pyworkflow.utils.path import cleanPath
from xmipp3.convert import writeSetOfParticles
import xmipp3

from validationReport import reportHistogram, readGuinier, reportMultiplePlots, reportPlot

def importMap(project, label, protImportMap):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportVolumes', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               filesPath=os.path.join(project.getPath(),protImportMap.outputVolume.getFileName()),
                               samplingRate=protImportMap.outputVolume.getSamplingRate())
    project.launchProtocol(prot, wait=True)
    return prot

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

def convertPDB(project, protImportMap, protAtom):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtConvertPdb', doRaise=True)
    protConvert = project.newProtocol(Prot,
                                      objLabel="Convert Pdb to map",
                                      inputPdbData=1,
                                      sampling=protImportMap.outputVolume.getSamplingRate(),
                                      vol=True)
    protConvert.pdbObj.set(protAtom.outputPdb)
    protConvert.volObj.set(protImportMap.outputVolume)
    project.launchProtocol(protConvert, wait=True)
    return protConvert

def fscq(project, report, protImportMap, protAtom, protConvert):
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

def guinierModel(project, report, protImportMap, protConvert, resolution):
    map = protImportMap.outputVolume
    Ts = map.getSamplingRate()

    fnAtom = protConvert.outputVolume.getFileName()
    fnOut = os.path.join(report.getReportDir(), "sharpenedModel.mrc")
    args = "-i %s -o %s --sampling %f --maxres %s --auto"%(fnAtom, fnOut, Ts, resolution)

    scipionHome = getScipionHome()
    scipion3 = os.path.join(scipionHome, 'scipion3')
    output = subprocess.check_output([scipion3, 'run xmipp_volume_correct_bfactor %s' % args])

    p = subprocess.Popen('%s run xmipp_volume_correct_bfactor %s' % (scipion3, args), shell=True,
                         stderr=subprocess.PIPE)

    dinv2, lnFMap, _ = readGuinier(os.path.join(report.getReportDir(), "sharpenedMap.mrc") + ".guinier")
    _, lnFAtom, _ = readGuinier(fnOut + ".guinier")
    lnFMapp = lnFMap+(np.mean(lnFAtom)-np.mean(lnFMap))

    R = np.corrcoef(lnFMapp, lnFAtom)
    R=R[0,1]

    fnPlot = os.path.join(report.getReportDir(), 'BfactorAtom.png')
    reportMultiplePlots(dinv2, [lnFAtom, lnFMapp], '1/Resolution^2 (1/A^2)', 'log Structure factor', fnPlot,
                        ['Atomic model', 'Experimental map'])

    secLabel = "sec:bfactorModel"
    msg = \
"""\\subsection{Level 6.d Map-Model Guinier analysis}
\\label{%s}
\\textbf{Explanation:}\\\\
We compared the Guinier plot \\cite{Rosenthal2003} of the atomic model and the experimental map. We made the mean
of both profiles to be equal (and equal to the mean of the atomic model) to make sure that they had comparable scales. 
\\\\
\\\\
\\textbf{Results:}\\\\
Fig. \\ref{fig:BfactorModel} shows the logarithm (in natural units) of the structure factor (the module squared of the
Fourier transform) of the atom model and the experimental map. The correlation between both profiles was %5.3f.

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Guinier plot of the atom model and experimental map. 
              The X-axis is the square of the inverse of the resolution in \\AA.}
    \\label{fig:BfactorModel}
\\end{figure}

""" % (secLabel, R, fnPlot)
    report.write(msg)

    warnings = []
    testWarnings = False
    if R<0.5 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The correlation is smaller than 0.5, it is %5.3f.}}"%R)
    report.writeWarningsAndSummary(warnings, "6.d Map-Model Guinier", secLabel)


def phenix(project, report, protImportMap, protAtom, resolution):
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
\\cite{Afonine2018}. There are several cross-correlations to assess the quality of the fitting:\\\\
\\begin{itemize}
    \\item CC (mask): Model map vs. experimental map correlation coefficient calculated considering map values inside 
a mask calculated around the macromolecule. 
    \\item CC (box): Model map vs. experimental map correlation coefficient calculated considering all grid points of the 
box.
    \\item CC (volume) and CC (peaks) compare only map regions with the highest density values and regions below a 
    certain contouring threshold level are ignored. CC (volume): The map region considered is defined by
    the N highest points inside the molecular mask. CC (peaks): In this case, calculations consider the union of 
    regions defined by the N highest peaks in the model-calculated map and the N highest peaks in the experimental map.
    \\item Local real-space correlation coefficients CC (main chain) and CC (side chain) involve main skeleton chain 
and lateral chains, respectively.
\\end{itemize}
There are also multiple ways of measuring the resolution:
\\begin{itemize}
    \\item d99: Resolution cutoff beyond which Fourier map coefficients are negligibly small. Calculated from the 
    full map.
    \\item Overall B-iso: Overall isotropic B-value.
    \\item d\_model: Resolution cutoff at which the model map is the most similar to the target (experimental)
 map. For d\_model to be meaningful, model is expected to fit the map as good as possible. d\_model (B\ factors = 0) 
 tries to avoid the blurring of the map.
    \\item d\_FSC\_model; Resolution cutoff up to which the model and map Fourier coefficients are similar at FSC values 
        of 0, 0.143, 0.5.
\\end{itemize}
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('phenix.protocols',
                                            'PhenixProtRunValidationCryoEM', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="6.e Phenix",
                               resolution=max(resolution,3.0))
    prot.inputVolume.set(protImportMap.outputVolume)
    prot.inputStructure.set(protAtom.outputPdb)
    project.launchProtocol(prot, wait=True)

    if prot.isFailed():
        report.writeSummary("6.e Phenix", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

    fnPkl = os.path.join(project.getPath(),prot._getExtraPath("validation_cryoem.pkl"))
    fnPklOut = os.path.join(report.getReportDir(),"validation_cryoem.pkl")
    phenixScript=\
"""import pickle

data=pickle.load(open("%s","r"))
dataOut = {}

# CC ----------------------------------------------------------------
dataOut["cc_mask"]=data.model_vs_data.cc.cc_mask
dataOut["cc_box"]=data.model_vs_data.cc.cc_box
dataOut["cc_peaks"]=data.model_vs_data.cc.cc_peaks
dataOut["cc_volume"]=data.model_vs_data.cc.cc_volume
dataOut["cc_main_chain"]=data.model_vs_data.cc.cc_main_chain.cc
try:
    dataOut["cc_side_chain"]=data.model_vs_data.cc.cc_side_chain.cc
except:
    pass

# CC per chain ------------------------------------------------------
dataOut['chain_list'] = []
chain_names = []
for item in data.model_vs_data.cc.cc_per_chain:
    if item.chain_id not in chain_names:
        chain_names.append(item.chain_id)
        dataOut['chain_list'].append((item.chain_id, item.cc))

dataOut['resseq_list'] = {}
for chain_name in chain_names:
    resseq_list = []
    residue_cc = []
    for item in data.model_vs_data.cc.cc_per_residue:
        if item.chain_id == chain_name:
            resseq_list.append(item.resseq)
            residue_cc.append(item.cc)
    dataOut['resseq_list'][chain_name]=(resseq_list,residue_cc)

# Resolutions ------------------------------------------------------
dataOut['*d99_full_masked']=data.data.masked.d99
dataOut["overall_b_iso_masked"]=data.data.masked.b_iso_overall
dataOut['*dmodel_masked']=data.data.masked.d_model
dataOut['d_model_b0_masked']=data.data.masked.d_model_b0
dataOut['*dFSCmodel_0_masked']=data.data.masked.d_fsc_model_0
dataOut['*dFSCmodel_0.143_masked']=data.data.masked.d_fsc_model_0143
dataOut['*dFSCmodel_0.5_masked']=data.data.masked.d_fsc_model_05
dataOut['mask_smoothing_radius']=data.data.masked.radius_smooth

dataOut['*d99_full_unmasked']=data.data.unmasked.d99
dataOut["overall_b_iso_unmasked"]=data.data.unmasked.b_iso_overall
dataOut['*dmodel_unmasked']=data.data.unmasked.d_model
dataOut['d_model_b0_unmasked']=data.data.unmasked.d_model_b0
dataOut['*dFSCmodel_0_unmasked']=data.data.unmasked.d_fsc_model_0
dataOut['*dFSCmodel_0.143_unmasked']=data.data.unmasked.d_fsc_model_0143
dataOut['*dFSCmodel_0.5_unmasked']=data.data.unmasked.d_fsc_model_05

# FSCs -------------------------------------------------------------
if data.data.masked.fsc_curve_model.fsc is not None:
    fsc_model_map_masked = []
    d_inv_model_map_masked = []
    for item in data.data.masked.fsc_curve_model.fsc:
        fsc_model_map_masked.append(item)
    dataOut['FSC_Model_Map_Masked'] = fsc_model_map_masked
    for item in data.data.masked.fsc_curve_model.d_inv:
        d_inv_model_map_masked.append(item)
    dataOut['d_inv_Model_Map_Masked'] = d_inv_model_map_masked
if data.data.unmasked.fsc_curve_model.fsc is not None:
    fsc_model_map_unmasked = []
    d_inv_model_map_unmasked = []
    for item in data.data.unmasked.fsc_curve_model.fsc:
        fsc_model_map_unmasked.append(item)
    dataOut['FSC_Model_Map_Unmasked'] = fsc_model_map_unmasked
    for item in data.data.unmasked.fsc_curve_model.d_inv:
        d_inv_model_map_unmasked.append(item)
    dataOut['d_inv_Model_Map_Unmasked'] = d_inv_model_map_unmasked
        
fh = open("%s",'wb')
pickle.dump(dataOut,fh)
fh.close()
"""%(fnPkl, fnPklOut)
    fnPhenixScript = os.path.join(report.getReportDir(),"validation_cryoem.py")
    fhPhenixScript = open(fnPhenixScript,"w")
    fhPhenixScript.write(phenixScript)
    fhPhenixScript.close()

    from phenix import Plugin
    Plugin.runPhenixProgram('',fnPhenixScript)

    data = pickle.load(open(fnPklOut, "rb"))

    # CC
    msg =\
"""To avoid ringing in Fourier space a smooth mask with a radius of %5.1f \\AA~has been applied.  \\\\
\\underline{Overall correlation coefficients}: \\\\
\\begin{tabular}{rc}
CC (mask) = & %5.3f\\\\
CC (box) = & %5.3f\\\\
CC (volume) = & %5.3f\\\\
CC (peaks) = & %5.3f\\\\
CC (main chain) = & %5.3f\\\\
"""%(data['mask_smoothing_radius'], data['cc_mask'],data['cc_box'],data['cc_volume'],data['cc_peaks'],
     data['cc_main_chain'])
    if 'cc_side_chain' in data:
        msg+="CC (side chain) = & %5.3f\\\\ \n"%data['cc_side_chain']
    msg+="\\end{tabular}\n\\\\\n"

    # CC per chain
    msg+=\
"""
\\underline{Correlation coefficients per chain}:\\\\
\\begin{tabular}{cc}
    \\textbf{Chain} & \\textbf{Cross-correlation} \\\\
"""
    for chain_id, cc in data['chain_list']:
        msg+="%s & %f\\\\ \n"%(chain_id, cc)
    msg+="\\end{tabular}\n\n"

    # CC per residues
    def plotCCResidue(chain_id, reportDir):
        fnPlot = os.path.join(reportDir, "ccresidue_%s.png"%chain_id)
        resseq_list, residue_cc = data['resseq_list'][chain_id]
        x = [x+1 for x in np.arange(0,len(residue_cc))]
        print('x',x)
        print('y',residue_cc)
        reportPlot(x, residue_cc, 'Aminoacid no.', 'Cross-correlation', fnPlot, addMean=True, title="Chain %s"%chain_id)
        return fnPlot

    msg+="""We now show the correlation profiles of the different chain per residue.\n"""
    for chain_id in data['resseq_list']:
        fnPlot = plotCCResidue(chain_id, report.getReportDir())
        msg+="""\\includegraphics[width=7cm]{%s}\n"""%fnPlot
        
    # Resolutions
    msg+=\
"""
\\underline{Resolutions estimated from the model}:\\\\
"""
    report.write(msg)

def emringer(project, report, protImportMap, protAtom):
    bblCitation = \
"""\\bibitem[Barad et~al., 2015]{Barad2015}
Barad, B.~A., Echols, N., Wang, R. Y.-R., Cheng, Y., DiMaio, F., Adams, P.~D.,
  and Fraser, J.~S. (2015).
\\newblock {EMR}inger: side chain-directed model and map validation for 3{D}
  cryo-electron microscopy.
\\newblock {\em Nature Methods}, 12(10):943--946.
"""
    report.addCitation("Barad2015", bblCitation)

    secLabel = "sec:emringer"
    msg = \
"""
\\subsection{Level 6.f EMRinger validation}
\\label{%s}
\\textbf{Explanation}:\\\\ 
EMringer \\cite{Barad2015} compares the side chains of the atomic model to the CryoEM map.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('phenix.protocols',
                                            'PhenixProtRunEMRinger', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="6.f EMRinger")
    prot.inputVolume.set(protImportMap.outputVolume)
    prot.inputStructure.set(protAtom.outputPdb)
    project.launchProtocol(prot, wait=True)

    if prot.isFailed():
        report.writeSummary("6.f EMRinger", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

def daq(project, report, protImportMap, protAtom):
    bblCitation = \
"""\\bibitem[Terashi et~al., 2022]{Terashi2022}
Terashi, G., Wang, X., Subramaniya, S.R.M.V., Tesmer, J.J.G. and Kihara, D. (2022).
\\newblock Residue-Wise Local Quality Estimation for Protein Models from {Cryo-EM} Maps.
\\newblock (submitted).
"""
    report.addCitation("Terashi2022", bblCitation)

    secLabel = "sec:daq"
    msg = \
"""
\\subsection{Level 6.g DAQ validation}
\\label{%s}
\\textbf{Explanation}:\\\\ 
DAQ \\cite{Terashi2022} DAQ is a computational tool using deep learning that can estimate the residue-wise local
quality for protein models from cryo-Electron Microscopy (EM) maps. The method calculates the likelihood that a given
density feature corresponds to an aminoacid, atom, and secondary structure. These likelihoods are combined into a
 score that ranges from -1 (bad quality) to 1 (good quality). \\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('kiharalab.protocols',
                                            'ProtDAQValidation', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="6.g DAQ")
    prot.inputVolume.set(protImportMap.outputVolume)
    prot.inputAtomStruct.set(protAtom.outputPdb)
    project.launchProtocol(prot, wait=True)

    if prot.isFailed():
        report.writeSummary("6.f EMRinger", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

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
    if protImportMap.outputVolume.hasHalfMaps():
        protImportMapWOHalves = importMap(project, "Import map2", protImportMap)
        protImportForPhenix = protImportMapWOHalves
    else:
        protImportForPhenix = protImportMap

    protAtom = importModel(project, "Import atomic", protImportMap, FNMODEL)

    reportInput(project, report, FNMODEL)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 6 analysis')
        protConvert = convertPDB(project, protImportMap, protAtom)
        # mapq(project, report, protImportMap, protAtom)
        # fscq(project, report, protImportMap, protAtom, protConvert)
        # if doMultimodel:
        #     multimodel(project, report, protImportMap, protAtom)
        # guinierModel(project, report, protImportMap, protConvert, resolution)
        phenix(project, report, protImportForPhenix, protAtom, resolution)
        # emringer(project, report, protImportForPhenix, protAtom)
        # daq(project, report, protImportMap, protAtom)
