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

import collections
import glob
import json
import numpy as np
import os
import pickle
import subprocess
import re
from datetime import datetime
from random import randint
from time import sleep

from scipion.utils import getScipionHome
import pyworkflow.plugin as pwplugin
from pyworkflow.utils.path import cleanPath, copyFile
from pwem.convert.atom_struct import AtomicStructHandler
from pwem.viewers.viewer_localres import replaceOcuppancyWithAttribute, makeResidueValuesDic
import pwem.convert.atom_struct
import xmipp3

from validationReport import reportHistogram, readGuinier, reportMultiplePlots, reportPlot, isHomogeneous
from resourceManager import waitOutput, sendToSlurm, waitUntilFinishes, createScriptForSlurm, checkIfJobFinished

import configparser

from tools.utils import saveIntermediateData, getFilename, getScoresFromWS, getFileFromWS
from tools.emv_utils import convert_2_json

from resources.constants import *

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'config.yaml'))
useSlurm = config['QUEUE'].getboolean('USE_SLURM')
chimeraProgram = config['MAPQ'].get('CHIMERA_PROGRAM_PATH')
mapq_path = config['MAPQ'].get('MAPQ_PATH')
validation_tools_path = config['EM-VALIDATION'].get('VALIDATION_TOOLS_PATH')
EMDB_entries_path = config['EMDB'].get('ENTRIES_PATH')


def importMap(project, label, protImportMap, mapCoordX, mapCoordY, mapCoordZ, priority=False):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportVolumes', doRaise=True)

    if mapCoordX is not None and mapCoordY is not None and mapCoordZ is not None:
        prot = project.newProtocol(Prot,
                                   objLabel=label,
                                   filesPath=os.path.join(project.getPath(),protImportMap.outputVolume.getFileName()),
                                   samplingRate=protImportMap.outputVolume.getSamplingRate(),
                                   setOrigCoord=True,
                                   x=mapCoordX,
                                   y=mapCoordY,
                                   z=mapCoordZ)
    else:
        prot = project.newProtocol(Prot,
                                   objLabel=label,
                                   filesPath=os.path.join(project.getPath(),protImportMap.outputVolume.getFileName()),
                                   samplingRate=protImportMap.outputVolume.getSamplingRate(),
                                   setOrigCoord=False)
    if useSlurm:
        sendToSlurm(prot, priority=True if priority else False)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputVolume')
    waitUntilFinishes(project, prot)
    return prot

def importModel(project, report, label, protImportMap, fnPdb, priority=False):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportPdb', doRaise=True)
    protImport = project.newProtocol(Prot,
                                     objLabel=label,
                                     inputPdbData=1,
                                     pdbFile=fnPdb)
    protImport.inputVolume.set(protImportMap.outputVolume)
    if useSlurm:
        sendToSlurm(protImport, priority=True if priority else False)
    project.launchProtocol(protImport)
    #waitOutput(project, protImport, 'outputPdb')
    waitUntilFinishes(project, protImport)
    if protImport.isFailed():
        raise Exception("Import atomic model did not work")
    if protImport.isAborted():
        raise Exception("Import atomic model was MANUALLY ABORTED")
    saveIntermediateData(report.fnReportDir, 'inputData', True, 'atomic model', os.path.join(project.getPath(), fnPdb.split('/')[-1].replace('.pdb', '.cif')), 'atomic model')

    return protImport

def moveOriginTo(newOrigin, handler):
    centerMass = handler.centerOfMass(geometric=True)
    for atom in handler.getStructure().get_atoms():
        coords = atom.get_coord()
        atom.coord = coords + np.asarray(newOrigin) - np.asarray(centerMass)

def phenixExecution(project, report, protImportMap, protAtom, resolution, label, priority=False):
    Prot = pwplugin.Domain.importFromPlugin('phenix.protocols',
                                            'PhenixProtRunValidationCryoEM', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               resolution=max(resolution,3.0))
    prot.inputVolume.set(protImportMap.outputVolume)
    prot.inputStructure.set(protAtom.outputPdb)
    if useSlurm:
        sendToSlurm(prot, priority=True if priority else False)
    project.launchProtocol(prot)
    waitUntilFinishes(project, prot)

    fnPkl = os.path.join(project.getPath(),prot._getExtraPath("validation_cryoem.pkl"))
    data = None

    if os.path.isfile(fnPkl):
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

        saveIntermediateData(report.getReportDir(), 'phenix', True, 'validation_cryoem.pkl', fnPkl, 'validation_cryoem.pkl file')
        saveIntermediateData(report.getReportDir(), 'phenix', True, 'validation_cryoem.py', fnPhenixScript, 'validation_cryoem.py file to get all phenix data from pickle')

        from phenix import Plugin
        Plugin.runPhenixProgram('',fnPhenixScript)

        data = pickle.load(open(fnPklOut, "rb"))
        saveIntermediateData(report.getReportDir(), 'phenix', False, 'dataDict', data, ['', 'phenix data dictionary containing all key params'])

    return prot, data


def phenixReporting(project, report, resolution, prot, data):

    secLabel = "sec:phenix"
    msg = \
"""
\\subsection{Level A.a Phenix validation}
\\label{%s}
\\textbf{Explanation}:\\\\ 
Phenix provides a number of tools to assess the agreement between the experimental map and its atomic model
(see this \\href{%s}{link} for more details). There are several cross-correlations to assess the quality of the fitting:\\\\
\\begin{itemize}
    \\item CC (mask): Model map vs. experimental map correlation coefficient calculated considering map values inside 
a mask calculated around the macromolecule. 
    \\item CC (box): Model map vs. experimental map correlation coefficient calculated considering all grid points of the 
box.
    \\item CC (volume) and CC (peaks) compare only map regions with the highest density values and regions below a 
    certain contouring threshold level are ignored. CC (volume): The map region considered is defined by
    the N highest points inside the molecular mask. CC (peaks): In this case, calculations consider the union of 
    regions defined by the N highest peaks in the model-calculated map and the N highest peaks in the experimental map.
    \\item Local real-space correlation coefficients CC (main chain) and CC (side chain) involve the main skeleton chain 
and side chains, respectively.
\\end{itemize}
There are also multiple ways of measuring the resolution:
\\begin{itemize}
    \\item d99: Resolution cutoff beyond which Fourier map coefficients are negligibly small. Calculated from the 
    full map.
    \\item d\_model: Resolution cutoff at which the model map is the most similar to the target (experimental)
 map. For d\_model to be meaningful, the model is expected to fit the map as well as possible. d\_model (B\ factors = 0) 
 tries to avoid the blurring of the map.
    \\item d\_FSC\_model; Resolution cutoff up to which the model and map Fourier coefficients are similar at FSC values 
        of 0, 0.143, 0.5.
\\end{itemize}
In addition to these resolution measurements the overall isotropic B factor is another indirect measure of the
quality of the map.
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, PHENIX_DOI)
    report.write(msg)

    if prot.isFailed():
        report.writeSummary("A.a Phenix", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED)
        phenixStdout = open(os.path.join(project.getPath(), prot.getStdoutLog()), "r").read()
        controlledErrors = ["Sorry: Input map is all zero after boxing", "Sorry: Map and model are not aligned", "Sorry: Fatal problems interpreting model file"]
        for error in controlledErrors:
            if error in phenixStdout:
                report.write("{\\color{red} \\textbf{REASON: %s.}}\\\\ \n" % error)
        report.write(STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_PHENIX)
        report.writeSummary("A.a Phenix", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot


    # CC
    msg =\
"""To avoid ringing in Fourier space a smooth mask with a radius of %5.1f \\AA~has been applied.  \\\\
\\underline{Overall correlation coefficients}: \\\\
\\\\
\\begin{center}
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
    msg+="\\end{center}\n\n"

    # CC per chain
    msg+=\
"""
\\underline{Correlation coefficients per chain}:\\\\
\\begin{center}
\\begin{tabular}{cc}
    \\textbf{Chain} & \\textbf{Cross-correlation} \\\\
"""
    for chain_id, cc in data['chain_list']:
        msg+="%s & %f\\\\ \n"%(chain_id, cc)
    msg+="\\end{tabular}\n\n"
    msg+="\\end{center}\n\n\n"

    # CC per residues
    allCCs = []
    def plotCCResidue(chain_id, reportDir, allCCs):
        fnPlot = os.path.join(reportDir, "ccresidue_%s.png"%chain_id)
        resseq_list, residue_cc = data['resseq_list'][chain_id]
        allCCs+=residue_cc
        x = [x+1 for x in np.arange(0,len(residue_cc))]
        reportPlot(x, residue_cc, 'Aminoacid no.', 'Cross-correlation', fnPlot, addMean=True, title="Chain %s"%chain_id)
        return fnPlot

    msg+="""We now show the correlation profiles of the different chain per residue.\n"""
    for chain_id in sorted(data['resseq_list']):
        fnPlot = plotCCResidue(chain_id, report.getReportDir(), allCCs)
        saveIntermediateData(report.fnReportDir, "phenix", True, "ccresidue_%s.png"%chain_id, fnPlot, 'Plot including the correlation profiles of the chain %s'%chain_id)
        msg+="""\\includegraphics[width=7cm]{%s}\n"""%fnPlot

    fnCCHist = os.path.join(report.getReportDir(),"ccModelHist.png")
    reportHistogram(allCCs, "Cross-correlation", fnCCHist)
    badResidues = np.sum(np.array(allCCs)<0.5)/len(allCCs)*100

    msg += \
"""

Fig. \\ref{fig:ccResidueHist} shows the histogram of all cross-correlations evaluated at the residues. The percentage
of residues whose correlation is below 0.5 is %4.1f \\%%.

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Histogram of the cross-correlation between the map and model evaluated for all residues.}
    \\label{fig:ccResidueHist}
\\end{figure}

"""%(badResidues, fnCCHist)

    saveIntermediateData(report.fnReportDir, "phenix", True, "ccModelHist.png", fnCCHist, 'Histogram of the cross-correlation between the map and model evaluated for all residues')
    saveIntermediateData(report.fnReportDir, "phenix", False, "percentageResidues05", badResidues, ['%', 'The percentage of residues whose correlation is below 0.5'])

    # Resolutions
    msg+=\
"""
\\underline{Resolutions estimated from the model}:\\\\
\\begin{center}
\\begin{tabular}{rcc}
    \\textbf{Resolution} (\\AA) & \\textbf{Masked} & \\textbf{Unmasked} \\\\
    d99 & %4.1f & %4.1f \\\\
    d\_model & %4.1f & %4.1f \\\\
    d\_model (B-factor=0) & %4.1f & %4.1f \\\\
    FSC\_model=0 & %4.1f & %4.1f \\\\
    FSC\_model=0.143 & %4.1f & %4.1f \\\\
    FSC\_model=0.5 & %4.1f & %4.1f \\\\
\\end{tabular}
\\end{center}

"""%(data['*d99_full_masked'],data['*d99_full_unmasked'],
     data['*dmodel_masked'],data['*dmodel_unmasked'],
     data['d_model_b0_masked'],data['d_model_b0_unmasked'],
     data['*dFSCmodel_0_masked'],data['*dFSCmodel_0_unmasked'],
     data['*dFSCmodel_0.143_masked'],data['*dFSCmodel_0.143_unmasked'],
     data['*dFSCmodel_0.5_masked'],data['*dFSCmodel_0.5_unmasked'])

    msg += \
"""
\\underline{Overall isotropic B factor}:\\\\
\\begin{center}
\\begin{tabular}{rcc}
    \\textbf{B factor} & \\textbf{Masked} & \\textbf{Unmasked} \\\\
    Overall B-iso & %4.1f & %4.1f \\\\
\\end{tabular}
\\end{center}

""" % (data["overall_b_iso_masked"], data["overall_b_iso_unmasked"])

    if 'FSC_Model_Map_Masked' in data and 'FSC_Model_Map_Unmasked' in data:
        fnFSCModel = os.path.join(report.getReportDir(),"fscModel.png")
        reportMultiplePlots(data['d_inv_Model_Map_Masked'],
                            [data['FSC_Model_Map_Masked'], data['FSC_Model_Map_Unmasked'],
                             0.5*np.ones(len(data['FSC_Model_Map_Masked']))],
                            "Resolution (A)", "FSC", fnFSCModel,
                            ['Masked','Unmasked','0.5 Threshold'], invertXLabels=True)
        msg+=\
"""Fig. \\ref{fig:fscModel} shows the FSC between the input map and the model.

\\begin{figure}[H]
    \centering
    \includegraphics[width=12cm]{%s}
    \\caption{FSC between the input map and model with and without a mask constructed from the model.
              The X-axis is the square of the inverse of the resolution in \\AA.}
    \\label{fig:fscModel}
\\end{figure}

"""%fnFSCModel
    report.write(msg)

    saveIntermediateData(report.fnReportDir, "phenix", True, "fscModel.png", fnFSCModel, 'Plot that shows FSC between the input map and model with and without a mask constructed from the model')

    warnings = []
    testWarnings = False
    if badResidues > 10 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of residues that have a cross-correlation below 0.5 " \
                        "is %4.1f, that is larger than 10\\%%}}" % badResidues)
    if resolution<0.8*data['*dFSCmodel_0.5_masked'] or testWarnings:
        warnings.append("{\\color{red} \\textbf{The resolution reported by the user, %4.1f \\AA, is significantly " \
                        "smaller than the resolution estimated between map and model (FSC=0.5), %4.1f \\AA}}" %\
                        (resolution,data['*dFSCmodel_0.5_masked']))
    report.addResolutionEstimate(data['*d99_full_masked'])
    report.addResolutionEstimate(data['*dmodel_masked'])
    report.addResolutionEstimate(data['d_model_b0_masked'])
    report.addResolutionEstimate(data['*dFSCmodel_0_masked'])
    report.addResolutionEstimate(data['*dFSCmodel_0.143_masked'])
    report.addResolutionEstimate(data['*dFSCmodel_0.5_masked'])

    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if 1) the percentage of residues whose correlation is
smaller than 0.5 is smaller than 10\\%, and 2) the resolution reported by the user is larger than 0.8 times the
resolution estimated between the map and model at FSC=0.5.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "A.a Phenix validation", secLabel)
    if len(warnings)>0:
        report.writeAbstract("According to phenix, it seems that there might be some mismatch between the map "\
                             "and its model (see Sec. \\ref{%s}). "%secLabel)

def phenix(project, report, protImportMap, protAtom, resolution, priority=False):
    label = "A.a Phenix"
    protPhenix, dataPhenix = phenixExecution(project, report, protImportMap, protAtom, resolution, label, priority)
    phenixReporting(project, report, resolution, protPhenix, dataPhenix)

def searchFitted(list_origins, original_cc_mask, project, report, protImportMap, FNMODEL, resolution, cc_mask_threshold, priority):
    h = AtomicStructHandler()
    h.read(FNMODEL)
    new_cc_mask_dict = {}
    for origin in list_origins:

        origin_dict = {}

        originStr = '-'.join(map(str, origin))
        modelName = os.path.splitext(os.path.basename(FNMODEL))[0]
        modelNewName = f"{modelName}_{originStr}.cif"
        newFnModel = os.path.join(report.getReportDir(), modelNewName)

        moveOriginTo(origin, h)
        h.writeAsCif(newFnModel)
        newProtAtom = importModel(project, report, f"Import atomic - origin {originStr}", protImportMap, newFnModel, priority=priority)
        protPhenix, dataPhenix = phenixExecution(project, report, protImportMap, newProtAtom, resolution, f"Phenix search fitted - origin {originStr}", priority)

        importProtId = newProtAtom.getObjId()
        phenixProtId = protPhenix.getObjId()

        origin_dict["cc_mask"] = dataPhenix["cc_mask"]
        origin_dict["import_prot_id"] = importProtId
        origin_dict["phenix_prot_id"] = phenixProtId

        new_cc_mask_dict[originStr] = origin_dict

    # Get greatest cc_mask that is also greater than cc_mask original and beyond cc_mask_threshold (that normally is 0.5)
    greatest_cc_mask = None
    greatest_cc_mask_key = None
    
    for key, sub_dict in new_cc_mask_dict.items():
        new_cc_mask = sub_dict.get("cc_mask", None)
        if new_cc_mask is not None and new_cc_mask > original_cc_mask and new_cc_mask > cc_mask_threshold:
            if greatest_cc_mask is None or new_cc_mask > greatest_cc_mask:
                greatest_cc_mask = new_cc_mask
                greatest_cc_mask_key = key

    if greatest_cc_mask:
        newOrigin = greatest_cc_mask_key.split("-")

        newImportProtId = new_cc_mask_dict[greatest_cc_mask_key]["import_prot_id"]
        newImportProt = project.getProtocol(int(newImportProtId), fromRuns=True)

        newPhenixProtId = new_cc_mask_dict[greatest_cc_mask_key]["phenix_prot_id"]
        newPhenixProt = project.getProtocol(int(newPhenixProtId), fromRuns=True)

        return newImportProt, newPhenixProt, dataPhenix, newOrigin
    else:
        return None, None, None, None

def checkFitted(project, report, EMDB_ID_NUM, section, secLabel, protImportMap, protAtom, FNMODEL, resolution, mapCoordX, mapCoordY, mapCoordZ, cc_mask_threshold, priority=False):
    label = "A.a Phenix"
    protPhenix, dataPhenix = phenixExecution(project, report, protImportMap, protAtom, resolution, label, priority)

    if protPhenix.isFailed():
        report.write(ERROR_MESSAGE_CHECK_FITTED_FAILED)
        return None, protPhenix, dataPhenix, protAtom

    if protPhenix.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_PHENIX)
        report.writeSummary(section, secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return None, protPhenix, dataPhenix, protAtom

    if dataPhenix["cc_mask"] > cc_mask_threshold:
        report.write(PROPERLY_FITTED)
        return True, protPhenix, dataPhenix, protAtom
       
    elif dataPhenix["cc_mask"] <= cc_mask_threshold:

        # Create list of origin of coordinate to test
        x, y, z = protImportMap.outputVolume.getDimensions()        

        list_origins = [[0,0,0], [x/2, y/2, z/2], [mapCoordX, mapCoordY, mapCoordZ]] #TODO: finish list of origin of coordinates to test (remaining: from header)
        
        protAtom, protPhenix, dataPhenix, newOrigin = searchFitted(list_origins, dataPhenix["cc_mask"], project, report, protImportMap, FNMODEL, resolution, cc_mask_threshold, priority)

        if protAtom is None and protPhenix is None:
            report.write(NOT_FOUND_NEW_FITTED)

            ######TODO: this part is for debugging. Remove when finishing debugging.######
            os.makedirs(EMDB_entries_path, exist_ok=True)
            with open(os.path.join(EMDB_entries_path, 'EMDB_fail_fitted.txt'), 'a') as output:
                if EMDB_ID_NUM:
                    output.writelines(f'{EMDB_ID_NUM}\n')
                else:
                    output.writelines('Unkown\n')
            ##############################################################################

            return False, protPhenix, dataPhenix, protAtom
        else:
            report.write(FITTED_NEW_ORIGIN % ', '.join(newOrigin))

            ######TODO: this part is for debugging. Remove when finishing debugging.######
            os.makedirs(EMDB_entries_path, exist_ok=True)
            with open(os.path.join(EMDB_entries_path, 'EMDB_new_fitted.txt'), 'a') as output:
                if EMDB_ID_NUM:
                    output.writelines(f'{EMDB_ID_NUM}\n')
                else:
                    output.writelines('Unkown\n')
            ##############################################################################

            return True, protPhenix, dataPhenix, protAtom


def convertPDB(project, report, protImportMap, protAtom, priority=False):

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtConvertPdb', doRaise=True)
    protConvert = project.newProtocol(Prot,
                                      objLabel="Convert Pdb to map",
                                      inputPdbData=1,
                                      sampling=protImportMap.outputVolume.getSamplingRate(),
                                      vol=True,
                                      centerPdb=False)
    protConvert.pdbObj.set(protAtom.outputPdb)
    protConvert.volObj.set(protImportMap.outputVolume)
    if useSlurm:
        sendToSlurm(protConvert, priority=True if priority else False)
    project.launchProtocol(protConvert)
    #waitOutput(project, protConvert, 'outputVolume')
    waitUntilFinishes(project, protConvert)
    secLabel = "sec:convertPdb2Map"
    msg = \
    """
    \\subsection{Level A.b Conversion PDB to map}
    \\label{%s}
    \\textbf{Explanation}:\\\\ 
    Convert a PDB file into a volume.\\\\
    \\\\
    \\textbf{Results:}\\\\
    \\\\
    """ % secLabel
    if protConvert.isFailed():     #TODO: check texts for ConverToPdb when fails
        report.writeSummary("A. Conversion PDB to map", secLabel, ERROR_CONVERT_MESSAGE)
        report.write(msg)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_CONVERT_MESSAGE)
        return None

    if protConvert.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_CONVERSION_PDB_MAP)
        report.writeSummary("A. Conversion PDB to map", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return protConvert

    # Check if volume is not empty
    volumeData = xmipp3.Image(protConvert.outputVolume.getFileName()).getData()
    if not np.sum(volumeData) > 0:
        report.writeSummary("A. Conversion PDB to map", secLabel, ERROR_CONVERT_MESSAGE)
        report.write(msg)
        report.write(ERROR_MESSAGE_EMPTY_VOL + STATUS_ERROR_MESSAGE)
        return None
    return protConvert

def fscq(project, report, protImportMap, protAtom, protConvert, protCreateSoftMask, fnMaskedMap, priority=False):
    if not protImportMap.outputVolume.hasHalfMaps():
        return

    secLabel = "sec:fscq"
    msg = \
"""
\\subsection{Level A.b FSC-Q}
\\label{%s}
\\textbf{Explanation}:\\\\ 
FSC-Q (see this \\href{%s}{link} for more details) compares the local FSC between the map and the atomic model to the local FSC of the two 
half maps. FSC-Qr is the normalized version of FSC-Q to facilitate comparisons. Typically, FSC-Qr should
take values between -1.5 and 1.5, being 0 an indicator of good matching between map and model.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, FSCQ_DOI)
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtValFit', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="A.b FSC-Q",
                               inputPDBObj=protAtom.outputPdb,
                               numberOfThreads=10)
    prot.inputVolume.set(protImportMap.outputVolume)
    prot.pdbMap.set(protConvert.outputVolume)
    prot.inputMask.set(protCreateSoftMask.outputMask)
    if useSlurm:
        sendToSlurm(prot, priority=True if priority else False)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputAtomStruct')
    waitUntilFinishes(project, prot)

    if prot.isFailed():
        report.writeSummary("A.b FSC-Q", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_FSCQ)
        report.writeSummary("A.b FSC-Q", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    V = xmipp3.Image(prot._getExtraPath("pdb_volume.map"))
    FSCQr = xmipp3.Image(prot._getExtraPath("diferencia_norm.map"))
    fscqr = FSCQr.getData()[V.getData()>0.97]
    avgFSCQr = np.mean(fscqr)
    ci = np.percentile(fscqr, [2.5, 97.5])
    f15=float(np.sum(abs(fscqr)>1.5))/fscqr.size*100

    fnHist = os.path.join(report.getReportDir(),"fscqrHist.png")
    reportHistogram(np.clip(fscqr, -1.5, 1.5), "FSC-Qr", fnHist)
    Bpercentiles = np.percentile(np.clip(fscqr, -1.5, 1.5), np.array([0.025, 0.25, 0.5, 0.75, 0.975])*100)


    msg =\
"""Fig. \\ref{fig:fscqHist} shows the histogram of FSC-Qr and Fig. 
\\ref{fig:fscq} the colored isosurface of the atomic model converted to map. The
average FSC-Qr is %5.2f, its 95\\%% confidence interval is [%5.2f,%5.2f]. The percentage of values
whose FSC-Qr absolute value is beyond 1.5 is %5.1f \\%%.

\\begin{figure}[H]
  \\centering
  \\includegraphics[width=8cm]{%s}
  \\caption{Histogram of the FSC-Qr limited to -1.5 and 1.5.}
  \\label{fig:fscqHist}
\\end{figure}

"""%(avgFSCQr, ci[0], ci[1], f15, fnHist)
    report.colorIsoSurfaces(msg, "Isosurface of the atomic model colored by FSC-Qr between -1.5 and 1.5",
                            "fig:fscq", project, "fscq", fnMaskedMap,
                            protImportMap.outputVolume.getSamplingRate(),
                            prot._getExtraPath("diferencia_norm.map"), -3.0, 3.0)

    saveIntermediateData(report.getReportDir(), 'FSCQ', False, 'averageFSQr', float(avgFSCQr), ['', 'average FSC-Qr'])
    saveIntermediateData(report.getReportDir(), 'FSCQ', False, 'llcinterval', float(ci[0]), ['', 'Lower limit 95% confidence interval'])
    saveIntermediateData(report.getReportDir(), 'FSCQ', False, 'ulcinterval', float(ci[1]), ['', 'Upper limit 95% confidence interval'])
    saveIntermediateData(report.getReportDir(), 'FSCQ', False, 'percentage15', f15, ['%', 'The percentage of values whose FSC-Qr absolute value is beyond 1.5'])
    saveIntermediateData(report.getReportDir(), 'FSCQ', False, 'fscqrList', np.clip(fscqr, -1.5, 1.5).tolist(), ['', 'List of FSCalues to create the histogram'])

    saveIntermediateData(report.getReportDir(), 'FSCQ', True, 'fscq_struct', os.path.join(project.getPath(), prot._getPath('fscq_struct.cif')), 'fscq_struct cif file')
    saveIntermediateData(report.getReportDir(), 'FSCQ', True, 'fscqHist', fnHist, 'FSCQ Histogram')
    saveIntermediateData(report.getReportDir(), 'FSCQ', True, 'fscqViewer',
                         [os.path.join(report.getReportDir(), 'fscq1.jpg'),
                          os.path.join(report.getReportDir(), 'fscq2.jpg'),
                          os.path.join(report.getReportDir(), 'fscq3.jpg')], 'FSCQ views')

    warnings = []
    testWarnings = False

    BperHomogeneous = isHomogeneous(Bpercentiles[0], Bpercentiles[-1], eps=0.1)

    if BperHomogeneous:
        warnings.append("{\\color{red} \\textbf{Program output seems to be too homogeneous. There might " \
                        "be some program issues analyzing the data.}}")

    if f15>10 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The percentage of voxels that have a FSC-Qr larger than 1.5 in "\
                        "absolute value is %5.1f, that is larger than 10\\%%}}"%f15)
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the percentage of residues whose FSC-Q is larger than 1.5
in absolute value is smaller than 10\\%.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "A.b FSC-Q", secLabel)
    if len(warnings)>0:
        report.writeAbstract("According to FSC-Q, it seems that there is a mismatch between the map and its model "\
                             "(see Sec. \\ref{%s}). "%secLabel)

def multimodel(project, report, protImportMap, protAtom, resolution, priority=False):

    secLabel = "sec:multimodel"
    msg = \
"""
\\subsection{Level A.c Multimodel stability}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method (see this \\href{%s}{link} for more details) estimates the ambiguity of the atomic model in each region of the CryoEM map due to 
the different local resolutions or local heterogeneity.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, MULTIMODEL_DOI)
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('rosetta.protocols',
                                            'ProtRosettaGenerateStructures', doRaise=True)
    prot1  = project.newProtocol(Prot,
                                 objLabel="A.c Multimodel ambiguity",
                                 inputStructure=protAtom.outputPdb,
                                 inputVolume=protImportMap.outputVolume,
                                 resolution=resolution,
                                 numMods=2)
    if useSlurm:
        sendToSlurm(prot1, GPU=True, priority=True if priority else False)
    project.launchProtocol(prot1)
    #waitOutput(project, prot1, 'outputAtomStructs')
    waitUntilFinishes(project, prot1)


    if prot1.isFailed() or not hasattr(prot1,"outputAtomStructs"):
        report.writeSummary("A.c Multimodel", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return

    if prot1.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_MULTIMODEL)
        report.writeSummary("A.c Multimodel", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot1

    Prot = pwplugin.Domain.importFromPlugin('atomstructutils.protocols',
                                            'ProtRMSDAtomStructs', doRaise=True)
    prot2 = project.newProtocol(Prot,
                                objLabel="A.c RMSD",
                                inputStructureSet=prot1.outputAtomStructs)
    if useSlurm:
        sendToSlurm(prot2, priority=True if priority else False)
    project.launchProtocol(prot2)
    #waitOutput(project, prot2, 'outputAtomStructs')
    waitUntilFinishes(project, prot2)
    if prot2.isFailed():
        report.writeSummary("A.c Multimodel", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return

    if prot2.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_MULTIMODEL)
        report.writeSummary("A.c Multimodel", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot2

    fnCifs = glob.glob(prot2._getPath('*.cif'))
    if len(fnCifs)==0:
        report.writeSummary("A.c Multimodel", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return

    fnCif = fnCifs[0]

    fnCifRMSD = os.path.join(report.getReportDir(),"atomicModelRMSD.cif")
    replaceOcuppancyWithAttribute(fnCif, "perResidueRMSD", fnCifRMSD)

    cifDic = AtomicStructHandler().readLowLevel(fnCifRMSD)
    rmsd = []
    for name, value in zip(cifDic['_scipion_attributes.name'],cifDic['_scipion_attributes.value']):
        if name=='perResidueRMSD':
            rmsd.append(float(value))
    fnRMSDHist = os.path.join(report.getReportDir(),"rmsdHist.png")
    reportHistogram(rmsd, "RMSD", fnRMSDHist)

    avgRMSD = np.mean(rmsd)
    msg =\
"""Fig. \\ref{fig:rmsdHist} shows the histogram of the RMSD of the different models. The average RMSD between models
is %4.2f \\AA. Fig. \\ref{fig:modelRMSD} shows the atomic model colored by RMSD.

\\begin{figure}[H]
  \\centering
  \\includegraphics[width=8cm]{%s}
  \\caption{Histogram of RMSD of the different atoms of the multiple models.}
  \\label{fig:rmsdHist}
\\end{figure}

"""%(avgRMSD, fnRMSDHist)

    report.atomicModel("modelRMSD", msg, "Atomic model colored by RMSD", fnCifRMSD, "fig:modelRMSD", occupancy=True)

    warnings=[]
    testWarnings = False
    if avgRMSD>2 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The average RMSD is too high, "\
                        "%4.1f\\%%}}"%avgRMSD)
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the average RMSD is smaller than 2\\AA.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "A.c Multimodel", secLabel)

    if len(warnings)>0:
        report.writeAbstract("It seems that the model is too ambiguous (see Sec. \\ref{%s}). "%\
                             secLabel)

def guinierModel(project, report, protImportMap, protConvert, resolution, priority=False):
    map = protImportMap.outputVolume
    Ts = map.getSamplingRate()

    fnAtom = protConvert.outputVolume.getFileName()
    fnOut = os.path.join(report.getReportDir(), "sharpenedModel.mrc")
    args = "-i %s -o %s --sampling %f --maxres %s --auto"%(fnAtom, fnOut, Ts, resolution)

    scipionHome = getScipionHome()
    scipion3 = os.path.join(scipionHome, 'scipion3')
    cmd = '%s run xmipp_volume_correct_bfactor %s' % (scipion3, args)

    if not useSlurm:
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    else:
        randomInt = int(datetime.now().timestamp()) + randint(0, 1000000)
        slurmScriptPath = createScriptForSlurm('xmipp_volume_correct_bfactor_levelA_' + str(randomInt), report.getReportDir(), cmd, priority=priority)
        # send job to queue
        subprocess.Popen('sbatch %s' % slurmScriptPath, shell=True)
        # check if job has finished
        while True:
            if checkIfJobFinished('xmipp_volume_correct_bfactor_levelA_' + str(randomInt)):
                break
        sleep(120)

    dinv2, lnFMap, _ = readGuinier(os.path.join(report.getReportDir() if not useSlurm else os.path.dirname(slurmScriptPath), 'sharpenedMap.mrc') + '.guinier')
    _, lnFAtom, _ = readGuinier(fnOut + ".guinier")
    lnFMapp = lnFMap+(np.mean(lnFAtom)-np.mean(lnFMap))

    R = np.corrcoef(lnFMapp, lnFAtom)
    R=R[0,1]

    fnPlot = os.path.join(report.getReportDir(), 'BfactorAtom.png')
    reportMultiplePlots(dinv2, [lnFAtom, lnFMapp], '1/Resolution^2 (1/A^2)', 'log Structure factor', fnPlot,
                        ['Atomic model', 'Experimental map'])

    secLabel = "sec:bfactorModel"
    msg = \
"""\\subsection{Level A.d Map-Model Guinier analysis}
\\label{%s}
\\textbf{Explanation:}\\\\
We compared the Guinier plot (see this \\href{%s}{link} for more details) of the atomic model and the experimental map. We made the mean
of both profiles to be equal (and equal to the mean of the atomic model) to make sure that they had comparable scales. 
\\\\
\\\\
\\textbf{Results:}\\\\
Fig. \\ref{fig:BfactorModel} shows the logarithm (in natural units) of the structure factor (the module squared of the
Fourier transform) of the atom model and the experimental map. The correlation between the two profiles was %5.3f.

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Guinier plot of the atom model and experimental map. 
              The X-axis is the square of the inverse of the resolution in \\AA.}
    \\label{fig:BfactorModel}
\\end{figure}

""" % (secLabel, BFACTOR_GUINIER_DOI, R, fnPlot)
    report.write(msg)

    warnings = []
    testWarnings = False
    if R<0.5 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The correlation is smaller than 0.5, it is %5.3f.}}"%R)
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the correlation between the two Guinier profiles is larger
than 0.5.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "A.d Map-Model Guinier", secLabel)
    if len(warnings)>0:
        report.writeAbstract("It seems that the Guinier plot of the map and its model do not match "\
                             "(see Sec. \\ref{%s}). "%secLabel)

    saveIntermediateData(report.getReportDir(), 'guinierModel', False, 'correlation', R, ['', 'The correlation between the structure factor of the atom model and the experimental map'])

    saveIntermediateData(report.getReportDir(), 'guinierModel', True, 'sharpenedModel.mrc.guinier', os.path.join(report.getReportDir(), 'sharpenedModel.mrc.guinier'), 'sharpenedModel.mrc.guinier file which contain the data to create the guinier plot')
    saveIntermediateData(report.getReportDir(), 'guinierModel', True, 'guinierPlot', fnPlot, 'guinier plot for Map-Model Guinier Analysis')

def mapq(project, report, protImportMap, protAtom, resolution, priority=False):

    secLabel = "sec:mapq"
    msg = \
"""
\\subsection{Level A.e MapQ}
\\label{%s}
\\textbf{Explanation}:\\\\ 
MapQ (see this \\href{%s}{link} for more details) computes the local correlation between the map and each one of its atoms assumed to
have a Gaussian shape.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, MAPQ_DOI)
    report.write(msg)

    # check if we have the precomputed data
    # https://3dbionotes.cnb.csic.es/bws/api/emv/7xzz/mapq/
    emdb_Id = getFilename(str(protImportMap.filesPath), withExt=False)
    pdbdb_Id = getFilename(str(protAtom.outputPdb._filename), withExt=False)
    print("Get MapQ scores from 3DBionotes-WS for %s" % pdbdb_Id)
    has_precalculated_data = False
    cif_data = getFileFromWS(pdbdb_Id, 'mapq')

    if cif_data:
        # save to report
        results_msg = \
            """
            \\Precalculated MapQ scores obtained from DB source via 3DBionotes-WS:
            \\\\
            \\url{https://3dbionotes.cnb.csic.es/bws/api/emv/%s/mapq/}
            \\\\
            """ % emdb_Id.lower().replace('_','-')
        report.write(results_msg)
        has_precalculated_data = True
    else:
        # if there is not precalculated data or failed to retrieve it
        print('- Could not get data for', pdbdb_Id)
        print('-- Proceed to calculate it localy')
        Prot = pwplugin.Domain.importFromPlugin('mapq.protocols',
                                                'ProtMapQ', doRaise=True)
        prot = project.newProtocol(Prot,
                                objLabel="A.e MapQ",
                                inputVol=protImportMap.outputVolume,
                                pdbs=[protAtom.outputPdb],
                                mapRes=resolution)
        if useSlurm:
            sendToSlurm(prot, priority=True if priority else False)
        project.launchProtocol(prot)
        #waitOutput(project, prot, 'scoredStructures')
        waitUntilFinishes(project, prot)
        if prot.isFailed():
            report.writeSummary("A.e MapQ", secLabel, ERROR_MESSAGE)
            report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
            return

        if prot.isAborted():
            print(PRINT_PROTOCOL_ABORTED + ": " + NAME_MAPQ)
            report.writeSummary("A.e MapQ", secLabel, ERROR_ABORTED_MESSAGE)
            report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
            return prot

        saveIntermediateData(report.getReportDir(), 'MapQ', True, 'cif', glob.glob(os.path.join(project.getPath(), prot._getExtraPath('*.cif')))[0], 'cif file')
        saveIntermediateData(report.getReportDir(), 'MapQ', True, 'Q__map_All', glob.glob(os.path.join(project.getPath(), prot._getExtraPath('*Q__map_All.txt')))[0], 'Q__map_All txt file')
        saveIntermediateData(report.getReportDir(), 'MapQ', True, 'Q__map.pdb', glob.glob(os.path.join(project.getPath(), prot._getExtraPath('*Q__map.pdb')))[0], 'Q__map pdb file')

        input_file = glob.glob(os.path.join(project.getPath(), prot._getExtraPath('*Q__map.pdb')))[0]
        # emd_26162_pdb_7txz_emv_mapq.json
        output_file = os.path.join(project.getPath(), prot._getExtraPath(), "%s_pdb_%s_emv_mapq.json" % (emdb_Id.lower().replace('-','_'), pdbdb_Id.lower()))
        json_file = convert_2_json(emdb_Id, pdbdb_Id, method='mapq', input_file=input_file, output_file=output_file)
        saveIntermediateData(report.getReportDir(), 'MapQ', True, 'EMV json file', json_file, 'MapQ scores in EMV json format')


    # get histogram
    mapq_scores = []
    if has_precalculated_data and cif_data:
        cifWSFilename = os.path.join(project.getPath(), project.getTmpPath(), pdbdb_Id + '_MapQFromWS.cif')
        pdbFromWS = open(cifWSFilename, 'w')
        pdbFromWS.write(cif_data)
        pdbFromWS.close()
        with open(cifWSFilename) as cif:
            lines = cif.readlines()
            for line in lines:
                if 'ATOM' in line:
                    line = re.sub(' +', ' ', line)
                    values = line.split(' ')
                    mapq_scores.append(float(values[15]))
    else:
        ASH = AtomicStructHandler()

        for struct in prot.scoredStructures:
            fileName = struct.getFileName()
            fields = ASH.readLowLevel(fileName)
            attributes = fields["_scipion_attributes.name"]
            values = fields["_scipion_attributes.value"]
            mapq_scores += [float(value) for attribute, value in zip(attributes, values) if attribute == "MapQ_Score"]

    fnHist = os.path.join(report.getReportDir(),"mapqHist.png")

    reportHistogram(mapq_scores, "MapQ score", fnHist)
    Bpercentiles = np.percentile(mapq_scores, np.array([0.025, 0.25, 0.5, 0.75, 0.975])*100)

    toWrite = \
"""
Fig. \\ref{fig:histMapQ} shows the histogram of the calculated Q-score. Some representative
percentiles are:

\\begin{center}
    \\begin{tabular}{|c|c|}
        \\hline
        \\textbf{Percentile} & \\textbf{MapQ score [-1,1]} \\\\
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

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Histogram of the per-atom Q-score.}
    \\label{fig:histMapQ}
\\end{figure}

""" % (Bpercentiles[0], Bpercentiles[1], Bpercentiles[2], Bpercentiles[3], Bpercentiles[4], fnHist)
    report.write(toWrite)

    if not has_precalculated_data:
        files = glob.glob(prot._getExtraPath("*All.txt"))
    else:
        cifFilename = os.path.join(report.getReportDir(), pdbdb_Id + '_MAPQFromWS.cif')
        cifFromWS = open(cifFilename, 'w')
        cifFromWS.write(getFileFromWS(pdbdb_Id, 'mapq'))
        cifFromWS.close()

        QStatsScript =\
"""
import sys
import chimera

from chimera import Molecule

mapq_path = "%s"
sys.path.append(mapq_path)

validation_tools_path = "%s"
sys.path.append(validation_tools_path)

import mmcif
import qscores
from mapq_utils import SaveQStats, ReadMol

print('\\nCreating MapQ Statistics...\\n')
mol = ReadMol("%s")

chimera.openModels.add([mol])

mols = chimera.openModels.list(modelTypes=[Molecule])
if len(mols) == 0:
    print(" - no molecules loaded")
else:
    print('\\nMolecules loaded correctly')
    for mi, mol in enumerate(mols):
        qscores.SetBBAts(mol)
        SaveQStats(mol, "All", 0.6, %d)

"""%(mapq_path, validation_tools_path, cifFilename, resolution)        
        fnQStatsScript = os.path.join(report.getReportDir(),"mapq_stats.py")
        fhQStatsScript = open(fnQStatsScript,"w")
        fhQStatsScript.write(QStatsScript)
        fhQStatsScript.close()

        args = "--nogui --script %s "%(fnQStatsScript)
        print("Running: %s %s" % (chimeraProgram, args))
        p = subprocess.Popen('%s %s' % (chimeraProgram, args), shell=True, stderr=subprocess.PIPE)
        p.wait()

        files = glob.glob(os.path.join(report.getReportDir(), "*All.txt"))

    fh = open(files[0])
    msg=\
""" The following table shows the average Q-score and estimated resolution for each chain.
\\begin{center}
    \\begin{tabular}{ccc}
        \\hline
        \\textbf{Chain} & \\textbf{Average Q-score [-1,-1]} & \\textbf{Estimated Resol. (\\AA)} \\\\
        \\hline
"""

    state = 0
    resolutions = []
    for line in fh.readlines():
        if state==0 and line.startswith('Chain'):
            state=1
        elif state==1:
            tokens = line.split()
            if len(tokens)>0:
                res = float(tokens[-1])
                msg+="      %s & %s & %4.1f \\\\ \n"%(tokens[0],tokens[3],res)
                resolutions.append(res)
            else:
                state=2
                break
    msg+=\
    """        \\hline
        \\end{tabular}
    \\end{center}

    """
    report.write(msg)
    fh.close()

    report.addResolutionEstimate(np.mean(resolutions))

    saveIntermediateData(report.getReportDir(), 'MapQ', False, 'estimatedResolution', np.mean(resolutions), ['\u212B', 'The estimated resolution (mean) in Angstroms obtained from MapQ'])

    # get colored models
    msg = "The atomic model colored by MapQ can be seen in Fig. \\ref{fig:mapq}.\n\n"
    if not has_precalculated_data:
        fnCifMapQ = os.path.join(project.getPath(), prot._getExtraPath("chimeraAttribute_MapQ_score.cif"))
        # make sure the output mapq file is correct
        with open(os.path.join(project.getPath(), prot._getExtraPath("%s.cif" % pdbdb_Id))) as cif:
            cifData = cif.read()
        cifData = cifData.replace('residues', 'atoms')
        with open(os.path.join(project.getPath(), prot._getExtraPath("%s.cif" % pdbdb_Id)), 'w') as cif:
            cif.write(cifData)

        replaceOcuppancyWithAttribute(os.path.join(project.getPath(), prot._getExtraPath("%s.cif" % pdbdb_Id)), "MapQ_Score",
                                      fnCifMapQ)
        report.atomicModel("mapqView", msg, "Atomic model colored by MapQ", fnCifMapQ, "fig:mapq", bfactor=False,
                           occupancy=True, rainbow=False, legendMin=-1, legendMax=1)

    else:
        # create .defattr file
        qscores = {}
        with open(cifWSFilename) as cif:
            lines = cif.readlines()
            for line in lines:
                if 'ATOM' in line:
                    line = re.sub(' +', ' ', line)
                    values = line.split(' ')
                    qscores[values[1]] = values[15]

        attributeFile = os.path.join(project.getPath(), project.getTmpPath(), pdbdb_Id + '_MapQFromWS.defattr')
        with open(attributeFile, 'a') as af:
            af.write('attribute: qscores\nrecipient: atoms\n')
            for atom in qscores:
                af.write('\t:%s\t%s\n' % (atom, qscores[atom]))

        report.atomicModel("mapqView", msg, "Atomic model colored by MapQ", cifWSFilename, "fig:mapq", bfactor=False,
                           occupancy=False, otherAttribute=[attributeFile, 'qscores'], rainbow=False, legendMin=-1, legendMax=1)

    saveIntermediateData(report.getReportDir(), 'MapQ', True, 'MapQView',
                         [os.path.join(report.getReportDir(), 'mapqView1.jpg'),
                          os.path.join(report.getReportDir(), 'mapqView1.jpg'),
                          os.path.join(report.getReportDir(), 'mapqView1.jpg')], 'MapQ views')

    # Warnings
    warnings=[]
    testWarnings = False

    BperHomogeneous = isHomogeneous(Bpercentiles[0], Bpercentiles[-1], eps=0.1)

    if BperHomogeneous:
        warnings.append("{\\color{red} \\textbf{Program output seems to be too homogeneous. There might " \
                        "be some program issues analyzing the data.}}")

    if Bpercentiles[2]<0.1 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The median Q-score is less than 0.1.}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the median Q-score is larger than 0.1.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "A.e MapQ", secLabel)
    if len(warnings)>0:
        report.writeAbstract("There seems to be a problem with its MapQ scores (see Sec. \\ref{%s}). "%secLabel)

def emringer(project, report, protImportMap, protAtom, priority=False):

    secLabel = "sec:emringer"
    msg = \
"""
\\subsection{Level A.f EMRinger validation}
\\label{%s}
\\textbf{Explanation}:\\\\ 
EMringer (see this \\href{%s}{link} for more details) compares the side chains of the atomic model to the CryoEM map. The following features are
reported:
\\begin{itemize}
    \\item Optimal Threshold: Electron potential map cutoff value at which the maximum EMRinger score was obtained.
    \\item Rotamer Ratio: Fraction of rotameric residues at the Optimal threshold value.
    \\item Max Zscore: Z-score computed to determine the significance of the distribution at the Optimal threshold value.
    \\item Model Length: Total of non-gamma-branched, non-proline aminoacids with a non-H gamma atom used in global EMRinger score computation.
    \\item EMRinger Score: Maximum EMRinger score calculated at the Optimal Threshold.
\\end{itemize}
A rotameric residue is one in which EMRinger peaks that fall within defined rotamers based on chi1, this often suggests 
a problem with the modelling of the backbone. In general, the user should look at the profiles and identify regions 
that may need improvement. 
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, EMRINGER_DOI)
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('phenix.protocols',
                                            'PhenixProtRunEMRinger', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="A.f EMRinger")
    prot.inputVolume.set(protImportMap.outputVolume)
    prot.inputStructure.set(protAtom.outputPdb)
    if useSlurm:
        sendToSlurm(prot, priority=True if priority else False)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'stringDataDict')
    #waitOutputFile(project, prot, '*_emringer_plots')
    waitUntilFinishes(project, prot)

    if prot.isFailed():
        report.writeSummary("A.f EMRinger", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED)
        emringerStdout = open(os.path.join(project.getPath(), prot.getStdoutLog()), "r").read()
        controlledErrors = ["Sorry: No residues could be scanned by EMRinger, so scores cannot be generated"]
        for error in controlledErrors:
            if error in emringerStdout:
                report.write("{\\color{red} \\textbf{REASON: %s.}}\\\\ \n" % error)
        report.write(STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_EMRINGER)
        report.writeSummary("A.f EMRinger", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    dataDict=json.loads(str(prot.stringDataDict), object_pairs_hook=collections.OrderedDict)

    maxScoreIndex = dataDict['_maxScoreIndex']
    optimalThreshold = str("%0.3f" % dataDict['_thresholds'][maxScoreIndex])

    fnScore = os.path.join(report.getReportDir(),"emringerThreshold_scan.png")
    copyFile(os.path.join(project.getPath(),
                          glob.glob(prot._getExtraPath("*_emringer_plots/Total.threshold_scan.png"))[0]),
             fnScore)
    fnResidueHist = os.path.join(report.getReportDir(),"residueHist.png")
    copyFile(os.path.join(project.getPath(),
                          glob.glob(prot._getExtraPath("*_emringer_plots/%s.histogram.png"%optimalThreshold))[0]),
             fnResidueHist)

    msg=\
"""\\underline{General results}:\\\\
\\begin{center}
\\begin{tabular}{rc}
    Optimal threshold & %f \\\\
    Rotamer ratio & %4.3f \\\\
    Max. Zscore & %5.2f \\\\
    Model length & %d \\\\
    EMRinger Score & %4.3f \\\\
\\end{tabular}
\\end{center}

Fig. \\ref{fig:emringerThreshold} shows the EMRinger score and fraction of rotameric residues as a function of
the map threshold. The optimal threshold was selected looking for the maximum EMRinger score in this plot.

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{EMRinger score and fraction of rotameric residues as a function of the map threshold.}
    \\label{fig:emringerThreshold}
\\end{figure}

Fig. \\ref{fig:emRingerPeaks} shows the histogram for rotameric (blue) and non-rotameric (red) residues at the 
optimal threshold.

\\begin{figure}[H]
    \centering
    \includegraphics[width=12cm]{%s}
    \\caption{Histogram for rotameric (blue) and non-rotameric (red) residues at the optimal threshold as a
    function of the angle Chi1.}
    \\label{fig:emRingerPeaks}
\\end{figure}

"""%(dataDict["Optimal Threshold"], dataDict["Rotamer-Ratio"], dataDict["Max Zscore"], dataDict["Model Length"],
     dataDict["EMRinger Score"], fnScore, fnResidueHist)
    
    saveIntermediateData(report.getReportDir(), 'EMRinger', False, 'dataDict', dataDict, ['', 'emringer data dictionary containing all key params'])

    saveIntermediateData(report.getReportDir(), 'EMRinger', True, 'emringerThreshold_scan.png', fnScore, 'emringer threshold scan plot showing the EMRinger score and fraction of rotameric residues as a function of the map threshold')
    saveIntermediateData(report.getReportDir(), 'EMRinger', True, 'residueHist.pngv', fnResidueHist, 'emringer histogram for rotameric (blue) and non-rotameric (red) residues at the optimal threshold')

    msg+=\
"""The following plots show the rolling window EMRinger analysis of the different chains to distinguish regions 
of improved model quality. This analysis was performed on rolling sliding 21-residue windows along the primary 
sequence of the protein chains.

"""
    for chain in sorted(dataDict['_chains']):
        fnPlot = os.path.join(project.getPath(),
                              glob.glob(prot._getExtraPath("*_emringer_plots/%s_rolling.png"%chain))[0])
        msg += "\\includegraphics[width=7cm]{%s}\n" % fnPlot
    msg+="\n"

    report.write(msg)

    warnings = []
    testWarnings = False
    if dataDict["EMRinger Score"] <1 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The EMRinger score is smaller than 1, it is %4.3f.}}"%\
                        dataDict["EMRinger Score"])
    if dataDict["Max Zscore"] <1 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The maximum Zscore is smaller than 1, it is %4.3f.}}"%\
                        dataDict["Max Zscore"])
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the EMRinger score and Max. Zscore are larger than 1.
\\\\

"""
    report.write(msg)

    report.writeWarningsAndSummary(warnings, "A.f EMRinger", secLabel)
    if len(warnings)>0:
        report.writeAbstract("The EMRinger score is negative, it seems that the model side chains do not match the "\
                            "map (see Sec. \\ref{%s}). "%secLabel)

    _emringer_plots = []
    for file in os.listdir(glob.glob(prot._getExtraPath('*_emringer_plots'))[0]):
        _emringer_plots.append(os.path.join(project.getPath(), prot._getExtraPath('*_emringer_plots'), file))

    saveIntermediateData(report.getReportDir(), 'EMRinger', True, 'emringer_csv', glob.glob(os.path.join(project.getPath(), prot._getExtraPath('*_emringer.csv')))[0], 'emringer_csv file')
    saveIntermediateData(report.getReportDir(), 'EMRinger', True, '7tmw_emringer.pkl', glob.glob(os.path.join(project.getPath(), prot._getExtraPath('*_emringer.pkl')))[0], 'emringer pickle file')
    saveIntermediateData(report.getReportDir(), 'EMRinger', True, 'emringer.map', glob.glob(os.path.join(project.getPath(), prot._getExtraPath('emringer.map')))[0], 'emringer.map file')

    saveIntermediateData(report.getReportDir(), 'EMRinger', True, '_emringer_plots', _emringer_plots, '_emringer_plots files')

def daq(project, report, protImportMap, protAtom, resolution, priority=False):

    secLabel = "sec:daq"
    msg = \
"""
\\subsection{Level A.g DAQ validation}
\\label{%s}
\\textbf{Explanation}:\\\\ 
DAQ (see this \\href{%s}{link} for more details) is a computational tool using deep learning that can estimate the residue-wise local
quality for protein models from cryo-Electron Microscopy maps. The method calculates the likelihood that a given
density feature corresponds to an aminoacid, atom, and secondary structure. These likelihoods are combined into a
 score that ranges from -1 (bad quality) to 1 (good quality). \\\\
\\\\
\\textbf{Results:}\\\\
""" % (secLabel, DAQ_DOI)
    report.write(msg)

    #TODO: API call
    # check if we have the precomputed data
    # https://3dbionotes.cnb.csic.es/bws/api/emv/7xzz/daq/
    emdb_Id = getFilename(str(protImportMap.filesPath), withExt=False)
    pdbdb_Id = getFilename(str(protAtom.outputPdb._filename), withExt=False)
    print("Get DAQ scores from 3DBionotes-WS for %s" % pdbdb_Id)
    has_precalculated_data = False
    json_data = getScoresFromWS(pdbdb_Id, 'daq')

    if json_data:
        # save to report
        results_msg = \
            """
            \\Precalculated DAQ scores obtained from DB source via 3DBionotes-WS:
            \\\\
            \\url{https://3dbionotes.cnb.csic.es/bws/api/emv/%s/daq/}
            \\\\
            """ % emdb_Id.lower().replace('_','-')
        report.write(results_msg)
        has_precalculated_data = True
    else:
        # if there is not precalculated data or failed to retrieve it
        print('- Could not get data for', pdbdb_Id)
        print('-- Proceed to calculate it localy')

        if resolution>5:
            report.writeSummary("A.g DAQ", secLabel, NOT_APPLY_MESSAGE)
            report.write(NOT_APPLY_WORSE_RESOLUTION % 5 + STATUS_NOT_APPLY)
            return None

        Prot = pwplugin.Domain.importFromPlugin('kiharalab.protocols',
                                                'ProtDAQValidation', doRaise=True)
        prot = project.newProtocol(Prot,
                                objLabel="A.g DAQ",
                                stride=3)
        prot.inputVolume.set(protImportMap.outputVolume)
        prot.inputAtomStruct.set(protAtom.outputPdb)
        if useSlurm:
            sendToSlurm(prot, GPU=True, priority=True if priority else False)
        project.launchProtocol(prot)
        #waitOutput(project, prot, 'outputAtomStruct')
        waitUntilFinishes(project, prot)

        if prot.isFailed():
            report.writeSummary("A.g DAQ", secLabel, ERROR_MESSAGE)
            report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
            return prot

        if prot.isAborted():
            print(PRINT_PROTOCOL_ABORTED + ": " + NAME_DAQ)
            report.writeSummary("A.g DAQ", secLabel, ERROR_ABORTED_MESSAGE)
            report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
            return prot
        
        saveIntermediateData(report.getReportDir(), 'DAQ', True, 'DAQcif', os.path.join(project.getPath(), prot._getPath('outputStructure.cif')), 'cif file containing DAQ scores')
        input_file = os.path.join(project.getPath(), prot._getPath('outputStructure.cif'))
        # emd_26162_pdb_7txz_emv_daq.json
        output_file = os.path.join(project.getPath(), prot._getExtraPath(), "%s_pdb_%s_emv_daq.json" % (emdb_Id.lower().replace('-','_'), pdbdb_Id.lower()))
        json_file = convert_2_json(emdb_Id, pdbdb_Id, method='daq', input_file=input_file, output_file=output_file)
        saveIntermediateData(report.getReportDir(), 'DAQ', True, 'EMV json file', json_file, 'DAQ scores in EMV json format')

    # get histogram
    try:            
        daqValues = []
        if has_precalculated_data and json_data:
            chain_data = json_data["chains"]
            for chain in chain_data:
                ch_seqData = chain["seqData"]
                for ch_residue in ch_seqData:
                    daqValues.append(float(ch_residue["scoreValue"]))
        else:
            # daqDic = prot.parseDAQScores(prot.outputAtomStruct.getFileName())
            # daqValues = [float(x) for x in list(daqDic.values())]
            from pwem.convert.atom_struct import AtomicStructHandler
            cifDic = AtomicStructHandler().readLowLevel(prot._getPath('outputStructure.cif'))
            for name, value in zip(cifDic['_scipion_attributes.name'],cifDic['_scipion_attributes.value']):
                if name=="DAQ_score":
                    daqValues.append(float(value))

        fnDAQHist = os.path.join(report.getReportDir(),"daqHist.png")
        reportHistogram(daqValues,"DAQ", fnDAQHist)
        Dpercentiles = np.percentile(daqValues, np.array([0.025, 0.25, 0.5, 0.75, 0.975])*100)
        saveIntermediateData(report.getReportDir(), 'DAQ', False, 'DAQHistData', daqValues, ['', 'DAQ values to create histogram'])
        saveIntermediateData(report.getReportDir(), 'DAQ', True, 'DAQHist', fnDAQHist, 'DAQ histogram')

        avgDaq = np.mean(daqValues)
        stdDaq = np.std(daqValues)

        msg =\
    """Fig. \\ref{fig:daqHist} shows the histogram of the DAQ values. The mean and standard deviation were %4.1f and
    %4.1f, respectively.
    
    \\begin{figure}[H]
        \centering
        \includegraphics[width=10cm]{%s}
        \\caption{Histogram of the per-residue DAQ values.}
        \\label{fig:daqHist}
    \\end{figure}
    """%(avgDaq, stdDaq, fnDAQHist)
        report.write(msg)

        saveIntermediateData(report.getReportDir(), 'DAQ', False, 'averageDAQ', avgDaq, ['', 'The mean of the DAQ values'])
        saveIntermediateData(report.getReportDir(), 'DAQ', False, 'stdDAQ', stdDaq, ['', 'The standard deviation'])

        # get colored models
        msg = "The atomic model colored by DAQ can be seen in Fig. \\ref{fig:daq}.\n\n"
        if has_precalculated_data:
            pdbFilename = os.path.join(project.getPath(), project.getTmpPath(), pdbdb_Id + '_DAQFromWS.pdb')
            pdbFromWS = open(pdbFilename, 'w')
            pdbFromWS.write(getFileFromWS(pdbdb_Id, 'daq'))
            pdbFromWS.close()
            report.atomicModel("daqView", msg, "Atomic model colored by DAQ", pdbFilename, "fig:daq", bfactor=True, occupancy=False, legendMin=-1, legendMax=1)

        else:
            fnCifDAQ = os.path.join(project.getPath(), prot._getExtraPath("chimeraAttribute_DAQ_score.cif"))
            replaceOcuppancyWithAttribute(os.path.join(project.getPath(),prot.outputAtomStruct.getFileName()), "DAQ_score", fnCifDAQ)
            report.atomicModel("daqView", msg, "Atomic model colored by DAQ", fnCifDAQ, "fig:daq", bfactor=False, occupancy=True, rainbow=False, legendMin=-1, legendMax=1)

        saveIntermediateData(report.getReportDir(), 'DAQ', True, 'DAQView',
                            [os.path.join(report.getReportDir(), 'daqView1.jpg'),
                            os.path.join(report.getReportDir(), 'daqView2.jpg'),
                            os.path.join(report.getReportDir(), 'daqView3.jpg')], 'DAQ views')

    except:
        report.writeSummary("A.g DAQ", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    warnings = []
    testWarnings = False

    DperHomogeneous = isHomogeneous(Dpercentiles[0], Dpercentiles[-1])

    if DperHomogeneous:
        warnings.append("{\\color{red} \\textbf{Program output seems to be too homogeneous. There might " \
                        "be some program issues analyzing the data.}}")

    if stdDaq==0 or testWarnings:
        warnings.append("{\\color{red} \\textbf{Could not be measured.}}")
    if avgDaq<0.5 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The average DAQ is smaller than 0.5.}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the average DAQ score is larger than 0.5.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "A.g DAQ", secLabel)
    if len(warnings)>0:
        report.writeAbstract("DAQ detects some mismatch between the map and its model (see Sec. \\ref{%s}). "%secLabel)

    if not has_precalculated_data:
        return prot
    return

def reportInput(project, report, FNMODEL, writeAtomicModelFailed=False):

    # Get file basename to write it in the report
    basenameFNMODEL = os.path.basename(FNMODEL)

    msg = \
"""
\\section{Atomic model}
\\label{sec:atomicModel}\n\n
Atomic model: %s \\\\
\\\\"""%basenameFNMODEL.replace('_','\_').replace('/','/\-')
    report.write(msg)

    if writeAtomicModelFailed:
        warnings = []
        warnings.append("{\\color{red} \\textbf{Atomic model file not valid. Some programs cannot handle it due to size or related issues.}}")
        report.writeWarningsAndSummary(warnings, "Atomic model", "sec:atomicModel")
        return True

    # try:
    #     h = AtomicStructHandler()
    #     h.read(FNMODEL)
    # except:
    #     warnings = []
    #     warnings.append("{\\color{red} \\textbf{Biopython cannot safely read this atomic model}}")
    #     report.writeWarningsAndSummary(warnings, "Atomic model", "sec:atomicModel")
    #     return True

    # try:
    #     fnPdb = os.path.join(report.getReportDir(),"tmp.pdb")
    #     h.writeAsPdb(fnPdb)
    #     cleanPath(fnPdb)
    # except:
    #     warnings = []
    #     warnings.append("{\\color{red} \\textbf{Biopython cannot safely write this PDB}}")
    #     report.writeWarningsAndSummary(warnings, "Atomic model", "sec:atomicModel")
    #     return True

    msg = "See Fig. \\ref{fig:modelInput}.\\\\"
    report.atomicModel("modelInput", msg, "Input atomic model", FNMODEL, "fig:modelInput")
    return False

def levelA(project, report, EMDB_ID_NUM, protImportMap, FNMODEL, fnPdb, writeAtomicModelFailed, resolution, doMultimodel, mapCoordX, mapCoordY, mapCoordZ, protCreateSoftMask, fnMaskedMapDict, skipAnalysis=False, priority=False):
    
    secLabel = "sec:AAnalysis"
    section = "Level A Analysis"

    if writeAtomicModelFailed:
        reportInput(project, report, FNMODEL, writeAtomicModelFailed)
        protAtom = None
        return protAtom
    else:
        if protImportMap.outputVolume.hasHalfMaps():
            protImportMapWOHalves = importMap(project, "Import map2", protImportMap, mapCoordX, mapCoordY, mapCoordZ, priority=priority)
            protImportForPhenix = protImportMapWOHalves
        else:
            protImportForPhenix = protImportMap

        protAtom = importModel(project, report, "Import atomic", protImportMap, fnPdb, priority=priority)

        skipAnalysis = skipAnalysis or reportInput(project, report, FNMODEL)

        # Quality Measures
        if not skipAnalysis:
            report.writeSection(section, secLabel)

            # Check if map and model are fitted with phenix
            cc_mask_threshold = 0.5
            fitted, protPhenix, dataPhenix, protAtom = checkFitted(project, report, EMDB_ID_NUM, section, secLabel, protImportMap, protAtom, FNMODEL, resolution, mapCoordX, mapCoordY, mapCoordZ, cc_mask_threshold, priority)

            if fitted is False: # Avoid execcuting level A if map and model are not fitted
                return protAtom
            else: # Continue executing level A
                protConvert = convertPDB(project, report, protImportMap, protAtom, priority=priority)
                if protConvert is not None:
                    phenixReporting(project, report, resolution, protPhenix, dataPhenix)
                    fscq(project, report, protImportMap, protAtom, protConvert, protCreateSoftMask, fnMaskedMapDict['fnSoftMaskedMap'], priority=priority)
                    if doMultimodel:
                        multimodel(project, report, protImportMap, protAtom, resolution, priority=priority)
                    guinierModel(project, report, protImportMap, protConvert, resolution, priority=priority)
                    mapq(project, report, protImportMap, protAtom, resolution, priority=priority)
                    emringer(project, report, protImportForPhenix, protAtom, priority=priority)
                    daq(project, report, protImportMap, protAtom, resolution, priority=priority)

    return protAtom
