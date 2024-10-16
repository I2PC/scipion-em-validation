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
from scipy.ndimage import gaussian_filter
import subprocess
from datetime import datetime
from random import randint
from time import sleep

from scipion.utils import getScipionHome
import pyworkflow.plugin as pwplugin
from pyworkflow.utils.path import cleanPattern, cleanPath
from pwem.emlib.image import ImageHandler
from validationReport import readMap, calculateSha256, CDFFromHistogram, CDFpercentile, reportPlot, \
    radialPlot, reportMultiplePlots, reportHistogram, isHomogeneous
import xmipp3

from resourceManager import waitOutput, waitOutputFile, sendToSlurm, waitUntilFinishes, createScriptForSlurm, checkIfJobFinished

import configparser

from tools.utils import saveIntermediateData

from resources.constants import *

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'config.yaml'))
useSlurm = config['QUEUE'].getboolean('USE_SLURM')
N_TASKS = config['QUEUE'].get('N_TASKS')
N_THREADS = config['SCIPION'].get('N_THREADS')

def importMap(project, label, fnMap, Ts, mapCoordX, mapCoordY, mapCoordZ, priority=False):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportVolumes', doRaise=True)
    fnDir, fnBase = os.path.split(fnMap)
    if mapCoordX is not None and mapCoordY is not None and mapCoordZ is not None:
        prot = project.newProtocol(Prot,
                                   objLabel=label,
                                   filesPath=fnDir,
                                   filesPattern=fnMap,
                                   samplingRate=Ts,
                                   setOrigCoord=True,
                                   x=mapCoordX,
                                   y=mapCoordY,
                                   z=mapCoordZ)
    else:
        prot = project.newProtocol(Prot,
                                   objLabel=label,
                                   filesPath=fnDir,
                                   filesPattern=fnMap,
                                   samplingRate=Ts,
                                   setOrigCoord=False)
    if useSlurm:
        sendToSlurm(prot, priority=True if priority else False)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputVolume')
    waitUntilFinishes(project, prot)

    return prot

def findFirstCross(x,y,y0,mode):
    if mode=='lesser':
        ycond = np.array(y)<y0
    else:
        ycond = np.array(y)>y0
    i0 = None
    for i in range(len(ycond)):
        if ycond[i]:
            i0=i
            break
    if i0 is not None and i0>1:
        f = scipy.interpolate.interp1d(y[(i0-1):(i0+1)], x[(i0-1):(i0+1)])
        return f(y0)
    else:
        return None

def globalResolution(project, report, label, protImportMap1, protImportMap2, resolution, priority=False):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtResolution3D', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label)
    prot.inputVolume.set(protImportMap1.outputVolume)
    prot.referenceVolume.set(protImportMap2.outputVolume)

    if useSlurm:
        sendToSlurm(prot, priority=True if priority else False)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputFSC')
    waitUntilFinishes(project, prot)

    secLabel = "sec:globalResolution"
    msg = \
"""
\\subsection{Level 1.a Global resolution}
\\label{%s}
\\textbf{Explanation}: The Fourier Shell Correlation (FSC) between the two half maps is the most standard 
method to determine the global resolution of a map. However, other measures exist such as the Spectral 
Signal-to-Noise Ratio and the Differential Phase Residual. There is a long debate about the right thresholds 
for these measures. Probably, the most clear threshold is the one of the SSNR (SSNR=1). For the DPR we have 
chosen 103.9$^\circ$ and for the FSC, the standard 0.143. For a deep discussion of all these thresholds, see
this \\href{%s}{link}. Note that these thresholds typically result in resolution values that are at the lower
extreme of the local resolution range, meaning that this resolution is normally in the first quarter.
It should not be understood as the average resolution of the map.\\\\
\\\\
Except for the noise, the FSC and DPR should be approximately monotonic. They should not have any ``coming back''
behavior. If they have, this is typically due to the presence of a mask in real space or non-linear processing.\\\\
\\\\
\\textbf{Results:} \\\\
""" % (secLabel, GLOBALRES_DOI)
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("1.a Global resolution", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_GLOBAL_RES)
        report.writeSummary("1.a Global resolution", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    saveIntermediateData(report.getReportDir(), 'globalResolution', True, 'fsc', os.path.join(project.getPath(), prot._getPath('fsc.xmd')), 'fsc xmd file containing data to create plots')
    saveIntermediateData(report.getReportDir(), 'globalResolution', True, 'structureFactor', os.path.join(project.getPath(), prot._getPath('structureFactor.xmd')), 'structureFactor xmd file containing data to create plots')

    md = xmipp3.MetaData()
    md.read(prot._getPath("fsc.xmd"))
    f = md.getColumnValues(xmipp3.MDL_RESOLUTION_FREQ)
    FSC = md.getColumnValues(xmipp3.MDL_RESOLUTION_FRC)
    DPR = md.getColumnValues(xmipp3.MDL_RESOLUTION_DPR)

    fFSC = findFirstCross(f,FSC,0.143,'lesser')
    fDPR = findFirstCross(f,DPR,102.9,'greater')
    if fFSC is None:
        strFSC = "The FSC does not cross the 0.143 threshold."
    else:
        strFSC = "The resolution according to the FSC is %5.2f\\AA."%(1/fFSC)
    if fDPR is None:
        strDPR = "The DPR does not cross the 103.9 threshold."
    else:
        strDPR = "The resolution according to the DPR is %5.2f\\AA."%(1/fDPR)

    def logistic(x, *args):
        a=args[0]
        b=args[1]
        x0=args[2]
        c=args[3]
        return 1-(a * np.reciprocal(1+np.exp(-b * (x-x0)) + c))
    f05 = findFirstCross(f,FSC,0.5,'lesser')
    yFitted=None
    if f05 is not None:
        opt, pcov = scipy.optimize.curve_fit(logistic, f, FSC, p0=[1, 1, f05, 0])
        yFitted = logistic(f, *opt)
        fsc90 = findFirstCross(f, yFitted, 0.9, 'lesser')
        if fsc90:
            strFSC+=" The map information is well preserved (FSC$>$0.9) up to %5.2f\\AA."%(1/fsc90)

    fnFSC = os.path.join(report.getReportDir(), "fsc.png")
    toPlot = [FSC, [0.143]*len(FSC)]
    legends = ['FSC','0.143']
    if yFitted is not None:
        toPlot.append(yFitted.tolist())
        legends.append("Fitted model")
    reportMultiplePlots(f, toPlot,
                        "Resolution (A)", "Fourier Shell Correlation", fnFSC, legends, invertXLabels=True)
    fnDPR = os.path.join(report.getReportDir(), "dpr.png")
    reportMultiplePlots(f[:-2], [DPR[:-2], [103.9]*len(DPR[:-2])],
                        "Resolution (A)", "Differential Phase Residual", fnDPR,
                        ['DPR','103.9'], invertXLabels=True)
    
    saveIntermediateData(report.getReportDir(), 'globalResolution', False, 'FSCresolution', 1/fFSC if fFSC else fFSC, ['\u212B', 'The resolution according to the FSC. FSCresolution = None means that FSC does not cross the 0.143 threshold'])
    saveIntermediateData(report.getReportDir(), 'globalResolution', False, 'DPRresolution', 1/fDPR if fDPR else fDPR, ['\u212B', 'The resolution according to the DPR. DPRresolution = None means that DPR does not cross the 103.9 threshold'])
    saveIntermediateData(report.getReportDir(), 'globalResolution', False, 'preservedFSC90resolution', 1/fsc90 if fsc90 else fsc90, ['\u212B', 'The resolution up to which the map information is well preserved (FSC>0.9)'])

    saveIntermediateData(report.getReportDir(), 'globalResolution', True, 'FSCPlot', fnFSC, 'Plot that shows the FSC and the 0.143 threshold')
    saveIntermediateData(report.getReportDir(), 'globalResolution', True, 'DPRPlot', fnDPR, 'Plot that shows the DPR and the 103.9 threshold')


    # SSNR
    V1 = xmipp3.Image(protImportMap1.outputVolume.getFileName())
    V2 = xmipp3.Image(protImportMap2.outputVolume.getFileName())
    VS = (V1.getData()+V2.getData())/2
    VN = (V1.getData()-V2.getData())/2

    VS2 = np.fft.fftshift(np.absolute(np.fft.fftn(VS)))
    VN2 = np.fft.fftshift(np.absolute(np.fft.fftn(VN)))
    VSSNR = np.divide(VS2,VN2)

    def radial_profile(V):
        z, y, x = np.indices((V.shape))
        center = [int(x/2) for x in V.shape]
        r = np.sqrt((x - center[2]) ** 2 + (y - center[1]) ** 2 + (z - center[0]) ** 2)
        r = r.astype(np.int64)

        tbin = np.bincount(r.ravel(), V.ravel())
        nr = np.bincount(r.ravel())
        radialprofile = tbin / nr
        return radialprofile
    radialSSNR=radial_profile(VSSNR)
    N = int(VSSNR.shape[0]/2)
    Ts = protImportMap1.outputVolume.getSamplingRate()
    f=np.arange(0,N)*2*Ts/VSSNR.shape[0]
    logRadialSSNR = np.log10(radialSSNR[0:N]-1)

    fSSNR=findFirstCross(f,logRadialSSNR,0,'lesser')
    if fSSNR is None:
        strSSNR= "The SSNR does not cross the 1 threshold."
    else:
        strSSNR = "The resolution according to the SSNR is %5.2f\\AA."%(1/fSSNR)

    fnSSNR = os.path.join(report.getReportDir(), "ssnr.png")
    reportMultiplePlots(f, [logRadialSSNR, [0]*len(f)],
                        "Resolution (A)", "log10(Radial SSNR)", fnSSNR,
                        ['log10(SSNR)','0'], invertXLabels=True)
    
    saveIntermediateData(report.getReportDir(), 'globalResolution', False, 'SSNRresolution', 1/fSSNR if fSSNR else fSSNR, ['\u212B', 'The resolution according to the SSNR. SSNRresolution = None means that SSNR does not cross the 1 threshold'])
    saveIntermediateData(report.getReportDir(), 'globalResolution', False, 'f', f.tolist(), ['\u212B\u207B\u00B9', 'frecuency data in SSNR plot'])
    saveIntermediateData(report.getReportDir(), 'globalResolution', False, 'logRadialSSNR', logRadialSSNR.tolist(), ['', 'logRadialSSNR data in SSNR plot'])
    saveIntermediateData(report.getReportDir(), 'globalResolution', True, 'SSNRPlot', fnSSNR, 'Plot that shows the SSNR and the 1 threshold')

    # Mean and uncertainty
    resolutionList = [1/value for value in [fFSC, fDPR, fSSNR] if value is not None]

    msg = \
"""Fig. \\ref{fig:FSC} shows the FSC and the 0.143 threshold. %s\\\\
Fig. \\ref{fig:DPR} shows the DPR and the 103.9$^\circ$ threshold. %s\\\\
Fig. \\ref{fig:SSNR} shows the SSNR and the SSNR=1 threshold. %s\\\\
""" % (strFSC, strDPR, strSSNR)
    report.write(msg)
    if len(resolutionList)>0:
        msg = \
"""The mean resolution between the three methods is %5.2f\AA~and its range is within the interval [%5.2f,%5.2f]\\AA."""  % (np.mean(resolutionList), np.min(resolutionList), np.max(resolutionList))
        report.write(msg)

    saveIntermediateData(report.getReportDir(), 'globalResolution', False, 'meanResolution3', np.mean(resolutionList), ['\u212B', 'The mean resolution between the three methods, FSC, DPR, SSNR'])
    saveIntermediateData(report.getReportDir(), 'globalResolution', False, 'llResolutionInterval', np.min(resolutionList), ['\u212B', 'Lower limit resolution interval for the three methods'])
    saveIntermediateData(report.getReportDir(), 'globalResolution', False, 'ulResolutionInterval', np.max(resolutionList), ['\u212B', 'Upper limit resolution interval for the three methods'])

    msg = \
"""\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Fourier Shell correlation between the two halves.}
    \\label{fig:FSC}
\\end{figure}

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Differential Phase Residual between the two halves.}
    \\label{fig:DPR}
\\end{figure}

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Spectral Signal-to-Noise Ratio estimated from the two halves.}
    \\label{fig:SSNR}
\\end{figure}
        """ % (fnFSC, fnDPR, fnSSNR)
    report.write(msg)

    # Warnings
    warnings=[]
    testWarnings = False
    if fFSC is not None:
        if resolution<0.8/fFSC or testWarnings:
            warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly high with respect "\
                            "to the resolution calculated by the FSC, %5.2f \\AA}}"%(resolution,1.0/fFSC))
    if fDPR is not None:
        if resolution<0.8/fDPR or testWarnings:
            warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly high with respect "\
                            "to the resolution calculated by the DPR, %5.2f\\AA.}}"%(resolution,1.0/fDPR))
    if fSSNR is not None:
        if resolution<0.8/fSSNR or testWarnings:
            warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly high with respect "\
                            "to the resolution calculated by the SSNR, %5.2f\\AA.}}"%(resolution,1.0/fSSNR))
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the user provided resolution is larger than 0.8 times the
resolution estimated by 1) FSC, 2) DPR, and 3) SSNR.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "1.a Global resolution", secLabel)
    if fFSC is not None:
        report.addResolutionEstimate(1 / fFSC)
    if fDPR is not None:
        report.addResolutionEstimate(1 / fDPR)
    if fSSNR is not None:
        report.addResolutionEstimate(1 / fSSNR)

    return prot

def fscPermutation(project, report, label, protImportMap1, protImportMap2, protMask, resolution, priority=False):

    secLabel = "sec:fscPermutation"
    msg = \
"""
\\subsection{Level 1.b FSC permutation}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method (see this \\href{%s}{link} for more details) calculates a global resolution by formulating a hypothesis test in which the 
distribution of the FSC of noise is calculated from the two maps.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, FSCPERMUTATION_DOI)
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('spoc.protocols',
                                            'ProtFscFdrControl', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label)
    prot.halfOne.set(protImportMap1.outputVolume)
    prot.halfTwo.set(protImportMap2.outputVolume)
    prot.mask.set(protMask.outputMask)

    if useSlurm:
        sendToSlurm(prot, priority=True if priority else False)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputFSC')
    waitUntilFinishes(project, prot)

    if prot.isFailed():
        report.writeSummary("1.b FSC permutation", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot
    
    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_FSC_PERMU)
        report.writeSummary("1.b FSC permutation", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    fh = open(prot._getPath("logs/run.stdout"))
    FDRResolution = None
    Bfactor = None
    for line in fh.readlines():
        if "FDR-FSC:" in line:
            tokens = line.split()
            FDRResolution = float(tokens[-2])
        if "B-factor" in line:
            tokens = line.split()
            Bfactor = float(tokens[4])
    fh.close()

    if FDRResolution is None:
        report.writeSummary("1.b FSC permutation", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    fscFile = np.loadtxt(prot._getExtraPath("FSC.txt"))
    fsc = fscFile[0]
    frecuency = fscFile[1]
    reportPlot(fsc, frecuency, "Frecuency (1/A)", "FSC", prot._getExtraPath("FSC.png"))
    saveIntermediateData(report.getReportDir(), 'FSCPermutation', True, 'FSC.txt', os.path.join(project.getPath(), prot._getExtraPath("FSC.txt")), 'fsc file')
    saveIntermediateData(report.getReportDir(), 'FSCPermutation', True, 'FSC.png', os.path.join(project.getPath(), prot._getExtraPath("FSC.png")), 'fsc plot file')

    msg=\
"""The resolution at 1\\%% of FDR was %4.1f. %sFig. \\ref{fig:fdrfsc} shows the
estimated FSC and resolution.

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{FSC and resolution estimated by a permutation test.}
    \\label{fig:fdrfsc}
\\end{figure}

"""%(FDRResolution, f'The estimated B-factor was {Bfactor} ' if Bfactor else ' ', os.path.join(project.getPath(),prot._getExtraPath("FSC.png")))
    report.write(msg)

    saveIntermediateData(report.getReportDir(), 'FSCPermutation', False, 'FDRResolution', FDRResolution, ['\u212B', 'The resolution at 1% of FDR'])
    saveIntermediateData(report.getReportDir(), 'FSCPermutation', False, 'Bfactor', Bfactor, ['\u212B\u207B\u00B2', 'Estimated B-factor'])

    warnings=[]
    testWarnings = False
    if resolution<0.8*FDRResolution or testWarnings:
        warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly high with respect "\
                        "to the resolution calculated by the FSC permutation, %5.2f \\AA}}"%(resolution,FDRResolution))
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the user provided resolution is larger than 0.8 times the
resolution estimated by FSC permutation.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "1.b FSC permutation", secLabel)
    report.addResolutionEstimate(FDRResolution)

    return prot

def blocres(project, report, label, protImportMap, protImportMap1, protImportMap2, protMask, resolution, fnMaskedMap, priority=False):

    secLabel = "sec:blocres"
    msg = \
"""
\\subsection{Level 1.c Local resolution with Blocres}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method (see this \\href{%s}{link} for more details) computes a local Fourier Shell Correlation (FSC) between the two half maps.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, BLOCRES_DOI)
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('bsoft.protocols',
                                            'BsoftProtBlocres', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               fill=0)
    prot.inputVolume.set(protImportMap1.outputVolume)
    prot.inputVolume2.set(protImportMap2.outputVolume)
    prot.mask.set(protMask.outputMask)

    if useSlurm:
        sendToSlurm(prot, nMPIs=10, priority=True if priority else False)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'resolution_Volume')
    waitUntilFinishes(project, prot)

    if prot.isFailed():
        report.writeSummary("1.c Blocres", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot
    
    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_BLOCRES)
        report.writeSummary("1.c Blocres", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    fnBlocres = os.path.join(project.getPath(), prot._getExtraPath("resolutionMap.map"))
    saveIntermediateData(report.getReportDir(), 'blocRes', True, 'resolutionMap.map', fnBlocres, 'Blocres output volume map')

    VblocRes = xmipp3.Image(prot._getExtraPath("resolutionMap.map:mrc")).getData()
    R = VblocRes[VblocRes>0]
    fnHist = os.path.join(report.getReportDir(),"blocRes.png")

    reportHistogram(R, "Local resolution (A)", fnHist)
    Rpercentiles = np.percentile(R, np.array([0.025, 0.25, 0.5, 0.75, 0.975])*100)
    resolutionP = np.sum(R<resolution)/R.size*100
    report.addResolutionEstimate(Rpercentiles[2])

    toWrite = \
"""
Fig. \\ref{fig:histBlocres} shows the histogram of the local resolution according to Blocres. Some representative
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

The reported resolution, %5.2f \AA, is at the percentile %4.1f. 
Fig. \\ref{fig:blocresColor} shows some representative views of the local resolution.

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Histogram of the local resolution according to blocres.}
    \\label{fig:histBlocres}
\\end{figure}

""" % (Rpercentiles[0], Rpercentiles[1], Rpercentiles[2], Rpercentiles[3], Rpercentiles[4], resolution,
       resolutionP, fnHist)
    report.write(toWrite)

    saveIntermediateData(report.getReportDir(), 'blocRes', False, 'resolutionPercentiles', Rpercentiles.tolist(), ['\u212B', 'List of local resolution in Angstroms at percentiles 2.5%, 25%, 50%, 75% and 97.5 %'])
    saveIntermediateData(report.getReportDir(), 'blocRes', False, 'resolutionPercentile', resolutionP, ['%', 'The percentile at which the reported resolution is'])
    saveIntermediateData(report.getReportDir(), 'blocRes', False, 'resolutionList', R.tolist(), ['\u212B', 'List of local resolution in Angstroms obtained from BlocRes to create the histogram'])
    saveIntermediateData(report.getReportDir(), 'blocRes', False, 'estimatedResolution', Rpercentiles[2], ['\u212B', 'The estimated resolution (median) in Angstroms obtained from BlocRes'])

    saveIntermediateData(report.getReportDir(), 'blocRes', True, 'blocResHist', fnHist, 'blocRes histogram')


    Ts = protImportMap.outputVolume.getSamplingRate()
    report.colorIsoSurfaces("", "Local resolution according to Blocres.", "fig:blocresColor",
                            project, "blocresViewer",
                            fnMaskedMap, Ts,
                            os.path.join(project.getPath(), prot._getExtraPath("resolutionMap.map")),
                            Rpercentiles[0], Rpercentiles[-1])
    saveIntermediateData(report.getReportDir(), 'deepRes', True, 'blocResViewer',
                         [os.path.join(report.getReportDir(), 'blocresViewer1.jpg'),
                          os.path.join(report.getReportDir(), 'blocresViewer2.jpg'),
                          os.path.join(report.getReportDir(), 'blocresViewer3.jpg')], 'blocRes views')

    # Warnings
    warnings = []
    testWarnings = False

    RperHomogeneous = isHomogeneous(Rpercentiles[0], Rpercentiles[-1])
    if RperHomogeneous:
        warnings.append("{\\color{red} \\textbf{Program output seems to be too homogeneous. There might " \
                        "be some program issues analyzing the data.}}")

    if resolutionP < 0.1 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly high with respect " \
                        "to the local resolution distribution. It occupies the %5.2f percentile}}" % \
                        (resolution, resolutionP))
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the percentile of the user provided resolution is larger than
0.1\\% of the percentile of the local resolution as estimated by BlocRes.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "1.c Blocres", secLabel)

def resmap(project, report, label,  protImportMap, protImportMap1, protImportMap2, protMask, resolution, fnMaskedMap, priority=False):

    secLabel = "sec:resmap"
    msg = \
"""
\\subsection{Level 1.d Local resolution with Resmap}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method (see this \\href{%s}{link} for more details) is based on a test hypothesis testing of the superiority of signal over noise at different frequencies.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, RESMAP_DOI)
    report.write(msg)

    fnVol1 = os.path.join(report.getReportDir(), "half1.mrc")
    fnVol2 = os.path.join(report.getReportDir(), "half2.mrc")
    ih = ImageHandler()
    ih.convert(protImportMap1.outputVolume, fnVol1)
    ih.convert(protImportMap2.outputVolume, fnVol2)
    Ts = protImportMap1.outputVolume.getSamplingRate()
    scipionHome = getScipionHome()
    resmap = None
    for file in  glob.glob(os.path.join(scipionHome, "software/em/resmap*/bin/ResMap*")):
        if not file.endswith(".so"):
            resmap = file
            break
    if resmap is None:
        report.writeSummary("1.d Resmap", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_NOT_BINARY + STATUS_ERROR_MESSAGE)
        return

    fnMask = os.path.join(project.getPath(),protMask.outputMask.getFileName())
    args = "--doBenchMarking --noguiSplit %s %s --vxSize=%f  --maskVol=%s"%(fnVol1, fnVol2, Ts, fnMask)
    print("Running: %s %s" % (resmap, args))
    cmd = '%s %s' % (resmap, args)

    if not useSlurm:
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        p.wait()
        sleep(120)
    else:
        randomInt = int(datetime.now().timestamp()) + randint(0, 1000000)
        slurmScriptPath = createScriptForSlurm('resmap_' + str(randomInt), report.getReportDir(), cmd, nTasks=int(N_TASKS), priority=priority)
        # send job to queue
        subprocess.Popen('sbatch %s' % slurmScriptPath, shell=True)
        # check if job has finished
        while True:
            if checkIfJobFinished('resmap_' + str(randomInt)):
                break
        sleep(120)

    fnResMap = os.path.join(report.getReportDir() if not useSlurm else os.path.dirname(slurmScriptPath), "half1_ori_resmap.mrc")
    if not os.path.exists(fnResMap):
        report.writeSummary("1.d Resmap", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return
    
    saveIntermediateData(report.getReportDir(), 'resMap', True, 'half1_ori_resmap.mrc', fnResMap, 'Resmap output volume map')

    Vres = xmipp3.Image(fnResMap+":mrc").getData()
    idx = Vres<100
    Vres[np.logical_not(idx)]=np.mean(Vres[idx])
    Vres = gaussian_filter(Vres,sigma=1.5)
    R = Vres[idx]
    if len(R) == 0: # resmap output is empty
        report.writeSummary("1.d Resmap", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_NO_RESULTS + STATUS_ERROR_MESSAGE)
        return
    else:
        fnHist = os.path.join(report.getReportDir(),"resmapHist.png")

        reportHistogram(R, "Local resolution (A)", fnHist)
        Rpercentiles = np.percentile(R, np.array([0.025, 0.25, 0.5, 0.75, 0.975])*100)
        resolutionP = np.sum(R<resolution)/R.size*100
        report.addResolutionEstimate(Rpercentiles[2])

        toWrite = \
    """
    Fig. \\ref{fig:histResmap} shows the histogram of the local resolution according to Resmap. Some representative
    percentiles are:
    
    \\begin{center}
        \\begin{tabular}{|c|c|}
            \\hline
            \\textbf{Percentile} & \\textbf{Resolution(\AA)} \\\\
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
    
    The reported resolution, %5.2f \AA, is at the percentile %4.0f. 
    Fig. \\ref{fig:resmapColor} shows some representative views of the local resolution.
    
    \\begin{figure}[H]
        \centering
        \includegraphics[width=10cm]{%s}
        \\caption{Histogram of the local resolution according to Resmap.}
        \\label{fig:histResmap}
    \\end{figure}
    
    """ % (Rpercentiles[0], Rpercentiles[1], Rpercentiles[2], Rpercentiles[3], Rpercentiles[4], resolution,
           resolutionP, fnHist)
        report.write(toWrite)

        saveIntermediateData(report.getReportDir(), 'resMap', False, 'resolutionPercentiles', Rpercentiles.tolist(), ['\u212B', 'List of local resolution in Angstroms at percentiles 2.5%, 25%, 50%, 75% and 97.5 %'])
        saveIntermediateData(report.getReportDir(), 'resMap', False, 'resolutionPercentile', resolutionP, ['%', 'The percentile at which the reported resolution is'])
        saveIntermediateData(report.getReportDir(), 'resMap', False, 'resolutionList', R.tolist(), ['\u212B', 'List of local resolution in Angstroms obtained from Resmap to create the histogram'])
        saveIntermediateData(report.getReportDir(), 'resMap', False, 'estimatedResolution', Rpercentiles[2], ['\u212B', 'The estimated resolution (median) in Angstroms obtained from Resmap'])

        saveIntermediateData(report.getReportDir(), 'resMap', True, 'resMapHist', fnHist, 'Resmap histogram')


        Ts = protImportMap.outputVolume.getSamplingRate()
        report.colorIsoSurfaces("", "Local resolution according to Resmap.", "fig:resmapColor",
                                project, "resmapViewer",
                                fnMaskedMap, Ts,
                                fnResMap, Rpercentiles[0], Rpercentiles[-1])
        saveIntermediateData(report.getReportDir(), 'resMap', True, 'resMapViewer',
                             [os.path.join(report.getReportDir(), 'resmapViewer1.jpg'),
                              os.path.join(report.getReportDir(), 'resmapViewer2.jpg'),
                              os.path.join(report.getReportDir(), 'resmapViewer3.jpg')], 'resMap views')

        # Warnings
        warnings = []
        testWarnings = False

        RperHomogeneous = isHomogeneous(Rpercentiles[0], Rpercentiles[-1])
        if RperHomogeneous:
            warnings.append("{\\color{red} \\textbf{Program output seems to be too homogeneous. There might " \
                            "be some program issues analyzing the data.}}")

        if resolutionP < 0.1 or testWarnings:
            warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly high with respect " \
                            "to the local resolution distribution. It occupies the %5.2f percentile}}" % \
                            (resolution, resolutionP))
        msg = \
    """\\textbf{Automatic criteria}: The validation is OK if the percentile of the user provided resolution is larger than
    0.1\\% of the percentile of the local resolution as estimated by Resmap.
    \\\\
    
    """
        report.write(msg)
        report.writeWarningsAndSummary(warnings, "1.d Resmap", secLabel)

        #cleanPath(fnResMap)
        cleanPath(fnVol1)
        cleanPath(fnVol2)

def monores(project, report, label, protImportMap, protCreateMask, resolution, fnMaskedMap, priority=False):
    Ts = protImportMap.outputVolume.getSamplingRate()

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtMonoRes', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               useHalfVolumes=True,
                               minRes=2*Ts,
                               maxRes=max(10,5*resolution),
                               numberOfThreads=N_THREADS)
    prot.associatedHalves.set(protImportMap.outputVolume)
    prot.mask.set(protCreateMask.outputMask)
    if useSlurm:
        sendToSlurm(prot, priority=True if priority else False)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'resolution_Volume')
    #waitOutputFile(project, prot, "hist.xmd")
    waitUntilFinishes(project, prot)

    secLabel = "sec:monores"
    msg = \
"""
\\subsection{Level 1.e Local resolution with MonoRes}
\\label{%s}
\\textbf{Explanation}:\\\\ 
MonoRes (see this \\href{%s}{link} for more details) evaluates the local energy of a point with respect to the distribution of 
energy in the noise. This comparison is performed at multiple frequencies and for each one, the monogenic 
transformation separates the amplitude and phase of the input map. Then the energy of the amplitude within the map
is compared to the amplitude distribution observed in the noise, and a hypothesis test is run for every voxel to check
if its energy is signficantly above the level of noise.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, MONORES_DOI)
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("1.e MonoRes", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_MONORES)
        report.writeSummary("1.e MonoRes", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    fnMonoRes = os.path.join(project.getPath(), prot._getExtraPath("monoresResolutionMap.mrc"))
    saveIntermediateData(report.getReportDir(), 'monoRes', True, 'monoresResolutionMap.mrc', fnMonoRes, 'Monores output volume map')

    md = xmipp3.MetaData()
    # md.read(prot._getExtraPath("hist.xmd"))
    md.read(os.path.join(project.getPath(), prot._getExtraPath("hist.xmd")))
    x_axis = md.getColumnValues(xmipp3.MDL_X)
    y_axis = md.getColumnValues(xmipp3.MDL_COUNT)

    fnHistMonoRes = os.path.join(report.getReportDir(), "histMonoRes.png")
    reportPlot(x_axis[:-2], y_axis[:-2], 'Resolution (A)', '# of voxels', fnHistMonoRes, plotType="bar",
               barWidth=(x_axis[-1] - x_axis[0]) / len(x_axis))

    R, RCDF=CDFFromHistogram(x_axis[:-2], y_axis[:-2])
    Rpercentiles = CDFpercentile(R, RCDF, Fp=[0.025, 0.25, 0.5, 0.75, 0.975])
    resolutionP = CDFpercentile(R, RCDF, xp=resolution)
    report.addResolutionEstimate(Rpercentiles[2])

    toWrite=\
"""
Fig. \\ref{fig:histMonores} shows the histogram of the local resolution according to MonoRes. Some representative
percentiles are:

\\begin{center}
    \\begin{tabular}{|c|c|}
        \\hline
        \\textbf{Percentile} & \\textbf{Resolution(\AA)} \\\\
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

The reported resolution, %5.2f \AA, is at the percentile %4.1f. 
Fig. \\ref{fig:monoresColor} shows some representative views of the local resolution

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Histogram of the local resolution according to MonoRes.}
    \\label{fig:histMonores}
\\end{figure}

"""%(Rpercentiles[0], Rpercentiles[1], Rpercentiles[2], Rpercentiles[3], Rpercentiles[4], resolution, resolutionP*100,
     fnHistMonoRes)
    report.write(toWrite)

    saveIntermediateData(report.getReportDir(), 'monoRes', False, 'resolutionPercentiles', [float(Rpercentil) for Rpercentil in Rpercentiles], ['\u212B', 'List of local resolution in Angstroms at percentiles 2.5%, 25%, 50%, 75% and 97.5 %'])
    saveIntermediateData(report.getReportDir(), 'monoRes', False, 'resolutionPercentile', resolutionP*100, ['%', 'The percentile at which the reported resolution is'])
    saveIntermediateData(report.getReportDir(), 'monoRes', False, 'estimatedResolution', [float(Rpercentil) for Rpercentil in Rpercentiles][2], ['\u212B', 'The estimated resolution (median) in Angstroms obtained from MonoRes'])


    saveIntermediateData(report.getReportDir(), 'monoRes', True, 'hist.xmd', os.path.join(project.getPath(), prot._getExtraPath("hist.xmd")), 'hist.xmd file containing data to create histogram')
    saveIntermediateData(report.getReportDir(), 'monoRes', True, 'monoResHist', fnHistMonoRes, 'monoRes histogram')


    report.colorIsoSurfaces("", "Local resolution according to Monores.", "fig:monoresColor",
                            project, "monoresViewer", fnMaskedMap,
                            Ts, prot._getExtraPath("monoresResolutionChimera.mrc"), Rpercentiles[0], Rpercentiles[-1])
    saveIntermediateData(report.getReportDir(), 'monoRes', True, 'monoResViewer',
                         [os.path.join(report.getReportDir(), 'monoresViewer1.jpg'),
                          os.path.join(report.getReportDir(), 'monoresViewer2.jpg'),
                          os.path.join(report.getReportDir(), 'monoresViewer3.jpg')], 'monoRes views')
    # Warnings
    warnings=[]
    testWarnings = False

    RperHomogeneous = isHomogeneous(Rpercentiles[0], Rpercentiles[-1])
    if RperHomogeneous:
        warnings.append("{\\color{red} \\textbf{Program output seems to be too homogeneous. There might " \
                        "be some program issues analyzing the data.}}")

    if resolutionP<0.001 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly high with respect "\
                        "to the local resolution distribution. It occupies the %5.2f percentile}}"%\
                        (resolution,resolutionP*100))
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the percentile of the user provided resolution is larger than
0.1\\% of the percentile of the local resolution as estimated by MonoRes.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "1.e MonoRes", secLabel)
    return prot

def monodir(project, report, label, protImportMap, protCreateMask, resolution, priority=False):

    secLabel = "sec:monodir"
    msg = \
"""
\\subsection{Level 1.f Local and directional resolution with MonoDir}
\\label{%s}
\\textbf{Explanation}:\\\\ 
MonoDir (see this \\href{%s}{link} for more details) extends the concept of local resolution to local and directional resolution by changing 
the shape of the filter applied to the input map. The directional analysis can reveal image alignment problems.\\n

The histogram of best resolution voxels per direction (Directional Histogram 1D) shows how many voxels in the
volume have their maximum resolution in that direction. Directions are arbitrarily numbered from 1 to N. This histogram
should be relatively flat. We perform a Kolmogorov-Smirnov test to check its uniformity. If the null hypothesis is
rejected, then the directional resolution is not uniform. It does not mean that it is wrong, and it could be caused
by several reasons: 1) the angular distribution is not uniform, 2) there are missing directions, 3) there is some
anisotropy in the data (including some preferential directional movement).

Ideally, the radial average of the minimum, maximum, and average resolution at each voxel (note that these are spatial
radial averages) should be flat and as low as possible. If they show some slope, this is associated with
inaccuracies in the angular assignment. These averages make sense when the shells are fully contained within the
protein. As the shells approach the outside of the protein, these radial averages make less sense.
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, MONODIR_DOI)
    report.write(msg)

    if resolution>10:
        report.writeSummary("1.f MonoDir", secLabel, NOT_APPLY_MESSAGE)
        report.write(NOT_APPLY_WORSE_RESOLUTION % 10 + STATUS_NOT_APPLY)
        return None

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtMonoDir', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               fast=True,
                               resstep=0.5,
                               numberOfThreads=N_THREADS)
    prot.inputVolumes.set(protImportMap.outputVolume)
    prot.Mask.set(protCreateMask.outputMask)
    if useSlurm:
        sendToSlurm(prot, priority=True if priority else False)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputVolume_doa')
    #waitOutput(project, prot, 'azimuthalVolume')
    #waitOutput(project, prot, 'radialVolume')
    waitUntilFinishes(project, prot)

    if prot.isFailed():
        report.writeSummary("1.f MonoDir", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot

    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_MONODIR)
        report.writeSummary("1.f MonoDir", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    # 1D Histogram
    fnHistDirMonoDir1 = os.path.join(report.getReportDir(), "histDirMonoDir1D.png")
    md = xmipp3.MetaData()
    md.read(prot._getExtraPath("hist_prefdir.xmd"))
    direction = md.getColumnValues(xmipp3.MDL_X)
    count = md.getColumnValues(xmipp3.MDL_COUNT)
    reportPlot(direction, count, 'Direction index', '# of voxels', fnHistDirMonoDir1, plotType="bar")

    # Test
    randomSample = np.random.choice([x for x in range(int(np.max(direction))+1)], size=50000,
                                    p=np.array(count,dtype=np.float64)/np.sum(count))
    D, p = scipy.stats.kstest(randomSample, scipy.stats.randint(0, int(np.max(direction))).cdf)

    # 2D Histogram
    fnHistDirMonoDir2 = os.path.join(report.getReportDir(), "histDirMonoDir2D.png")
    rot = md.getColumnValues(xmipp3.MDL_ANGLE_ROT)
    tilt = md.getColumnValues(xmipp3.MDL_ANGLE_TILT)
    radialPlot(rot, tilt, count, fnHistDirMonoDir2)

    # Radial averages
    fnMonodirRadial = os.path.join(report.getReportDir(), "monoDirRadial.png")
    md = xmipp3.MetaData()
    md.read(prot._getExtraPath("Radial_averages.xmd"))
    x = md.getColumnValues(xmipp3.MDL_IDX)
    minResolution = md.getColumnValues(xmipp3.MDL_VOLUME_SCORE3)
    maxResolution = md.getColumnValues(xmipp3.MDL_VOLUME_SCORE4)
    avgResolution = md.getColumnValues(xmipp3.MDL_AVG)
    reportMultiplePlots(x, [minResolution, maxResolution, avgResolution], "Radius (voxels)", "Resolution (A)",
                        fnMonodirRadial, ['Min. Resolution', 'Max. Resolution', 'Average Resolution'])
    avgDirResolution = np.mean(avgResolution)

    msg=\
"""Fig. \\ref{fig:histDirMonoDir1} shows the 1D directional histogram and Fig. \\ref{fig:histDirMonoDir2} the 2D
directional histogram. We compared the 1D directional histogram to a uniform distribution using a Kolmogorov-Smirnov
test. The D statistic was %f, and the p-value of the null hypothesis %f.

The radial average of the minimum, maximum and average resolution at each voxel is shown in
Fig. \\ref{fig:monoDirRadial}. The overall mean of the directional resolution is %5.2f

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Histogram 1D of the best direction at each voxel.}
    \\label{fig:histDirMonoDir1}
\\end{figure}
\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Histogram 2D of the best direction at each voxel. The azimuthal rotation is circular, while the tilt
    angle is the radius. The size of the point is proportional to the number of voxels whose maximum resolution
    is in that direction (this count can be seen in Fig. \\ref{fig:histDirMonoDir1}. }
    \\label{fig:histDirMonoDir2}
\\end{figure}

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Radial averages (in space) of the minimum, maximum and average resolution at each voxel.}
    \\label{fig:monoDirRadial}
\\end{figure}
"""%(D, p, avgDirResolution, fnHistDirMonoDir1, fnHistDirMonoDir2, fnMonodirRadial)
    report.write(msg)

    saveIntermediateData(report.getReportDir(), 'monoDir', False, 'D-statistic', D, ['', 'D-statistic after applying Kolmogorov-Smirnov test to compare the 1D directional histogram to a uniform distribution'])
    saveIntermediateData(report.getReportDir(), 'monoDir', False, 'p-value', p, ['', 'p-value of the null hypothesis after applying Kolmogorov-Smirnov test to compare the 1D directional histogram to a uniform distribution'])
    saveIntermediateData(report.getReportDir(), 'monoDir', False, 'avgResolution', avgDirResolution, ['\u212B', 'The overall mean of the directional resolution'])

    saveIntermediateData(report.getReportDir(), 'monoDir', True, 'hist_prefdir', os.path.join(project.getPath(), prot._getExtraPath("hist_prefdir.xmd")), 'hist_prefdir xmd file')
    saveIntermediateData(report.getReportDir(), 'monoDir', True, 'thresholds', os.path.join(project.getPath(), prot._getExtraPath("thresholds.xmd")), 'thresholds xmd file')
    saveIntermediateData(report.getReportDir(), 'monoDir', True, 'Radial_averages', os.path.join(project.getPath(), prot._getExtraPath("Radial_averages.xmd")), 'Radial_averages xmd file')
    saveIntermediateData(report.getReportDir(), 'monoDir', True, 'hist_DoA', os.path.join(project.getPath(), prot._getExtraPath("hist_DoA.xmd")), 'hist_DoA xmd file')
    saveIntermediateData(report.getReportDir(), 'monoDir', True, 'hist_DoA2', os.path.join(project.getPath(), prot._getExtraPath("hist_DoA2.xmd")), 'hist_DoA2 xmd file')
    saveIntermediateData(report.getReportDir(), 'monoDir', True, 'hist_radial', os.path.join(project.getPath(), prot._getExtraPath("hist_radial.xmd")), 'hist_radial xmd file')
    saveIntermediateData(report.getReportDir(), 'monoDir', True, 'hist_azimuthal', os.path.join(project.getPath(), prot._getExtraPath("hist_azimuthal.xmd")), 'hist_azimuthal xmd file')

    saveIntermediateData(report.getReportDir(), 'monoDir', True, 'monoDir1DHist', fnHistDirMonoDir1, 'Histogram 1D of the best direction at each voxel')
    saveIntermediateData(report.getReportDir(), 'monoDir', True, 'monoDir2DHist', fnHistDirMonoDir2, 'Histogram 2D of the best direction at each voxel')
    saveIntermediateData(report.getReportDir(), 'monoDir', True, 'monoDirRadialPlot', fnMonodirRadial, 'Plot of Radial averages (in space) of the minimum, maximum and average resolution at each voxel.')

    # Warnings
    warnings=[]
    testWarnings = False
    if p<0.001 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The distribution of best resolution is not uniform in all directions. "\
                        "The associated p-value is %f.}}"%p)
        report.writeAbstract("The resolution does not seem to be uniform in all directions (see Sec. \\ref{%s}). "%\
                             secLabel)
    if resolution<0.8*avgDirResolution or testWarnings:
        warnings.append("{\\color{red} \\textbf{The resolution reported by the user, %5.2f\\AA, is at least 80\\%% "\
                        "smaller than the average directional resolution, %5.2f \\AA.}}" % (resolution, avgDirResolution))
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if 1) the null hypothesis that the directional resolution is not
uniform is not rejected with a threshold of 0.001 for the p-value, and 2) the resolution provided by the user is not 
smaller than 0.8 times the average directional resolution.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "1.f MonoDir", secLabel)
    report.addResolutionEstimate(avgDirResolution)

def fso(project, report, label, protImportMap, protMask, resolution, priority=False):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtFSO', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               halfVolumesFile=True,
                               estimate3DFSC=False)
    prot.inputHalves.set(protImportMap.outputVolume)
    prot.mask.set(protMask.outputMask)

    if useSlurm:
        sendToSlurm(prot, priority=True if priority else False)
    project.launchProtocol(prot)
    waitUntilFinishes(project, prot)

    secLabel = "sec:fso"
    msg = \
"""
\\subsection{Level 1.g Fourier Shell Occupancy}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method (see this \\href{%s}{link} for more details) calculates the anisotropy of the energy distribution in Fourier shells. This is an indirect measure of
anisotropy of the angular distribution or the presence of heterogeneity. A natural threshold for this measure is 0.5.
However, 0.9 and 0.1 are also interesting values that define the frequency at which the occupancy is 90\\%% and 10\\%%,
respectively. This region is shaded in the plot.
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, FSO_DOI)
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("1.g FSO", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        return prot
    
    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_FSO)
        report.writeSummary("1.g FSO", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    md = xmipp3.MetaData(prot._getExtraPath("fso.xmd"))
    f = md.getColumnValues(xmipp3.MDL_RESOLUTION_FREQ)
    fso = md.getColumnValues(xmipp3.MDL_RESOLUTION_FSO)
    anisotropy = md.getColumnValues(xmipp3.MDL_RESOLUTION_ANISOTROPY)

    f05 = findFirstCross(f,fso,0.5,'lesser')
    f09 = findFirstCross(f,fso,0.9,'lesser')
    f01 = findFirstCross(f,fso,0.1,'lesser')
    if f05 is None:
        strFSO = "The FSO does not cross the 0.5 threshold."
    else:
        strFSO = "The resolution according to the FSO is %5.2f\\AA."%(1/f05)
    if f01 is not None:
        strFSO += " Fourier shells are occupied at between 90 and than 10\\%% in the range [%5.2f,%5.2f]\\AA."%\
                  (1/f09, 1/f01)

    fnFSO = os.path.join(report.getReportDir(), "fso.png")
    reportMultiplePlots(f, [fso, anisotropy, [0.5]*len(f)], "Resolution (A)", "",
                        fnFSO, ['FSO', "Anisotropy", '0.5 threshold'], invertXLabels=True,
                        xshade0=f09, xshadeF=f01)

    md = xmipp3.MetaData(prot._getExtraPath('Resolution_Distribution.xmd'))
    rot = md.getColumnValues(xmipp3.MDL_ANGLE_ROT)
    tilt = md.getColumnValues(xmipp3.MDL_ANGLE_TILT)
    counts = md.getColumnValues(xmipp3.MDL_RESOLUTION_FRC)
    fnContour = os.path.join(report.getReportDir(), "fsoDirectional.png")
    try:
        radialPlot(rot, tilt, counts, fnContour, plotType="contour")
    except:
        pass

    msg = \
        """Fig. \\ref{fig:fso} shows the Fourier Shell Occupancy and its anisotropy. The directional resolution is shown in
        Fig. \\ref{fig:fsoContour}. %s
    
        \\begin{figure}[H]
            \centering
            \includegraphics[width=9cm]{%s}
            \\caption{FSO and anisotropy.}
            \\label{fig:fso}
        \\end{figure}
        """ % (strFSO, fnFSO)
    if os.path.exists(fnContour):
        msg += \
            """
            \\begin{figure}[H]
                \centering
                \includegraphics[width=9cm]{%s}
                \\caption{Directional resolution in the projection sphere.}
                \\label{fig:fsoContour}
            \\end{figure}
            """ % (fnContour)
    report.write(msg)

    saveIntermediateData(report.getReportDir(), 'FSO', False, 'FSOresolution', 1/f05 if f05 else f05, ['\u212B', 'The resolution according to the FSO. FSCresolution = None means that FSC does not cross the 0.5 threshold'])
    saveIntermediateData(report.getReportDir(), 'FSO', False, 'llFSO', 1/f09 if f01 else f01, ['\u212B', 'Lower limit of the range at which the fourier shells are occupied between 90 and than 10%'])
    saveIntermediateData(report.getReportDir(), 'FSO', False, 'ulFSO', 1/f01 if f01 else f01, ['\u212B', 'Upper limit of the range at which the fourier shells are occupied between 90 and than 10%'])

    saveIntermediateData(report.getReportDir(), 'FSO', True, 'fso', os.path.join(project.getPath(), prot._getExtraPath("fso.xmd")), 'fso xmd file')
    saveIntermediateData(report.getReportDir(), 'FSO', True, 'GlobalFSC', os.path.join(project.getPath(), prot._getExtraPath("GlobalFSC.xmd")), 'GlobalFSC xmd file')
    saveIntermediateData(report.getReportDir(), 'FSO', True, 'Resolution_Distribution', os.path.join(project.getPath(), prot._getExtraPath("Resolution_Distribution.xmd")), 'Resolution_Distribution xmd file')

    saveIntermediateData(report.getReportDir(), 'FSO', True, 'fso.png', fnFSO, 'FSO and anisotropy plot')
    saveIntermediateData(report.getReportDir(), 'FSO', True, 'fsoDirectional.png', fnContour, 'Directional resolution in the projection sphere.')

    # Warnings
    warnings=[]
    testWarnings = False
    if (f05 is not None and resolution<0.8/f05) or testWarnings:
        warnings.append("{\\color{red} \\textbf{The resolution reported by the user, %5.2f\\AA, is at least 80\\%% "\
                        "smaller than the resolution estimated by FSO, %5.2f \\AA.}}" % (resolution, 1/f05))
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the resolution provided by the user is not 
smaller than 0.8 times the resolution estimated by the first cross of FSO below 0.5.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "1.g FSO", secLabel)

    cleanPattern(os.path.join(project.getPath(),"fscDirection*.xmd"))
    if f05 is not None:
        report.addResolutionEstimate(1/f05)

    return prot

def resizeMapToTargetResolution(project, map, TsTarget, priority=False):
    Xdim = map.getDim()[0]
    Ts = map.getSamplingRate()
    AMap = Xdim * Ts

    Xdimp = AMap/TsTarget
    Xdimp = int(2*math.floor(Xdimp/2))

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCropResizeVolumes', doRaise=True)
    protResizeMap = project.newProtocol(Prot,
                                        objLabel="Resize Volume Ts=%2.1f"%TsTarget,
                                        doResize=True,
                                        resizeSamplingRate=TsTarget,
                                        doWindow=True,
                                        windowOperation=1,
                                        windowSize=Xdimp)
    protResizeMap.inputVolumes.set(map)
    if useSlurm:
        sendToSlurm(protResizeMap, priority=True if priority else False)
    project.launchProtocol(protResizeMap)
    waitUntilFinishes(project, protResizeMap)

    if protResizeMap.isFailed():
        return None

    return protResizeMap


def fsc3d(project, report, label, protImportMap, protImportMap1, protImportMap2, protCreateSoftMask, resolution, priority=False):
    Xdim = protImportMap.outputVolume.getDim()[0]

    secLabel = "sec:fsc3d"
    msg = \
"""
\\subsection{Level 1.h Fourier Shell Correlation 3D}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method (see this \\href{%s}{link} for more details) analyzes the FSC in different directions and evaluates its homogeneity.
\\\\
\\textbf{Results:}\\\\
\\\\
""" % (secLabel, FSC3D_DOI)
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('fsc3d.protocols',
                                            'Prot3DFSC', doRaise=True)
    prot = project.newProtocol(Prot,
                                objLabel=label,
                                provideHalfMaps=True,
                                applyMask=True,
                                useGpu=True,
                                hpFilter=200,
                                numThr=0)
    prot.inputVolume.set(protImportMap.outputVolume)
    prot.volumeHalf1.set(protImportMap1.outputVolume)
    prot.volumeHalf2.set(protImportMap2.outputVolume)
    prot.maskVolume.set(protCreateSoftMask.outputMask)

    if useSlurm:
        sendToSlurm(prot, GPU=True, priority=True if priority else False)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputVolume')
    waitUntilFinishes(project, prot)

    if prot.isFailed():
        report.writeSummary("1.h FSC3D", secLabel, ERROR_MESSAGE)
        report.write(ERROR_MESSAGE_PROTOCOL_FAILED + STATUS_ERROR_MESSAGE)
        fsc3dStdout = open(os.path.join(project.getPath(), prot.getStderrLog()), "r").read()
        if "ValueError: could not broadcast input array from shape" in fsc3dStdout:
            report.write("{\\color{red} \\textbf{REASON: FSC3D software experienced an internal error.}}\\\\ \n")
        if "UnboundLocalError: cannot access local variable 'highpassfilter_fouriershell' where it is not associated with a value" in fsc3dStdout:
            report.write("{\\color{red} \\textbf{REASON: FSC3D software experienced an internal error.}}\\\\ \n")
        return prot
    
    if prot.isAborted():
        print(PRINT_PROTOCOL_ABORTED + ": " + NAME_FSC3D)
        report.writeSummary("1.h FSC3D", secLabel, ERROR_ABORTED_MESSAGE)
        report.write(ERROR_MESSAGE_ABORTED + STATUS_ERROR_ABORTED_MESSAGE)
        return prot

    Ts = protImportMap.outputVolume.getSamplingRate()
    md=np.genfromtxt(prot._getExtraPath(os.path.join('Results_vol','Plotsvol.csv')), delimiter=' ')
    N=md.shape[0]
    f = np.arange(0,N)*2*Ts/Xdim
    fscx = md[:,0].tolist()
    fscy = md[:,1].tolist()
    fscz = md[:,2].tolist()
    fscg = md[:,4].tolist()

    fx = findFirstCross(f,fscx,0.143,'lesser')
    fy = findFirstCross(f,fscy,0.143,'lesser')
    fz = findFirstCross(f,fscz,0.143,'lesser')
    fg = findFirstCross(f,fscg,0.143,'lesser')
    if fx is None or fy is None or fz is None or fg is None:
        strFSC3D = "The FSC 3D did not cross the 0.143 threshold in at least one direction."
        saveIntermediateData(report.getReportDir(), 'FSC3D', False, 'llResolutionRange', None, ['\u212B', 'The FSC 3D resolution range lower limit. llResolutionRange = None means that FSC 3D does not cross the 0.143 threshold in at least one direction'])
        saveIntermediateData(report.getReportDir(), 'FSC3D', False, 'ulResolutionRange', None, ['\u212B', 'The FSC 3D resolution range upper limit. ulResolutionRange = None means that FSC 3D does not cross the 0.143 threshold in at least one direction'])

    else:
        fList = [1/fx, 1/fy, 1/fz, 1/fg]
        strFSC3D = "The FSC 3D resolutions at a 0.143 threshold in X, Y, and Z are %5.2f, %5.2f, and %5.2f \AA, "\
                    "respectively. The global resolution at the same threshold is %5.2f \AA. The resolution range is "\
                    "[%5.2f,%5.2f]\AA."%(1/fx, 1/fy, 1/fz, 1/fg, np.min(fList), np.max(fList))
        saveIntermediateData(report.getReportDir(), 'FSC3D', False, 'llResolutionRange', np.min(fList), ['\u212B', 'The FSC 3D resolution range lower limit. llResolutionRange = None means that FSC 3D does not cross the 0.143 threshold in at least one direction'])
        saveIntermediateData(report.getReportDir(), 'FSC3D', False, 'ulResolutionRange', np.max(fList), ['\u212B', 'The FSC 3D resolution range upper limit. ulResolutionRange = None means that FSC 3D does not cross the 0.143 threshold in at least one direction'])

    fnDir = os.path.join(project.getPath(),prot._getExtraPath('Results_vol','Plotsvol.jpg'))
    fnHist = os.path.join(project.getPath(),prot._getExtraPath('Results_vol','histogram.png'))
    fnPower = os.path.join(project.getPath(),prot._getExtraPath('Results_vol','FTPlotvol.jpg'))

    msg = \
"""Fig. \\ref{fig:fsc3DDir} shows the FSCs in X, Y, Z, and the global FSC. Fig. \\ref{fig:fsc3DHist} shows the global
FSC and the histogram of the directional FSC. Finally, Fig. \\ref{fig:FSC3DFTPower} shows the rotational average of
the map power in Fourier space. %s

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{FSC in  X, Y, Z, the global FSC, and the Average Cosine Phase.}
    \\label{fig:fsc3DDir}
\\end{figure}

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Global FSC and histogram of the directional FSC.}
    \\label{fig:fsc3DHist}
\\end{figure}

\\begin{figure}[H]
    \centering
    \includegraphics[width=9cm]{%s}
    \\caption{Logarithm of the radial average of the input map power in Fourier space.}
    \\label{fig:FSC3DFTPower}
\\end{figure}

""" % (strFSC3D, fnDir, fnHist, fnPower)
    report.write(msg)

    saveIntermediateData(report.getReportDir(), 'FSC3D', False, 'FSC3DresolutionX', 1/fx if fx else fx, ['?', 'The FSC 3D resolutions at a 0.143 threshold in X. FSC3DresolutionX = None means that FSC 3D does not cross the 0.143 threshold in X'])
    saveIntermediateData(report.getReportDir(), 'FSC3D', False, 'FSC3DresolutionY', 1/fy if fy else fy, ['?', 'The FSC 3D resolutions at a 0.143 threshold in Y. FSC3DresolutionY = None means that FSC 3D does not cross the 0.143 threshold in Y'])
    saveIntermediateData(report.getReportDir(), 'FSC3D', False, 'FSC3DresolutionZ', 1/fz if fz else fz, ['?', 'The FSC 3D resolutions at a 0.143 threshold in Z. FSC3DresolutionZ = None means that FSC 3D does not cross the 0.143 threshold in Z'])
    saveIntermediateData(report.getReportDir(), 'FSC3D', False, 'FSC3DresolutionGlobal', 1/fg if fg else fg, ['', 'The estimated global resolution according to the FSC 3D. FSC3DresolutionGlobal = None means that FSC 3D does not cross the 0.143 threshold'])

    saveIntermediateData(report.getReportDir(), 'FSC3D', True, 'Plotsvol', fnDir, 'Plotsvol image file')
    saveIntermediateData(report.getReportDir(), 'FSC3D', True, 'histogram', fnHist, 'FSC3D histogram')
    saveIntermediateData(report.getReportDir(), 'FSC3D', True, 'FTPlotvol', fnPower, 'FTPlotvol image file')
    saveIntermediateData(report.getReportDir(), 'FSC3D', True, 'histogram_raw',
                            os.path.join(project.getPath(), prot._getExtraPath('Results_vol', 'histogram_raw.csv')), 'histogram_raw csv file')
    saveIntermediateData(report.getReportDir(), 'FSC3D', True, 'Plotsvol',
                            os.path.join(project.getPath(), prot._getExtraPath('Results_vol', 'Plotsvol.csv')), 'Plotsvol csv file')
    saveIntermediateData(report.getReportDir(), 'FSC3D', True, 'ResEMvolOutglobalFSC',
                            os.path.join(project.getPath(), prot._getExtraPath('Results_vol', 'ResEMvolOutglobalFSC.csv')),
                            'ResEMvolOutglobalFSC csv file')
    saveIntermediateData(report.getReportDir(), 'FSC3D', True, 'histogram_values',
                            os.path.join(project.getPath(), prot._getExtraPath('Results_vol', 'histogram_values.lst')), 'histogram_values lst file')

    # Warnings
    warnings=[]
    testWarnings = False
    if (fg is not None and resolution<0.8/fg) or testWarnings:
        warnings.append("{\\color{red} \\textbf{The resolution reported by the user, %5.2f\\AA, is at least 80\\%% "\
                        "smaller than the resolution estimated by FSC3D, %5.2f \\AA.}}" % (resolution, 1/fg))
    if fg is not None or testWarnings:
        warnings.append("{\\color{red} \\textbf{We could not estimate the global FSC3D}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the resolution provided by the user is not 
smaller than 0.8 the resolution estimated by the first cross of the global directional FSC below 0.143.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "1.h FSC3D", secLabel)
    if fg is not None:
        report.addResolutionEstimate(1/fg)


def reportInput(project, report, fnMap1, fnMap2, protImportMap1, protImportMap2):

    # Get file basenames to write them in the report
    basenameFnMap1 = os.path.basename(fnMap1)
    basenameFnMap2 = os.path.basename(fnMap2)

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
"""%(basenameFnMap1.replace('_','\_').replace('/','/\-'), calculateSha256(fnMap1),\
     basenameFnMap2.replace('_','\_').replace('/','/\-'), calculateSha256(fnMap2))
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

def level1(project, report, fnMap1, fnMap2, Ts, resolution, mapCoordX, mapCoordY, mapCoordZ, protImportMap, protCreateHardMask, protCreateSoftMask, fnMaskedMapDict, skipAnalysis = False, priority=False):
    # Import maps
    protImportMap1 = importMap(project, "import half1", fnMap1, Ts, mapCoordX, mapCoordY, mapCoordZ, priority=priority)
    if protImportMap1.isFailed():
        raise Exception("Import map did not work")
    if protImportMap1.isAborted():
        raise Exception("Import map was MANUALLY ABORTED")
    protImportMap2 = importMap(project, "import half2", fnMap2, Ts, mapCoordX, mapCoordY, mapCoordZ, priority=priority)
    if protImportMap2.isFailed():
        raise Exception("Import map did not work")
    if protImportMap2.isAborted():
        raise Exception("Import map was MANUALLY ABORTED")
    reportInput(project, report, fnMap1, fnMap2, protImportMap1, protImportMap2)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 1 analysis')
        globalResolution(project, report, "1.a Global", protImportMap1, protImportMap2, resolution, priority=priority)
        fscPermutation(project, report, "1.b FSC permutation", protImportMap1, protImportMap2, protCreateSoftMask,
                       resolution, priority=priority)
        blocres(project, report, "1.c Blocres", protImportMap, protImportMap1, protImportMap2, protCreateSoftMask, resolution, fnMaskedMapDict['fnHardMaskedMap'], priority=priority)
        resmap(project, report, "1.d Resmap", protImportMap, protImportMap1, protImportMap2, protCreateHardMask, resolution, fnMaskedMapDict['fnHardMaskedMap'], priority=priority)
        monores(project, report, "1.e MonoRes", protImportMap, protCreateHardMask, resolution, fnMaskedMapDict['fnHardMaskedMap'], priority=priority)
        monodir(project, report, "1.f MonoDir", protImportMap, protCreateHardMask, resolution, priority=priority)
        fso(project, report, "1.g FSO", protImportMap, protCreateSoftMask, resolution, priority=priority)
        fsc3d(project, report, "1.h FSC3D", protImportMap, protImportMap1, protImportMap2,
              protCreateSoftMask, resolution, priority=priority)

    return protImportMap1, protImportMap2
