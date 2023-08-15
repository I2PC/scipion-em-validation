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
import numpy as np
import os
import scipy
import subprocess

from scipion.utils import getScipionHome
import pyworkflow.plugin as pwplugin
from validationReport import readMap, readGuinier, latexEnumerate, calculateSha256, reportPlot, reportMultiplePlots,\
    reportHistogram

import xmipp3

from resourceManager import sendToSlurm, waitOutput, skipSlurm, waitOutputFile, waitUntilFinishes

import configparser

from tools.utils import saveIntermediateData

# used by the ProtImportVolumes protocol, volumes will be downloaded from EMDB
IMPORT_FROM_EMDB = 1

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'config.yaml'))
useSlurm = config['QUEUE'].getboolean('USE_SLURM')
gpuIdSkipSlurm = config['QUEUE'].getint('GPU_ID_SKIP_SLURM')

def importMap(project, report, label, fnMap, fnMap1, fnMap2, Ts, mapCoordX, mapCoordY, mapCoordZ):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportVolumes', doRaise=True)

    fnDir, fnBase = os.path.split(fnMap)
    if mapCoordX is not None and mapCoordY is not None and mapCoordZ is not None:
        prot = project.newProtocol(Prot,
                                objLabel=label,
                                filesPath=os.path.join(fnDir,fnMap),
                                samplingRate=Ts,
                                setOrigCoord=True,
                                x=mapCoordX,
                                y=mapCoordY,
                                z=mapCoordZ)
    else:
        prot = project.newProtocol(Prot,
                                objLabel=label,
                                filesPath=os.path.join(fnDir,fnMap),
                                samplingRate=Ts,
                                setOrigCoord=False)
    if fnMap1 is not None and fnMap2 is not None:
        prot.setHalfMaps.set(True)
        prot.half1map.set(fnMap1)
        prot.half2map.set(fnMap2)

    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputVolume')
    waitUntilFinishes(project, prot)
    saveIntermediateData(report.fnReportDir, 'inputData', True, 'map', str(prot.filesPath), 'map from EMDB')
    if fnMap1 is not None and fnMap2 is not None:
        saveIntermediateData(report.fnReportDir, 'inputData', True, 'map', str(prot.half1map), 'halfmap1 from EMDB')
        saveIntermediateData(report.fnReportDir, "inputData", True, "map", str(prot.half2map), 'halfmap2 from EMDB')
    return prot

def createMask(project, label, map, Ts, threshold):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols.protocol_preprocess',
                                            'XmippProtCreateMask3D', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               inputVolume=map,
                               threshold=threshold,
                               doBig=True,
                               doMorphological=True,
                               doSmooth=True,
                               sigmaConvolution=2.0,
                               elementSize=math.ceil(2/Ts)) # Dilation by 2A
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputMask')
    waitUntilFinishes(project, prot)
    return prot

def resizeProject(project, protMap, protMask, resolution):
    Xdim = protMap.outputVolume.getDim()[0]
    Ts = protMap.outputVolume.getSamplingRate()
    AMap = Xdim * Ts

    TsTarget = resolution/2
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
    protResizeMap.inputVolumes.set(protMap.outputVolume)
    if useSlurm:
        sendToSlurm(protResizeMap)
    project.launchProtocol(protResizeMap)
    #waitOutput(project, protResizeMap, 'outputVol')
    waitUntilFinishes(project, protResizeMap)

    projectPath = os.path.join(project.getPath())
    subprocess.call(['chmod', '-R', 'o+r', projectPath])

    protResizeMask = project.newProtocol(Prot,
                                         objLabel="Resize Mask Ts=%2.1f"%TsTarget,
                                         doResize=True,
                                         resizeSamplingRate=TsTarget,
                                         doWindow=True,
                                         windowOperation=1,
                                         windowSize=Xdimp)
    protResizeMask.inputVolumes.set(protMask.outputMask)
    if useSlurm:
        sendToSlurm(protResizeMask)
    project.launchProtocol(protResizeMask)
    #waitOutput(project, protResizeMask, 'outputVol')
    waitUntilFinishes(project, protResizeMask)

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtPreprocessVolumes', doRaise=True)
    protPreprocessMask = project.newProtocol(Prot,
                                             objLabel="Binarize",
                                             doThreshold=True,
                                             thresholdType=1,
                                             threshold=0.5,
                                             fillType=1)
    protPreprocessMask.inputVolumes.set(protResizeMask.outputVol)
    if useSlurm:
        sendToSlurm(protPreprocessMask)
    project.launchProtocol(protPreprocessMask)
    #waitOutput(project, protPreprocessMask, 'outputVol')
    waitUntilFinishes(project, protPreprocessMask)

    return protResizeMap, protPreprocessMask

def properMask(mask):
    M = readMap(mask.getFileName()).getData()
    totalMass = np.sum(M)
    if totalMass > 0:
        return True
    else:
        return False

def massAnalysis(report, volume, mask, Ts):
    V = readMap(volume.getFileName()).getData()
    M = readMap(mask.getFileName()).getData()

    X,Y,Z = M.shape

    ix = np.where(np.sum(M,axis=(1,2))>0)[0]
    iy = np.where(np.sum(M,axis=(0,2))>0)[0]
    iz = np.where(np.sum(M,axis=(0,1))>0)[0]

    x0 = Ts*(ix[0]) # Left space
    xF = Ts*(X-ix[-1]) # Right space
    y0 = Ts*(iy[0])
    yF = Ts*(Y-iy[-1])
    z0 = Ts*(iz[0])
    zF = Ts*(Z-iz[-1])

    dx = abs(xF-x0)/(Ts*X)*100
    dy = abs(yF-y0)/(Ts*Y)*100
    dz = abs(zF-z0)/(Ts*Z)*100

    secLabel = "sec:massAnalysis"
    toWrite=\
"""
\\subsection{Level 0.a Mass analysis}
\\label{%s}
\\textbf{Explanation:}\\\\
The reconstructed map must be relatively well centered in the box, and there should be at least 30\AA~(the exact size 
depends on the CTF) on each side to make sure that the CTF can be appropriately corrected.
\\\\
\\\\
"""%(secLabel)

    #test
    totalMass = np.sum(M)

    if totalMass > 0:
        ix = np.where(np.sum(M,axis=(1,2))>0)[0]
        iy = np.where(np.sum(M,axis=(0,2))>0)[0]
        iz = np.where(np.sum(M,axis=(0,1))>0)[0]

        x0 = Ts*(ix[0]) # Left space
        xF = Ts*(X-ix[-1]) # Right space
        y0 = Ts*(iy[0])
        yF = Ts*(Y-iy[-1])
        z0 = Ts*(iz[0])
        zF = Ts*(Z-iz[-1])

        dx = abs(xF-x0)/(Ts*X)*100
        dy = abs(yF-y0)/(Ts*Y)*100
        dz = abs(zF-z0)/(Ts*Z)*100

        toWrite+= \
    """
    \\textbf{Results:}\\\\
    The space from the left and right in X are %6.2f and %6.2f \AA, respectively. 
    There is a decentering ratio (abs(Right-Left)/Size)\\%% of %5.2f\\%%\\\\
    \\\\
    The space from the left and right in Y are %6.2f and %6.2f \AA, respectively.
    There is a decentering ratio (abs(Right-Left)/Size)\\%% of %5.2f\\%%\\\\
    \\\\
    The space from the left and right in Z are %6.2f and %6.2f \AA, respectively.
    There is a decentering ratio (abs(Right-Left)/Size)\\%% of %5.2f\\%%\\\\
    \\\\
    """%(x0, xF, dx, y0, yF, dy, z0, zF, dz)
        

        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "leftSpaceX", x0, ['\u212B', 'Space to the left of the reconstructed map in the box on the x-axis in Angstroms'])
        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "rightSpaceX", xF, ['\u212B', 'Space to the right of the reconstructed map in the box on the x-axis in Angstroms'])
        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "decenteringRatioX", dx, ['%', '(abs(Right-Left)/Size) %']) # 'details': [units, description]

        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "leftSpaceY", y0, ['\u212B', 'Space to the left of the reconstructed map in the box on the y-axis in Angstroms'])
        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "rightSpaceY", yF, ['\u212B', 'Space to the right of the reconstructed map in the box on the y-axis in Angstroms'])
        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "decenteringRatioY", dy, ['%', '(abs(Right-Left)/Size) %'])

        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "leftSpaceZ", z0, ['\u212B', 'Space to the left of the reconstructed map in the box on the z-axis in Angstroms'])
        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "rightSpaceZ", zF, ['\u212B', 'Space to the right of the reconstructed map in the box on the z-axis in Angstroms'])
        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "decenteringRatioZ", dz, ['%', '(abs(Right-Left)/Size) %'])

        # Analysis of Center of mass
        cx, cy, cz = scipy.ndimage.measurements.center_of_mass(V)
        dcx = abs(cx-X/2)/X*100
        dcy = abs(cy-Y/2)/Y*100
        dcz = abs(cz-Z/2)/Z*100

        toWrite += \
    """
    The center of mass is at (x,y,z)=(%6.2f,%6.2f,%6.2f). The decentering of the center of mass (abs(Center)/Size)\\%% is
    %5.2f, %5.2f, and %5.2f, respectively.\\\\
    
    """%(cx,cy,cz,dcx,dcy,dcz)
                
        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "centerOfMass", [cx,cy,cz], ['\u212B', '(x,y,z)'])
        saveIntermediateData(report.fnReportDir, "massAnalysis", False, "decenteringCenterOfMass", [dcx,dcy,dcz], ['%', '(abs(Right-Left)/Size) %'])

        warnings=[]
        testWarnings = False
        if dx>20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{The volume might be significantly decentered in X.}}")
        if dy>20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{The volume might be significantly decentered in Y.}}")
        if dz>20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{The volume might be significantly decentered in Z.}}")
        if x0<20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{There could be little space from X left to effectively correct for the CTF.}}")
        if y0<20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{There could be little space from Y left to effectively correct for the CTF.}}")
        if z0<20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{There could be little space from Z left to effectively correct for the CTF.}}")
        if xF<20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{There could be little space from X right to effectively correct for the CTF.}}")
        if yF<20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{There could be little space from Y right to effectively correct for the CTF.}}")
        if zF<20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{There could be little space from Z right to effectively correct for the CTF.}}")
        if dcx>20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{The center of mass in X may be significantly shifted.}}")
        if dcy>20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{The center of mass in Y may be significantly shifted.}}")
        if dcz>20 or testWarnings:
            warnings.append("{\\color{red} \\textbf{The center of mass in Z may be significantly shifted.}}")

    else:
        warnings = []
        warnings.append("{\\color{red} \\textbf{Threshold parameter is too high, try to lower it.}}")

    report.write(toWrite)

    msg=\
"""\\textbf{Automatic criteria}: The validation is OK if 1) the decentering and center of mass less than 20\\% of the map 
dimensions in all directions, and 2) the extra space on each direction is more than 20\\% of the map dimensions. For local
and focused refinement, or similar, warnings are expected. 
\\\\
    
    """
    report.write(msg)

    report.writeWarningsAndSummary(warnings, "0.a Mass analysis", secLabel)
    if len(warnings)==0:
        report.writeAbstract("The map seems to be well centered. ")
    else:
        report.writeAbstract("The map seems to have some problem in its centering or extra space (see Sec. "\
                             "\\ref{%s}). "%secLabel)

def maskAnalysis(report, volume, mask, Ts, threshold):
    V = readMap(volume.getFileName()).getData()
    Ts3 = math.pow(Ts,3)

    # Analysis of the raw mask
    rawM = np.where(V>=threshold,1,0)

    # Connected components
    structure = np.ones((3, 3, 3), dtype=np.int)
    labeled, ncomponents = scipy.ndimage.measurements.label(rawM, structure)
    sumRawM=np.sum(rawM)
    secLabel="sec:maskAnalysis"
    toWrite=\
"""
\\subsection{Level 0.b Mask analysis}
\\label{%s}
\\textbf{Explanation:}\\\\
The map at the suggested threshold should have most of its mass concentrated in a single connected component.
It is normal that after thresholding there are a few thousands of very small, disconnected noise blobs. However,
there total mass should not exceed 10\\%%. 
The raw mask (just thresholding) and the mask constructed for the analysis (thresholding +
largest connected component + dilation) should significantly overlap. Overlap is defined by the overlapping coefficient
(size(Raw AND Constructed)/size(Raw)) that is a number between 0 and 1, the closer to 1, the more they
agree.
\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
\\underline{Raw mask}: At threshold %f, there are %d connected components with a total number of voxels of %d and
a volume of %5.2f \\AA$^3$ (see Fig. \\ref{fig:rawMask}). 
The size and percentage of the total number of voxels for the raw mask are listed below (up to 95\\%% of the mass or
the first 100 clusters, whatever happens first),
the list contains (No. voxels (volume in \AA$^3$), percentage, cumulated percentage):\\\\
\\\\
"""%(secLabel,threshold, ncomponents, sumRawM, sumRawM*Ts3)
    
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "connectedComponents", ncomponents, ['', ''])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "sizeCC", int(sumRawM), ['voxels', 'Size in terms of number of voxels of connected components']) # convert to int() since int64 is not json serializable
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "volumeCC", sumRawM*Ts3, ['\u212B\u00B3', 'Connected components volume in cubic Angstroms'])

    # individualMass = [np.sum(labeled==i) for i in range(1,ncomponents+1)]
    individualMass = np.zeros(ncomponents+1)
    for l in np.nditer(labeled):
        individualMass[l]+=1
    idx = np.argsort(-np.asarray(individualMass)) # Minus is for sorting i descending order
    cumulatedMass = 0
    i = 0
    toWrite2 = ""
    nvoxels95 = []
    volume95 = []
    percentage95 = []
    cumulatedPercentage = []
    while cumulatedMass/sumRawM<0.95:
        if idx[i]>0:
            massi = individualMass[idx[i]]
            cumulatedMass += massi
            if i<100:
                if len(toWrite2)>0:
                    toWrite2 += ", "
                toWrite+="(%d (%5.2f), %5.2f, %5.2f)"%(massi, massi*Ts3, 100.0*massi/sumRawM, 100.0*cumulatedMass/sumRawM)
                nvoxels95.append(massi)
                volume95.append(massi*Ts3)
                percentage95.append(100.0*massi/sumRawM)
                cumulatedPercentage.append(100.0*cumulatedMass/sumRawM)
        i+=1
    toWrite += toWrite2
    ncomponents95 = i
    toWrite+="\\\\ \\\\Number of components to reach 95\\%% of the mass: %d\\\\ \\\\"%ncomponents95
    ncomponentsRemaining = ncomponents-ncomponents95
    voxelsRemaining = sumRawM-cumulatedMass
    avgVolumeRemaining = voxelsRemaining/voxelsRemaining*Ts3
    maxVolumeRemaining = individualMass[idx[ncomponents95]]*Ts3
    minVolumeRemaining = individualMass[idx[-1]]*Ts3
    toWrite+="The average size of the remaining %d components is %5.2f voxels (%5.2f \AA$^3$). "\
              "Their size go from %d voxels (%5.2f \AA$^3$) to %d voxels (%5.2f \AA$^3$). \\\\ \\\\"%\
             (ncomponentsRemaining, voxelsRemaining/ncomponentsRemaining, avgVolumeRemaining,
              int(individualMass[idx[ncomponents95]]), maxVolumeRemaining,
              int(individualMass[idx[-1]]), minVolumeRemaining)
    report.write(toWrite)

    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "size95", nvoxels95, ['voxels', 'Sizes in terms of number of voxels for the list components to reach the 95% of the mass or the first 100 clusters (whatever happens first)'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "volume95", volume95, ['\u212B\u00B3', 'Volumes in cubic Angstroms for the list components to reach the 95% of the mass or the first 100 clusters (whatever happens first)'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "percentage95", percentage95, ['%', 'Percentages for the list components to reach the 95% of the mass or the first 100 clusters (whatever happens first)'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "cumulatedPercentage", cumulatedPercentage, ['%', 'Cumulated percentages for the list components to reach the 95% of the mass or the first 100 clusters (whatever happens first)'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "connectedComponents95", ncomponents95, ['', 'Number of components to reach 95% of the mass'])

    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "ncomponentsRemaining", ncomponentsRemaining, ['', 'Number of remaining components to reach 100% of the mass'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "averageSizeComponentsRemaining", voxelsRemaining/ncomponentsRemaining, ['voxels', 'Average size in terms of number of voxels for the remaining components'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "volumeComponentsRemaining", avgVolumeRemaining, ['\u212B\u00B3', 'Average volume in cubic Angstroms of the remaining components'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "maxSizeComponentsRemaining", int(individualMass[idx[ncomponents95]]), ['voxels', 'Max size in terms of number of voxels of the remaining components'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "maxVolumeComponentsRemaining", maxVolumeRemaining, ['\u212B\u00B3', 'Max volume in cubic Angstroms of the remaining components'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "minSizeComponentsRemaining", int(individualMass[idx[-1]]), ['voxels', 'Min size in terms of number of voxels of the remaining components'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "minVolumeComponentsRemaining", minVolumeRemaining, ['\u212B\u00B3', 'Min volume in cubic Angstroms of the remaining components'])
   
    msg = "The slices of the raw mask can be seen in Fig. \\ref{fig:rawMask}.\\\\"
    report.orthogonalSlices("rawMask", msg, "Maximum variance slices in the three dimensions of the raw mask", rawM,
                            "fig:rawMask", maxVar=True)

    # Analysis of different thresholds
    maxVal = np.max(V)
    stepGray = maxVal/25
    g = np.arange(stepGray, maxVal, stepGray)
    w = np.zeros(g.size)
    toWrite="\nThe following table shows the variation of the mass enclosed at different thresholds "\
            "(see Fig. \\ref{fig:mass}):\n\n"
    toWrite+="\\begin{small}\n"
    toWrite+="\\begin{center}\n"
    toWrite+="\\begin{tabular}{|c|c|c|c|}\n"
    toWrite+="\\hline\n"
    toWrite+="\\textbf{Threshold} & \\textbf{Voxel mass} & \\textbf{Molecular mass(kDa)} & \\textbf{\\# Aminoacids}\\\\ \n"
    toWrite+="\\hline\n"

    mm = []
    aa = []
    for i in range(g.size):
        w[i] = np.sum(V>g[i])
        toWrite+="%5.4f & %5.2f & %5.2f & %5.2f \\\\ \n"%(g[i], w[i], (w[i]*Ts3/(1.207*1000)), w[i]*Ts3/(110*1.207))
        mm.append(w[i]*Ts3/(1.207*1000))
        aa.append(w[i]*Ts3/(110*1.207))

    toWrite+="\\hline\n"
    toWrite += "\\end{tabular}\n"
    toWrite += "\\end{center}\n\n"
    toWrite += "\\end{small}\n"
    fnFigMass = os.path.join(report.getReportDir(), "mass.png")
    toWrite +=\
"""
\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Voxel mass as a function of the gray level.}
    \\label{fig:mass}
\\end{figure}

"""%fnFigMass
    report.write(toWrite)

    reportPlot(g,w, 'Gray level', 'Voxel mass', fnFigMass, yscale="log")

    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "grayLevel", g.tolist(), ['', 'List of thresholds in table and plot'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "voxelMass", w.tolist(), ['voxels', 'List of voxel mass in table and plot'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "molecularMass", mm, ['kDa', ''])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "Aminoacids", aa, ['', ''])

    saveIntermediateData(report.fnReportDir, "maskAnalysis", True, "massPlot", fnFigMass, 'Voxel mass as a function of the gray level')
    
    # Constructed mask
    M = readMap(mask.getFileName()).getData()
    sumM = np.sum(M)
    overlap = np.sum(np.multiply(M,rawM))/sumRawM
    toWrite=\
"""
\\underline{Constructed mask}: After keeping the largest component of the previous mask and dilating it by 2\AA,
there is a total number of voxels of %d and a volume of %5.2f \\AA$^3$. The overlap between the
raw and constructed mask is %5.2f.\\\\

"""%(sumM, sumM*Ts3, overlap)
    report.write(toWrite)

    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "sizeConstructedMask", float(sumM), ['voxels', 'Size in terms of number of voxels of constructued mask'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "volumeConstructedMask", sumM*Ts3, ['\u212B\u00B3', 'Volume in cubic Angstroms of the constructed mask'])
    saveIntermediateData(report.fnReportDir, "maskAnalysis", False, "overlap", overlap, ['', 'Overlap between raw and constructed mask'])

    # Warnings
    warnings=[]
    testWarnings = False
    if ncomponents95>5 or testWarnings:
        warnings.append("{\\color{red} \\textbf{There might be a problem of connectivity at this threshold because "\
                        "more than 5 connected components are needed to reach 95\\% of the total mask. Probably a "\
                        "smaller threshold will not cause this issue.}}")
    if avgVolumeRemaining>5 or testWarnings:
        warnings.append("{\\color{red} \\textbf{There might be a problem with noise and artifacts, because the "\
                        "average noise blob has a volume of %f \AA$^3$.}}"%avgVolumeRemaining)
    if overlap<0.75 or testWarnings:
        warnings.append("{\\color{red} \\textbf{There might be a problem in the construction of the mask, because the "\
                        "overlap is smaller than 0.75. A common reason is that the suggested threshold causes too many "\
                        "disconnected components.}}")

    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if 1) to keep 95\\% of the mass we need to keep at most 5
connected components; and 2) the average volume of the blobs outside the given threshold has a size smaller than
5\\AA$^3$; and 3) the overlap between the raw mask and the mask constructed for the analysis is larger than 75\\%.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "0.b Mask analysis", secLabel)
    if len(warnings)==0:
        report.writeAbstract("There is no problem with the suggested threshold. ")
    else:
        report.writeAbstract("There seems to be a problem with the suggested threshold (see Sec. "\
                             "\\ref{%s}). "%secLabel)

def backgroundAnalysis(report, volume, mask):
    V = readMap(volume.getFileName()).getData()

    secLabel = "sec:bgAnalysis"
    toWrite=\
"""

\\subsection{Level 0.c Background analysis}
\\label{%s}
\\textbf{Explanation:}\\\\
Background is defined as the region outside the macromolecule mask. The background mean should be zero, and the number
of voxels with a very low or very high value (below 5 standard deviations of the noise) should be very small and 
they should be randomly distributed without any specific structure. Sometimes, you can see some structure due to
the symmetry of the structure.
\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
"""%secLabel

    M = 1-readMap(mask.getFileName()).getData() # Background mask
    Vbg = V[M>0]
    [t,p] = scipy.stats.ttest_1samp(Vbg,0)

    meanBg = np.mean(Vbg)
    stdBg = np.std(Vbg)
    fractionLarge = np.sum(np.abs(Vbg)>5*stdBg)/Vbg.size

    cdf5 = 2*scipy.stats.norm.cdf(-5)
    cdf5Ratio = fractionLarge/cdf5

    toWrite+="The null hypothesis that the background mean is 0 was tested with a one-sample Student's t-test. The "\
             "resulting t-statistic and p-value were %5.2f and %f, respectively.\\\\ \n"\
             "\\\\ \n"\
             "The mean and standard deviation (sigma) of the background were %f and %f. "\
             "The percentage of background voxels whose absolute value is larger than 5 times the standard "\
             "deviation is %5.2f \\%% (see Fig. \\ref{fig:sigma5}). The same percentage from a Gaussian would be "\
             "%f\\%% (ratio between the two percentages: %f).\\\\ \n\\\\ \n"%\
             (t,p,meanBg, stdBg,fractionLarge*100,cdf5*100,cdf5Ratio)
    
    saveIntermediateData(report.fnReportDir, "backgroundAnalysis", False, "t-statistic", t, ['', "t-statistics after testing with a a one-sample Student's t-test the null hypothesis that the background mean is 0"])
    saveIntermediateData(report.fnReportDir, "backgroundAnalysis", False, "p-value", p, ['', "p-value after testing with a a one-sample Student's t-test the null hypothesis that the background mean is 0"])
    saveIntermediateData(report.fnReportDir, "backgroundAnalysis", False, "mean", float(meanBg), ['', 'The mean of the background (gray level)'])
    saveIntermediateData(report.fnReportDir, "backgroundAnalysis", False, "standardDeviation", float(stdBg), ['', 'The standard deviation (sigma) of the background (gray level)'])
    saveIntermediateData(report.fnReportDir, "backgroundAnalysis", False, "percentageVoxelsLarger5", fractionLarge*100, ['%', 'The percentage of background voxels whose absolute value is larger than 5 times the standard deviation'])
    saveIntermediateData(report.fnReportDir, "backgroundAnalysis", False, "percentageGaussian", cdf5*100, ['%', 'The same percentage from a Gaussian'])
    saveIntermediateData(report.fnReportDir, "backgroundAnalysis", False, "percentageRatio", cdf5Ratio, ['', 'Ration between the two percentages'])

    report.write(toWrite)

    Vshooting = np.where(np.logical_and(M, np.abs(V)>5*stdBg),V,0)
    msg = "Slices of the background beyond 5*sigma can be seen in Fig. \\ref{fig:sigma5}.\\\\"
    report.orthogonalSlices("sigma5", msg, "Maximum variance slices in the three dimensions of the parts of the "\
                            "background beyond 5*sigma", Vshooting, "fig:sigma5", maxVar=True)

    # Warnings
    warnings=[]
    testWarnings = False
    if p<0.001 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The null hypothesis that the background mean is 0 has been rejected "\
                        "because the p-value of the comparison is smaller than 0.001}}")
    if cdf5Ratio>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{There is a significant proportion of outlier values in the background "\
                        "(cdf5 ratio=%5.2f})}"%cdf5Ratio)
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if 1) the p-value of the null hypothesis that the background
has 0 mean is larger than 0.001; and 2) the number of voxels above or below 5 sigma is smaller than 20 times the 
amount expected for a Gaussian with the same standard deviation whose mean is 0.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "0.c Background analysis", secLabel)
    if len(warnings)==0:
        report.writeAbstract("There is no problem with its background. ")
    else:
        report.writeAbstract("There seems to be a problem with the map's background (see Sec. "\
                             "\\ref{%s}). "%secLabel)


def bFactorAnalysis(report, map, resolution):
    fnIn = map.getFileName()
    if fnIn.endswith(".mrc"):
        fnIn+=":mrc"
    fnOut = os.path.join(report.getReportDir(), "sharpenedMap.mrc")
    Ts = map.getSamplingRate()
    args = "-i %s -o %s --sampling %f --maxres %s --auto"%(fnIn, fnOut, Ts, resolution)

    scipionHome = getScipionHome()
    scipion3 = os.path.join(scipionHome,'scipion3')
    output = subprocess.check_output([scipion3, 'run xmipp_volume_correct_bfactor %s'%args])

    p = subprocess.Popen('%s run xmipp_volume_correct_bfactor %s'%(scipion3, args), shell=True, stderr=subprocess.PIPE)
    outputLines = p.stderr.read().decode('utf-8').split('\n')
    tokens = outputLines[0].split()
    a = float(tokens[2])
    b = float(tokens[4])

    tokens = outputLines[1].split()
    bfactor = float(tokens[3])

    dinv2, lnF, lnFc = readGuinier(fnOut+".guinier")
    fitted = a*dinv2 + b
    fnPlot = os.path.join(report.getReportDir(),'Bfactor.png')
    reportMultiplePlots(dinv2, [lnF, fitted, lnFc], '1/Resolution^2 (1/A^2)', 'log Structure factor', fnPlot,
                        ['Experimental', 'Fitted', 'Corrected'])

    bblCitation = \
        """\\bibitem[Rosenthal and Henderson, 2003]{Rosenthal2003}
        Rosenthal, P.~B. and Henderson, R. (2003).
        \\newblock Optimal determination of particle orientation, absolute hand, and
          contrast loss in single particle electron-cryomicroscopy.
        \\newblock {\em J. {M}olecular {B}iology}, 333:721--745."""
    report.addCitation("Rosenthal2003", bblCitation)

    secLabel = "sec:bfactor"
    msg=\
"""\\subsection{Level 0.d B-factor analysis}
\\label{%s}
\\textbf{Explanation:}\\\\
The B-factor line \\cite{Rosenthal2003} fitted between 15\AA and the resolution reported should have a slope that 
is between 0 and 300 \AA$^2$.
\\\\
\\\\
\\textbf{Results:}\\\\
Fig. \\ref{fig:Bfactor} shows the logarithm (in natural units) of the structure factor (the module squared of the
Fourier transform) of the experimental map, its fitted line, and the corrected map. The estimated B-factor was
%5.1f. The fitted line was $\\log(|F|^2)=%4.1f/R^2 + (%4.1f)$. 

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Guinier plot. The X-axis is the square of the inverse of the resolution in \\AA.}
    \\label{fig:Bfactor}
\\end{figure}

"""%(secLabel, bfactor, a, b, fnPlot)
    report.write(msg)

    saveIntermediateData(report.fnReportDir, "bFactorAnalysis", False, "bfactor", bfactor, ['\u212B\u207B\u00B2', 'The estimated B-factor'])
    saveIntermediateData(report.fnReportDir, "bFactorAnalysis", False, "a", a, ['', ''])
    saveIntermediateData(report.fnReportDir, "bFactorAnalysis", False, "b", b, ['', ''])

    saveIntermediateData(report.getReportDir(), 'bFactorAnalysis', True, 'sharpenedMap.mrc.guinier', os.path.join(report.getReportDir(), 'sharpenedMap.mrc.guinier'), 'sharpenedMap.mrc.guinier file which contain the data to create the guinier plot')
    saveIntermediateData(report.getReportDir(), 'bFactorAnalysis', True, 'guinierPlot', fnPlot, 'guinier plot for B-Factor Analysis')

    msg = "\\underline{\\textbf{Orthogonal slices of maximum variance of the B-factor corrected map}}\\\\"\
          "\\textbf{Results}:\\\\"\
          "See Fig. \\ref{fig:maxVarBfactor}.\\\\"
    report.orthogonalSlices("maxVarSlicesBfactor", "", "Slices of maximum variation in the three dimensions of the "\
                            "B-factor corrected map", fnOut, "fig:maxVarBfactor", maxVar=True)

    # Warnings
    warnings=[]
    testWarnings = False
    if bfactor<-300 or bfactor>0 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The B-factor is out of the interval [-300,0]}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the B-factor is in the range [-300,0].
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "0.d B-factor analysis", secLabel)
    if len(warnings)>0:
        report.writeAbstract("There seems to be a problem with its B-factor (see Sec. \\ref{%s}). "%secLabel)

    return bfactor

def xmippDeepRes(project, report, label, map, mask, resolution):
    bblCitation= \
"""\\bibitem[Ram\\'{\i}rez-Aportela et~al., 2019]{Ramirez2019}
Ram\\'{\i}rez-Aportela, E., Mota, J., Conesa, P., Carazo, J.~M., and Sorzano, C.
 O.~S. (2019).
\\newblock {DeepRes}: a new deep-learning- and aspect-based local resolution
 method for electron-microscopy maps.
\\newblock {\em IUCRj}, 6:1054--1063."""
    report.addCitation("Ramirez2019",bblCitation)

    secLabel = "sec:deepres"
    msg = \
"""
\\subsection{Level 0.e Local resolution with DeepRes}
\\label{%s}
\\textbf{Explanation}:\\\\ 
DeepRes \\cite{Ramirez2019} measures the local resolution using a neural network that has been trained on 
the appearance of atomic structures at different resolutions. Then, by comparing the local appearance of the
input map to the appearance of the atomic structures a local resolution label can be assigned.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
"""%secLabel
    report.write(msg)

    if resolution<2:
        report.writeSummary("0.e DeepRes", secLabel, "{\\color{brown} Does not apply}")
        report.write("This method cannot be applied to maps with a resolution better than 2\\AA.\\\\ \n")
        return None

    if resolution>13:
        report.writeSummary("0.e DeepRes", secLabel, "{\\color{brown} Does not apply}")
        report.write("This method cannot be applied to maps with a resolution worse than 13\\AA.\\\\ \n")
        return None

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtDeepRes', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               inputVolume=map,
                               Mask=mask)
    if useSlurm:
        skipSlurm(prot, gpuIdSkipSlurm)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'resolution_Volume')
    waitUntilFinishes(project, prot)

    if prot.isFailed():
        report.writeSummary("0.e DeepRes", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        deepresStderr = open(os.path.join(project.getPath(), prot.getStderrLog()), "r").read()
        if "ran out of memory trying to allocate" in deepresStderr:
            report.write("{\\color{red} \\textbf{REASON: %s.}}\\\\ \n" % "System ran out of memory. Try to launch it again.")
        return prot

    fnRes = os.path.join(project.getPath(), prot._getExtraPath("deepRes_resolution.vol"))
    if not os.path.exists(fnRes):
        report.writeSummary("0.e DeepRes", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return
    
    #TODO: launch convertion from deepres_resolution_origainlSize.vol to .pdb using EMV scripts and save it using saveIntermediateData

    Vres = xmipp3.Image(fnRes).getData()
    R = Vres[Vres >0]
    fnHist = os.path.join(report.getReportDir(), "deepresHist.png")

    reportHistogram(R, "Local resolution (A)", fnHist)
    Rpercentiles = np.percentile(R, np.array([0.025, 0.25, 0.5, 0.75, 0.975])*100)
    resolutionP = np.sum(R < resolution) / R.size * 100
    report.addResolutionEstimate(Rpercentiles[2])

    toWrite = \
"""
Fig. \\ref{fig:histDeepres} shows the histogram of the local resolution according to DeepRes. Some representative
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
Fig. \\ref{fig:deepresColor} shows some representative views of the local resolution.

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Histogram of the local resolution according to deepres.}
    \\label{fig:histDeepres}
\\end{figure}

""" % (Rpercentiles[0], Rpercentiles[1], Rpercentiles[2], Rpercentiles[3], Rpercentiles[4], resolution,
       resolutionP, fnHist)
    report.write(toWrite)

    saveIntermediateData(report.getReportDir(), 'deepRes', False, 'resolutionPercentiles', Rpercentiles.tolist(), ['\u212B', 'List of local resolution in Angstroms at percentiles 2.5%, 25%, 50%, 75% and 97.5 %'])
    saveIntermediateData(report.getReportDir(), 'deepRes', False, 'resolutionPercentile', resolutionP, ['%', 'The percentile at which the reported resolution is'])
    saveIntermediateData(report.getReportDir(), 'deepRes', False, 'resolutionList', R.tolist(), ['\u212B', 'List of local resolution in Angstroms obtained from DeepRes to create the histogram'])
    saveIntermediateData(report.getReportDir(), 'deepRes', False, 'estimatedResolution', Rpercentiles[2], ['\u212B', 'The estimated resolution (median) in Angstroms obtained from DeepRes'])


    saveIntermediateData(report.getReportDir(), 'deepRes', True, 'deepResHist', fnHist, 'deepRes histogram')

    Ts = map.getSamplingRate()
    report.colorIsoSurfaces("", "Local resolution according to DeepRes.", "fig:deepresColor",
                            project, "deepresViewer",
                            os.path.join(project.getPath(), prot._getExtraPath("originalVolume.vol")), Ts,
                            fnRes, Rpercentiles[0], Rpercentiles[-1])
    saveIntermediateData(report.getReportDir(), 'deepRes', True, 'deepResViewer',
                         [os.path.join(report.getReportDir(), 'deepresViewer1.jpg'),
                          os.path.join(report.getReportDir(), 'deepresViewer2.jpg'),
                          os.path.join(report.getReportDir(), 'deepresViewer3.jpg')], 'deepRes views')

    # Warnings
    warnings = []
    testWarnings = False
    if resolutionP < 0.1 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly high with respect " \
                        "to the local resolution distribution. It occupies the %5.2f percentile}}" % \
                        (resolution, resolutionP))
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the percentile of the user provided resolution is larger than
0.1\\% of the percentile of the local resolution as estimated by DeepRes.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "0.e DeepRes", secLabel)
    return prot

def locBfactor(project, report, label, map, mask, resolution):
    bblCitation = \
"""\\bibitem[Kaur et~al., 2021]{Kaur2021}
Kaur, S., Gomez-Blanco, J., Khalifa, A.~A., Adinarayanan, S., Sanchez-Garcia,
  R., Wrapp, D., McLellan, J.~S., Bui, K.~H., and Vargas, J. (2021).
\\newblock Local computational methods to improve the interpretability and
  analysis of cryo-{EM} maps.
\\newblock {\em Nature Communications}, 12(1):1--12."""
    report.addCitation("Kaur2021", bblCitation)

    secLabel = "sec:locbfactor"
    msg = \
"""
\\subsection{Level 0.f Local B-factor}
\\label{%s}
\\textbf{Explanation}:\\\\ 
LocBfactor \\cite{Kaur2021} estimates a local resolution B-factor by decomposing the input map into a 
local magnitude and phase term using the spiral transform.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('ucm.protocols',
                                            'ProtLocBFactor', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               vol=map,
                               mask_in_molecule=mask,
                               max_res=resolution,
                               numberOfThreads=1)
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'bmap')
    #waitOutputFile(project, prot, "bmap.mrc")
    waitUntilFinishes(project, prot)

    fnBfactor = prot._getExtraPath("bmap.mrc")
    if prot.isFailed() or not os.path.exists(fnBfactor):
        report.writeSummary("0.f LocBfactor", secLabel, "{\\color{brown} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot
    
    #TODO: launch convertion from bmap.mrc to .pdb using EMV scripts and save it using saveIntermediateData


    V = xmipp3.Image(fnBfactor+":mrc").getData()
    M = xmipp3.Image(mask.getFileName()).getData()
    B = V[M>0.5]
    fnHist = os.path.join(report.getReportDir(),"locBfactorHist.png")

    reportHistogram(B, "Local B-factor (A^-2)", fnHist)
    Bpercentiles = np.percentile(B, np.array([0.025, 0.25, 0.5, 0.75, 0.975])*100)

    toWrite = \
"""
Fig. \\ref{fig:histLocBfactor} shows the histogram of the local B-factor according to LocBfactor. Some representative
percentiles are:

\\begin{center}
    \\begin{tabular}{|c|c|}
        \\hline
        \\textbf{Percentile} & \\textbf{Local B-factor (\AA$^{-2}$)} \\\\
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

Fig. \\ref{fig:locBfactorColor} shows some representative views of the local B-factor.

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Histogram of the local B-factor according to LocBfactor.}
    \\label{fig:histLocBfactor}
\\end{figure}

""" % (Bpercentiles[0], Bpercentiles[1], Bpercentiles[2], Bpercentiles[3], Bpercentiles[4], fnHist)
    report.write(toWrite)


    saveIntermediateData(report.getReportDir(), 'locBfactor', False, 'bfactorPercentiles', Bpercentiles.tolist(), ['\u212B\u207B\u00B2', 'List of local resolution B-factor in Angstroms^-2 at percentiles 2.5%, 25%, 50%, 75% and 97.5 %'])
    saveIntermediateData(report.getReportDir(), 'locBfactor', False, 'bfactorList', B.tolist(), ['\u212B\u207B\u00B2', 'List of local resolution B-factor in Angstroms^-2 to create the histogram'])

    saveIntermediateData(report.getReportDir(), 'locBfactor', True, 'locBfactorHist', fnHist, 'locBfactor histogram')

    Ts = map.getSamplingRate()
    report.colorIsoSurfaces("", "Local B-factor according to LocBfactor.", "fig:locBfactorColor",
                            project, "locBfactorViewer",
                            os.path.join(project.getPath(), map.getFileName()), Ts,
                            fnBfactor, Bpercentiles[0], Bpercentiles[-1])
    saveIntermediateData(report.getReportDir(), 'locBfactor', True, 'locBfactorViewer',
                         [os.path.join(report.getReportDir(), 'locBfactorViewer1.jpg'),
                          os.path.join(report.getReportDir(), 'locBfactorViewer2.jpg'),
                          os.path.join(report.getReportDir(), 'locBfactorViewer3.jpg')], 'locBfactor views')

    # Warnings
    warnings=[]
    testWarnings = False
    if Bpercentiles[2]<-300 or Bpercentiles[2]>0 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The median B-factor is out of the interval [-300,0]}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the median B-factor is in the range [-300,0].
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "0.f LocBfactor", secLabel)
    if len(warnings)>0:
        report.writeAbstract("There seems to be a problem with its local B-factor (see Sec. \\ref{%s}). "%secLabel)


def locOccupancy(project, report, label, map, mask, resolution):
    bblCitation = \
"""\bibitem[Kaur et~al., 2021]{Kaur2021}
Kaur, S., Gomez-Blanco, J., Khalifa, A.~A., Adinarayanan, S., Sanchez-Garcia,
  R., Wrapp, D., McLellan, J.~S., Bui, K.~H., and Vargas, J. (2021).
\newblock Local computational methods to improve the interpretability and
  analysis of cryo-{EM} maps.
\newblock {\em Nature Communications}, 12(1):1--12."""
    report.addCitation("Kaur2021", bblCitation)

    secLabel = "sec:locOccupancy"
    msg = \
"""
\\subsection{Level 0.g Local Occupancy}
\\label{%s}
\\textbf{Explanation}:\\\\ 
LocOccupancy \\cite{Kaur2021} estimates the occupancy of a voxel by the macromolecule.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('ucm.protocols',
                                            'ProtLocOccupancy', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               vol=map,
                               mask_in_molecule=mask,
                               max_res=resolution,
                               numberOfThreads=1)
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'omap')
    waitUntilFinishes(project, prot)

    fnOccupancy = prot._getExtraPath("omap.mrc")
    if prot.isFailed() or not os.path.exists(fnOccupancy):
        report.writeSummary("0.g LocOccupancy", secLabel, "{\\color{brown} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot
    
    #TODO: launch convertion from bmap.mrc to .pdb using EMV scripts and save it using saveIntermediateData

    V = xmipp3.Image(fnOccupancy+":mrc").getData()
    M = xmipp3.Image(mask.getFileName()).getData()
    B = V[M>0.5]
    fnHist = os.path.join(report.getReportDir(),"locOccupancyHist.png")

    reportHistogram(B, "Local occupancy", fnHist)
    Bpercentiles = np.percentile(B, np.array([0.025, 0.25, 0.5, 0.75, 0.975])*100)

    toWrite = \
"""
Fig. \\ref{fig:histLocOccupancy} shows the histogram of the local occupancy according to LocOccupancy. Some representative
percentiles are:

\\begin{center}
    \\begin{tabular}{|c|c|}
        \\hline
        \\textbf{Percentile} & \\textbf{Local Occupancy [0-1]} \\\\
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

Fig. \\ref{fig:locOccupancyColor} shows some representative views of the local occupancy.

\\begin{figure}[H]
    \centering
    \includegraphics[width=10cm]{%s}
    \\caption{Histogram of the local occupancy according to LocOccupancy.}
    \\label{fig:histLocOccupancy}
\\end{figure}

""" % (Bpercentiles[0], Bpercentiles[1], Bpercentiles[2], Bpercentiles[3], Bpercentiles[4], fnHist)
    report.write(toWrite)

    saveIntermediateData(report.getReportDir(), 'locOccupancy', False, 'locOccupancyPercentiles', Bpercentiles.tolist(), ['', 'List of local occupancy at percentiles 2.5%, 25%, 50%, 75% and 97.5 %'])
    saveIntermediateData(report.getReportDir(), 'locOccupancy', False, 'locOccupancyList', B.tolist(), ['', 'List of local occupancy to create the histogram'])

    saveIntermediateData(report.getReportDir(), 'locOccupancy', True, 'locOccupancyHist', fnHist, 'locOccupancy histogram')

    Ts = map.getSamplingRate()
    report.colorIsoSurfaces("", "Local occupancy according to LocOccupancy.", "fig:locOccupancyColor",
                            project, "locOccupancyViewer",
                            os.path.join(project.getPath(), map.getFileName()), Ts,
                            fnOccupancy, Bpercentiles[0], Bpercentiles[-1])
    saveIntermediateData(report.getReportDir(), 'locOccupancy', True, 'locOccupancyViewer',
                         [os.path.join(report.getReportDir(), 'locOccupancyViewer1.jpg'),
                          os.path.join(report.getReportDir(), 'locOccupancyViewer2.jpg'),
                          os.path.join(report.getReportDir(), 'locOccupancyViewer3.jpg')], 'locOccupancy views')

    # Warnings
    warnings=[]
    testWarnings = False
    if Bpercentiles[2]<0.5 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The median occupancy is less than 50\\%}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the median occupancy is larger than 50\\%.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "0.g LocOccupancy", secLabel)
    if len(warnings)>0:
        report.writeAbstract("There seems to be a problem with its local occupancy (see Sec. \\ref{%s}). "%secLabel)


def deepHand(project, report, label, resolution, map, threshold):
    secLabel = "sec:deepHand"
    msg = \
"""
\\subsection{Level 0.h Hand correction}
\\label{%s}
\\textbf{Explanation}:\\\\ 
Deep Hand determines the correction of the hand for those maps with a resolution smaller than 5\\AA. The method
calculates a value between 0 (correct hand) and 1 (incorrect hand) using a neural network to assign its hand.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    if resolution>5:
        toWrite="This method cannot be applied to maps whose resolution is worse than 5\AA.\\\\"\
                "\\textbf{STATUS}: {\\color{blue} OK}\\\\ \n"
        report.write(toWrite)
        report.writeSummary("0.h Deep hand", secLabel, "{\\color{blue} OK}")
        return

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtDeepHand', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               inputVolume=map,
                               threshold=threshold)
    if useSlurm:
        sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'outputHand')
    #waitOutput(project, prot, 'outputVol')
    waitUntilFinishes(project, prot)

    if prot.isFailed():
        report.writeSummary("0.h DeepHand", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

    hand = prot.outputHand.get()
    msg="Deep hand assigns a score of %4.3f to the input volume.\\\\ \n"%hand
    report.write(msg)

    saveIntermediateData(report.getReportDir(), 'deepHand', False, 'score', hand, ['', 'Value between 0 and 1 that determines the correction of the hand (0 correct hand, 1 incorrect hand)'])


    # Warnings
    warnings=[]
    testWarnings = False
    if hand>0.5 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The volume seems to be flipped.}}")
    if hand>0.4 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The orientation of the volume is uncertain.}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the deep hand score is smaller than 0.5.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "0.h DeepHand", secLabel)
    if len(warnings)>0:
        report.writeAbstract("There seems to be a problem with the map hand (see Sec. \\ref{%s}). "%secLabel)

def reportInput(project, report, fnMap, Ts, threshold, resolution, protImportMap, protCreateMask):

    # Get file basename to write it in the report
    basenameFnMap = os.path.basename(fnMap)

    toWrite=\
"""
\\section{Input data}
Input map: %s \\\\
SHA256 hash: %s \\\\ 
Voxel size: %f (\AA) \\\\
Visualization threshold: %f \\\\
Resolution estimated by user: %f \\\\

"""%(basenameFnMap.replace('_','\_').replace('/','/\-'), calculateSha256(fnMap), Ts, threshold, resolution)
    report.write(toWrite)

    fnImportMap = os.path.join(project.getPath(),protImportMap.outputVolume.getFileName())
    msg = "\\underline{\\textbf{Orthogonal slices of the input map}}\\\\"\
          "\\textbf{Explanation}:\\\\ In the orthogonal slices of the map, the noise outside the protein should not "\
          "have any structure (stripes going out, small blobs, particularly high or low densities, ...)\\\\ \\\\"\
          "\\textbf{Results}:\\\\"\
          "See Fig. \\ref{fig:centralInput}.\\\\"
    report.orthogonalSlices("centralSlicesInputMap", msg, "Central slices of the input map in the three dimensions",
                            fnImportMap, "fig:centralInput")

    msg = "\\underline{\\textbf{Orthogonal slices of maximum variance of the input map}}\\\\"\
          "\\textbf{Results}:\\\\"\
          "See Fig. \\ref{fig:maxVarInput}.\\\\"
    report.orthogonalSlices("maxVarSlicesInputMap", msg, "Slices of maximum variation in the three dimensions",
                            fnImportMap, "fig:maxVarInput", maxVar=True)

    msg = "\\underline{\\textbf{Orthogonal projections of the input map}}\\\\" \
          "\\textbf{Explanation}:\\\\ In the projections there should not be stripes (this is an indication of " \
          "directional overweighting, or angular attraction), and there should not be a dark halo around or inside the structure" \
          " (this is an indication of incorrect CTF correction or the reconstruction of a biased map).\\\\ \\\\" \
          "\\textbf{Results}:\\\\"\
          "See Fig. \\ref{fig:projInput}.\\\\"
    report.orthogonalProjections("projInput", msg, "Projections in the three dimensions",
                                fnImportMap, "fig:projInput")

    msg = "\\underline{\\textbf{Isosurface views of the input map}}\\\\" \
          "\\textbf{Explanation}:\\\\ An isosurface is the surface of all points that have the same gray value. "\
          "In these views there should not be many artifacts or noise blobs around the map.\\\\ \\\\" \
          "\\textbf{Results}:\\\\"\
          "See Fig. \\ref{fig:isoInput}.\\\\"
    report.isoSurfaces("isoInput", msg, "Isosurface at threshold=%f."%threshold,
                       fnImportMap, threshold, "fig:isoInput")

    fnMask = os.path.join(project.getPath(),protCreateMask.outputMask.getFileName())
    msg = "\\underline{\\textbf{Orthogonal slices of maximum variance of the mask}}\\\\"\
          "\\textbf{Explanation}:\\\\ The mask has been calculated at the suggested threshold %f, "\
          "the largest connected component was selected, and then dilated by 2\AA.\\\\ \\\\"\
          "\\textbf{Results}:\\\\"\
          "See Fig. \\ref{fig:maxVarMask}.\\\\"%threshold
    report.orthogonalSlices("maxVarMask", msg, "Slices of maximum variation in the three dimensions of the mask",
                            fnMask, "fig:maxVarMask")

def level0(project, report, fnMap, fnMap1, fnMap2, Ts, threshold, resolution, mapCoordX, mapCoordY, mapCoordZ, skipAnalysis = False):
    # Import map
    protImportMap = importMap(project, report, "import map", fnMap, fnMap1, fnMap2, Ts, mapCoordX, mapCoordY, mapCoordZ)
    if protImportMap.isFailed():
        raise Exception("Import map did not work")
    saveIntermediateData(report.fnReportDir, 'inputData', False, 'sampling rate', Ts, ['\u212B', 'Sampling rate from EMDB map in Angstroms'])
    saveIntermediateData(report.fnReportDir, 'inputData', False, 'threshold', threshold, ['', 'Threshold from EMDB map'])
    saveIntermediateData(report.fnReportDir, 'inputData', False, 'resolution', resolution, ['\u212B', 'Resolution from EMDB map'])

    protCreateMask = createMask(project, "create mask", protImportMap.outputVolume, Ts, threshold)
    if protCreateMask.isFailed():
        raise Exception("Create mask did not work")
    reportInput(project, report, fnMap, Ts, threshold, resolution, protImportMap, protCreateMask)
    # properMask = properMask(protCreateMask.outputMask)

    # Resize to the given resolution
    protResizeMap, protResizeMask = resizeProject(project, protImportMap, protCreateMask, resolution)

    # Quality Measures
    report.writeSection('Level 0 analysis')
    massAnalysis(report, protImportMap.outputVolume, protCreateMask.outputMask, Ts)
    maskAnalysis(report, protImportMap.outputVolume, protCreateMask.outputMask, Ts, threshold)
    backgroundAnalysis(report, protImportMap.outputVolume, protCreateMask.outputMask)
    bfactor=bFactorAnalysis(report, protImportMap.outputVolume, resolution)

    if not skipAnalysis:
        xmippDeepRes(project, report, "0.e deepRes", protImportMap.outputVolume, protCreateMask.outputMask, resolution)
        locBfactor(project, report, "0.f locBfactor", protResizeMap.outputVol, protResizeMask.outputVol, resolution)
        locOccupancy(project, report, "0.g locOccupancy", protResizeMap.outputVol, protResizeMask.outputVol, resolution)
        deepHand(project, report, "0.h deepHand", resolution, protImportMap.outputVolume, threshold)
    return protImportMap, protCreateMask, bfactor, protResizeMap, protResizeMask
