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
from scipy import ndimage

import pyworkflow.plugin as pwplugin
from validationReport import readMap

def importMap(project, report, label, fnMap, Ts):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportVolumes', doRaise=True)
    fnDir, fnBase = os.path.split(fnMap)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               filesPath=fnDir,
                               filesPattern=fnMap,
                               samplingRate=Ts)
    project.launchProtocol(prot, wait=True)

    subsection = "Orthogonal slices of the input map"
    msg = "\\textbf{Explanation}:\\\\ In the orthogonal slices of the map, the noise outside the protein should not "\
          "have any structure (stripes going out, small blobs, particularly high or low densities, ...)\\\\ \\\\"\
          "\\textbf{Results}:\\\\"
    if prot.isFailed():
        report.writeFailedSection(subsection, msg)
    else:
        report.orthogonalSlices(subsection, msg, "Central slices in the three dimensions", fnMap)

    subsection = "Orthogonal slices of maximum variance of the input map"
    if prot.isFailed():
        report.writeFailedSubsection(subsection, msg)
    else:
        report.orthogonalSlices(subsection, msg, "Central slices in the three dimensions", fnMap, maxVar=True)

    subsection = "Orthogonal projections of the input map"
    msg = "\\textbf{Explanation}:\\\\ In the projections there should not be stripes (this is an indication of "\
          "directional overweighting, or angular attraction), and there should not be a dark halo around or inside the structure"\
          " (this is an indication of incorrect CTF correction or the reconstruction of a biased map).\\\\ \\\\"\
          "\\textbf{Results}:\\\\"
    if prot.isFailed():
        report.writeFailedSubsection(subsection, msg)
    else:
        report.orthogonalProjections(subsection, msg, "Projections in the three dimensions", fnMap)

    return prot

def createMask(project, report, label, map, Ts, threshold):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols.protocol_preprocess',
                                            'XmippProtCreateMask3D', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               inputVolume=map,
                               threshold=threshold,
                               doBig=True,
                               doMorphological=True,
                               elementSize=math.ceil(2/Ts)) # Dilation by 2A
    project.launchProtocol(prot, wait=True)

    subsection = "Orthogonal slices of mask"
    msg = "The mask has been calculated at the suggested threshold %f, and then dilated by 2\AA."%threshold
    if prot.isFailed():
        report.writeFailedSection(subsection, msg)
    else:
        report.orthogonalSlices(subsection, msg, "Central slices in the three dimensions", prot.outputMask.getFileName())

    subsection = "Orthogonal slices of maximum variance of mask"
    if prot.isFailed():
        report.writeFailedSubsection(subsection, msg)
    else:
        report.orthogonalSlices(subsection, msg, "Central slices in the three dimensions", prot.outputMask.getFileName(),
                                maxVar=True)

    return prot

def massAnalysis(project, report, label, volume, mask, Ts):
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

    toWrite=\
"""
\\subsection{Mass analysis}
\\textbf{Explanation:}\\\\
The reconstructed map must be relatively well centered in the box, and there should be at least 30\AA~(the exact size 
depends on the CTF) on each side to make sure that the CTF can be appropriately corrected.
\\\\
\\\\
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

    # Analysis of Center of mass
    cx, cy, cz = ndimage.measurements.center_of_mass(V)
    dcx = abs(cx-X/2)/X*100
    dcy = abs(cy-Y/2)/Y*100
    dcz = abs(cz-Z/2)/Z*100

    toWrite += \
"""

The center of mass is at (x,y,z)=(%6.2f,%6.2f,%6.2f). The decentering of the center of mass (abs(Center)/Size)\\%% is
%5.2f, %5.2f, and %5.2f, respectively.\\%%\\\\
\\\\
"""%(cx,cy,cz,dcx,dcy,dcz)

    warnings=""
    if dx>20:
        warnings+="\\textbf{The volume is significantly decentered in X.\\\\ \n"
    if dy>20:
        warnings+="\\textbf{The volume is significantly decentered in Y.\\\\ \n"
    if dz>20:
        warnings+="\\textbf{The volume is significantly decentered in Z.\\\\ \n"
    if x0<20:
        warnings+="\\textbf{There could be little space from X left to effectively correct for the CTF.\\\\ \n"
    if y0<20:
        warnings+="\\textbf{There could be little space from Y left to effectively correct for the CTF.\\\\ \n"
    if z0<20:
        warnings+="\\textbf{There could be little space from Z left to effectively correct for the CTF.\\\\ \n"
    if xF<20:
        warnings+="\\textbf{There could be little space from X right to effectively correct for the CTF.\\\\ \n"
    if yF<20:
        warnings+="\\textbf{There could be little space from Y right to effectively correct for the CTF.\\\\ \n"
    if zF<20:
        warnings+="\\textbf{There could be little space from Z right to effectively correct for the CTF.\\\\ \n"
    if dcx>20:
        warnings+="\\textbf{The center of mass in X may be significantly shifted.}\\\\ \n"
    if dcy>20:
        warnings+="\\textbf{The center of mass in Y may be significantly shifted.}\\\\ \n"
    if dcz>20:
        warnings+="\\textbf{The center of mass in Z may be significantly shifted.}\\\\ \n"
    if warnings!="":
        toWrite+="\\textbf{WARNINGS}:\\\\"+warnings

    report.write(toWrite)

def xmippDeepRes(project, report, label, map, mask):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtDeepRes', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               inputVolume=map,
                               Mask=mask)
    project.launchProtocol(prot, wait=True)

    subsection = "DeepRes"
    msg = \
"""
Deepres \\cite{Ramirez2020} measures the local resolution using a neural network that has been trained on 
the appearance of atomic structures at different resolutions. Then, by comparing the local appearance of the
input map to the appearance of the atomic structures a local resolution label can be assigned.
"""
    report.addReference("Ram")
    return prot

def reportInput(report, fnMap, Ts, threshold):
    toWrite=\
"""
\\section{Input data}
Input map: %s \\\\
Voxel size: %f (\AA) \\\\
Visualization threshold: %f \\\\

"""%(fnMap.replace('_','\_').replace('/','/\-'), Ts, threshold)
    report.write(toWrite)

def level0(project, report, fnMap, Ts, threshold):
    reportInput(report, fnMap, Ts, threshold)

    # Import map
    protImportMap = importMap(project, report, "import map", fnMap, Ts)
    protCreateMask = createMask(project, report, "create mask", protImportMap.outputVolume, Ts, threshold)

    # Quality Measures
    report.writeSection('Level 0 analysis')
    massAnalysis(project, report, "0.a Mass analysis", protImportMap.outputVolume, protCreateMask.outputMask, Ts)
    # xmippDeepRes(project, report, "0.c deepRes", protImportMap.outputVolume, protCreateMask.outputMask)
    return protImportMap, protCreateMask