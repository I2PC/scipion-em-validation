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

import hashlib
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy

import pyworkflow.plugin as pwplugin
from validationReport import readMap, latexEnumerate

def importMap(project, report, label, fnMap, Ts, threshold):
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
        return
    msg+=" See Fig. \\ref{fig:centralInput}.\\\\"
    report.orthogonalSlices(subsection, msg, "Central slices in the three dimensions", fnMap, "fig:centralInput")

    subsection = "Orthogonal slices of maximum variance of the input map"
    msg=""
    msg+=" See Fig. \\ref{fig:maxVarInput}.\\\\"
    report.orthogonalSlices(subsection, msg, "Slices of maximum variation in the three dimensions", fnMap,
                            "fig:maxVarInput", maxVar=True)

    subsection = "Orthogonal projections of the input map"
    msg = "\\textbf{Explanation}:\\\\ In the projections there should not be stripes (this is an indication of "\
          "directional overweighting, or angular attraction), and there should not be a dark halo around or inside the structure"\
          " (this is an indication of incorrect CTF correction or the reconstruction of a biased map).\\\\ \\\\"\
          "\\textbf{Results}:\\\\"
    msg+= " See Fig. \\ref{fig:projInput}.\\\\"
    report.orthogonalProjections(subsection, msg, "Projections in the three dimensions", fnMap, "fig:projInput")

    subsection = "Isosurface views"
    msg = "\\textbf{Explanation}:\\\\ An isosurface is the surface of all points that have the same gray value. \\\\ \\\\"\
          "\\textbf{Results}:\\\\"
    msg+= " See Fig. \\ref{fig:isoInput}.\\\\"
    report.isoSurfaces(subsection, msg, "Isosurface at threshold=%f"%threshold, fnMap, threshold, "fig:isoInput")

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
        msg+=" See. Fig. \\ref{fig:centralMask}\\\\"
        report.orthogonalSlices(subsection, msg, "Central slices in the three dimensions",
                                prot.outputMask.getFileName(), "fig:centralMask")

    subsection = "Orthogonal slices of maximum variance of mask"
    msg = ""
    if prot.isFailed():
        report.writeFailedSubsection(subsection, msg)
    else:
        msg+=" See. Fig. \\ref{fig:maxVarMask}\\\\"
        report.orthogonalSlices(subsection, msg, "Slices of maximum variation in the three dimensions",
                                prot.outputMask.getFileName(), "fig:maxVarMask", maxVar=True)

    return prot

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

    toWrite=\
"""
\\subsection{Level 0.a Mass analysis}
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
    cx, cy, cz = scipy.ndimage.measurements.center_of_mass(V)
    dcx = abs(cx-X/2)/X*100
    dcy = abs(cy-Y/2)/Y*100
    dcz = abs(cz-Z/2)/Z*100

    toWrite += \
"""
The center of mass is at (x,y,z)=(%6.2f,%6.2f,%6.2f). The decentering of the center of mass (abs(Center)/Size)\\%% is
%5.2f, %5.2f, and %5.2f, respectively.\\%%\\\\

"""%(cx,cy,cz,dcx,dcy,dcz)

    warnings=[]
    testWarnings = False
    countWarnings=0
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
    if len(warnings)>0:
        countWarnings=len(warnings)
        toWrite+="\\textbf{WARNINGS}: %d warnings\\\\ \n%s\n"%(countWarnings,latexEnumerate(warnings))
        report.writeSummary("0.a Mass analysis", "{\\color{red} %d warnings}"%countWarnings)
    else:
        toWrite += "\\textbf{STATUS}: {\\color{blue} OK}\\\\ \n"
        report.writeSummary("0.a Mass analysis", "{\\color{blue} OK}")

    report.write(toWrite)

def maskAnalysis(report, volume, mask, Ts, threshold):
    V = readMap(volume.getFileName()).getData()
    Ts3 = math.pow(Ts,3)

    # Analysis of the raw mask
    rawM = np.where(V>threshold,1,0)

    # Connected components
    structure = np.ones((3, 3, 3), dtype=np.int)
    labeled, ncomponents = scipy.ndimage.measurements.label(rawM, structure)
    sumRawM=np.sum(rawM)
    toWrite=\
"""
\\subsection{Level 0.b Mask analysis}
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
The size and percentage of the total number of voxels for the raw mask are listed below (up to 95\\%% of the mass),
the list contains (No. voxels (volume in \AA$^3$), percentage, cumulatedPercentage):\\\\
\\\\
"""%(threshold, ncomponents, sumRawM, sumRawM*Ts3)

    individualMass = [np.sum(labeled==i) for i in range(1,ncomponents+1)]
    idx = np.argsort(-np.asarray(individualMass)) # Minus is for sorting i descending order
    cumulatedMass = 0
    i = 0
    while cumulatedMass/sumRawM<0.95:
        massi = individualMass[idx[i]]
        cumulatedMass += massi
        if i>0:
            toWrite+=", "
        toWrite+="(%d (%5.2f), %5.2f, %5.2f)"%(massi, massi*Ts3, 100.0*massi/sumRawM, 100.0*cumulatedMass/sumRawM)
        i+=1
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

    report.orthogonalSlices("", "", "Maximum variance slices in the three dimensions of the raw mask", rawM,
                            "fig:rawMask", maxVar=True, fnRoot="rawMask")

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
    for i in range(g.size):
        w[i] = np.sum(V>g[i])
        toWrite+="%5.4f & %5.2f & %5.2f & %5.2f \\\\ \n"%(g[i], w[i], (w[i]*Ts3/(1.207*1000)), w[i]*Ts3/(110*1.207))
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

    matplotlib.use('Agg')
    plt.plot(g,w)
    plt.yscale('log')
    plt.grid(True)
    plt.xlabel('Gray level')
    plt.ylabel('Voxel mass')
    plt.savefig(fnFigMass)

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

    # Warnings
    warnings=[]
    testWarnings = False
    countWarnings=0
    if ncomponents95>5 or testWarnings:
        warnings.append("{\\color{red} \\textbf{There might be a problem of connectivity at this threshold because "\
                        "more than 5 connected components are needed to reach 95\\% of the total mask.}}")
    if avgVolumeRemaining>5 or testWarnings:
        warnings.append("{\\color{red} \\textbf{There might be a problem with noise and artifacts, because the "\
                        "average noise blob has a volume of %f \AA$^3$.}}"%avgVolumeRemaining)
    if overlap<0.75 or testWarnings:
        warnings.append("{\\color{red} \\textbf{There might be a problem in the construction of the mask, because the "\
                        "overlap is smaller than 0.75. A common reason is that the suggested threshold causes too many "\
                        "disconnected components.}}")
    if len(warnings)>0:
        countWarnings=len(warnings)
        toWrite+="\\textbf{WARNINGS}: %d warnings\\\\ \n%s\n"%(countWarnings,latexEnumerate(warnings))
        report.writeSummary("0.b Mask analysis", "{\\color{red} %d warnings}"%countWarnings)
    else:
        toWrite += "\\textbf{STATUS}: {\\color{blue} OK}\\\\ \n"
        report.writeSummary("0.b Mask analysis", "{\\color{blue} OK}")

    report.write(toWrite)

def backgroundAnalysis(report, volume, mask):
    V = readMap(volume.getFileName()).getData()

    toWrite=\
"""

\\subsection{Level 0.c Background analysis}
\\textbf{Explanation:}\\\\
Background is defined as the region outside the macromolecule mask. The background mean should be zero, and the number
of voxels with a very low or very high value (below 5 standard deviations of the noise) should be very small and 
they should be randomly distributed without any specific structure. Sometimes, you can see some structure due to
the symmetry of the structure.
\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
"""

    M = 1-readMap(mask.getFileName()).getData() # Background mask
    Vbg = V[M>0]
    [t,p] = scipy.stats.ttest_1samp(Vbg,0)

    meanBg = np.mean(Vbg)
    stdBg = np.std(Vbg)
    fractionLarge = np.sum(np.abs(Vbg)>5*stdBg)/Vbg.size

    cdf5 = 2*scipy.stats.norm.cdf(-5)
    cdf5Ratio = fractionLarge/cdf5

    toWrite+="The null hypothesis that the background mean is 0 was tested with a one-sample Student's t-test. The "\
             "resulting t-statistic and p-value where %5.2f and %f, respectively.\\\\ \n"\
             "\\\\ \n"\
             "The mean and standard deviation of the background where %f and %f. "\
             "The percentage of background voxels whose absolute value is larger than 5 times the standard "\
             "deviation is %5.2f \\%% (see Fig. \\ref{fig:sigma5}). The same percentage from a Gaussian would be "\
             "%f\\%% (ratio between the two percentages: %f).\\\\ \n\\\\ \n"%\
             (t,p,meanBg, stdBg,fractionLarge*100,cdf5*100,cdf5Ratio)

    Vshooting = np.where(np.logical_and(M, np.abs(V)>5*stdBg),V,0)
    report.orthogonalSlices("", "", "Maximum variance slices in the three dimensions of the parts of the "\
                            "background beyond 5*sigma", Vshooting, "fig:sigma5",
                            maxVar=True, fnRoot="sigma5")
    # Warnings
    warnings=[]
    testWarnings = False
    if p<0.05 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The null hypothesis that the background mean is 0 has been rejected "\
                        "because the p-value of the comparison is smaller than 0.05}}")
    if cdf5Ratio>20 or testWarnings:
        warnings.append("{\\color{red} \\textbf{There is a significant proportion of outlier values in the background "\
                        "(cdf5 ratio=%5.2f})}"%cdf5Ratio)
    if len(warnings)>0:
        countWarnings=len(warnings)
        toWrite+="\\textbf{WARNINGS}: %d warnings\\\\ \n%s\n"%(countWarnings,latexEnumerate(warnings))
        report.writeSummary("0.c Background analysis", "{\\color{red} %d warnings}"%countWarnings)
    else:
        toWrite += "\\textbf{STATUS}: {\\color{blue} OK}\\\\ \n"
        report.writeSummary("0.c Background analysis", "{\\color{blue} OK}")

    report.write(toWrite)

def xmippDeepRes(project, report, label, map, mask):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtDeepRes', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               inputVolume=map,
                               Mask=mask)
    project.launchProtocol(prot, wait=True)

    bblCitation= \
"""\\bibitem[Ram\\'{\i}rez-Aportela et~al., 2019]{Ramirez2019}
Ram\\'{\i}rez-Aportela, E., Mota, J., Conesa, P., Carazo, J.~M., and Sorzano, C.
 O.~S. (2019).
\\newblock {DeepRes}: a new deep-learning- and aspect-based local resolution
 method for electron-microscopy maps.
\\newblock {\em IUCRj}, 6:1054--1063."""
    report.addCitation("Ramirez2019",bblCitation)

    msg = \
"""
\\textbf{Explanation}:\\\\ 
DeepRes \\cite{Ramirez2019} measures the local resolution using a neural network that has been trained on 
the appearance of atomic structures at different resolutions. Then, by comparing the local appearance of the
input map to the appearance of the atomic structures a local resolution label can be assigned.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
"""
    report.writeSubsection("0.d DeepRes", msg)

    if prot.isFailed():
        report.writeSummary("0.d Deepres", "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

    return prot

def reportInput(report, fnMap, Ts, threshold):
    sha256_hash = hashlib.sha256()
    with open(fnMap, "rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    toWrite=\
"""
\\section{Input data}
Input map: %s \\\\
SHA256 hash: %s \\\\ 
Voxel size: %f (\AA) \\\\
Visualization threshold: %f \\\\

"""%(fnMap.replace('_','\_').replace('/','/\-'), sha256_hash.hexdigest(), Ts, threshold)
    report.write(toWrite)

def level0(project, report, fnMap, Ts, threshold, skipAnalysis = False):
    reportInput(report, fnMap, Ts, threshold)

    # Import map
    protImportMap = importMap(project, report, "import map", fnMap, Ts, threshold)
    if protImportMap.isFailed():
        raise Exception("Import map did not work")
    protCreateMask = createMask(project, report, "create mask", protImportMap.outputVolume, Ts, threshold)
    if protCreateMask.isFailed():
        raise Exception("Create mask did not work")

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 0 analysis')
        massAnalysis(report, protImportMap.outputVolume, protCreateMask.outputMask, Ts)
        maskAnalysis(report, protImportMap.outputVolume, protCreateMask.outputMask, Ts, threshold)
        backgroundAnalysis(report, protImportMap.outputVolume, protCreateMask.outputMask)
        xmippDeepRes(project, report, "0.d deepRes", protImportMap.outputVolume, protCreateMask.outputMask)
    return protImportMap, protCreateMask