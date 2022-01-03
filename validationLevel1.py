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
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy

import pyworkflow.plugin as pwplugin
from validationReport import readMap, latexEnumerate, calculateSha256, CDFFromHistogram, CDFpercentile, reportPlot, \
    radialPlot, reportMultiplePlots
import xmipp3

def importMap(project, label, fnMap, Ts):
    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportVolumes', doRaise=True)
    fnDir, fnBase = os.path.split(fnMap)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               filesPath=fnDir,
                               filesPattern=fnMap,
                               samplingRate=Ts)
    project.launchProtocol(prot, wait=True)
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

def globalResolution(project, report, label, protImportMap1, protImportMap2, resolution):
    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtResolution3D', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label)
    prot.inputVolume.set(protImportMap1.outputVolume)
    prot.referenceVolume.set(protImportMap2.outputVolume)

    project.launchProtocol(prot, wait=True)

    bblCitation = \
"""\\bibitem[Sorzano et~al., 2017]{Sorzano2017}
Sorzano, C. O.~S., Vargas, J., Oton, J., Abrishami, V., de~la Rosa-Trevin,
  J.~M., Gomez-Blanco, J., Vilas, J.~L., Marabini, R., and Carazo, J.~M.
  (2017).
\\newblock A review of resolution measures and related aspects in {3D} electron
  microscopy.
\\newblock {\em Progress in biophysics and molecular biology}, 124:1--30."""
    report.addCitation("Sorzano2017", bblCitation)

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
\\cite{Sorzano2017}. Note that these thresholds typically result in resolution values that are at the lower
extreme of the local resolution range, meaning that this resolution is normally in the first quarter.
It should not be understood as the average resolution of the map.\\\\
\\\\
Except for the noise, the FSC and DPR should be approximately monotonic. They should not have any ``coming back''
behavior. If they have, this is typically due to the presence of a mask in real space or non-linear processing.\\\\
\\textbf{Results:} \\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("1.a Global resolution", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

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
        r = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2)
        r = r.astype(np.int)

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
    if fDPR is None:
        strSSNR= "The SSNR does not cross the 1 threshold."
    else:
        strSSNR = "The resolution according to the SSNR is %5.2f\\AA."%(1/fSSNR)

    fnSSNR = os.path.join(report.getReportDir(), "ssnr.png")
    reportMultiplePlots(f, [logRadialSSNR, [0]*len(f)],
                        "Resolution (A)", "log10(Radial SSNR)", fnSSNR,
                        ['log10(SSNR)','0'], invertXLabels=True)

    # Mean and uncertainty
    resolutionList = [1/fFSC, 1/fDPR, 1/fSSNR]

    msg = \
"""Fig. \\ref{fig:FSC} shows the FSC and the 0.143 threshold. %s\\\\
Fig. \\ref{fig:DPR} shows the DPR and the 103.9$^\circ$ threshold. %s\\\\
Fig. \\ref{fig:SSNR} shows the SSNR and the SSNR=1 threshold. %s\\\\
The mean resolution between the three methods is %5.2f\AA~and its range is within the interval [%5.2f,%5.2f]\\AA.

\\begin{figure}[H]
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
        """ % (strFSC, strDPR, strSSNR, np.mean(resolutionList), np.min(resolutionList), np.max(resolutionList),
               fnFSC, fnDPR, fnSSNR)
    report.write(msg)

    # Warnings
    warnings=[]
    testWarnings = False
    if resolution<0.8/fFSC or testWarnings:
        warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly with respect "\
                        "to the resolution calculated by the FSC, %5.2f \\AA}}"%(resolution,1.0/fFSC))
    if resolution<0.8/fDPR or testWarnings:
        warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly with respect "\
                        "to the resolution calculated by the DPR, %5.2f\\AA.}}"%(resolution,1.0/fDPR))
    if resolution<0.8/fSSNR or testWarnings:
        warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly with respect "\
                        "to the resolution calculated by the SSNR, %5.2f\\AA.}}"%(resolution,1.0/fSSNR))
    report.writeWarningsAndSummary(warnings, "1.a Global resolution", secLabel)

    return prot

def blocres(project, report, label, map, mask):
    bblCitation = \
"""\\bibitem[Cardone et~al., 2013]{Cardone2013}
Cardone, G., Heymann, J.~B., and Steven, A.~C. (2013).
\\newblock One number does not fit all: Mapping local variations in resolution
  in cryo-em reconstructions.
\\newblock {\em J. Structural Biology}, 184:226--236."""
    report.addCitation("Cardone2013", bblCitation)

    secLabel = "sec:blocres"
    msg = \
"""
\\subsection{Level 1.b Local resolution with Blocres}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method \\cite{Cardone2013} computes a local Fourier Shell Correlation (FSC) between the two half maps.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    report.writeSummary("1.b Blocres", secLabel, "{\\color{red} Binary installation fails}")
    report.write("{\\color{red} \\textbf{ERROR: Binary installation fails.}}\\\\ \n")

def resmap(project, report, label, map, mask):
    bblCitation = \
"""\\bibitem[Kucukelbir et~al., 2014]{Kucukelbir2014}
Kucukelbir, A., Sigworth, F.~J., and Tagare, H.~D. (2014).
\\newblock Quantifying the local resolution of cryo-{EM} density maps.
\\newblock {\em Nature Methods}, 11:63--65."""
    report.addCitation("Kucukelbir2014", bblCitation)

    secLabel = "sec:resmap"
    msg = \
"""
\\subsection{Level 1.c Local resolution with Resmap}
\\label{%s}
\\textbf{Explanation}:\\\\ 
This method \\cite{Kucukelbir2014} is based on a test hypothesis testing of the superiority of signal over noise at different frequencies.\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)

    report.writeSummary("1.c Resmap", secLabel, "{\\color{red} Not fully automatic}")
    report.write("{\\color{red} \\textbf{ERROR: Not fully automatic.}}\\\\ \n")

def monores(project, report, label, protImportMap, protCreateMask, resolution):
    Ts = protImportMap.outputVolume.getSamplingRate()

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtMonoRes', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               useHalfVolumes=True,
                               minRes=2*Ts,
                               maxRes=max(10,5*resolution))
    prot.associatedHalves.set(protImportMap.outputVolume)
    prot.mask.set(protCreateMask.outputMask)
    project.launchProtocol(prot, wait=True)

    bblCitation = \
"""\\bibitem[Vilas et~al., 2018]{Vilas2018}
Vilas, J.~L., G{\\'o}mez-Blanco, J., Conesa, P., Melero, R., de~la
  Rosa~Trev\\'{\i}n, J.~M., Ot{\\'o}n, J., Cuenca, J., Marabini, R., Carazo,
  J.~M., Vargas, J., and Sorzano, C. O.~S. (2018).
\\newblock {MonoRes}: automatic and unbiased estimation of local resolution for
  electron microscopy maps.
\\newblock {\em Structure}, 26:337--344."""
    report.addCitation("Vilas2018", bblCitation)

    secLabel = "sec:monores"
    msg = \
"""
\\subsection{Level 1.d Local resolution with MonoRes}
\\label{%s}
\\textbf{Explanation}:\\\\ 
MonoRes \\cite{Vilas2018} evaluates the local energy of a point with respect to the distribution of 
energy in the noise. This comparison is performed at multiple frequencies and for each one, the monogenic 
transformation separates the amplitude and phase of the input map.\\\\
\\\\
\\textbf{Results:}\\\\
\\\\
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("1.d MonoRes", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

    md = xmipp3.MetaData()
    md.read(prot._getExtraPath("hist.xmd"))
    x_axis = md.getColumnValues(xmipp3.MDL_X)
    y_axis = md.getColumnValues(xmipp3.MDL_COUNT)

    fnHistMonoRes = os.path.join(report.getReportDir(), "histMonoRes.png")
    reportPlot(x_axis[:-2], y_axis[:-2], 'Resolution (A)', '# of voxels', fnHistMonoRes, plotType="bar",
               barWidth=(x_axis[-1] - x_axis[0]) / len(x_axis))

    R, RCDF=CDFFromHistogram(x_axis[:-2], y_axis[:-2])
    Rpercentiles = CDFpercentile(R, RCDF, Fp=[0.025, 0.25, 0.5, 0.75, 0.975])
    resolutionP = CDFpercentile(R, RCDF, xp=resolution)

    toWrite=\
"""
Fig. \\ref{fig:histMonores} shows the histogram of the local resolution according to MonoRes. Some representative
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

The reported resolution, %5.2f \AA, is at the percentile %5.2f. 
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

    report.colorIsoSurfaces("", "Local resolution according to Monores.", "fig:monoresColor",
                            project, "monoresViewer", protImportMap.outputVolume.getFileName(),
                            Ts, prot._getExtraPath("monoresResolutionChimera.mrc"), Rpercentiles[0], Rpercentiles[-1])

    # Warnings
    warnings=[]
    testWarnings = False
    if resolutionP<0.05 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The reported resolution, %5.2f \\AA, is particularly with respect "\
                        "to the local resolution distribution. It occupies the %5.2f percentile}}"%\
                        (resolution,resolutionP*100))
    report.writeWarningsAndSummary(warnings, "1.d MonoRes", secLabel)
    return prot

def monodir(project, report, label, protImportMap, protCreateMask, resolution):
    Ts = protImportMap.outputVolume.getSamplingRate()

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtMonoDir', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel=label,
                               fast=True,
                               resstep=resolution/3)
    prot.inputVolumes.set(protImportMap.outputVolume)
    prot.Mask.set(protCreateMask.outputMask)
    project.launchProtocol(prot, wait=True)

    bblCitation = \
"""\\bibitem[Vilas et~al., 2020]{Vilas2020}
Vilas, J.~L., Tagare, H.~D., Vargas, J., Carazo, J.~M., and Sorzano, C. O.~S.
  (2020).
\\newblock Measuring local-directional resolution and local anisotropy in
  cryo-{EM} maps.
\\newblock {\em Nature communications}, 11:55.
"""
    report.addCitation("Vilas2020", bblCitation)

    secLabel = "sec:monodir"
    msg = \
"""
\\subsection{Level 1.e Local and directional resolution with MonoDir}
\\label{%s}
\\textbf{Explanation}:\\\\ 
MonoDir \\cite{Vilas2020} extends the concept of local resolution to local and directional resolution by changing 
the shape of the filter applied to the input map. The directional analysis can reveal image alignment problems.

The histogram of best resolution voxels per direction (Directional Histogram 1D) shows how many voxels in the
volume has their maximum resolution in that direction. Directions are arbitrarily numbered from 1 to N. This histogram
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
""" % secLabel
    report.write(msg)
    if prot.isFailed():
        report.writeSummary("1.e MonoDir", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return prot

    # 1D Histogram
    fnHistDirMonoDir1 = os.path.join(report.getReportDir(), "histDirMonoDir1D.png")
    md = xmipp3.MetaData()
    md.read(prot._getExtraPath("hist_prefdir.xmd"))
    direction = md.getColumnValues(xmipp3.MDL_X)
    count = md.getColumnValues(xmipp3.MDL_COUNT)
    reportPlot(direction, count, 'Resolution (A)', '# of voxels', fnHistDirMonoDir1, plotType="bar")

    # Test
    randomSample = np.random.choice([x for x in range(int(np.max(direction))+1)], size=50000,
                                    p=np.array(count,dtype=np.float)/np.sum(count))
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
"""
Fig. \\ref{fig:histDirMonoDir1} shows the 1D directional histogram and Fig. \\ref{fig:histDirMonoDir2} the 2D
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

    # Warnings
    warnings=[]
    testWarnings = False
    if p<0.05 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The distribution of best resolution is not uniform in all directions. "\
                        "The associated p-value is %f.}}"%p)
    if resolution<0.8*avgDirResolution or testWarnings:
        warnings.append("{\\color{red} \\textbf{The resolution reported by the user, %5.2f\\AA, is at least 80\\%% "\
                        "smaller than the average directional resolution, %5.2f \\AA.}}" % (resolution, avgDirResolution))
    report.writeWarningsAndSummary(warnings, "1.e MonoDir", secLabel)


def reportInput(project, report, fnMap1, fnMap2, protImportMap1, protImportMap2):
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
"""%(fnMap1.replace('_','\_').replace('/','/\-'), calculateSha256(fnMap1),\
     fnMap2.replace('_','\_').replace('/','/\-'), calculateSha256(fnMap2))
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

def level1(project, report, fnMap1, fnMap2, Ts, resolution, protImportMap, protCreateMask, skipAnalysis = False):
    # Import maps
    protImportMap1 = importMap(project, "import half1", fnMap1, Ts)
    if protImportMap1.isFailed():
        raise Exception("Import map did not work")
    protImportMap2 = importMap(project, "import half2", fnMap2, Ts)
    if protImportMap2.isFailed():
        raise Exception("Import map did not work")
    reportInput(project, report, fnMap1, fnMap2, protImportMap1, protImportMap2)

    # Quality Measures
    if not skipAnalysis:
        report.writeSection('Level 1 analysis')
        globalResolution(project, report, "1.a Global", protImportMap1, protImportMap2, resolution)
        # blocres(project, report, "1.b Blocres", protImportMap, protCreateMask)
        # resmap(project, report, "1.c Resmap", protImportMap, protCreateMask)
        # monores(project, report, "1.d MonoRes", protImportMap, protCreateMask, resolution)
        # monodir(project, report, "1.e MonoDir", protImportMap, protCreateMask, resolution)

    return protImportMap1, protImportMap2