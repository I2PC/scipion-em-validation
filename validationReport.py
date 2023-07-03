import hashlib
from math import radians
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import subprocess
import PIL

from pyworkflow.protocol import StringParam
from pyworkflow.utils.path import makePath, cleanPath
from pwem.viewers import LocalResolutionViewer
from pwem.emlib.metadata import iterRows

import xmipp3

def readMap(fnMap):
    return xmipp3.Image(fnMap)

def readStack(fnXmd):
    md = xmipp3.MetaData(fnXmd)
    S = None
    i=0
    for objId in md:
        I = xmipp3.Image(md.getValue(xmipp3.MDL_IMAGE, objId))
        Xdim, Ydim, Zdim, Ndim = I.getDimensions()
        if S is None:
            S=np.zeros((md.size(), Ydim, Xdim))
        S[i,:,:]=I.getData()
        i+=1
    return S

def readGuinier(fnGuinier):
    fh = open(fnGuinier)
    lineNo = 0
    content = []
    for line in fh.readlines():
        if lineNo > 0:
            content.append([float(x) for x in line.split()])
        lineNo += 1
    X = np.array(content)
    fh.close()

    dinv2 = X[:, 0]
    lnF = X[:, 1]
    lnFc = X[:, 3]
    return dinv2, lnF, lnFc

def writeImage(I, fnOut, scale=True):
    if scale:
        I = I.astype(np.float)
        I -= np.min(I)
        I *= 255/np.max(I)
    Iout = PIL.Image.fromarray(I.astype(np.uint8))
    Iout.save(fnOut)

def latexItemize(itemList):
    toWrite="\\begin{itemize}\n"
    for item in itemList:
        toWrite+="   \item %s\n"%item
    toWrite+="\\end{itemize}\n\n"
    return toWrite

def latexEnumerate(itemList):
    toWrite="\\begin{enumerate}\n"
    for item in itemList:
        toWrite+="   \item %s\n"%item
    toWrite+="\\end{enumerate}\n\n"
    return toWrite

def generateChimeraView(fnWorkingDir, fnMap, fnView, isMap=True, threshold=0, angX=0, angY=0, angZ=0, bfactor=False,\
                        occupancy=False):
    chimeraScript=\
"""
set bgColor white
open %s
"""%fnMap
    if isMap:
        chimeraScript+=\
"""volume #1 level %f
volume #1 color #4e9a06
lighting soft
"""%threshold
    else:
        chimeraScript+=\
"""hide atoms
show cartoons
"""
        if bfactor:
            chimeraScript+="color bfactor\n"
        if occupancy:
            chimeraScript += "color byattribute occupancy palette rainbow\n"
    chimeraScript+=\
"""turn x %f
turn y %f
turn z %f
save %s
exit
"""%(angX, angY, angZ, fnView)
    fnTmp = os.path.join(fnWorkingDir, "chimeraScript.cxc")
    fh = open(fnTmp,"w")
    fh.write(chimeraScript)
    fh.close()

    from chimera import Plugin
    args = "chimeraScript.cxc"
    Plugin.runChimeraProgram(Plugin.getProgram(), args, cwd=fnWorkingDir)
    cleanPath(fnTmp)

def generateChimeraColorView(fnWorkingDir, project, fnRoot, fnMap, Ts, fnColor, minVal, maxVal):
    viewer = LocalResolutionViewer(project=project)
    viewer.colorMap = StringParam()
    viewer.colorMap.set("jet")
    cmdFile = os.path.join(fnWorkingDir, fnRoot + ".py")
    viewer.createChimeraScript(cmdFile, fnColor, fnMap, Ts,
                               numColors=11,
                               lowResLimit=maxVal,
                               highResLimit=minVal)

    fn1 = os.path.join(fnWorkingDir, fnRoot + "1.jpg")
    fn2 = os.path.join(fnWorkingDir, fnRoot + "2.jpg")
    fn3 = os.path.join(fnWorkingDir, fnRoot + "3.jpg")
    fhCmd = open(cmdFile, "a")
    toWrite = \
"""
run(session, 'save %s')
run(session, 'turn x 90')
run(session, 'save %s')
run(session, 'turn x -90')
run(session, 'turn y 90')
run(session, 'save %s')
run(session, 'exit')
""" % (fn1, fn2, fn3)
    fhCmd.write(toWrite)
    fhCmd.close()

    from chimera import Plugin
    args = "--script %s"%cmdFile
    Plugin.runChimeraProgram(Plugin.getProgram(), args, cwd=fnWorkingDir)

def formatInv(value, pos):
    """ Format function for Matplotlib formatter. """
    inv = 999.
    if value:
        inv = 1 / value
    return "1/%0.2f" % inv

def reportPlot(x,y, xlabel, ylabel, fnOut, yscale="linear",grid=True, plotType="plot", barWidth=1,
               invertXLabels=False, addMean=False, title=""):
    matplotlib.use('Agg')
    plt.figure()
    if plotType=="plot":
        plt.plot(x,y)
        if addMean:
            ymean = np.mean(y)
            plt.plot(x,ymean*np.ones(len(y)),'--')
    elif plotType=="bar":
        plt.bar(x,y,barWidth)
    elif plotType=="scatter":
        plt.scatter(x,y)
    plt.yscale(yscale)
    if invertXLabels:
        plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(formatInv))
    plt.grid(grid)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if title!="":
        plt.title(title)
    plt.savefig(fnOut, bbox_inches='tight')
    plt.close('all')

def reportHistogram(y, ylabel, fnOut):
    matplotlib.use('Agg')
    plt.figure()
    plt.hist(y, bins=25)
    plt.grid(True)
    plt.xlabel(ylabel)
    plt.ylabel("Count")
    plt.savefig(fnOut, bbox_inches='tight')
    plt.close('all')

def reportMultiplePlots(x,yList, xlabel, ylabel, fnOut, legends, invertXLabels=False, xshade0=None, xshadeF=None):
    matplotlib.use('Agg')
    plt.figure()
    for i in range(len(yList)):
        plt.plot(x,yList[i], label=legends[i])
    plt.legend()
    if invertXLabels:
        plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(formatInv))
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if xshade0 is not None and xshadeF is not None:
        plt.gca().axes.axvspan(xshade0, xshadeF, alpha=0.3, color='green')
    plt.savefig(fnOut, bbox_inches='tight')
    plt.close('all')

def radialPlot(thetas, radii, weights, fnOut, plotType="points"):
    max_w = max(weights)
    min_w = min(weights)
    min_p = 5
    max_p = 40

    matplotlib.use('Agg')
    plt.figure()

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    if plotType=="points":
        for theta, radius, w in zip(thetas, radii, weights):
            pointsize = int((w - min_w) / (max_w - min_w + 0.001) * (max_p - min_p) + min_p)
            ax.plot(radians(theta), radius, markerfacecolor='blue', marker='.', markersize=pointsize)
    elif plotType=="contour":
        thetasp = np.radians(np.linspace(0, 360, 360))
        radiip = np.arange(0, np.max(radii)+1, 1)
        rmesh, thetamesh = np.meshgrid(radiip, thetasp)
        values = np.zeros((len(thetasp), len(radiip)))

        for i in range(0, len(thetas)):
            values[int(thetas[i]), int(radii[i])] = weights[i]
        stp = 0.1
        lowlim = max(0.0, values.min())
        highlim = values.max() + stp
        pc = plt.contourf(thetamesh, rmesh, values, np.arange(lowlim, highlim, stp))
        plt.colorbar(pc)
    ax.set_rmax(np.max(radii))
    ax.set_rticks([15, 30, 45, 60, 75, 90])  # Less radial ticks
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.grid(True)

    plt.savefig(fnOut, bbox_inches='tight')
    plt.close('all')

def calculateSha256(fn):
    sha256_hash = hashlib.sha256()
    with open(fn, "rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def CDFFromHistogram(x,count):
    xWidth=x[1]-x[0]
    newX = np.array(x+[x[-1]+xWidth])
    newX -= xWidth/2
    newY = np.array([0]+count)
    F = np.cumsum(newY)
    F = F/F[-1]
    return newX, F

def CDFpercentile(x,F,xp=None,Fp=None):
    if xp is not None:
        f = scipy.interpolate.interp1d(x, F, bounds_error=False, fill_value=(0.0,1.0))
        if type(xp)==list:
            return [f(xpp) for xpp in xp]
        else:
            return f(xp)
    if Fp is not None:
        f = scipy.interpolate.interp1d(F, x, bounds_error=False, fill_value=(np.min(x),np.max(x)))
        if type(Fp)==list:
            return [f(Fpp) for Fpp in Fp]
        return f(Fp)
    return None

def plotMicrograph(fnMic, fnOut, coords=None, boxSize=0, Ts=0):
    I = xmipp3.Image()
    I.read(fnMic)
    Xdim, Ydim, _, _ = I.getDimensions()
    I.readPreviewSmooth(fnMic, 800)
    Xdimp, Ydimp, _, _ = I.getDimensions()
    boxSizep = int(boxSize*Ts/float(Xdim)*Xdimp)
    matplotlib.use('Agg')
    plt.figure()
    plt.imshow(I.getData(), cmap='gray')
    plt.axis('off')
    Kx = float(Xdimp)/Xdim
    Ky = float(Ydimp)/Ydim
    if coords is not None:
        for x,y in coords:
            xp = int(x*Kx)
            yp = int(y*Ky)
            circle = plt.Circle((xp, yp), boxSizep, color='green', fill=False)
            plt.gca().add_artist(circle)
    plt.savefig(fnOut, bbox_inches='tight')
    plt.close('all')

class ValidationReport:

    def __init__(self, fnDir, levels):
        self.fnProjectDir = fnDir
        self.fnReportDir = os.path.join(fnDir,"validationReport")
        makePath(self.fnReportDir)
        self.fnReport = os.path.join(self.fnReportDir,"report.tex")
        self.fh = open(self.fnReport,"w")
        self.fnContext = os.path.join(self.fnReportDir, "context.tex")
        self.fnContext = open(self.fnContext, "w")
        self.fnContext.write("\\begin{center}\\textbf{Context}\\end{center}  \n"\
                             "Cryo-electron microscopy is currently one of the most active techniques in Structural Biology. The number of maps deposited at the \\href{https://www.ebi.ac.uk/emdb/}{Electron Microscopy Data Bank} is rapidly growing every year and keeping the quality of the submitted maps is essential to maintain the scientific quality of the field. \\\\ \n" \
                             "The ultimate quality measure is the consistency of the map and an atomic model. However, this is only possible for high resolution maps. Over the years there have been many suggestions about validation measures of 3DEM maps. Unfortunately, most of these measures are not currently in use for their spread in multiple software tools and the associated difficulty to access them. To alleviate this problem, we made available a validation grading system that evaluate the information provided to assess the map. \\\\ \n" \
                             "This system grades a map from 0 to 5 depending on the amount of information available. In this way, a map could be validated at Level 0 (the deposited map), 1 (two half maps), 2 (2D classes), 3 (particles), 4 (... + angular assignment), 5 (... + micrographs and coordinates). In addition, we can have three optional qualifiers: A (... + atomic model), W (... + image processing workflow), and O (... + other techniques). \\\\ \n\n" \
                             "This Validation Report Service is explained in more detail in the paper:\\\\ \n\n" \
                             "C.O.S. Sorzano, J.L. Vilas, E. Ramírez-Aportela, D. del Hoyo, D. Herreros, E. Fernandez-Giménez, D. Marchán, F. de Isidro Gómez, J.R. Macías, I. Sánchez, L. del Caño, Y. Fonseca-Reyna, P. Conesa, A. García-Mena, J. Burguet, J. García Condado, J. Méndez García, M. Martínez, A. Muñoz Barrutia, R. Marabini, J. Vargas, J.M. Carazo. \\textbf{Image processing tools for the validation of CryoEM maps.} Faraday Discuss., 2022,240, 210-227.\n" \
                             "Open access (peer-reviewed journals) \\href{http://creativecommons.org/licenses/by-nc/3.0/}{Creative Commons BY-NC license} DOI: \\href{https://doi.org/10.1039/D2FD00059H}{10.1039/D2FD00059H}")
        self.fnContext.close()
        self.fnAbstract = os.path.join(self.fnReportDir,"abstract.tex")
        self.fhAbstract = open(self.fnAbstract,"w")
        self.fnSummary = os.path.join(self.fnReportDir,"summary.tex")
        self.fhSummary = open(self.fnSummary,"w")
        self.fnSummaryWarnings = os.path.join(self.fnReportDir,"summaryWarnings.tex")
        self.fhSummaryWarnings = open(self.fnSummaryWarnings,"w")
        self.fhSummaryWarnings.write("\\textbf{\\underline{Summary of the warnings across sections.}}\\\\ \n\n\n")
        self.citations = {}
        self.writePreamble(levels)
        self.resolutionEstimates = []
        self.score = 0
        self.scoreN = 0

        self.abstract = ""

    def getReportDir(self):
        return self.fnReportDir

    def writePreamble(self, levels):
        toWrite = \
"""
\\documentclass[12pt, letterpaper]{article}
\\usepackage[utf8]{inputenc}
\\usepackage{graphicx}
\\usepackage{float}
\\usepackage{subfig}
\\usepackage{hyperref}
\\usepackage{xcolor}
\\usepackage[us,12hr]{datetime}
\\usepackage{longtable}
\\usepackage{enumitem}
\\usepackage{tikz}
\\usepackage{fancyhdr}
\\setlist{nosep}

% Define own colors
\\definecolor{mygreen}{RGB}{116,183,46}
\\definecolor{myblue}{RGB}{37,150,190} 

% Set hyperlinks configuration
\\hypersetup{
    colorlinks=true,
    linkcolor=red,
    filecolor=magenta,      
    urlcolor=blue,
    pdftitle={Validation Report Service},
    pdfpagemode=FullScreen,
    citecolor=mygreen,
}

\\begin{document}

\\input{frontpage.tex}

\\input{context.tex}
\\clearpage

\\renewcommand{\\abstractname}{Summarized overall quality}
\\begin{abstract}
\\input{abstract.tex}
\\clearpage
\\input{summary.tex}
\\clearpage
\\input{summaryWarnings.tex}
\\end{abstract}

\\clearpage

\\tableofcontents

\\clearpage

"""%", ".join(levels)
        self.fh.write(toWrite)

        toWrite = \
"""
\\begin{tabular}{lrr}
"""
        self.fhSummary.write(toWrite)

    def write(self,msg):
        self.fh.write(msg)

    def writeWarningsAndSummary(self, warnings, section, secLabel):
        if warnings is None:
            toWrite = "\\textbf{STATUS}: {\\color{brown} Cannot be automatically evaluated}\\\\ \n"
            self.writeSummary(section, secLabel, "{\\color{brown} Cannot be automated}")
        else:
            self.scoreN+=1
            if len(warnings) > 0:
                countWarnings = len(warnings)
                toWrite = "\\textbf{WARNINGS}: %d warnings\\\\ \n%s\n" % (countWarnings, latexEnumerate(warnings))
                self.writeSummary(section, secLabel, "{\\color{red} %d warnings}" % countWarnings)

                self.fhSummaryWarnings.write("\\underline{Section \\ref{%s} (%s)}\n" % (secLabel, section))
                self.fhSummaryWarnings.write("\\begin{enumerate}\n")
                for warning in warnings:
                    self.fhSummaryWarnings.write("\\item %s\n"%warning)
                self.fhSummaryWarnings.write("\\end{enumerate}\n\n\n")
            else:
                toWrite = "\\textbf{STATUS}: {\\color{blue} OK}\\\\ \n"
                self.writeSummary(section, secLabel, "{\\color{blue} OK}")
                self.score += 1
        self.write(toWrite)

    def writeSummary(self, test, sec, msg):
        toWrite = "%s & Sec. \\ref{%s} & %s \\\\ \n"%(test,sec,msg)
        self.fhSummary.write(toWrite)

    def writeAbstract(self, msg):
        self.fhAbstract.write(msg)

    def addResolutionEstimate(self, R):
        self.resolutionEstimates.append(R)

    def abstractResolution(self, resolution):
        if len(self.resolutionEstimates)>0:
            msg="\n\n\\vspace{0.5cm}The average resolution of the map estimated by various methods goes from %4.1f\\AA~to %4.1f\\AA~ with an "\
                "average of %4.1f\\AA. The resolution provided by the user was %4.1f\\AA."%\
                (np.min(self.resolutionEstimates), np.max(self.resolutionEstimates), np.mean(self.resolutionEstimates),
                 resolution)
            if resolution<0.8*np.mean(self.resolutionEstimates):
                msg+=" The resolution reported by the user may be overestimated."
            msg+="\n\n\\vspace{0.5cm}"
            self.writeAbstract(msg)

    def addCitation(self, key, bblEntry):
        # Bibliographystyle: apalike
        # \begin{thebibliography}{}
        #
        # \bibitem[Greenwade, 1993] {greenwade93}
        # Greenwade, G. ~D.(1993).
        # \newblock The {C}omprehensive {T}ex {A}rchive {N}etwork({CTAN}).
        # \newblock {\em TUGBoat}, 14(3): 342 - -351.
        #
        # \end{thebibliography}
        if not key in self.citations:
            self.citations[key] = bblEntry

    def writeFailedSubsection(self, subsection, msg):
        toWrite=\
"""
\\subsection{%s}

%s

\\textbf{Failed, it could not be computed}

"""%(subsection, msg)
        self.fh.write(toWrite)

    def writeSection(self, section):
        toWrite=\
"""
\\section{%s}

"""%section
        self.fh.write(toWrite)

    def orthogonalProjections(self, fnRoot, msg, caption, fnMap, label):
        V = readMap(fnMap).getData()

        fnZ = os.path.join(self.fnReportDir,fnRoot+"_Z.jpg")
        writeImage(np.sum(V,axis=2), fnZ)
        fnY = os.path.join(self.fnReportDir,fnRoot+"_Y.jpg")
        writeImage(np.sum(V,axis=1), fnY)
        fnX = os.path.join(self.fnReportDir,fnRoot+"_X.jpg")
        writeImage(np.sum(V,axis=0), fnX)

        toWrite="""
%s

\\begin{figure}[H]
  \\centering
  \\subfloat[X Projection]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[Y Projection]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[Z Projection]{\includegraphics[width=4cm]{%s}}
  \\caption{%s}
  \\label{%s}
\\end{figure}

"""%(msg, fnX, fnY, fnZ, caption, label)
        self.fh.write(toWrite)

    def isoSurfaces(self, fnRoot, msg, caption, fnMap, threshold, label):
        generateChimeraView(self.fnReportDir, fnMap, fnRoot + "1.jpg", True, threshold, 0, 0, 0)
        generateChimeraView(self.fnReportDir, fnMap, fnRoot + "2.jpg", True, threshold, 90, 0, 0)
        generateChimeraView(self.fnReportDir, fnMap, fnRoot + "3.jpg", True, threshold, 0, 90, 0)

        caption +=  " Views generated by ChimeraX at a the following "\
                    "X, Y, Z angles: View 1 (0,0,0), View 2 (90, 0, 0), View 3 (0, 90, 0)."

        fn1 = os.path.join(self.fnReportDir,fnRoot + "1.jpg")
        fn2 = os.path.join(self.fnReportDir,fnRoot + "2.jpg")
        fn3 = os.path.join(self.fnReportDir,fnRoot + "3.jpg")
        toWrite = \
"""
%s

\\begin{figure}[H]
  \\centering
  \\subfloat[View 1]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[View 2]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[View 3]{\includegraphics[width=4cm]{%s}}
  \\caption{%s}
  \\label{%s}
\\end{figure}

""" % (msg, fn1, fn2, fn3, caption, label)
        self.fh.write(toWrite)

    def colorIsoSurfaces(self, msg, caption, label, project, fnRoot, fnMap, Ts, fnColor, minVal, maxVal):
        generateChimeraColorView(self.getReportDir(), project, fnRoot, fnMap, Ts, fnColor, minVal, maxVal)

        caption +=  " Views generated by ChimeraX at a the following "\
                    "X, Y, Z angles: View 1 (0,0,0), View 2 (90, 0, 0), View 3 (0, 90, 0)."

        fn1 = os.path.join(self.fnReportDir, fnRoot + "1.jpg")
        fn2 = os.path.join(self.fnReportDir, fnRoot + "2.jpg")
        fn3 = os.path.join(self.fnReportDir, fnRoot + "3.jpg")
        toWrite = \
"""
%s

\\begin{figure}[H]
  \\centering
  \\subfloat[View 1]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[View 2]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[View 3]{\includegraphics[width=4cm]{%s}}
  \\caption{%s}
  \\label{%s}
\\end{figure}

""" % (msg, fn1, fn2, fn3, caption, label)
        self.fh.write(toWrite)

    def atomicModel(self, fnRoot, msg, caption, fnModel, label, bfactor=False, occupancy=False):
        generateChimeraView(self.fnReportDir, fnModel, fnRoot + "1.jpg", False, 0, 0, 0, 0, bfactor=bfactor,
                            occupancy=occupancy)
        generateChimeraView(self.fnReportDir, fnModel, fnRoot + "2.jpg", False, 0, 90, 0, 0, bfactor=bfactor,
                            occupancy=occupancy)
        generateChimeraView(self.fnReportDir, fnModel, fnRoot + "3.jpg", False, 0, 0, 90, 0, bfactor=bfactor,
                            occupancy=occupancy)

        caption += " Views generated by ChimeraX at a the following " \
                   "X, Y, Z angles: View 1 (0,0,0), View 2 (90, 0, 0), View 3 (0, 90, 0)."

        fn1 = os.path.join(self.fnReportDir, fnRoot + "1.jpg")
        fn2 = os.path.join(self.fnReportDir, fnRoot + "2.jpg")
        fn3 = os.path.join(self.fnReportDir, fnRoot + "3.jpg")
        toWrite = \
            """
            %s
        
            \\begin{figure}[H]
              \\centering
              \\subfloat[View 1]{\includegraphics[width=4cm]{%s}}
              \\hspace{0.2cm}
              \\subfloat[View 2]{\includegraphics[width=4cm]{%s}}
              \\hspace{0.2cm}
              \\subfloat[View 3]{\includegraphics[width=4cm]{%s}}
              \\caption{%s}
              \\label{%s}
            \\end{figure}
        
            """ % (msg, fn1, fn2, fn3, caption, label)
        self.fh.write(toWrite)

    def orthogonalSlices(self, fnRoot, msg, caption, map, label, maxVar=False):
        if type(map) is str:
            V = readMap(map)
            Xdim, Ydim, Zdim, _ = V.getDimensions()
            mV = V.getData()
        else:
            mV = map
            Xdim, Ydim, Zdim = mV.shape

        if maxVar:
            ix = np.argmax([np.var(mV[i, :, :]) for i in range(Xdim)])
            iy = np.argmax([np.var(mV[:, i, :]) for i in range(Ydim)])
            iz = np.argmax([np.var(mV[:, :, i]) for i in range(Zdim)])
        else:
            ix = int(Xdim / 2)
            iy = int(Ydim / 2)
            iz = int(Zdim / 2)

        fnZ = os.path.join(self.fnReportDir, fnRoot + "_Z.jpg")
        writeImage(mV[:, :, iz], fnZ)
        fnY = os.path.join(self.fnReportDir, fnRoot + "_Y.jpg")
        writeImage(mV[:, iy, :], fnY)
        fnX = os.path.join(self.fnReportDir, fnRoot + "_X.jpg")
        writeImage(mV[ix, :, :], fnX)

        toWrite = msg

        toWrite+=\
"""\\begin{figure}[H]
  \\centering
  \\subfloat[X Slice %d]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[Y Slice %d]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[Z Slice %d]{\includegraphics[width=4cm]{%s}}
  \\caption{%s}
  \\label{%s}
\\end{figure}

""" % (ix, fnX, iy, fnY, iz, fnZ, caption, label)
        self.fh.write(toWrite)

    def setOfImages(self, fnImgs, mdlLabel, caption, label, fnRoot, width, Ximgs, imgMax=1e6):
        toWrite = \
"""\\begin{figure}[H]
  \\centering
"""
        self.write(toWrite)

        md = xmipp3.MetaData(fnImgs)
        imgNames = md.getColumnValues(mdlLabel)
        toWrite = ""
        for i in range(len(imgNames)):
            I = xmipp3.Image(imgNames[i])
            fnOut = "%s_%04d.jpg"%(fnRoot,i)
            I.write(fnOut)
            toWrite+=\
"""\\includegraphics[width=%s]{%s}
"""%(width,fnOut)
            if (i+1)%Ximgs==0:
                toWrite+="\\\\ \n"
            if i==imgMax:
                break
        toWrite+=\
"""\\caption{%s}
  \\label{%s}
\\end{figure}

"""%(caption, label)
        self.write(toWrite)

    def showj(self, md, mdLabelList, renderList, formatList, headers, fnRoot, width):
        toWrite = \
"""\\begin{longtable}{%s}
  \\centering
"""%('c'*len(headers))
        i = 0
        for header in headers:
            if i>0:
                toWrite+=" & "
            toWrite+="\\textbf{%s}"%header
            i+=1
        toWrite+="\\\\ \n"
        self.write(toWrite)

        toWrite = ""
        idx=0
        for row in iterRows(md):
            i=0
            for label, render, format in zip(mdLabelList, renderList, formatList):
                content = row.getValue(label)
                if i > 0:
                    toWrite += " & "
                if render:
                    I = xmipp3.Image(content)
                    fnOut = "%s_%05d.jpg" % (fnRoot, idx)
                    I.write(fnOut)
                    idx+=1
                    toWrite += "\\includegraphics[width=%s]{%s} " % (width, fnOut)
                else:
                    toWrite+=format%content
                i+=1
            toWrite+="\\\\ \n"

        toWrite += \
"""\\end{longtable}

"""
        self.write(toWrite)

    def closeReport(self):
        toWrite = "\n\n\\begin{thebibliography}{}\n\n"
        for key in self.citations:
            toWrite +="%s\n\n"%self.citations[key]
        toWrite += "\\end{thebibliography}\n\n"
        toWrite += "\\end{document}\n"
        self.fh.write(toWrite)
        self.fh.close()

        self.fhAbstract.write('\n\n')
        self.fhAbstract.write('\\textbf{The overall score (passing tests) of this report is %d out of %d evaluable '\
                              'items.}\n\n'%(self.score,self.scoreN))
        self.fhAbstract.close()

        self.fhSummaryWarnings.write('\n\n')
        self.fhSummaryWarnings.close()

        toWrite = \
"""
\\end{tabular}
"""
        self.fhSummary.write(toWrite)
        self.fhSummary.close()

        os.chdir(self.fnReportDir)
        subprocess.run(["pdflatex","-interaction=nonstopmode","report.tex"],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.STDOUT)
        subprocess.run(["pdflatex","-interaction=nonstopmode","report.tex"],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.STDOUT)
        os.chdir(self.fnProjectDir)
