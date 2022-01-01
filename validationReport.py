import hashlib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import scipy
import subprocess
import PIL

from pyworkflow.protocol import StringParam
from pyworkflow.utils.path import makePath, cleanPath
from pwem.viewers import LocalResolutionViewer

import xmipp3

def readMap(fnMap):
    return xmipp3.Image(fnMap)

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

def generateChimeraView(fnWorkingDir, fnMap, threshold, fnView, angX=0, angY=0, angZ=0):
    chimeraScript=\
"""
set bgColor white
open %s
volume #1 level %f
volume #1 color #4e9a06
turn x %f
turn y %f
turn z %f
save %s
exit
"""%(fnMap,threshold, angX, angY, angZ, fnView)
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

def reportPlot(x,y, xlabel, ylabel, fnOut, yscale="linear", grid=True, plotType="plot", barWidth=1):
    matplotlib.use('Agg')
    plt.figure()
    if plotType=="plot":
        plt.plot(x,y)
    elif plotType=="bar":
        plt.bar(x,y,barWidth)
    plt.yscale(yscale)
    plt.grid(grid)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(fnOut)

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

class ValidationReport:

    def __init__(self, fnDir, levels):
        self.fnProjectDir = fnDir
        self.fnReportDir = os.path.join(fnDir,"validationReport")
        makePath(self.fnReportDir)
        self.fnReport = os.path.join(self.fnReportDir,"report.tex")
        self.fh = open(self.fnReport,"w")
        self.fnSummary = os.path.join(self.fnReportDir,"summary.tex")
        self.fhSummary = open(self.fnSummary,"w")
        self.citations = {}
        self.writePreamble(levels)

    def getReportDir(self):
        return self.fnReportDir

    def writePreamble(self, levels):
        maxLevel = np.max(levels)

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

\\title{Validation report of Level %d}
\\author{I$^2$PC Validation server}
\\date{\\today \\\\ \\currenttime}

\\begin{document}

\\begin{titlepage}
\\maketitle
\\end{titlepage}

\\begin{abstract}
\\input{summary.tex}
\\end{abstract}

\\clearpage

\\tableofcontents

\\clearpage

"""%maxLevel
        self.fh.write(toWrite)

        toWrite = \
"""
\\begin{tabular}{lrr}
"""
        self.fhSummary.write(toWrite)

    def write(self,msg):
        self.fh.write(msg)

    def writeWarningsAndSummary(self, warnings, section, secLabel):
        if len(warnings) > 0:
            countWarnings = len(warnings)
            toWrite = "\\textbf{WARNINGS}: %d warnings\\\\ \n%s\n" % (countWarnings, latexEnumerate(warnings))
            self.writeSummary(section, secLabel, "{\\color{red} %d warnings}" % countWarnings)
        else:
            toWrite = "\\textbf{STATUS}: {\\color{blue} OK}\\\\ \n"
            self.writeSummary(section, secLabel, "{\\color{blue} OK}")
        self.write(toWrite)

    def writeSummary(self, test, sec, msg):
        toWrite = "%s & Sec. \\ref{%s} & %s \\\\ \n"%(test,sec,msg)
        self.fhSummary.write(toWrite)

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
        generateChimeraView(self.fnReportDir, fnMap, threshold, fnRoot + "1.jpg", 0, 0, 0)
        generateChimeraView(self.fnReportDir, fnMap, threshold, fnRoot + "2.jpg", 90, 0, 0)
        generateChimeraView(self.fnReportDir, fnMap, threshold, fnRoot + "3.jpg", 0, 90, 0)

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

    def closeReport(self):
        toWrite = "\n\n\\begin{thebibliography}{}\n\n"
        for key in self.citations:
            toWrite +="%s\n\n"%self.citations[key]
        toWrite += "\\end{thebibliography}\n\n"
        toWrite += "\\end{document}\n"
        self.fh.write(toWrite)
        self.fh.close()

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
