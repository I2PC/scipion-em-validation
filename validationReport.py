import numpy as np
import os
import subprocess
import PIL

from pyworkflow.utils.path import makePath
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

class ValidationReport:

    def __init__(self, fnDir):
        self.fnProjectDir = fnDir
        self.fnReportDir = os.path.join(fnDir,"validationReport")
        makePath(self.fnReportDir)
        self.fnReport = os.path.join(self.fnReportDir,"report.tex")
        self.fh = open(self.fnReport,"w")
        self.writePreamble()

    def writePreamble(self):
        toWrite = \
"""
\\documentclass[12pt, letterpaper]{article}
\\usepackage[utf8]{inputenc}
\\usepackage{graphicx}
\\usepackage{float}
\\usepackage{subfig}
\\usepackage[us,12hr]{datetime}

\\title{Validation report}
\\author{I$^2$PC Validation server}
\\date{\\today \\\\ \\currenttime}

\\begin{document}

\\begin{titlepage}
\\maketitle
\\end{titlepage}

"""
        self.fh.write(toWrite)

    def write(self,msg):
        self.fh.write(msg)

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

    def writeSubsection(self, section, msg):
        toWrite = \
"""
\\subsection{%s}

%s

""" %(section,msg)
        self.fh.write(toWrite)

    def orthogonalProjections(self, subsection, msg, caption, fnMap):
        V = readMap(fnMap).getData()

        fnRoot = subsection.replace(' ','_')
        fnZ = os.path.join(self.fnReportDir,fnRoot+"_Z.jpg")
        writeImage(np.sum(V,axis=2), fnZ)
        fnY = os.path.join(self.fnReportDir,fnRoot+"_Y.jpg")
        writeImage(np.sum(V,axis=1), fnY)
        fnX = os.path.join(self.fnReportDir,fnRoot+"_X.jpg")
        writeImage(np.sum(V,axis=0), fnX)

        toWrite="""
\\subsection{%s}

%s

\\begin{figure}[H]
  \\centering
  \\subfloat[X Projection]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[Y Projection]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[Z Projection]{\includegraphics[width=4cm]{%s}}
  \\caption{%s}
\\end{figure}

"""%(subsection, msg, fnX, fnY, fnZ, caption)
        self.fh.write(toWrite)

    def orthogonalSlices(self, subsection, msg, caption, map, maxVar=False, fnRoot=""):
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

        if fnRoot=="":
            fnRoot = subsection.replace(' ', '_')
        fnZ = os.path.join(self.fnReportDir, fnRoot + "_Z.jpg")
        writeImage(mV[:, :, iz], fnZ)
        fnY = os.path.join(self.fnReportDir, fnRoot + "_Y.jpg")
        writeImage(mV[:, iy, :], fnY)
        fnX = os.path.join(self.fnReportDir, fnRoot + "_X.jpg")
        writeImage(mV[ix, :, :], fnX)

        toWrite = ""
        if subsection!="":
            toWrite+=\
"""

\\subsection{%s}

%s
"""%(subsection, msg)

        toWrite+=\
"""\\begin{figure}[H]
  \\centering
  \\subfloat[X Slice %d]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[Y Slice %d]{\includegraphics[width=4cm]{%s}}
  \\hspace{0.2cm}
  \\subfloat[Z Slice %d]{\includegraphics[width=4cm]{%s}}
  \\caption{%s}
\\end{figure}

""" % (ix, fnX, iy, fnY, iz, fnZ, caption)
        self.fh.write(toWrite)

    def closeReport(self):
        toWrite = \
"""
\\end{document}
"""
        self.fh.write(toWrite)
        self.fh.close()

        os.chdir(self.fnReportDir)
        subprocess.run(["pdflatex","-interaction=nonstopmode","report.tex"],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.STDOUT)
        os.chdir(self.fnProjectDir)
