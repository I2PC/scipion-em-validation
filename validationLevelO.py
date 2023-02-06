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
import json
import numpy as np
import os

import pyworkflow.plugin as pwplugin

from validationReport import calculateSha256, reportMultiplePlots, radialPlot, reportHistogram

from resourceManager import waitOutput, sendToSlurm, waitOutputFile, waitUntilFinishes


def xlmValidation(project, report, protAtom, XLM):
    bblCitation = \
"""\\bibitem[Sinnott et~al., 2020]{Sinnott2020}
Sinnott, M., Malhotra, S., Madhusudhan, M.~S., Thalassinos, K., and Topf, M.
(2020).
\\newblock Combining information from crosslinks and monolinks in the modeling
of protein structures.
\\newblock {\em Structure}, 28:1061--1070.e3.

"""
    report.addCitation("Sinnott2020", bblCitation)

    secLabel = "sec:xlm"
    msg = \
"""
\\subsection{O.a Mass-spectroscopy}
\\label{%s}

\\textbf{Explanation}:\\\\ 
The method in \\cite{Sinnott2020} uses information from cross- and mono-links to validate the atomic model. In a 
cross-link, the linker is bound to two residues, while in a mono-link, it is bound to only one residue and it can
be thought as a measure of the residue surface exposure.\\\\
\\\\
\\textbf{Results:}\\\\
""" % (secLabel)
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('xlmtools.protocols',
                                            'ProtWLM', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="O.a XLM",
                               xlList=XLM)
    prot.pdbs.set([protAtom.outputPdb])
    sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutput(project, prot, 'crosslinkStruct_1')
    waitUntilFinishes(project, prot)

    if prot.isFailed():
        report.writeSummary("O.a XLM", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return None

    fileList = glob.glob(prot._getExtraPath("Jwalk_results/*.txt"))
    if len(fileList)==0:
        report.writeSummary("O.a XLM", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol did not produce any result.}}\\\\ \n")
        return None

    msg=\
"""The user provided the following cross-link or mono-link constraints: 
\\begin{enumerate}"""
    fh = open(XLM)
    lineNo = 0
    for line in fh.readlines():
        msg+="\\item "+line+"\n"
        lineNo+=1
    fh.close()
    msg+="\\end{enumerate}\n\n"
    Nconstraints = lineNo

    msg+="From these constraints, the program has validated the following ones:\n"
    fnResults = fileList[0]
    fh = open(fnResults)
    lineNo = 0
    validated = []
    for line in fh.readlines():
        if lineNo>0:
            tokens = line.strip().split()
            validated.append((tokens[2], tokens[3], float(tokens[4]), float(tokens[5])))
        lineNo+=1
    fh.close()
    Nvalidated=len(validated)

    if Nvalidated>0:
        msg+="\\begin{center} \n \\begin{tabular}{cccc} \n"
        msg+="\\textbf{Atom1} & \\textbf{Atom2} & \\textbf{SASD} & \\textbf{Distance (\\AA)} \\\\ \n"
        for atom1, atom2, sasd, dist in validated:
            msg+="   %s & %s & %4.1f & %4.1f \\\\ \n"%(atom1, atom2, sasd, dist)
        msg+="\\end{tabular} \n \\end{center}\n\n"
    else:
        msg+="None of them\n"

    report.write(msg)
    warnings=[]
    testWarnings = False
    if Nvalidated<0.5*Nconstraints:
        warnings.append("{\\color{red} \\textbf{Less than half of the cross/mono-links were fulfilled by the atomic "\
                        "model. Precisely, %d out of %d}}"%(Nvalidated,Nconstraints))
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if more than 50\\% of the constraints are validated by the
program.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "O.a Mass spectroscopy", secLabel)
    if len(warnings)>0:
        report.writeAbstract("It seems that many of the mass spectroscopy constraints are not met (see Sec. \\ref{%s}). "%secLabel)

def saxsValidation(project, report, protMap, protMask, SAXS):
    bblCitation = \
"""\\bibitem[Jim\\'enez et~al., 2019]{Jimenez2019}
Jim\\'enez, A., Jonic, S., Majtner, T., O\\'on, J., Vilas, J.~L., Maluenda, D.,
  Mota, J., Ram\\'irez-Aportela, E., Mart\\'inez, M., Rancel, Y., Segura, J.,
  S\\'anchez-Garc\\'ia, R., Melero, R., {Del Ca\~no}, L., Conesa, P., Skjaerven,
  L., Marabini, R., Carazo, J.~M., and Sorzano, C. O.~S. (2019).
\\newblock Validation of electron microscopy initial models via small angle
  {X}-ray scattering curves.
\\newblock {\em Bioinformatics}, 35:2427--2433.
"""
    report.addCitation("Jimenez2019", bblCitation)

    secLabel = "sec:saxs"
    msg = \
"""
\\subsection{O.b SAXS}
\\label{%s}
SAXS file: %s \\\\
SHA256 hash: %s \\\\ 

\\textbf{Explanation}:\\\\ 
The method in \\cite{Jimenez2019} compares the expected energy profile from the reconstructed map to the one 
obtained by a SAXS experiment. \\\\
\\\\
\\textbf{Results:}\\\\
""" % (secLabel, SAXS.replace('_', '\_').replace('/', '/\-'), calculateSha256(SAXS))
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('continuousflex.protocols',
                                            'FlexProtConvertToPseudoAtoms', doRaise=True)
    protPseudo = project.newProtocol(Prot,
                                     objLabel="O.b Convert Map to Pseudo",
                                     maskMode=2,
                                     pseudoAtomRadius=1.5)
    protPseudo.inputStructure.set(protMap.outputVolume)
    protPseudo.volumeMask.set(protMask.outputMask)
    sendToSlurm(protPseudo)
    project.launchProtocol(protPseudo)
    #waitOutput(project, protPseudo, 'outputVolume')
    #waitOutput(project, protPseudo, 'outputPdb')
    waitUntilFinishes(project, protPseudo)
    if protPseudo.isFailed():
        report.writeSummary("O.b SAXS", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return None

    Prot = pwplugin.Domain.importFromPlugin('atsas.protocols',
                                            'AtsasProtConvertPdbToSAXS', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="O.b SAXS",
                               experimentalSAXS=SAXS)
    prot.inputStructure.set(protPseudo.outputPdb)
    sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitOutputFile(project, prot, "crysol_summary.txt")
    waitUntilFinishes(project, prot)
    if prot.isFailed():
        report.writeSummary("O.b SAXS", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return None

    fnSummary = prot._getExtraPath("crysol_summary.txt")
    fh = open(fnSummary)
    for line in fh.readlines():
        tokens = line.strip().split()
        if len(tokens)>0 and tokens[0]=="Model:":
            Rg = float(tokens[3])
            chi2 = float(tokens[7])

            fnSaxs = os.path.join(report.getReportDir(),"saxs.png")
            fnResults = prot._getPath("pseudoatoms00.fit")
            X = np.loadtxt(fnResults, skiprows=1)
            reportMultiplePlots(X[:,0],[np.log10(X[:,1]), np.log10(X[:,3])], 'Frequency (A^-1)', 'log10(SAXS)',
                                fnSaxs,['Simulated', 'Experimental'])

            msg=\
"""The radius of gyration was %5.1f \\AA. The $\\chi^2$ between the simulated curve and the experimental one was %4.1f.
Fig. \\ref{fig:saxs} shows the two SAXS profiles for comparison.

\\begin{figure}[H]
  \\centering
  \\includegraphics[width=10cm]{%s}
  \\caption{Simulated and experimental SAXS curves.}
  \\label{fig:saxs}
\\end{figure}

"""%(Rg, chi2, fnSaxs)
            report.write(msg)
            break
    fh.close()
    warnings=None
    report.writeWarningsAndSummary(warnings, "O.b SAXS", secLabel)

def tiltPairValidation(project, report, protMap, UNTILTEDMIC, TILTEDMIC, TILTKV, TILTCS, TILTQ0, TILTTS, TILTANGLE,
                       UNTILTEDCOORDS, TILTEDCOORDS, SYM):
    bblCitation = \
"""\\bibitem[Henderson et~al., 2011]{Henderson2011}
Henderson, R., Chen, S., Chen, J.~Z., Grigorieff, N., Passmore, L.~A.,
  Ciccarelli, L., Rubinstein, J.~L., Crowther, R.~A., Stewart, P.~L., and
  Rosenthal, P.~B. (2011).
\\newblock Tilt-pair analysis of images from a range of different specimens in
  single-particle electron cryomicroscopy.
\\newblock {\em J. {M}olecular {B}iology}, 413(5):1028--1046.
"""
    report.addCitation("Henderson2011", bblCitation)

    secLabel = "sec:tiltpair"
    msg = \
"""
\\subsection{O.c Tilt pair}
\\label{%s}
\\textbf{Explanation}:\\\\ 
this method is capable of experimentally validating the hand of the reconstructed map by comparing the angular 
assignment of two sets of particles related by a single-axis tilt \\cite{Henderson2011}.\\\\
\\\\
\\textbf{Results:}\\\\
""" % (secLabel)
    report.write(msg)

    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportMicrographsTiltPairs', doRaise=True)
    protImport = project.newProtocol(Prot,
                               objLabel="O.c Import tilt pairs",
                               patternUntilted=UNTILTEDMIC,
                               patternTilted=TILTEDMIC,
                               voltage=TILTKV,
                               ampContrast=TILTQ0,
                               sphericalAberration=TILTCS,
                               samplingRate=TILTTS)
    sendToSlurm(protImport)
    project.launchProtocol(protImport)
    #waitOutput(project, protImport, 'outputMicrographsTiltPair')
    waitUntilFinishes(project, protImport)
    if protImport.isFailed():
        report.writeSummary("O.c Tilt pair", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return None

    Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportCoordinatesPairs', doRaise=True)
    x,y,z = protMap.outputVolume.getDimensions()
    Ts = protMap.outputVolume.getSamplingRate()
    dMap = x*Ts
    boxSize = int(dMap/TILTTS)

    protCoords = project.newProtocol(Prot,
                                     objLabel="O.c Import paired coordinates",
                                     patternUntilted=UNTILTEDCOORDS,
                                     patternTilted=TILTEDCOORDS,
                                     boxSize=boxSize)
    if UNTILTEDCOORDS.endswith('.json'):
        protCoords.importFrom.set(1)
    protCoords.inputMicrographsTiltedPair.set(protImport.outputMicrographsTiltPair)
    sendToSlurm(protCoords)
    project.launchProtocol(protCoords)
    #waitOutput(project, protCoords, 'outputCoordinatesTiltPair')
    waitUntilFinishes(project, protCoords)
    if protCoords.isFailed():
        report.writeSummary("O.c Tilt pair", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return None

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtExtractParticlesPairs', doRaise=True)
    protExtract = project.newProtocol(Prot,
                                      objLabel="O.c Extract pairs",
                                      boxSize=boxSize,
                                      doInvert=True)
    protExtract.inputCoordinatesTiltedPairs.set(protCoords.outputCoordinatesTiltPair)
    sendToSlurm(protExtract)
    project.launchProtocol(protExtract)
    #waitOutput(project, protExtract, 'outputParticlesTiltPair')
    waitUntilFinishes(project, protExtract)
    if protExtract.isFailed():
        report.writeSummary("O.c Tilt pair", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return None

    Prot = pwplugin.Domain.importFromPlugin('xmipp3.protocols',
                                            'XmippProtCropResizeVolumes', doRaise=True)
    protResize = project.newProtocol(Prot,
                                      objLabel="Resize and resample Pair",
                                      doResize=True,
                                      resizeSamplingRate=TILTTS,
                                      doWindow=True,
                                      windowOperation=1,
                                      windowSize=boxSize)
    protResize.inputVolumes.set(protMap.outputVolume)
    sendToSlurm(protResize)
    project.launchProtocol(protResize)
    #waitOutput(project, protResize, 'outputVol')
    waitUntilFinishes(project, protResize)

    Prot = pwplugin.Domain.importFromPlugin('eman2.protocols',
                                            'EmanProtTiltValidate', doRaise=True)
    prot = project.newProtocol(Prot,
                               objLabel="O.c Tilt pair validation",
                               symmetry=SYM,
                               maxtilt=60,
                               delta=5,
                               doContourPlot=True)
    if SYM=="o":
        prot.symmetry.set("oct")
    elif SYM=="i1":
        prot.symmetry.set("icos")
    prot.inputVolume.set(protResize.outputVol)
    prot.inputTiltPair.set(protExtract.outputParticlesTiltPair)
    sendToSlurm(prot)
    project.launchProtocol(prot)
    #waitUntilFinishes(project, prot)
    waitUntilFinishes(project, prot)
    if prot.isFailed():
        report.writeSummary("O.c Tilt pair", secLabel, "{\\color{red} Could not be measured}")
        report.write("{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n")
        return None
    fnAngles = os.path.join(project.getPath(),prot._getExtraPath("TiltValidate_01/perparticletilts.json"))
    with open(fnAngles) as jsonFile:
        jsonDict = json.load(jsonFile)
    tiltpairs = jsonDict["particletilt_list"]

    r = []
    theta = []
    for tp in tiltpairs:
        r.append(tp[1])
        theta.append(tp[2])
    fnTiltPair = os.path.join(report.getReportDir(),"tiltpair.png")
    radialPlot(theta, r, [1]*len(theta), fnTiltPair)

    fnTiltHist = os.path.join(report.getReportDir(),"tiltpairHist.png")
    reportHistogram(theta, "Tilt angle", fnTiltHist)

    meanTilt = np.mean(r)
    diff = abs(meanTilt-TILTANGLE)

    msg=\
"""The average tilt of the dataset was %4.1f degrees. The experimental tilt angle of the validation was %4.1f. The
difference is %4.1f. Fig. \\ref{fig:tiltpair} shows the radial plot of the angular assignment and the histogram
of the tilts.

\\begin{figure}[H]
  \\centering
  \\includegraphics[width=10cm]{%s}\\\\
  \\includegraphics[width=10cm]{%s}\\\\
  \\caption{Radial plot of the angular assignment and the histogram of the tilt angles.}
  \\label{fig:tiltpair}
\\end{figure}
"""%(meanTilt, TILTANGLE, diff, fnTiltPair, fnTiltHist)
    report.write(msg)
    warnings=[]
    testWarnings = False
    if diff>15 or testWarnings:
        warnings.append("{\\color{red} \\textbf{The difference between the estimated mean tilt angle and "\
                        "the experimental one is larger than 15 degrees}}")
    msg = \
"""\\textbf{Automatic criteria}: The validation is OK if the difference between the estimated mean tilt angle and
the experimental one is smaller than 15 degrees.
\\\\

"""
    report.write(msg)
    report.writeWarningsAndSummary(warnings, "O.c Tilt pair", secLabel)
    if len(warnings)>0:
        report.writeAbstract("It seems that the map does not fulfill the tilt pair validation (see Sec. \\ref{%s}). "%\
                             secLabel)

def levelO(project, report, protMap, protMask, protAtom, XLM, SAXS,
           UNTILTEDMIC, TILTEDMIC, TILTKV, TILTCS, TILTQ0, TILTTS, TILTANGLE, UNTILTEDCOORDS, TILTEDCOORDS, SYM,
           skipAnalysis=False):
    checkXlm = XLM is not None and protAtom is not None
    checkSaxs = SAXS is not None
    checkPair = UNTILTEDMIC is not None

    if not skipAnalysis and (checkXlm or checkSaxs or checkPair):
        msg = "\\section{Other experimental techniques}\n\n"
        report.write(msg)

        if checkXlm: 
            xlmValidation(project, report, protAtom, XLM)
        if checkSaxs:
            saxsValidation(project, report, protMap, protMask, SAXS)
        if checkPair:
            tiltPairValidation(project, report, protMap, UNTILTEDMIC, TILTEDMIC, TILTKV, TILTCS, TILTQ0, TILTTS, TILTANGLE,
                               UNTILTEDCOORDS, TILTEDCOORDS, SYM)
