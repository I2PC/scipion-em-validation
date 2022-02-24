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

import graphviz
import json
import os
import urllib.request

from validationReport import calculateSha256

def levelW(project, report, WORKFLOW, skipAnalysis=False):
    fnWorkflow = WORKFLOW
    if WORKFLOW.startswith("http"):
        fnWorkflow = os.path.join(report.getReportDir(),"workflow.json")
        urlWorkflow = WORKFLOW
        if not urlWorkflow.endswith('.json'):
            urlWorkflow="http://nolan.cnb.csic.es/cryoemworkflowviewer/static/uploadedFiles/%s/workflow.json"%\
                        WORKFLOW.split('/')[-1]
            urllib.request.urlretrieve(urlWorkflow, fnWorkflow)

    fh = open(fnWorkflow)
    workflow = json.load(fh)
    fh.close()

    secLabel = "sec:workflow"
    if WORKFLOW.startswith("http"):
        msgWorkflow = "\\url{%s}" % WORKFLOW
    else:
        msgWorkflow = WORKFLOW.replace('_', '\_').replace('/', '/\-')
    msg = \
"""
\\section{Workflow}
\\label{%s}
Workflow file: %s \\\\
SHA256 hash: %s \\\\ 
\\\\
""" % (secLabel, msgWorkflow, calculateSha256(fnWorkflow))
    report.write(msg)

    dot = graphviz.Digraph()

    # create nodes
    for protocol in workflow:
        dot.node(protocol['object.id'], protocol['object.label'])

    # create edges
    for protocol in workflow:
        for key in protocol:
            if 'input' in key:
                source = protocol[key]
                if isinstance(source, str):
                    source = source.split('.')[0]
                    dot.edge(source, protocol['object.id'])
                if isinstance(source, list):
                    for s in source:
                        source = s.split('.')[0]
                        dot.edge(source, protocol['object.id'])

    fnGraph = os.path.join(report.getReportDir(),"workflow.png")
    dot.render(outfile=fnGraph, view=False)

    msg=\
"""Fig. \\ref{fig:workflow} shows the image processing workflow followed in Scipion to achieve these results.

\\begin{figure}[H]
    \centering
    \includegraphics[width=15cm]{%s}
    \\caption{Image processing workflow in Scipion to achieve these results.}
    \\label{fig:workflow}
\\end{figure}
"""%fnGraph

    report.write(msg)
    report.writeWarningsAndSummary(None, "W Workflow", secLabel)

