######################### METHOD NAME CONSTANTS #########################
NAME_DEEPRES = "DeepRes"
NAME_LOCBFACTOR = "LocBfactor"
NAME_LOCOCCUPANCY = "LocOccupancy"
NAME_DEEPHAND = "DeepHand"
NAME_GLOBAL_RES = "Global resolution"
NAME_FSC_PERMU = "FSC permutation"
NAME_BLOCRES = "Blocres"
NAME_MONORES = "MonoRes"
NAME_MONODIR = "MonoDir"
NAME_FSO = "FSO"
NAME_FSC3D = "FSC3D"
NAME_COMPARE_PROJECTIONS = "Compare reprojections"
NAME_OUTLIER_DETECTION = "Outlier detection"
NAME_CORE_ANALYSIS = "Core analysis"
NAME_CLASSIFICATION_EXT_CONSISTENCY = "Classification external consistency"
NAME_SIMILARITY_CRITERIA = "Similarity criteria"
NAME_ALIGN_SMOOTHNESS =  "Alignability smoothness"
NAME_ALIGNABILITY = "Alignability"
NAME_RELION_ALIGN = "Relion alignment"
NAME_CRYOSPARC_ALIGN = "CryoSparc alignment"
NAME_RELION_CRYOSPARC_ALIGN = "Relion/CryoSparc alignments"
NAME_CLASSIFICATION_WITHOUT_ALIGN = "Classification without alignment"
NAME_OVERFITTING_DETECT = "Overfitting detection"
NAME_ANGULAR_DISTR_EFFICIENCY = "Angular distribution efficiency"
NAME_CTF_STABILITY = "CTF stability"
NAME_MICROGRAPH_CLEANER = "Micrograph cleaner"
NAME_MAPQ = "MapQ"
NAME_CONVERSION_PDB_MAP = "Conversion PDB to map"
NAME_FSCQ = "FSC-Q"
NAME_MULTIMODEL = "Multimodel"
NAME_PHENIX = "Phenix"
NAME_EMRINGER = "EMRinger"
NAME_DAQ = "DAQ"
NAME_XLM = "XLM"
NAME_SAXS = "SAXS"
NAME_TILT_PAIR = "Tilt pair"

####################### TERMINAL MESSAGES CONSTANTS #######################
PRINT_PROTOCOL_ABORTED = "PROTOCOL MANUALLY ABORTED"

####################### REPORT MESSAGES CONSTANTS #######################
# Summary
OK_MESSAGE = "{\\color{blue} OK}"
WARNINGS_MESSAGE = "{\\color{red} %d warnings}"
NOT_APPLY_MESSAGE = "{\\color{brown} Does not apply}"
ERROR_MESSAGE = "{\\color{red} Could not be measured}"
ERROR_CONVERT_MESSAGE = "{\\color{red} Could not be converted}"
ERROR_ABORTED_MESSAGE = "{\\color{red} MANUALLY ABORTED}"

# STATUS MESSAGES
STATUS = "\\textbf{STATUS}: "
STATUS_OK = STATUS + OK_MESSAGE + "\\\\ \n\n"
STATUS_NOT_APPLY = STATUS + NOT_APPLY_MESSAGE + "\\\\ \n\n"
STATUS_ERROR_MESSAGE = STATUS + ERROR_MESSAGE + "\\\\ \n\n"
STATUS_ERROR_CONVERT_MESSAGE = STATUS + ERROR_CONVERT_MESSAGE + "\\\\ \n\n"
STATUS_ERROR_ABORTED_MESSAGE = STATUS + ERROR_ABORTED_MESSAGE + "\\\\ \n\n"

# RESULT SECTION MESSAGES
ERROR_MESSAGE_PROTOCOL_FAILED = "{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n\n"
ERROR_MESSAGE_PROTOCOL_X_FAILED = "{\\color{red} \\textbf{ERROR: The protocol %s failed.}}\\\\ \n\n"
ERROR_MESSAGE_NOT_BINARY = "{\\color{red} \\textbf{ERROR: Cannot find the binary}}\\\\ \n\n"
ERROR_MESSAGE_HALF1_NOT_RESIZED = "{\\color{red} \\textbf{ERROR: Half map 1 couldn't be resized}}\\\\ \n\n"
ERROR_MESSAGE_HALF2_NOT_RESIZED = "{\\color{red} \\textbf{ERROR: Half map 2 couldn't be resized}}\\\\ \n\n"
ERROR_MESSAGE_NOT_CLASSES = "{\\color{red} \\textbf{ERROR: Cannot find the output classes.}}\\\\ \n\n"
ERROR_MESSAGE_EMPTY_VOL = "{\\color{red} \\textbf{ERROR: The volume is empty.}}\\\\ \n\n"
ERROR_MESSAGE_NO_RESULTS = "{\\color{red} \\textbf{ERROR: The protocol did not produce any result.}}\\\\ \n\n"
ERROR_MESSAGE_ABORTED = "{\\color{red} \\textbf{ERROR: The protocol could not be finished because it was manually aborted.}}\\\\ \n\n"
ERROR_MESSAGE_CHECK_FITTED_FAILED = "{\\color{red} \\textbf{The protocol to check whether the model and map are aligned has failed. Some problems or bad results may occur when running the Level A analysis.}}\\\\ \n\n"

NOT_APPY_NO_RESOLUTION = "{\\color{brown} This method cannot be applied to maps with no resolution reported.}\\\\ \n\n"
NOT_APPLY_BETTER_RESOLUTION = "{\\color{brown} This method cannot be applied to maps with a resolution better than %d\\AA.}\\\\ \n\n"
NOT_APPLY_WORSE_RESOLUTION = "{\\color{brown} This method cannot be applied to maps with a resolution worse than %d\\AA.}\\\\ \n\n"


SUMMARY_WARNINGS_TITLE = "\\textbf{\\underline{Summary of the warnings across sections.}}\\\\ \n\n\n"




