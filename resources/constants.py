####################### REPORT MESSAGES CONSTANTS #######################
# Summary
OK_MESSAGE = "{\\color{blue} OK}"
WARNINGS_MESSAGE = "{\\color{red} %d warnings}"
NOT_APPLY_MESSAGE = "{\\color{brown} Does not apply}"
ERROR_MESSAGE = "{\\color{red} Could not be measured}"
ERROR_CONVERT_MESSAGE = "{\\color{red} Could not be converted}"

# STATUS MESSAGES
STATUS = "\\textbf{STATUS}: "
STATUS_OK = STATUS + OK_MESSAGE + "\\\\ \n"
STATUS_NOT_APPLY = STATUS + NOT_APPLY_MESSAGE + "\\\\ \n"
STATUS_ERROR_MESSAGE = STATUS + ERROR_MESSAGE + "\\\\ \n"
STATUS_ERROR_CONVERT_MESSAGE = STATUS + ERROR_CONVERT_MESSAGE + "\\\\ \n"

# RESULT SECTION MESSAGES
ERROR_MESSAGE_PROTOCOL_FAILED = "{\\color{red} \\textbf{ERROR: The protocol failed.}}\\\\ \n"
ERROR_MESSAGE_NOT_BINARY = "{\\color{red} \\textbf{ERROR: Cannot find the binary}}\\\\ \n"
ERROR_MESSAGE_NOT_RESIZED = "{\\color{red} \\textbf{ERROR: Half maps couldn't be resized}}\\\\ \n"
ERROR_MESSAGE_NOT_CLASSES = "{\\color{red} \\textbf{ERROR: Cannot find the output classes.}}\\\\ \n"
ERROR_MESSAGE_EMPTY_VOL = "{\\color{red} \\textbf{ERROR: The volume is empty.}}\\\\ \n"
ERROR_MESSAGE_NO_RESULTS = "{\\color{red} \\textbf{ERROR: The protocol did not produce any result.}}\\\\ \n"


SUMMARY_WARNINGS_TITLE = "\\textbf{\\underline{Summary of the warnings across sections.}}\\\\ \n\n\n"

