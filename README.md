# scipion-em-validation
Core code of the cryo-EM Validation Report Service presented in https://doi.org/10.1039/D2FD00059H

This scheme has been implemented in a website https://biocomp.cnb.csic.es/EMValidationService/ to which reconstructed maps and their ESI can be uploaded.

#### Commands to launch validations over EMDB
##### Single entry
- scipion3 python validationLevels.py EMDBid=EMD-29734 doLevels=0,1 > "/home/EMDB/logLauncher/$(date +'%Y-%m-%d_%H-%M').txt" 2>&1
- scipion3 python validationLevels.py EMDBid=EMD-29734 doLevels=0,1 --isTest > "/home/EMDB/logLauncher/$(date +'%Y-%m-%d_%H-%M').txt" 2>&1

IMPORTANT: --isTest flag is meant to avoid storing intermediate data in /data/ScipionUserDataToRetrieve

##### All entries for a specific level
- scipion3 python tools/EMDBlauncher.py --launchAll --level=0 -n=5744 --startEntry=1 > "/home/EMDB/logLauncher/$(date +'%Y-%m-%d_%H-%M').txt" 2>&1
- scipion3 python tools/EMDBlauncher.py --launchAll --level=0 -n=5744 --startEntry=1 --isTest > "/home/EMDB/logLauncher/$(date +'%Y-%m-%d_%H-%M').txt" 2>&1

IMPORTANT: --isTest flag is meant to avoid storing intermediate data in /data/ScipionUserDataToRetrieve

##### Random set of entries for a specific level
- scipion3 python tools/EMDBlauncher.py --launchAll --level=0,1,A --random 10 > "/home/EMDB/logLauncher/$(date +'%Y-%m-%d_%H-%M').txt" 2>&1
- scipion3 python tools/EMDBlauncher.py --launchAll --level=0,1,A --random 10 --isTest > "/home/EMDB/logLauncher/$(date +'%Y-%m-%d_%H-%M').txt" 2>&1

IMPORTANT: --isTest flag is meant to avoid storing intermediate data in /data/ScipionUserDataToRetrieve

##### Entries from a list of a specific level
- scipion3 python tools/EMDBlauncher.py --launchList --inputList /home/EMDB/entries/level0+1_n50_first_batch.txt level=0,1 > "/home/EMDB/logLauncher/$(date +'%Y-%m-%d_%H-%M').txt" 2>&1
- scipion3 python tools/EMDBlauncher.py --launchList --inputList /home/EMDB/entries/level0+1_n50_first_batch.txt level=0,1 --isTest > "/home/EMDB/logLauncher/$(date +'%Y-%m-%d_%H-%M').txt" 2>&1

IMPORTANT: --isTest flag is meant to avoid storing intermediate data in /data/ScipionUserDataToRetrieve
