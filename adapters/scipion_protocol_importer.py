from pyworkflow.protocol import Protocol
import pyworkflow.plugin as pwplugin

class ScipionProtocolImporter:
    def import_protocol() -> Protocol:

        Prot = pwplugin.Domain.importFromPlugin('pwem.protocols',
                                            'ProtImportVolumes', doRaise=True)
        
        return {}
