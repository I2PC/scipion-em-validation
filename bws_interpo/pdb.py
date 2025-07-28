import os
from Bio.PDB import MMCIFParser, PDBIO
from io import StringIO

def load_cif_as_pdb(cif_path: os.PathLike):
    parser = MMCIFParser(QUIET=True)
    
    # Parse structure from the CIF file
    structure = parser.get_structure("structure", cif_path)
    
    # Create PDBIO object and save as PDB
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)

    pdb_buffer = StringIO()
    pdb_io.save(pdb_buffer)

    # Get PDB contents as string
    pdb_contents = pdb_buffer.getvalue()
    pdb_buffer.close()
    
    return pdb_contents
