import os
import re
import logging
from pathlib import Path
from .find_files import _find_file, find_dependency_filenames
from functools import partial

from .pdb import load_cif_as_pdb
from .bws import save_for_bws
# from .scipion_bridge.environment import configure_default_env
# from .ffi.scipion import xmipp_pdb_label_from_volume
# from .utils.bws import save_for_bws
from .download import download_emdb_metadata
from .external_call import foreign_function, Domain

import argparse
from collections import namedtuple
from tempfile import NamedTemporaryFile
from typing import Optional

InputFiles = namedtuple("InputFile", ["volume", "mask", "structure"])
xmipp_func = partial(foreign_function, domain=Domain(
    "XMIPP", ["scipion", "run"])
)


def fetch_files(protocol: str, *, project_root: os.PathLike, volume: str):
    project_root = Path(project_root)

    suffix = "".join(Path(volume).suffixes)[1:]

    vol_path = _find_file(
        project_root / "Runs" / f"*_{protocol}" / "extra",
        suffix=suffix,
        pattern=f"(.*){volume}",
        label="Volume",
    )

    (mask_filename,) = find_dependency_filenames(
        project_root,
        protocol=protocol,
        query=["XmippProtCreateMask3D"],
    )

    mask_path = _find_file(
        project_root / "Runs" / mask_filename,
        suffix="mrc",
        pattern="(.*).mrc",
        label="DeepRes mask",
    )

    structure_path = _find_file(project_root, suffix="cif", label="CIF file")

    return InputFiles(volume=vol_path, mask=mask_path, structure=structure_path)


def find_emdb_identifier(project_root: os.PathLike):
    project_root = str(project_root)  # type: ignore

    pattern = re.compile("EMD-[0-9]+")
    matches = re.findall(pattern, project_root)  # type: ignore

    if len(matches) > 1:
        logging.warning(
            "Multiple possible EMDB entries found; behavior is undefined")

    if not matches:
        raise ValueError("No EMDB entry found")

    return int(matches[0][4:])


@partial(
    xmipp_func,
    args_map={
        "outputs": "o",
        "volume": "vol",
    },
)
def xmipp_pdb_label_from_volume(
    outputs: str,
    *,
    pdb: str,
    volume: str,
    mask: str,
    sampling: str,
    origin: str,
):
    pass


def convert(
    protocol: str, emb_entry: Optional[str] = None, *, project_root: os.PathLike, volume: str, **kwargs
):

    inputs = fetch_files(protocol, project_root=project_root, volume=volume)

    structure = load_cif_as_pdb(inputs.structure)
    # structure = TempFileProxy.proxy_for_string(structure, file_ext="pdb")

    if emb_entry is None:
        emb_id = find_emdb_identifier(project_root)
        emb_entry = f"EMD-{emb_id}"

    pdb_entry = os.path.split(inputs.structure)[-1][:-4]
    metadata = download_emdb_metadata(emb_entry)  # type: ignore

    pdb_path = Path(project_root) / f"{protocol}.pdb"

    with open(pdb_path, mode="x+") as f:
        f.write(structure)

        path_atomic_model = Path(project_root) / f"{protocol}.atom.pdb"
        path_bws = Path(project_root) / \
            f"{protocol}_{emb_entry}_{pdb_entry}.json"

        xmipp_pdb_label_from_volume(
            outputs=path_atomic_model,
            pdb=str(pdb_path),
            volume=inputs.volume,
            mask=inputs.mask,
            sampling=metadata.sampling,
            origin="%f %f %f" % (
                metadata.org_x, metadata.org_y, metadata.org_z
            ),
        )

        save_for_bws(path_atomic_model, path_bws,
                     emb_entry, pdb_entry, title=protocol
                     )

        return path_bws

# Example Usage:
# DeepRes:
# scipion3 python -m emv_tools.convert_eval_results
# -project '/home/max/Documents/val-server/data/val-report-service/EMD-41510'
# -o /home/max/Documents/val-server/EMV-Script-fork/emv-tools/data/converted.json
# -n XmippProtDeepRes
# -vol deepRes_resolution_originalSize.vol

# MonoRes:
# scipion3 python -m emv_tools.convert_eval_results
# -project '/home/max/Documents/val-server/data/val-report-service/EMD-41510'
# -o /home/max/Documents/val-server/EMV-Script-fork/emv-tools/data/converted.json
# -n XmippProtMonoRes
# -vol monoresResolutionMap.mrc

# BlocRes:
# scipion3 python -m emv_tools.convert_eval_results
# -project '/home/max/Documents/val-server/data/val-report-service/EMD-41510'
# -o /home/max/Documents/val-server/EMV-Script-fork/emv-tools/data/converted.json
# -n BsoftProtBlocres
# -vol resolutionMap.map

# FSC-Q
# scipion3 python -m emv_tools.convert_eval_results
# -project '/home/max/Documents/val-server/data/val-report-service/EMD-41510'
# -o /home/max/Documents/val-server/EMV-Script-fork/emv-tools/data/converted.json
# -n XmippProtValFit
# -vol diferencia.map
