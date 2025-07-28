import os
from pathlib import Path
import requests
import json
import gzip
import shutil

from collections import namedtuple

EMDB_EBI_REPOSITORY = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures"
EMDB_EBI_JSON_REPOSITORY = "https://www.ebi.ac.uk/emdb/api/entry"
PDB_EBI_REPOSITORY = "https://www.ebi.ac.uk/pdbe/entry-files/download"

EMDBMetadata = namedtuple("EMDBMetadata", ("pdb_id", "resolution",
                          "sampling", "size", "org_x", "org_y", "org_z"))


def download_emdb_metadata(entry: str) -> EMDBMetadata:
    """
    Downloads metadata for an emdb entry

    Args:
        - entry: The EMDB ID to download. Use only the numbers after the 'EMD-' prefix
    """

    url = f"{EMDB_EBI_JSON_REPOSITORY}/{entry}"

    response = requests.get(url)
    response.raise_for_status()

    raw_data = json.loads(response.content)

    map_info = raw_data["map"]
    resolution = float(
        raw_data["structure_determination_list"]["structure_determination"][0][
            "image_processing"][0]["final_reconstruction"]["resolution"]["valueOf_"]
    )

    try:
        pdb_id = raw_data["structure_determination_list"]["structure_determination"][
            0]["image_processing"][0]["startup_model"][0]["pdb_model"]["pdb_id"]
        pdb_id = str(pdb_id).lower()
    except KeyError:
        pdb_id = None

    metadata = EMDBMetadata(
        pdb_id=pdb_id,
        resolution=resolution,
        sampling=float(map_info["pixel_spacing"]["y"]["valueOf_"]),
        size=int(map_info["dimensions"]["col"]),
        org_x=-(map_info["origin"]["col"]),
        org_y=-(map_info["origin"]["sec"]),
        org_z=-(map_info["origin"]["row"]),
    )

    return metadata


def _download(url, location: os.PathLike):
    location = Path(location)

    response = requests.get(url)
    response.raise_for_status()

    with open(location, mode="wb") as file:
        file.write(response.content)

    return os.path.abspath(location)


def download_emdb_map(map_id, path: os.PathLike):
    path = Path(path)

    filename = f"emd_{map_id}.map"
    gz_filename = f"{filename}.gz"

    url = f"{EMDB_EBI_REPOSITORY}/EMD-{map_id}/map/{gz_filename}"
    dest_path_compressed = path / gz_filename
    dest_path = path / filename

    _ = _download(url, dest_path_compressed)

    with gzip.open(dest_path_compressed, 'rb') as f_in:
        with open(dest_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return dest_path


def download_pdb_model(pdb_id, path: os.PathLike):
    path = Path(path)

    filename = f"pdb{pdb_id}.ent"

    url = f"{PDB_EBI_REPOSITORY}/pdb{pdb_id}.ent"
    dest_path = path / filename

    return _download(url, dest_path)
