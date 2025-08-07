import os
import pathlib
import logging
import re
from collections import namedtuple
import glob
import sqlite3

from typing import List

InputFiles = namedtuple("InputFile", ["emdb_map", "deepres_vol", "structure", "mask"])
ScipionProtocolRun = namedtuple("ScipionProtocolRun", ["identifier", "classname"])


def _glob_re(pattern, strings):
    pattern = re.compile(pattern)
    return filter(lambda x: pattern.match(str(x)), strings)


def _find_file(path, *, suffix, label, pattern=".*"):
    paths = glob.glob(str(path))
    if len(paths) > 1:
        cands_str = "\n".join(map(str, paths))
        logging.warning(
            f"More than one candidate found when expanding paths, the behavior is undefined:\n{cands_str}"
        )

    if len(paths) == 0:
        return None

    path = paths[0]

    cands = list(_glob_re(pattern, pathlib.Path(path).glob(f"*.{suffix}")))
    if len(cands) > 1:
        cands_str = "\n".join(map(str, cands))
        logging.warning(
            f"More than one candidate found for {label}, the behavior is undefined:\n{cands_str}"
        )
    elif len(cands) == 0:
        return None

    return cands[0]


def find_files(*, scipion_project_root: os.PathLike):
    scipion_project_root = pathlib.Path(scipion_project_root)

    # EMDB map
    emdb_map = _find_file(
        scipion_project_root / "Runs" / "000002_ProtImportVolumes" / "extra",
        suffix="map",
        pattern="(.*)emd_([0-9]+).map",
        label="map file",
    )

    deepres_vol = _find_file(
        scipion_project_root / "Runs" / "*_XmippProtDeepRes" / "extra",
        suffix="vol",
        pattern="(.*)deepRes_resolution_originalSize.vol",
        label="DeepRes volume",
    )

    # 000289_XmippProtCreateMask3D
    deepres_mask = _find_file(
        scipion_project_root / "Runs" / "000289_XmippProtCreateMask3D",
        suffix="mrc",
        pattern="(.*).mrc",
        label="DeepRes mask",
    )

    structure = _find_file(scipion_project_root, suffix="cif", label="CIF file")

    return InputFiles(
        str(emdb_map), str(deepres_vol), str(structure), str(deepres_mask)
    )


def find_dependencies(scipion_project_root, *, protocol: str):
    def _fetch_protocol_for_output(identifier, *, cursor: sqlite3.Cursor):
        res = cursor.execute(f"SELECT parent_id FROM Objects WHERE id={identifier}")
        if res.arraysize > 1:
            logging.warning(
                "Multiple parents found upstream protocol; behavior is undefined."
            )

        if res.arraysize == 0:
            raise ValueError(f"No instance found for upstream protocol")

        return res.fetchone()[0]

    def _fetch_protocol_name(identifier, *, cursor: sqlite3.Cursor):
        res = cursor.execute(f"SELECT classname FROM Objects WHERE id={identifier}")
        return res.fetchone()[0]

    def _is_digit(value):
        try:
            return value.isdigit()
        except:
            return False

    db_path = pathlib.Path(scipion_project_root) / "project.sqlite"
    db = sqlite3.connect(db_path)

    cursor = db.cursor()
    res = cursor.execute(f"SELECT id FROM Objects WHERE classname='{protocol}'")

    protocol_id = res.fetchall()
    if len(protocol_id) > 1:
        logging.warning(
            "Multiple instances found for target protocol; behavior is undefined."
        )

    if len(protocol_id) == 0:
        raise ValueError("No instance found for target protocol")

    (protocol_id,) = protocol_id[0]

    res = cursor.execute(
        f"SELECT value FROM Objects where parent_id == {protocol_id} AND classname=='Pointer'"
    )

    # print(res.fetchall())

    output_ids = [
        int(identifier) for identifier, in res.fetchall() if _is_digit(identifier)
    ]
    upstream_ids = [_fetch_protocol_for_output(i, cursor=cursor) for i in output_ids]

    upstream = [
        ScipionProtocolRun(i, _fetch_protocol_name(i, cursor=cursor))
        for i in upstream_ids
    ]

    return upstream


def find_dependency_filenames(scipion_project_root, *, protocol: str, query: List[str]):
    upstream = find_dependencies(scipion_project_root, protocol=protocol)
    upstream = {u.classname: u for u in upstream}

    try:
        file_info = [upstream[n] for n in query]
    except KeyError as e:
        raise KeyError(f"Could not find protocol '{e.args[0]}' in project")

    filenames = [f"{str(f.identifier).zfill(6)}_{f.classname}" for f in file_info]
    return filenames


if __name__ == "__main__":
    find_dependency_filenames(
        "/home/max/Documents/val-server/data/val-report-service/EMD-41510",
        protocol="XmippProtDeepRes",
        query=["ProtImportVolumes", "XmippProtCreateMask3D"],
    )
