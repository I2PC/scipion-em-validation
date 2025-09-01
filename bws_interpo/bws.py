import json
import datetime
import itertools
import statistics

from collections import namedtuple

Atom = namedtuple("Atom", ["name", "res_seq_name",
                  "res_seq_number", "chain", "value"])


def save_for_bws(input_path, output_path, emd_entry, pdb_entry, *, title):
    with open(input_path) as file:

        chains_dict = {}
        for line in file.readlines():
            if not line.startswith("ATOM"):
                continue

            # _validate_atom_line(line)

            atom_name = line[12:16].strip()
            res_seq_name = line[17:20]
            res_seq_number = int(line[22:26])
            chain_id = line[21]
            value = float(line[54:60].strip())

            atom = Atom(
                atom_name,
                res_seq_name=res_seq_name,
                res_seq_number=res_seq_number,
                chain=chain_id,
                value=value,
            )
            chains_dict[chain_id] = (
                chains_dict[chain_id] +
                [atom] if chain_id in chains_dict else [atom]
            )

        chains = []
        for name, values in chains_dict.items():
            sequence = [
                {
                    "resSeqName": v.res_seq_name,
                    "resSeqNumber": v.res_seq_number,
                    "scoreValue": v.value,
                }
                for v in values
            ]

            grouped_seq = []
            for _, g in itertools.groupby(sequence, key=lambda x: x["resSeqNumber"]):
                g = list(g)
                avg = statistics.median([e["scoreValue"] for e in g])
                grouped_seq.append(
                    {
                        "resSeqName": g[0]["resSeqName"],
                        "resSeqNumber": g[0]["resSeqNumber"],
                        "scoreValue": avg,
                    }
                )

            chain = {"name": name, "seqData": grouped_seq}

            chains.append(chain)

        outputs = {
            "resource": f"{title} (VRS)",
            "entry": {
                "volume_map": emd_entry,
                "atomic_model": pdb_entry,
                "date": datetime.date.today().isoformat(),
            },
            "chains": chains,
        }

        with open(output_path, "x") as f:
            json.dump(outputs, f)
