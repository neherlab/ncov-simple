import argparse
import json


def extract_spike_mutations(node_data):
    data = {}
    for name, node in node_data["nodes"].items():
        smuts = node.get("aa_muts", {}).get("S", [])
        if smuts:
            data[name] = ", ".join(smuts)
    return data


def extract_insertions(node_data):
    data = {}
    for name, node in node_data["nodes"].items():
        insertions = node["muts"]
        if insertions:
            data[name] = insertions
    return data


def extract_clade_labels(node_data):
    data = {}
    for name, node in node_data["nodes"].items():
        if "clade_annotation" in node:
            data[name] = node["clade_annotation"]
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove extraneous colorings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input", type=str, metavar="JSON", required=True, help="input Auspice JSON"
    )
    parser.add_argument("--mutations", type=str, required=False, help="mutations node data file")
    parser.add_argument("--insertions", type=str, required=False, help="insertions node data file")
    parser.add_argument(
        "--output", type=str, metavar="JSON", required=True, help="output Auspice JSON"
    )
    args = parser.parse_args()

    with open(args.input, "r") as f:
        auspice_json = json.load(f)

    if args.mutations:
        with open(args.mutations, "r") as f:
            spike_mutations = extract_spike_mutations(json.load(f))
    else:
        spike_mutations = {}

    if args.insertions:
        with open(args.insertions, "r") as f:
            insertions = extract_insertions(json.load(f))
    else:
        insertions = {}

    def attach_labels(n):  # closure
        if n["name"] in spike_mutations:
            if "branch_attrs" not in n:
                n["branch_attrs"] = {}
            if "labels" not in n["branch_attrs"]:
                n["branch_attrs"]["labels"] = {}
            n["branch_attrs"]["labels"]["spike_mutations"] = spike_mutations[n["name"]]

        if n["name"] in insertions:
            if "branch_attrs" not in n:
                n["branch_attrs"] = {}
            if "labels" not in n["branch_attrs"]:
                n["branch_attrs"]["labels"] = {}
            if "aa" in n["branch_attrs"]["labels"]:
                n["branch_attrs"]["labels"]["aa"] += "; "
            else:
                n["branch_attrs"]["labels"]["aa"] = ""
            insertion_string = "; ".join(
                [
                    f'{insertion_type}: {", ".join(insertions[n["name"]][insertion_type])}'
                    for insertion_type in insertions[n["name"]]
                ]
            )
            n["branch_attrs"]["labels"]["aa"] += insertion_string

            n["branch_attrs"]["labels"]["insertions"] = insertion_string
            if n["name"].startswith("NODE_"):
                n["branch_attrs"]["labels"]["insertions (internal)"] = insertion_string

        if "children" in n:
            for c in n["children"]:
                attach_labels(c)

    attach_labels(auspice_json["tree"])

    with open(args.output, "w") as f:
        json.dump(auspice_json, f, indent=2)
