from pathlib import Path
import argparse
import csv


COHORT_COL = "Cohort"
SAMPLE_COL = "Sample Name"
FEATURE_COL = "Feature ID"
CLASSIFICATION_COL = "Classification"
CANCER_TYPE_COL = "Cancer Type"
GRAPH_PATH_COL = "Graph file path / notes"


def clean_graph_path(path_text, project_dir):
    path_text = path_text.strip().strip('"').strip("'")
    path_text = path_text.replace("\\", "/")
    path_text = path_text.replace("/Project/", "", 1)

    return project_dir / Path(path_text)


def parse_position(position_text):
    chromosome, position_with_strand = position_text.split(":", 1)

    strand = position_with_strand[-1]
    position = int(position_with_strand[:-1])

    return chromosome, position, strand


def parse_edge(edge_text):
    start_position, end_position = edge_text.split("->", 1)
    
    return start_position, end_position


def make_sv_key(chr1, pos1, strand1, chr2, pos2, strand2):
    endpoint1 = (chr1, pos1, strand1)
    endpoint2 = (chr2, pos2, strand2)

    first, second = sorted([endpoint1, endpoint2])

    return (
        f"{first[0]}:{first[1]}:{first[2]}|"
        f"{second[0]}:{second[1]}:{second[2]}"
    )


def get_sv_type(chr1, chr2):
    if chr1 == chr2:
        return "intrachromosomal"

    return "interchromosomal"


def extract_discordant_from_graph(graph_path, sample_info):
    sv_rows = []

    with open(graph_path, "r", encoding="utf-8", errors="replace") as graph_file:
        for line_number, line in enumerate(graph_file, start=1):
            line = line.strip()

            if line == "":
                continue

            parts = line.split()

            if parts[0].lower() != "discordant":
                continue

            start_position, end_position = parse_edge(parts[1])

            chr1, pos1, strand1 = parse_position(start_position)
            chr2, pos2, strand2 = parse_position(end_position)

            copy_count = parts[2]
            read_pairs = parts[3]

            sv_type = get_sv_type(chr1, chr2)
            sv_key = make_sv_key(chr1, pos1, strand1, chr2, pos2, strand2)

            sv_rows.append({
                "cohort": sample_info["cohort"],
                "sample_id": sample_info["sample_id"],
                "feature_id": sample_info["feature_id"],
                "classification": sample_info["classification"],
                "cancer_type": sample_info["cancer_type"],
                "graph_file": sample_info["graph_file"],
                "line_number": line_number,
                "raw_line": line,
                "chr1": chr1,
                "pos1": pos1,
                "strand1": strand1,
                "chr2": chr2,
                "pos2": pos2,
                "strand2": strand2,
                "copy_count": copy_count,
                "read_pairs": read_pairs,
                "sv_type": sv_type,
                "sv_key": sv_key,
            })

    return sv_rows


def main():
    parser = argparse.ArgumentParser(
        description="Extract discordant edges from AA graph files."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input selected sample TSV file."
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output CSV file for candidate SVs."
    )

    args = parser.parse_args()

    project_dir = Path.cwd()
    input_tsv = Path(args.input)
    output_csv = Path(args.output)

    output_csv.parent.mkdir(exist_ok=True)

    all_svs = []

    with open(input_tsv, "r", encoding="utf-8-sig", newline="") as sample_file:
        reader = csv.DictReader(sample_file, delimiter="\t")

        for row in reader:
            graph_path = clean_graph_path(row[GRAPH_PATH_COL], project_dir)

            sample_info = {
                "cohort": row[COHORT_COL],
                "sample_id": row[SAMPLE_COL],
                "feature_id": row[FEATURE_COL],
                "classification": row[CLASSIFICATION_COL],
                "cancer_type": row[CANCER_TYPE_COL],
                "graph_file": str(graph_path),
            }

            sample_svs = extract_discordant_from_graph(graph_path, sample_info)
            all_svs.extend(sample_svs)

            print(f"{row[SAMPLE_COL]}: {len(sample_svs)} discordant SVs")

    columns = [
        "cohort",
        "sample_id",
        "feature_id",
        "classification",
        "cancer_type",
        "graph_file",
        "line_number",
        "raw_line",
        "chr1",
        "pos1",
        "strand1",
        "chr2",
        "pos2",
        "strand2",
        "copy_count",
        "read_pairs",
        "sv_type",
        "sv_key",
    ]

    with open(output_csv, "w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(all_svs)

    print()
    print(f"Total discordant SVs: {len(all_svs)}")
    print(f"Wrote: {output_csv}")

main()