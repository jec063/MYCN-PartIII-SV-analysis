from pathlib import Path
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        description="Count recurrent SVs from candidate SV CSV."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input candidate SV CSV file."
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output recurrent SV CSV file."
    )

    parser.add_argument(
        "-n",
        "--total-samples",
        required=True,
        type=int,
        help="Total number of selected samples."
    )

    args = parser.parse_args()

    input_csv = Path(args.input)
    output_csv = Path(args.output)
    total_samples = args.total_samples

    output_csv.parent.mkdir(exist_ok=True)

    sv_groups = {}

    with open(input_csv, "r", encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file)

        for row in reader:
            sv_key = row["sv_key"]

            if sv_key not in sv_groups:
                sv_groups[sv_key] = {
                    "chr1": row["chr1"],
                    "pos1": row["pos1"],
                    "strand1": row["strand1"],
                    "chr2": row["chr2"],
                    "pos2": row["pos2"],
                    "strand2": row["strand2"],
                    "sv_type": row["sv_type"],
                    "samples": set(),
                }

            sv_groups[sv_key]["samples"].add(row["sample_id"])

    recurrent_rows = []

    for sv_key, sv in sv_groups.items():
        samples = sorted(sv["samples"])
        frequency_count = len(samples)
        frequency_fraction = frequency_count / total_samples

        recurrent_rows.append({
            "chr1": sv["chr1"],
            "pos1": sv["pos1"],
            "strand1": sv["strand1"],
            "chr2": sv["chr2"],
            "pos2": sv["pos2"],
            "strand2": sv["strand2"],
            "sv_type": sv["sv_type"],
            "frequency_count": frequency_count,
            "total_samples": total_samples,
            "frequency_fraction": round(frequency_fraction, 4),
            "samples": ";".join(samples),
        })

    recurrent_rows.sort(
        key=lambda row: row["frequency_count"],
        reverse=True
    )

    for index, row in enumerate(recurrent_rows, start=1):
        row["SV_ID"] = f"SV_{index:03d}"

    columns = [
        "SV_ID",
        "chr1",
        "pos1",
        "strand1",
        "chr2",
        "pos2",
        "strand2",
        "sv_type",
        "frequency_count",
        "total_samples",
        "frequency_fraction",
        "samples",
    ]

    with open(output_csv, "w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(recurrent_rows)

    print(f"Unique SVs: {len(recurrent_rows)}")
    print(f"Wrote: {output_csv}")

main()