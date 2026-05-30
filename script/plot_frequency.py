from pathlib import Path
import argparse
import csv
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description="Plot frequent recurrent SVs."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input recurrent SV CSV file."
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output plot file, such as PNG or PDF."
    )

    parser.add_argument(
        "-t",
        "--top",
        type=int,
        default=10,
        help="Number of top SVs to plot."
    )

    args = parser.parse_args()

    input_csv = Path(args.input)
    output_plot = Path(args.output)
    top_n = args.top

    output_plot.parent.mkdir(exist_ok=True)

    rows = []

    with open(input_csv, "r", encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file)

        for row in reader:
            rows.append(row)

    rows = rows[:top_n]

    sv_ids = []
    frequencies = []

    for row in rows:
        sv_ids.append(row["SV_ID"])
        frequencies.append(int(row["frequency_count"]))

    plt.figure(figsize=(10, 5))
    plt.bar(sv_ids, frequencies)

    plt.xlabel("Structural variant")
    plt.ylabel("Number of samples")
    plt.title("Top recurrent structural variants in MYCN-amplified samples")

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    plt.savefig(output_plot, dpi=300)
    plt.close()

    print(f"Wrote: {output_plot}")

main()