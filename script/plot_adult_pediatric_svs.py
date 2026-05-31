from pathlib import Path
import argparse
import csv
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Plot adult vs pediatric counts for top recurrent SVs."
    )

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-t", "--top", type=int, default=15)

    args = parser.parse_args()

    input_csv = Path(args.input)
    output_png = Path(args.output)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    rows = []

    with open(input_csv, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    rows = rows[:args.top]

    sv_ids = [row["SV_ID"] for row in rows]
    adult_counts = [int(row["adult_count"]) for row in rows]
    pediatric_counts = [int(row["pediatric_count"]) for row in rows]

    plt.figure(figsize=(12, 5))

    plt.bar(sv_ids, adult_counts, label="PCAWG adult")
    plt.bar(sv_ids, pediatric_counts, bottom=adult_counts, label="Pediatric")

    plt.xlabel("Structural variant")
    plt.ylabel("Number of samples")
    plt.title("Adult vs pediatric recurrence of top MYCN-associated SVs")
    plt.xticks(rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_png, dpi=300)
    plt.close()

    print(f"Wrote: {output_png}")


if __name__ == "__main__":
    main()
