from pathlib import Path
import argparse
import csv


MYCN_CHR = "chr2"
MYCN_START = 16_225_000
MYCN_END = 16_227_000


def norm_chr(chrom):
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        return chrom
    return "chr" + chrom


def distance_to_mycn(chrom, pos):
    chrom = norm_chr(chrom)

    if chrom != MYCN_CHR:
        return ""

    if MYCN_START <= pos <= MYCN_END:
        return 0

    if pos < MYCN_START:
        return MYCN_START - pos

    return pos - MYCN_END


def get_near_mycn(chrom1, pos1, chrom2, pos2, window):
    d1 = distance_to_mycn(chrom1, pos1)
    d2 = distance_to_mycn(chrom2, pos2)

    distances = []
    for d in [d1, d2]:
        if d != "":
            distances.append(d)

    if not distances:
        return "no", ""

    min_distance = min(distances)

    if min_distance <= window:
        return "yes", min_distance

    return "no", min_distance


def get_cohort_pattern(samples):
    adult = 0
    pediatric = 0

    for sample in samples:
        if sample.startswith("DO"):
            adult += 1
        elif sample.startswith("BS"):
            pediatric += 1

    if adult > 0 and pediatric == 0:
        pattern = "PCAWG/adult only"
    elif pediatric > 0 and adult == 0:
        pattern = "Pediatric only"
    elif adult > 0 and pediatric > 0:
        pattern = "Shared adult and pediatric"
    else:
        pattern = "Unknown"

    return adult, pediatric, pattern


def make_interpretation(near_mycn, cohort_pattern, sv_type):
    if near_mycn == "yes":
        if cohort_pattern == "Pediatric only":
            return f"Pediatric-specific recurrent {sv_type} breakpoint near MYCN"
        if cohort_pattern == "PCAWG/adult only":
            return f"Adult-specific recurrent {sv_type} breakpoint near MYCN"
        if cohort_pattern == "Shared adult and pediatric":
            return f"Shared recurrent {sv_type} breakpoint near MYCN"
        return f"Recurrent {sv_type} breakpoint near MYCN"

    return f"Recurrent {sv_type} breakpoint outside the MYCN window"


def main():
    parser = argparse.ArgumentParser(
        description="Create annotation table for top recurrent MYCN SVs."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input recurrent SV CSV file."
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output annotated top SV CSV file."
    )

    parser.add_argument(
        "-t", "--top",
        type=int,
        default=15,
        help="Number of top SVs to annotate."
    )

    parser.add_argument(
        "--window",
        type=int,
        default=1_000_000,
        help="Distance window around MYCN."
    )

    args = parser.parse_args()

    input_csv = Path(args.input)
    output_csv = Path(args.output)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    rows = []

    with open(input_csv, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            rows.append(row)

    rows = rows[:args.top]

    annotated = []

    for row in rows:
        chr1 = norm_chr(row["chr1"])
        chr2 = norm_chr(row["chr2"])
        pos1 = int(row["pos1"])
        pos2 = int(row["pos2"])

        samples = row["samples"].split(";") if row.get("samples") else []
        adult_count, pediatric_count, cohort_pattern = get_cohort_pattern(samples)

        near_mycn, distance = get_near_mycn(chr1, pos1, chr2, pos2, args.window)

        interpretation = make_interpretation(
            near_mycn,
            cohort_pattern,
            row["sv_type"]
        )

        annotated.append({
            "SV_ID": row["SV_ID"],
            "breakpoint_1": f"{chr1}:{pos1}{row['strand1']}",
            "breakpoint_2": f"{chr2}:{pos2}{row['strand2']}",
            "sv_type": row["sv_type"],
            "frequency_count": row["frequency_count"],
            "total_samples": row["total_samples"],
            "frequency_fraction": row["frequency_fraction"],
            "adult_count": adult_count,
            "pediatric_count": pediatric_count,
            "cohort_pattern": cohort_pattern,
            "near_MYCN_1Mb": near_mycn,
            "distance_to_MYCN_bp": distance,
            "samples": row["samples"],
            "interpretation": interpretation,
        })

    columns = [
        "SV_ID",
        "breakpoint_1",
        "breakpoint_2",
        "sv_type",
        "frequency_count",
        "total_samples",
        "frequency_fraction",
        "adult_count",
        "pediatric_count",
        "cohort_pattern",
        "near_MYCN_1Mb",
        "distance_to_MYCN_bp",
        "samples",
        "interpretation",
    ]

    with open(output_csv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(annotated)

    print(f"Wrote: {output_csv}")


if __name__ == "__main__":
    main()