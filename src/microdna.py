import pysam
from util import compare
from util import norm_score
import argparse


def get_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True, type=str, help="file to read")
    parser.add_argument(
        "-t",
        "--threshold",
        type=int,
        default=10,
        help="minimum observations of a read to continue",
    )
    parser.add_argument(
        "-c", "--cutoff", type=int, default=1.5, help="score reporting cutoff"
    )
    return parser.parse_args(args)


def main(args=None):
    args = get_args()
    # Path to your BAM file "data/SRR413984.sorted.NC_000001.10.bam"
    bam_file = args.file

    threshold = args.threshold
    excluded_ops = [1, 2, 3, 5, 6, 7, 8, 9]

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        window = []
        sc_start_length = 0
        sc_end_length = 0
        start_homo_penalty = 0
        end_homo_penalty = 0

        potential_starts = []
        potential_ends = []

        for read in bam:
            # exclude reads w/ operations other that match/soft clip
            if any(t[0] in excluded_ops for t in read.cigartuples):
                continue

            # check if the read has beginning soft clips in the CIGAR string
            if read.cigartuples[0][0] == 4:
                if len(window) == 0:
                    window.append([read.pos, read.pos])
                    sc_start_length = read.cigartuples[0][1]
                    if len(read.cigartuples) > 1:
                        for t in read.cigartuples[1:]:
                            if t[0] == 4:
                                start_homo_penalty += t[1]
                else:
                    if read.pos == window[0][0]:
                        window.append([read.pos, read.pos])
                        if read.cigartuples[0][1] > sc_start_length:
                            sc_start_length = read.cigartuples[0][1]
                        if len(read.cigartuples) > 1:
                            for t in read.cigartuples[1:]:
                                if t[0] == 4:
                                    start_homo_penalty += t[1]
                    else:
                        if len(window) >= threshold:
                            potential_starts.append(
                                (
                                    window[0][0],
                                    len(window),
                                    sc_start_length,
                                    start_homo_penalty,
                                )
                            )
                        window = [[read.pos, read.pos]]
                        sc_start_length = read.cigartuples[0][1]
                        start_homo_penalty = 0
                        if len(read.cigartuples) > 1:
                            for t in read.cigartuples[1:]:
                                if t[0] == 4:
                                    start_homo_penalty += t[1]

            # Check if the read has ending soft clips in the CIGAR string
            if read.cigartuples[-1][0] == 4:
                # pos of the last-aligned base
                rel_pos = read.pos + read.cigartuples[0][1] - 1
                true_pos = read.pos
                if len(window) == 0:
                    window.append([true_pos, rel_pos])
                    sc_end_length = read.cigartuples[-1][1]
                    if len(read.cigartuples) > 1:
                        for t in read.cigartuples[:-1]:
                            if t[0] == 4:
                                end_homo_penalty += t[1]
                else:
                    if rel_pos == window[0][0]:
                        window.append([true_pos, rel_pos])
                        if read.cigartuples[-1][1] > sc_end_length:
                            sc_end_length = read.cigartuples[-1][1]
                        if len(read.cigartuples) > 1:
                            for t in read.cigartuples[:-1]:
                                if t[0] == 4:
                                    end_homo_penalty += t[1]
                    else:
                        if len(window) >= threshold:
                            potential_ends.append(
                                (
                                    window[0][0],
                                    window[0][1],
                                    len(window),
                                    sc_end_length,
                                    end_homo_penalty,
                                )
                            )
                        window = [[true_pos, rel_pos]]
                        sc_end_length = read.cigartuples[-1][1]
                        end_homo_penalty = 0
                        if len(read.cigartuples) > 1:
                            for t in read.cigartuples[:-1]:
                                if t[0] == 4:
                                    end_homo_penalty += t[1]
        # find potential (start, end) matches
        potential_combos = [
            (start, end)
            for start in potential_starts
            for end in potential_ends
            if 200 <= end[1] - start[0] <= 400
        ]

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        read_pos_map = {}
        for read in bam:
            if any(t[0] in excluded_ops for t in read.cigartuples):
                continue
            if read.cigartuples[0][0] == 4 or read.cigartuples[-1][0] == 4:
                read_pos_map[read.pos] = read

        combo_reports = []
        for combo in potential_combos:
            combo_reports.append(compare(combo, read_pos_map))

    best_results = norm_score(combo_reports)

    sorted_results = sorted(
        best_results.values(), key=lambda x: x["score"], reverse=True
    )

    for i, entry in enumerate(sorted_results):
        if entry["score"] > args.cutoff:
            print(
                f"\n{i + 1}. {entry['start']} -- {entry['end']}: Score {entry['score']:.4f}"
            )


if __name__ == "__main__":
    main()
