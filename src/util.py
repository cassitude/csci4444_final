import pysam
from sw import sw


def compare(pair, read_pos_map):
    starter = read_pos_map.get(pair[0][0])
    ender = read_pos_map.get(pair[1][0])

    start_sc = starter.query_sequence[0 : starter.cigartuples[0][1]]
    start_seq = starter.query_alignment_sequence
    end_sc = ender.query_sequence[-ender.cigartuples[-1][0] :]
    end_seq = ender.query_alignment_sequence

    if len(start_sc) < len(end_sc):
        sc_sw, align_A, align_B, match = sw(end_sc, start_sc, -2, -1, 1)
    else:
        sc_sw, align_A, align_B, match = sw(start_sc, end_sc, -2, -1, 1)

    if len(start_seq) < len(end_seq):
        seq_sw, align_A, align_B, match = sw(end_seq, start_seq, -2, -1, 1)
    else:
        seq_sw, align_A, align_B, match = sw(start_seq, end_seq, -2, -1, 1)

    if abs(pair[0][1] - pair[1][2]) > 150:
        read_evidence = min(pair[0][1], pair[1][2])
    else:
        read_evidence = (pair[0][1] + pair[1][2]) / 2

    if abs(pair[0][2] - pair[1][3]) > 5:
        max_sc_length = min(pair[0][2], pair[1][3])
    else:
        max_sc_length = (pair[0][2] + pair[1][3]) / 2

    comb_micro_penalty = (pair[0][3] + pair[1][4]) / 2

    return (
        pair[0][0],
        pair[1][1],
        read_evidence,
        max_sc_length,
        comb_micro_penalty,
        sc_sw,
        seq_sw,
    )


def norm_score(combo_reports):
    identifiers = [(report[0], report[1]) for report in combo_reports]
    data_to_norm = [report[2:] for report in combo_reports]

    # get min and max for each attribute to normalize
    mins = [float("inf")] * len(data_to_norm[0])
    maxs = [float("-inf")] * len(data_to_norm[0])

    for row in data_to_norm:
        for i, val in enumerate(row):
            mins[i] = min(mins[i], val)
            maxs[i] = max(maxs[i], val)

    normalized_data = []
    for row in data_to_norm:
        normalized_row = []
        for i, val in enumerate(row):
            if mins[i] == maxs[i]:
                normalized_val = 1.0
            else:
                normalized_val = (val - mins[i]) / (maxs[i] - mins[i])
            normalized_row.append(normalized_val)
        normalized_data.append(normalized_row)

    scores = []
    for row in normalized_data:
        read_evidence = row[0]
        max_sc_length = row[1]
        comb_micro_penalty = row[2]
        sc_sw = row[3]
        eq_sw = row[4]

        score = (
            0.75 * read_evidence
            + 2.5 * max_sc_length
            - 1.5 * comb_micro_penalty
            + 1.75 * sc_sw
            + eq_sw
        )
        scores.append(score)

    combined_data = []
    for i in range(len(combo_reports)):
        combined_data.append(
            {"start": identifiers[i][0], "end": identifiers[i][1], "score": scores[i]}
        )

    best_entries = {}
    for entry in combined_data:
        i = entry["start"]
        if i not in best_entries or entry["score"] > best_entries[i]["score"]:
            best_entries[i] = entry

    return best_entries
