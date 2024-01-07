import requests
from typing import List
from io import StringIO
from Bio import SeqIO, SeqRecord
import time


GENOME_NAME = "DQ011153.1"  # Monkeypox, 0.44 sec
# GENOME_NAME = "NC_060948.1"  # Human, 125 sec

# algorithm parameters
MIN_WINDOW_SIZE = 200
MIN_GC_PERCENTAGE = 0.5
MIN_OBSERVED_TO_EXPECTED_CPG = 0.6
MIN_ISLAND_MERGE_GAP = 100


class Island:
    def __init__(
        self,
        interval_begin: int,
        interval_end: int,
        gc_perc: float,
        obs_exp: float,
        win_length: str,
    ):
        self.interval_begin = interval_begin
        self.interval_end = interval_end
        self.gc_perc = gc_perc
        self.obs_exp = obs_exp
        self.win_length = win_length

    def __str__(self):
        return f"Island at ({self.interval_begin}, {self.interval_end})\tgcper:{self.gc_perc}\tobs_exp:{self.obs_exp}\twin length:{self.win_length}"


def get_sequence_from_ncbi(genome_id: str) -> str:
    response = requests.get(
        "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=" + genome_id
    )
    genome = response.content.decode()
    if genome.startswith("Error"):
        raise Exception("Error, genome not found")
    print("genome found")

    fasta_io = StringIO(genome)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    return records[0].seq


def find_islands_simple(record: SeqRecord) -> List[Island]:
    record_length = len(record)
    if record_length < MIN_WINDOW_SIZE:
        return []
    found_islands = []
    window_position = 0

    window = record[window_position : window_position + MIN_WINDOW_SIZE]
    gc_count = window.count("G") + window.count("C")
    obs_cpg = window.count("CG")
    while True:
        # gc_count = window.count("G") + window.count("C")
        # obs_cpg = window.count("CG")
        gc_perc = gc_count / MIN_WINDOW_SIZE
        exp_cpg = ((gc_count / 2) ** 2) / MIN_WINDOW_SIZE
        try:
            obs_exp = obs_cpg / exp_cpg
        except ZeroDivisionError:
            obs_exp = 0.0
        if gc_perc > MIN_GC_PERCENTAGE and obs_exp > MIN_OBSERVED_TO_EXPECTED_CPG:
            found_islands.append(
                Island(
                    interval_begin=window_position,
                    interval_end=window_position + MIN_WINDOW_SIZE,
                    gc_perc=gc_perc,
                    obs_exp=obs_exp,
                    win_length=MIN_WINDOW_SIZE,
                )
            )

        window_position += 1
        if window_position + MIN_WINDOW_SIZE > record_length:
            break

        window = record[window_position : window_position + MIN_WINDOW_SIZE]
        # update CGs
        if record[window_position - 1] in ["C", "G"]:
            gc_count -= 1
        if window[MIN_WINDOW_SIZE - 1] in ["C", "G"]:
            gc_count += 1
        # update CpGs
        if record[window_position - 1] == "C" and window[0] == "G":
            obs_cpg -= 1
        if window[MIN_WINDOW_SIZE - 2] == "C" and window[MIN_WINDOW_SIZE - 1] == "G":
            obs_cpg += 1

    return found_islands


def find_island_1bp_shifts(
    record_position: int, record: SeqRecord, min_window_size: str, record_length: int
) -> int:
    window = record[record_position : record_position + min_window_size]
    gc_count, obs_cpg = new_window_gcs(window)
    while True:
        gc_perc, obs_exp, is_island = island_evaluation(gc_count, obs_cpg, min_window_size)
        if is_island:
            print(f"island found {record_position, record_position + min_window_size}")
            return record_position
            # found_islands.append(Island(begin=record_position,
            #                             end=window_position + MIN_WINDOW_SIZE,
            #                             gc_perc=gc_perc,
            #                             obs_exp=obs_exp,
            #                             win_length=MIN_WINDOW_SIZE))
        if record_position + min_window_size + 1 > record_length:
            print("reached the end of record")
            return -1
        gc_count, obs_cpg, record_position = shift_window_1bp(
            gc_count, obs_cpg, record, record_position
        )


def new_window_gcs(window):
    gc_count = window.count("G") + window.count("C")
    obs_cpg = window.count("CG")
    return gc_count, obs_cpg


def find_islands(record: SeqRecord) -> List[Island]:
    record_length = len(record)
    if record_length < MIN_WINDOW_SIZE:
        return []
    found_islands = []

    record_position = 0

    while True:
        # step 1 find window that is an island
        print("- Stage 1 1bp shifts")
        record_position = find_island_1bp_shifts(
            record_position, record, MIN_WINDOW_SIZE, record_length
        )
        if record_position == -1:
            break

        # step 2 window-length shifts
        print("- Stage 2 window-length shifts")
        island_start = record_position
        record_position += MIN_WINDOW_SIZE  # perform a shift and evaluate the following windows
        island_end = extend_island_window_shifts(record, record_length, record_position)
        print(f"island: {island_start}, {island_end}")

        # step 3 shift the last window by 1bp until it meets the criteria
        print("- Stage 3 1bp rollback")
        window = record[island_end - MIN_WINDOW_SIZE : island_end]
        gc_count, obs_cpg = new_window_gcs(window)

        for i in range(1, MIN_WINDOW_SIZE):
            if record[island_end - i] in ["C", "G"]:
                gc_count -= 1
            if record[island_end - MIN_WINDOW_SIZE - i] in ["C", "G"]:
                gc_count += 1
            if record[island_end - i - 2 : island_end - i] == "CG":
                obs_cpg -= 1
            if (
                record[island_end - MIN_WINDOW_SIZE - i - 2 : island_end - MIN_WINDOW_SIZE - i]
                == "CG"
            ):
                obs_cpg += 1
            gc_perc, obs_exp, is_island = island_evaluation(gc_count, obs_cpg, MIN_WINDOW_SIZE)
            if is_island:
                print(f"island ending shrunk by {i}")
                break
        island_end = island_end - i
        print(f"Updated island after 1bp shrinking {island_start}, {island_end}")

        # step 4 shrink the whole island by 1bp until it meets the criteria
        print("- Stage 4 1bp two-fold shrinking")
        window = record[island_start:island_end]
        gc_count, obs_cpg = new_window_gcs(window)

        while True:
            gc_perc, obs_exp, is_island = island_evaluation(
                gc_count, obs_cpg, island_end - island_start
            )
            if is_island:
                break
            if record[island_start] in ["C", "G"]:
                gc_count -= 1
            if record[island_end] in ["C", "G"]:
                gc_count -= 1
            if record[island_start : island_start + 1] == ["CG"]:
                gc_count -= 1
            island_start += 1
            island_end -= 1
        print(f"updated island after shrinking: {island_start}, {island_end}")

        found_islands.append(
            Island(
                interval_begin=island_start,
                interval_end=island_end,
                gc_perc=gc_perc,
                obs_exp=obs_exp,
                win_length=island_end - island_start,
            )
        )
        print("AAA")
        record_position = island_end

    print("-- INITIAL SCAN COMPLETE")

    # step 5 at the end merge islands that are at least 100bp apart
    print("- Stage 5 merging islands")
    merged_islands = []

    found_islands.append(
        Island(
            interval_begin=189260, interval_end=189480, gc_perc=0.51, obs_exp=0.9, win_length=200
        )
    )
    found_islands.append(
        Island(
            interval_begin=189570, interval_end=189770, gc_perc=0.51, obs_exp=0.9, win_length=200
        )
    )
    found_islands.append(
        Island(
            interval_begin=189990, interval_end=190300, gc_perc=0.51, obs_exp=0.9, win_length=200
        )
    )

    if len(found_islands) >= 2:
        merged_islands.append(found_islands[0])

        for found_island in found_islands:
            if found_island.interval_begin - merged_islands[-1].interval_end < MIN_ISLAND_MERGE_GAP:
                print(
                    f"Merged Islands at {merged_islands[-1].interval_end}, {found_island.interval_end}"
                )
                merged_islands[-1].interval_end = found_island.interval_end
            else:
                merged_islands.append(found_island)
    return merged_islands


def extend_island_window_shifts(record, record_length, record_position):
    while True:
        real_window_length = min(MIN_WINDOW_SIZE, record_length - 1 - record_position)
        if real_window_length < 2:
            print("reached end of sequence in big jump")
            island_end = record_position
            break
        window = record[record_position:real_window_length]
        gc_count, obs_cpg = new_window_gcs(window)
        gc_perc, obs_exp, is_island = island_evaluation(gc_count, obs_cpg, real_window_length)
        record_position += real_window_length
        if not is_island:
            print("window is no longer an island")
            island_end = record_position
            break
    return island_end


def island_evaluation(gc_count, obs_cpg, window_size) -> (float, float, bool):
    is_island = False
    gc_perc = gc_count / window_size
    exp_cpg = ((gc_count / 2) ** 2) / window_size
    try:
        obs_exp = obs_cpg / exp_cpg
    except ZeroDivisionError:
        obs_exp = 0.0
    if gc_perc > MIN_GC_PERCENTAGE and obs_exp > MIN_OBSERVED_TO_EXPECTED_CPG:
        is_island = True
    return gc_perc, obs_exp, is_island


def shift_window_1bp(
    gc_count: int, obs_cpg: int, record: SeqRecord, window_position: int
) -> (int, int, int):
    window_position += 1
    window = record[window_position : window_position + MIN_WINDOW_SIZE]
    # update CGs
    if record[window_position - 1] in ["C", "G"]:
        gc_count -= 1
    if window[MIN_WINDOW_SIZE - 1] in ["C", "G"]:
        gc_count += 1
    # update CpGs
    if record[window_position - 1] == "C" and window[0] == "G":
        obs_cpg -= 1
    if window[MIN_WINDOW_SIZE - 2] == "C" and window[MIN_WINDOW_SIZE - 1] == "G":
        obs_cpg += 1
    return gc_count, obs_cpg, window_position


def main():
    record = get_sequence_from_ncbi(GENOME_NAME)
    print(f"record len: {len(record)}")
    # assert record_length >= MIN_WINDOW_SIZE  # TODO check this or maybe not needed

    start = time.time()
    islands = find_islands(record)  # time: 0.3477661609649658  ->
    end = time.time()

    for island in islands:
        print(island)

    print(f"Time: {end - start}")
    print("Done")
    # record = records[0].seq
    # print(record)

    # fasta_sequences = SeqIO.parse(response.content, 'fasta')


if __name__ == "__main__":
    main()
