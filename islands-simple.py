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


class Island:
    def __init__(
        self,
        interval_start: int,
        interval_end: int,
        gc_perc: float,
        obs_exp: float,
        win_length: int,
    ):
        self.interval_start = interval_start
        self.end = interval_end
        self.gc_perc = gc_perc
        self.obs_exp = obs_exp
        self.win_length = win_length

    def __str__(self):
        return f"Island at ({self.interval_start}, {self.end})\tgcper:{self.gc_perc}\tobs_exp:{self.obs_exp}\twin length:{self.win_length}"


def get_sequence_from_ncbi(genome_id: str) -> str:
    response = requests.get(
        "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=" + genome_id
    )
    genome = response.content.decode()
    if genome.startswith("Error"):
        raise Exception("Error, genome not found")

    fasta_io = StringIO(genome)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    return records[0].seq


def find_islands_simple_algorithm(record: SeqRecord) -> List[Island]:
    record_length = len(record)
    if record_length < MIN_WINDOW_SIZE:
        return []
    found_islands = []
    window_position = 0

    window = record[window_position : window_position + MIN_WINDOW_SIZE]
    gc_count = window.count("G") + window.count("C")
    obs_cpg = window.count("CG")
    while True:
        gc_perc = gc_count / MIN_WINDOW_SIZE
        exp_cpg = ((gc_count / 2) ** 2) / MIN_WINDOW_SIZE
        try:
            obs_exp = obs_cpg / exp_cpg
        except ZeroDivisionError:
            obs_exp = 0.0
        if gc_perc > MIN_GC_PERCENTAGE and obs_exp > MIN_OBSERVED_TO_EXPECTED_CPG:
            found_islands.append(
                Island(
                    interval_start=window_position,
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


record = get_sequence_from_ncbi(GENOME_NAME)
print(f"record len: {len(record)}")

start = time.time()
islands = find_islands_simple_algorithm(record)  # time: 0.43636059761047363
# islands = find_islands_algorithm(record)  # time: 0.3477661609649658  ->
end = time.time()

for island in islands:
    print(island)


print(f"Time: {end - start}")
print("Done")
# record = records[0].seq
# print(record)


# fasta_sequences = SeqIO.parse(response.content, 'fasta')
