import requests
from typing import List
from io import StringIO
from Bio import SeqIO, SeqRecord, Seq
import time
import logging

NCBI_LINK = "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id="

# genome name to download and analyze from NCBI
GENOME_NAME = "DQ011153.1"  # Monkeypox, 0.44 sec
# GENOME_NAME = "NC_060948.1"  # Human, 125 sec

# algorithm parameters explained in README.md
MIN_WINDOW_SIZE = 200
MIN_GC_PERCENTAGE = 0.5
MIN_OBSERVED_TO_EXPECTED_CPG = 0.6
MIN_ISLAND_MERGE_GAP = 100
FIRST_NUCLEOTIDE = "C"
SECOND_NUCLEOTIDE = "G"


class Window:
    """Sliding window of the CpG island search algorithm."""

    def __init__(
        self,
        record: SeqRecord,
        window_begin: int,
        window_size: int,
    ):
        self.record = record
        self.record_len = len(record)
        self.window_begin = window_begin
        self.window_size = window_size
        self.window_end = window_begin + window_size
        self.gc_count = 0
        self.obs_cpg = 0
        self.gc_perc = 0
        self.obs_exp = 0
        self.update_gc_count()

    def __str__(self):
        return f"Window at ({self.window_begin}, {self.window_end})\tgcper:{self.gc_perc}\tobs_exp:{self.obs_exp}\twin_length:{self.window_size}"

    def is_island(self) -> bool:
        """Evaluate the current window and return True if it fits the definition of an island."""
        self._evaluate()
        return self.gc_perc > MIN_GC_PERCENTAGE and self.obs_exp > MIN_OBSERVED_TO_EXPECTED_CPG

    def _evaluate(self):
        self.gc_perc = self.gc_count / self.window_size
        exp_cpg = ((self.gc_count / 2) ** 2) / self.window_size
        try:
            self.obs_exp = self.obs_cpg / exp_cpg
        except ZeroDivisionError:
            self.obs_exp = 0.0

    def update_gc_count(self):
        """Count the nucleotide content in the current window and reevaluate relevant statistics."""
        record_cut = self.record[self.window_begin : self.window_end]
        self.gc_count = record_cut.count(FIRST_NUCLEOTIDE) + record_cut.count(SECOND_NUCLEOTIDE)
        self.obs_cpg = record_cut.count(f"{FIRST_NUCLEOTIDE}{SECOND_NUCLEOTIDE}")
        self._evaluate()

    def shrink_left(self):
        """Shrinks window by 1 bp from the 5' end."""
        if self.record[self.window_begin] in [FIRST_NUCLEOTIDE, SECOND_NUCLEOTIDE]:
            self.gc_count -= 1
        if (
            self.record[self.window_begin] == FIRST_NUCLEOTIDE
            and self.record[self.window_begin + 1] == SECOND_NUCLEOTIDE
        ):
            self.obs_cpg -= 1
        self.window_begin += 1

    def expand_right(self) -> bool:
        """Expands window by 1 bp to 3' end,

        :return: False if the record end is reached.
        """
        if self.window_end >= self.record_len:
            return False
        if self.record[self.window_end] in [FIRST_NUCLEOTIDE, SECOND_NUCLEOTIDE]:
            self.gc_count += 1
        if (
            self.record[self.window_end - 1] == FIRST_NUCLEOTIDE
            and self.record[self.window_end] == SECOND_NUCLEOTIDE
        ):
            self.obs_cpg += 1
        self.window_end += 1
        return True

    def shrink_right(self):
        """Shrinks window by 1 bp from the 3' end."""
        if self.record[self.window_end - 1] in [FIRST_NUCLEOTIDE, SECOND_NUCLEOTIDE]:
            self.gc_count -= 1
        if (
            self.record[self.window_end - 2] == FIRST_NUCLEOTIDE
            and self.record[self.window_end - 1] == SECOND_NUCLEOTIDE
        ):
            self.obs_cpg -= 1
        self.window_end -= 1

    def expand_left(self):
        """Expands window by 1 bp to 5' end."""
        if self.record[self.window_begin - 1] in [FIRST_NUCLEOTIDE, SECOND_NUCLEOTIDE]:
            self.gc_count += 1
        if (
            self.record[self.window_begin - 1] == FIRST_NUCLEOTIDE
            and self.record[self.window_begin] == SECOND_NUCLEOTIDE
        ):
            self.obs_cpg += 1
        self.window_begin -= 1

    def expand_right_island(self) -> bool:
        """Expands window by minimal window size to the 3' end. It creates new window, checks if
        it's an island and joins it to this island if so.

        :return:  True if this window was extended, false otherwise or at the end of record.
        """
        if not (remaining_window_length := self.record_len - self.window_end):
            return False
        if remaining_window_length < MIN_WINDOW_SIZE:
            logging.debug("reached end of sequence in big jump, be careful")
            new_window_begin = self.record_len - MIN_WINDOW_SIZE
        else:
            new_window_begin = self.window_end
        next_window = Window(self.record, new_window_begin, MIN_WINDOW_SIZE)
        self.window_end = new_window_begin + MIN_WINDOW_SIZE
        self.window_size = self.window_end - self.window_begin
        return next_window.is_island()

    def rollback_until_island(self):
        """Rollback the last sub-window of minimal window size to the 5' end until the sub-window
        is an island.
        """
        last_window = Window(self.record, self.window_end - MIN_WINDOW_SIZE, MIN_WINDOW_SIZE)
        while not last_window.is_island():
            last_window.expand_left()
            last_window.shrink_right()
        self.join_window(last_window)

    def join_window(self, other_window):
        """Expands this window so that its end is and the end of the other_window. Recalculates the
        nucleotide content of the merged window.
        """
        self.window_end = other_window.window_end
        self.window_size = self.window_end - self.window_begin
        self.update_gc_count()


def get_sequence_from_ncbi(genome_id: str) -> Seq:
    """Retrieves record from NCBI file server given its name.

    :raises: Exception if the genome is not found on the server.
    """
    print("Pulling data from NCBI")
    response = requests.get(NCBI_LINK + genome_id)
    genome = response.content.decode()
    if genome.startswith("Error"):
        raise Exception("Error, genome not found")

    fasta_io = StringIO(genome)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    return records[0].seq


def find_island_1bp_shifts(window: Window) -> bool:
    """Searches the record in 1bp increments until the first island is found.

    :return: False if end of record is reached before island is found, True otherwise.
    """
    while True:
        if window.is_island():
            logging.debug(f"island found {window.window_begin, window.window_end}")
            return True
        if not window.expand_right():
            logging.debug("reached the end of record")
            return False
        window.shrink_left()


def two_fold_shrinking(window: Window) -> bool:
    """Shrinks big windows from both sides until they meet island criteria.

    :return: True if shrunk window is an island, False if shrinking results in a window smaller
    than minimal window size.
    """
    while not window.is_island():
        window.shrink_right()
        window.shrink_left()
        window.window_size -= 2
        if window.window_size < MIN_WINDOW_SIZE:
            logging.debug("! Shrinking resulted in window smaller than minimal window size.")
            return False
    return True


def merge_islands(found_islands: List[Window]) -> List[Window]:
    """Merge islands in the :param: found_islands that are close to each other."""
    merged_islands = [found_islands[0]]
    if len(found_islands) >= 2:
        for found_island in found_islands[1:]:
            if found_island.window_begin - merged_islands[-1].window_end < MIN_ISLAND_MERGE_GAP:
                logging.debug(
                    f"Merged Islands at {merged_islands[-1].window_begin}, "
                    f"{merged_islands[-1].window_end}, {found_island.window_begin}, "
                    f"{found_island.window_end}"
                )
                merged_islands[-1].join_window(found_island)
            else:
                merged_islands.append(found_island)
    return merged_islands


def extend_island_window_shifts(window: Window):
    """Expands window by minimal window size increments to the 3' end as long as the added windows
    are islands.
    """
    while True:
        extended = window.expand_right_island()
        if not extended:
            break

    window.update_gc_count()
    logging.debug(f"island: {window.window_begin}, {window.window_end}")


def find_islands(record: SeqRecord) -> List[Window]:
    """Finds CpG islands in the given record.

    :return: List of found islands
    """
    found_islands = []
    record_position = 0
    while True:
        window = Window(record, record_position, MIN_WINDOW_SIZE)

        # step 1 find window that is an island
        logging.debug("- Stage 1: 1bp shifts")
        record_continues = find_island_1bp_shifts(window)
        if not record_continues:
            break

        # step 2 window-length shifts
        logging.debug("- Stage 2: window-length shifts")
        extend_island_window_shifts(window)

        # step 3 shift the last window by 1bp until it meets the criteria
        logging.debug("- Stage 3: 1bp rollback")
        window.rollback_until_island()
        logging.debug(
            f"Updated island after 1bp shrinking {window.window_begin}, {window.window_end}"
        )

        # step 4 shrink the whole island by 1bp until it meets the criteria
        logging.debug("- Stage 4 1bp two-fold shrinking")
        if not two_fold_shrinking(window):
            record_position = window.window_begin + 1
            continue
        logging.debug(f"Updated island after shrinking: {window.window_begin}, {window.window_end}")

        found_islands.append(window)
        record_position = window.window_end

    # step 5 at the end merge islands that are at least 100bp apart
    logging.debug("- Stage 5 merging islands")
    if not found_islands:
        return []

    merged_islands = merge_islands(found_islands)
    return merged_islands


def main():
    # logging.debug = print  # uncomment for debug logs
    record = get_sequence_from_ncbi(GENOME_NAME)

    print(f"record len: {len(record)}")

    start = time.time()
    islands = find_islands(record)
    end = time.time()

    for island in islands:
        print(island)

    logging.debug(f"Time: {end - start}")
    print("Done")


if __name__ == "__main__":
    main()
