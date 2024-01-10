# CpG Islands

This is a prototype implementation of the [Takai and Jones’ algorithm](https://www.pnas.org/doi/10.1073/pnas.052410099) to find CpG islands. 

## Installation and start
Install `Python 3.10` or higher and run the following.

```
pip install -r requirements.txt
python islands.py
```

## Usage

The script in `islands.py` is controlled through constants at its beginning. The algorithm itself comes fully parametrized, including the nucleotide content of the islands (normally C and G in this order). The algorithm defaults are taken from its original specification.

Currently, the script only prints its output to the console, to store it on Linux distributions, run the program with the following command:
```
python islands.py > output.txt
```

#### Program parameters
* `NCBI_LINK` - Link to the NCBI file server
* `GENOME_NAME` - Name of the genome on the NCBI server, which will get downloaded and analyzed

#### Algorithm parameters

* `MIN_WINDOW_SIZE` - Minimal bp size of the window that can be considered an island (default `200`)
* `MIN_GC_PERCENTAGE` - Minimal required nucleotide content of preferred (default `0.5`)
* `MIN_OBSERVED_TO_EXPECTED_CPG` - Minimal required value of observed to expected CpGs (default `0.6`)
* `MIN_ISLAND_MERGE_GAP` - Minimal bp gap between two islands which will cause them to merge into one (default `100`)
* `FIRST_NUCLEOTIDE` - Fist nucleotide of the island (default `"C"`)
* `SECOND_NUCLEOTIDE` - Second nucleotide of the island (default `"G"`)

For comparison, there is also `islands-simple.py`, which contains the naive algorithm and is used more or less the same way that the main script. It only searches for windows of specified size and does not attempt to merge them. It is also slower (always O(n)), whereas the original, while also in O(n) runs a lot faster whenever there is an island larger than the minimal specified size.

## Notes

_The original specification of the algorithm is rather ambiguous, but this implementation attempts to stick to the specification as closely as possible. Many of the following notes are related to both this implementation and the original specification._

- If there are multiple records returned form the NCBI, the computation is only done for the first one.
- The algorithm, as per the specification drops small islands near the minimal window size in search of bigger ones. It can also shrink islands slightly. Consequently, its results are not comparable to the [naive implementation](https://www.bioinformatics.org/sms2/cpg_islands.html), which just finds islands of length 200, even if we merge these islands in case they are next to each other.
- Merging windows with small gaps in between may result in the larger windows not meeting the simple window criteria, but being windows because of the condition that connecting windows breached by a reasonable gap results in a window.
- Since this is a prototype, the output of the tool is not in any way "pretty" or to be directly used as a csv. 