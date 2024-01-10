# CpG Islands

This is a prototype implementation of the [Takai and Jones’ algorithm](https://www.pnas.org/doi/10.1073/pnas.052410099) to find CpG islands. 

## Usage

The scripts is controlled through constants at its beginning. The algorithm itself comes fully parametrized, including the nucleotide content of the islands (normally C and G in this order). The algorithm defaults are taken from its original specification.

#### Program parameters
* `NCBI_LINK` - Link to the NCBI file server
* `GENOME_NAME` - Name of the genome on the NCBI server, which will get downloaded and analyzed

#### Algorithm parameters

* `MIN_WINDOW_SIZE` - Minimal bp size of the window that can be considered an island (default `200`)
* `MIN_GC_PERCENTAGE` - Minimal required nucleotide content of prefered (default `0.5`)
* `MIN_OBSERVED_TO_EXPECTED_CPG` - Minimal required value of observed to expected CpGs (default `0.6`)
* `MIN_ISLAND_MERGE_GAP` - Minimal bp gap between two islands which will cause them to merge into one (default `100`)
* `FIRST_NUCLEOTIDE` - Fist nucleotide of the island (default `"C"`)
* `SECOND_NUCLEOTIDE` - Second nucleotide of the island (default `"G"`)

## Notes

- The original specification of the algorithm is rather ambiguous, but this implementation sticks to the specification in as literal sense as possible and all the following notes are related to the for this implementation and the original specification.
- If there are multiple records returned form the NCBI, the computation is only done for the first one.
- The algorithm, as per the specification drops small islands near the minimal window size in search of bigger ones. It can also shrink islands sligthly. Consequently, its results are not comparable to the [naive implementation](https://www.bioinformatics.org/sms2/cpg_islands.html), which just finds islands of length 200, even if we merge these islands in case they are next to each other.
- Merging windows with small gaps in between may result in the larger windows not meeting the simple window criteria, but being windows because of the condition that connecting windows breached by a reasonable gap results in a window. 