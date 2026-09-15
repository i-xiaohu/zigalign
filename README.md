# ZigAlign
This project explores the model of sequence alignment with tandem duplications, or DSI (Duplications, Substitutions, Indels) model. 
Tandem repeats are enriched in complex genomes, highly variable, and confuse downstream analysis.
Incorporating tandem copy events into sequence model will be instrumental for genomic studies.

`zigalign` 1) identify tandem repeats within a sequence through self-alignment dynamic programming and 2) compare two sequences
with channels between breakpoints where tandem duplications happen.

Run the command to compile it.
```
mkdir build; cd build
cmake ..; make -j8
```

Its usage is described below. The parameters should be tuned for sequences of different mutation/variation rate.
```
Usage: zigalign [options] seq1.fa seq2.fa > aln.paf
  Common options:
    -t [INT]  number of threads
    -v [STR]  intermediate results prefix
  Scoring parameters for pairwise alignment:
    -A [INT]  match score [1]
    -B [INT]  mismatch penalty [-120]
    -O [INT]  open gap(indel) penalty [-180]
    -E [INT]  extend gap penalty [-30]
    -D [INT]  repeat unit deletion penalty [-20]
  Scoring options for self-alignment:
    -u [INT]  minimum repeat unit size [100]
    -d [INT]  open tandem repeat penalty [-2]
    -p [INT]  close tandem repeat penalty [-6]
    -a [INT]  match score [2]
    -b [INT]  mismatch penalty [-3]
    -o [INT]  open gap(indel) penalty [-3]
    -e [INT]  extend gap penalty [-1]
Note: self-alignment scoring matrix must reward more and/or 
  penalize less than regular matrix to discover tandem repeats.
  Pairwise scoring matrix must penalize discrepancies harder as DSI models do.
```

# Update
`ZigAlign` now can align extreme-long duplication sequences, such as centromeres. It outputs alignment results in PAF
format. The CIGAR is extended with duplication information, where repeat units are given in square brackets. See the 
illustration below.

```
cg:Z:[339D][340M][339D][340D][340D][340D][339D][340D][339D][340D][340D][340D][339D][340D][340D]
[339D]: a repeat unit size of 339 is deleted
[340M]: a repeat unit size of 340 is matched
[340I]: a repeat unit size of 340 is inserted
```

The square brackets do not affect the synergy between aligned sequences, and it is totally ok to remove them when drawing
co-linear plots.