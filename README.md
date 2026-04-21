# Zigalign
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
Usage: zigalign [options] seq1.fa seq2.fa
  Regular Scoring options:
    -A [INT]  match score [1]
    -B [INT]  mismatch penalty [-4]
    -O [INT]  open gap(indel) penalty [-6]
    -E [INT]  extend gap penalty [-1]
    -D [INT]  delete duplication penalty [-5]
  Scoring options for self-alignment:
    -u [INT]  minimum repeat unit size [5]
    -d [INT]  open tandem repeat penalty [-2]
    -p [INT]  close tandem repeat penalty [-6]
    -a [INT]  match score [2]
    -b [INT]  mismatch penalty [-3]
    -o [INT]  open gap(indel) penalty [-3]
    -e [INT]  extend gap penalty [-1]
Note: self-alignment scoring matrix must reward more and/or
  penalize less than regular matrix to discover tandem repeats.
```

# Work in progress
This program is still in development stage. I plan to add these features in future:
- Reduce memory consumption
- CIGAR interface
- Visualization
- Banded alignment
- SIMD acceleration