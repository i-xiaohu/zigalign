unialigner="/mnt/d/ClionProjects/software/unialigner-0.1/tandem_aligner/build/bin/tandem_aligner"
zigalign="/mnt/d/ClionProjects/zigalign/cmake-build-debug/zigalign"

cd /mnt/d/ClionProjects/zigalign/pairwise

for i in {1..22} X; do {
   pstat ${zigalign} -t8 -v ours/cen${i}/cen${i} centromeres/chm13_cen${i}.fa centromeres/hg002_Mcen${i}.fa # \
      # 2> ours/cen${i}/cen${i}.log
  break

   # Unialigner
#  ${unialigner} --first centromeres/chm13_cen${i}.fa --second centromeres/hg002_Mcen${i}.fa -o unialigner/cen${i}
#  mv unialigner/cen${i}/cigar.txt unialigner/cigar/cen${i}_cigar.txt
#  rm -r unialigner/cen${i}
#  python cigar2plot.py unialigner/cigar/cen${i}_cigar.txt unialigner/plots/cen${i}_plot.png

  # Zigalign
#  pstat ${zigalign} -t8 -v ours/cen${i}/cen${i} centromeres/chm13_cen${i}.fa centromeres/hg002_Mcen${i}.fa \
#       2> ours/cen${i}/cen${i}.log 1> ours/cen${i}/cen${i}.paf
#   python cigar2plot.py ours/cen${i}/cen${i}.paf ours/plots/cen${i}_plot.png

  # edlib
#  edlib-aligner -p -f CIG_EXT centromeres/chm13_cen${i}.fa centromeres/hg002_Mcen${i}.fa > edlib/cigar/cen${i}.log
#  python edlib/edlib_plot.py edlib/cigar/cen${i}.log edlib/plots/cen${i}.png

  # minimap2
#  minimap2 -a -asm -x asm20 -K 8G --eqx -r 500000,500000 \
#      centromeres/chm13_cen${i}.fa centromeres/hg002_Mcen${i}.fa > minimap2/cigar/cen${i}.sam
#  samtools view -F 260 -o minimap2/cigar/cen${i}_flt.sam minimap2/cigar/cen${i}.sam
#  mv minimap2/cigar/cen${i}_flt.sam minimap2/cigar/cen${i}.sam
#  python minimap2/mini_plot.py minimap2/cigar/cen${i}.sam minimap2/plots/cen${i}.png

  # merged plot
#  python sum_plot.py edlib/plots/cen${i}.png minimap2/plots/cen${i}.png \
#    unialigner/plots/cen${i}_plot.png ours/plots/cen${i}_plot.png \
#    -t centromere${i} --subtitles edlib minimap2 unialigner zigalign -o merged_plots/cen${i}.png
} done