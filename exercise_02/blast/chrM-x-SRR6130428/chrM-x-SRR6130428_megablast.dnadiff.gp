set terminal png small
set output "/Users/test/Desktop/UB/Advanced_Bioinformatics/exercise_02/blast/chrM-x-SRR6130428/chrM-x-SRR6130428_megablast.dnadiff.png"
set ytics ( \
 "7984" 1, \
 "*6342" 3453, \
 "8566" 3755, \
 "*6606" 10232, \
 "6882" 10602, \
 "*7358" 11118, \
 "9370" 12396, \
 "*7700" 29858, \
 "*7106" 32142, \
 "6994" 32953, \
 "*7714" 33616, \
 "8350" 35921, \
 "*5950" 41028, \
 "*8054" 41257, \
 "6240" 44943, \
 "*7312" 45222, \
 "8170" 46379, \
 "*6642" 50573, \
 "*8516" 50959, \
 "8502" 57088, \
 "8530" 63160, \
 "*6812" 69403, \
 "*8848" 69872, \
 "7596" 78856, \
 "6756" 80778, \
 "*7570" 81222, \
 "7736" 83074, \
 "7420" 85475, \
 "3047" 86917, \
 "6376" 87005, \
 "9644" 87312, \
 "6650" 117850, \
 "2717" 118242, \
 "2569" 118323, \
 "1401" 118402, \
 "3567" 118469, \
 "5652" 118572, \
 "727" 118769, \
 "2265" 118833, \
 "1879" 118908, \
 "2349" 118979, \
 "2551" 119055, \
 "187" 119134, \
 "2377" 119197, \
 "6628" 119273, \
 "" 119699 \
)
set size 1,1
set grid
set nokey
set border 10
set tics scale 0
set xlabel "Scer_chrM"
set ylabel "QRY"
set format "%.0f"
set xrange [1:85779]
set yrange [1:119699]
set linestyle 1  lt 1 lw 3 pt 6 ps 1
set linestyle 2  lt 3 lw 3 pt 6 ps 1
set linestyle 3  lt 2 lw 3 pt 6 ps 1
plot \
 "/Users/test/Desktop/UB/Advanced_Bioinformatics/exercise_02/blast/chrM-x-SRR6130428/chrM-x-SRR6130428_megablast.dnadiff.fplot" title "FWD" w lp ls 1, \
 "/Users/test/Desktop/UB/Advanced_Bioinformatics/exercise_02/blast/chrM-x-SRR6130428/chrM-x-SRR6130428_megablast.dnadiff.rplot" title "REV" w lp ls 2
