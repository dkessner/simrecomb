set title "Hotspots and Observed Switches"
set xlabel "chromosome position"
set ylabel "# switches (per 2kb)"

plot [27200000:27700000] 'hist_2kb_0.txt' using 1:3 with lines t "observed switches",\
'hist_2kb_0.txt' using 1:5 with lines t "unique switches",\
'../genetic_map_chr21_b36.txt' using 1:($2*4) with lines t "HapMap genetic map (arbitrarily scaled)"
