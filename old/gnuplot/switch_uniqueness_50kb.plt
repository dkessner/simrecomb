set title "Uniqueness of observed switches"
set xlabel "# observed switches (per 50kb)"
set ylabel "# unique switches"
set label "Population 100000" at 22, 10
set label "Subsample 1000" at 22, 8

plot 'hist_50kb_0.txt' using 3:5 notitle

