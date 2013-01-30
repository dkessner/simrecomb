set title "Uniqueness of observed switches"
set xlabel "# observed switches (per 2kb)"
set ylabel "# unique switches"
set label "Population 100000" at 80, 10

plot 'hist_2kb_0.txt' using 3:5 notitle

