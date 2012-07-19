set title "Uniqueness of observed switches (relative)"
set xlabel "# observed switches (per 2kb)"
set ylabel "fraction of unique switches"
set label "Population 100000" at 500, .1

plot 'hist_2kb_0.txt' using 3:($5/$3) notitle, .5 notitle
