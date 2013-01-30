set title "Uniqueness of observed recombination events (relative)"
set xlabel "# observed recombination events (per 2kb)"
set ylabel "fraction of unique events"
set label "Population 100000" at 2000, .1

plot 'hist_2kb_0.txt' using 2:($4/$2) notitle, .5 notitle
