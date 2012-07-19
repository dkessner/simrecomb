set title "Histogram of observed switches"
set xlabel "# observed switches (per 50kb window)"
set ylabel "# 50kb windows"
set label "Population 100000" at 22, 15
set label "Subsample 1000" at 22, 11 

plot 'switch_hist_50kb' using 1:2 with impulses notitle

