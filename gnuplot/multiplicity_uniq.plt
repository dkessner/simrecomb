set title "Multiple observations of identical-by-descent recombination events"
set xlabel "Number of observations"
set ylabel "Fraction of unique recombination events"
set label "Population 100000" at 12, .3
set label "Subsample 1000" at 12, .25
plot [0:15] 'events_0.txt' using 1:3 with impulses t "ch0",\
 'events_0.txt' using ($1+.1):7 with impulses t "ch0 switch",\
 'events_1.txt' using ($1+.2):3 with impulses t "ch1",\
 'events_1.txt' using ($1+.3):7 with impulses t "ch1 switch",\
 'events_2.txt' using ($1+.4):3 with impulses t "ch2",\
 'events_2.txt' using ($1+.5):7 with impulses t "ch2 switch"
