set title "Multiple observations of identical-by-descent recombination events"
set xlabel "Number of observations"
set ylabel "Fraction of observed recombination events"
set label "Population 100000" at 12, .2
set label "Subsample 3500" at 12, .15
plot [0:15] 'events_0.txt' using 1:5 with impulses t "ch0",\
 'events_0.txt' using ($1+.1):9 with impulses t "ch0 switch",\
 'events_1.txt' using ($1+.2):5 with impulses t "ch1",\
 'events_1.txt' using ($1+.3):9 with impulses t "ch1 switch",\
 'events_2.txt' using ($1+.4):5 with impulses t "ch2",\
 'events_2.txt' using ($1+.5):9 with impulses t "ch2 switch"
