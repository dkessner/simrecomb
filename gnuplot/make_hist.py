#!/usr/bin/env python



event_hist = {}
switch_hist = {}


for line in open('hist_50kb_0.txt'):
    tokens = line.split()
    if tokens[0] == '#': continue
    position = tokens[0]
    events = int(tokens[1])
    switches = int(tokens[2])
    events_uniq = tokens[3]
    switches_uniq = tokens[4]

    if events in event_hist:
        event_hist[events] = event_hist[events] + 1
    else:
        event_hist[events] = 1

    if switches in switch_hist:
        switch_hist[switches] = switch_hist[switches] + 1
    else:
        switch_hist[switches] = 1


for key in switch_hist:
    print key, switch_hist[key]

