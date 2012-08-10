#!/bin/bash

datestamp=$(date +%y%m%d)

pushd docs
pdflatex simrecomb_docs.tex > /dev/null
popd

bjam

zip -j simrecomb_$datestamp.zip bin/simrecomb bin/subsample_population bin/recombine_data docs/simrecomb_docs.pdf genetic_map_chr21_b36.txt

