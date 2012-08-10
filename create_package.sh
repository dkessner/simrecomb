#!/bin/bash

datestamp=$(date +%y%m%d)

if [ $(uname) == "Darwin" ]; then
    platform=osx
else
    platform=linux
fi

filename_archive=simrecomb_${platform}_${datestamp}.zip

echo Building docs.
pushd docs
pdflatex simrecomb_docs.tex > /dev/null
popd

echo Building executables.
bjam

echo Creating $filename_archive.
zip -j $filename_archive bin/simrecomb bin/subsample_population bin/recombine_data docs/simrecomb_docs.pdf genetic_map_chr21_b36.txt

