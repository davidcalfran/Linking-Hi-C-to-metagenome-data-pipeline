#!/bin/bash

out=merged_tophits_mash.tab
# delete the file prior to doing concatenation
# or if ran twice it would be counted in the input files!
rm -f "$out"

for f in *.fna.tab.tab.tab
do
   if [ -s "$f" ] ; then
      #cat "$f" | sed 's/^/$f,/'  # cat+sed is too much here
      sed "s/^/$f,/" "$f"
   else
      echo "$f,"
   fi
done > $out