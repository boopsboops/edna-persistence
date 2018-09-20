#!/usr/bin/env sh
# designing primers with Primer3

# installing
sudo apt-get install primer3

# ask for help
man primer3_core
# also 
# http://primer3.sourceforge.net/primer3_manual.htm


# to run primer3
# need to make a iput file first though (example below)
# with readable output on screen
primer3_core -format_output ../data/shanny.p3
primer3_core -format_output ../data/crab.p3
# with standard output to disk
primer3_core -output="../temp/shanny.p3.so.out" ../data/shanny.p3
primer3_core -output="../temp/crab.p3.so.out" ../data/crab.p3

# run a for loop to convert primer3 results to MFEprimer input
# need the curly braces for the var, and double quotes for the grep to not interpret $ as a regex
# 0..4..1 # from 0 to 4 by 1
for i in {0..11..1}
    do
        grep -e "PRIMER_LEFT_${i}_SEQUENCE=" -e "PRIMER_RIGHT_${i}_SEQUENCE=" ../temp/shanny.p3.so.out > ../temp/shanny.p$i.fas
        sed -i -e 's/^/>/g' -e 's/=/\n/g' ../temp/shanny.p$i.fas
    done

# crabs ...
for i in {0..11..1}
    do
        grep -e "PRIMER_LEFT_${i}_SEQUENCE=" -e "PRIMER_RIGHT_${i}_SEQUENCE=" ../temp/crab.p3.so.out > ../temp/crab.p$i.fas
        sed -i -e 's/^/>/g' -e 's/=/\n/g' ../temp/crab.p$i.fas
    done