#!/usr/bin/env sh
# MFE primer checks for primer specificity against a user supplied database

## installing MFEprimer
# install dependancy
sudo apt install python-psutil
sudo pip install psutil --upgrade
# cd and download git repo
cd Software
git clone https://github.com/quwubin/MFEprimer.git
cd MFEprimer
# add to path
export PATH=~/Software/MFEprimer:$PATH
#echo 'export PATH=~/Software/MFEprimer:$PATH' >> ~/.bashrc


## RUNNING MFEprimer

# make the db - liberal setting
mv ../data/qpcr_specificity.fas ../temp/qpcr_specificity.fas
IndexDb.sh ../temp/qpcr_specificity.fas 5

# crab loop
for i in {0..11..1}
    do
        MFEprimer.py -k 5 -i ../temp/crab.p$i.fas -d ../temp/qpcr_specificity.fas > ../temp/crab.p$i.fas.results
    done
# shanny loop
for i in {0..11..1}
    do
        MFEprimer.py -k 5 -i ../temp/shanny.p$i.fas -d ../temp/qpcr_specificity.fas > ../temp/shanny.p$i.fas.results
    done

# clean results
grep -o "==> .*" ../temp/*.results | uniq | sed 's/==> //g' > ../temp/shanny-crab.MFE.results.csv
