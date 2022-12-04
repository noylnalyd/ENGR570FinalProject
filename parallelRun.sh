#!/bin/bash

nProc=$1
nSamples=$2
aFile="a.txt"
bFile="b.txt"
cFile="c.txt"
partitionFile="p.txt"
FILE=../$partitionFile

g++ -o generateSamples.exe generateSamplesExample.cpp -O3
g++ -o computeSamples.exe computeSamplesExample.cpp -O3

./generateSamples.exe $nSamples $aFile $bFile $cFile $nProc $partitionFile
rm -rf runs
mkdir -p runs
cd runs
pID=0
while read -r line; do
    ((pID = pID + 1))
    mkdir -p $pID 
    cp ../computeSamples.exe ./$pID
    cd $pID
    echo $line
    ./computeSamples.exe $line ../../$aFile ../../$bFile $pID > textout.txt &
    cd ..
done < "$FILE"
wait

cd ..

python3 mcAnalysis.py $nProc $nSamples 4