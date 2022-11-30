#!/bin/bash

nProc=$1
nSamples=$2
aFile="a.txt"
bFile="b.txt"
cFile="c.txt"
partitionFile="p.txt"
FILE=$partitionFile

g++ -o generateSamples.exe generateSamples.cpp -O3
g++ -o computeSamples.exe computeSamples.cpp -O3

./generateSamples.exe $nSamples $aFile $bFile $cFile $nProc $partitionFile
for proc in $(seq 1 $nProc); do
    mkdir -p "pID${proc}" 
done
while read -r line; do
    echo "$line"
    ./computeSamples.exe $line $aFile $bFile &
done < "$FILE"
wait