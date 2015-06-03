#!/bin/bash

oIFS=$IFS
declare -a FOLDERS
FOLDERS=`find . -mindepth 1 -maxdepth 1 -type d`
tempo_max=0
for folder in $FOLDERS
do
 declare -a opic=`tail -2 $folder/opic.txt |head -1`
 IFS=' '
 read -ra tempi <<< "$opic"
 tempo=${tempi[4]}
 IFS=$oIFS
 if [[ $tempo -gt $tempo_max ]]
 then 
  tempo_max=$tempo
 fi
 echo $folder : $tempo
done
echo
tempo_max_h=`echo " $tempo_max / 3600 " | bc -l`
printf "Slowest required %d seconds, equivalent to %3.1f hours\n" "$tempo_max" "$tempo_max_h"


