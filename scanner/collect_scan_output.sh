#! /bin/bash
module load gnu/4.9.2
SIM_HEADER="pre_"

EXP_FIT_SOFTWARE=$HOME/bin/exponential_fit

#il tool seguente e` scan-columns in tools-Propaga
#serve per splittare i risultati collezionati dallo script su diversi files,
# in funzione di un parametro, per produrre piu` plots
SCANNER=$HOME/bin/scan-columns
DO_SCAN=true

#per il seguente, copiarsi dal prepare_scan_input la riga che genera tutti i valori (in questo caso di bulk lengths scannerizzate)
columns_values=$(awk 'BEGIN{for(i=2.0;i<=4.0;i+=1.0)print i}')
column=4
# 1:{PREPLASMA_LENGTH} 2:{DENSITY} 3:{RAMP_LENGTH} 4:{BULK_LENGTH} 5:{CONT_LENGTH} 
# 6:{PROTON_MAX_ENERGY} 7:{PROTON_TOT_ENERGY} 8:{PROTON_AVE_ENERGY} 9:{PROTON_TOT_NUMBER}" 
# 10:{FRONT_PROTON_AVE_ENERGY} 11:{FRONT_PROTON_TOT_NUMBER} 12:{REAR_PROTON_AVE_ENERGY} 13:{REAR_PROTON_TOT_NUMBER}

SPEC_TIME_TO_BE_READ=50
OUTPUT_FILE="energy_scan.txt"


SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))

rm -f ${OUTPUT_FILE}
touch ${OUTPUT_FILE}

printf "# PreplasmaLength \t Density \t RampLength \t BulkLength \t ContLength \t ${PARTICLE_TYPE}MaxEnergy \t       ${PARTICLE_TYPE}TotEnergy \t       ${PARTICLE_TYPE}AveEnergy \t       ${PARTICLE_TYPE}TotNumber \t   Front${PARTICLE_TYPE}AveEnergy \t   Front${PARTICLE_TYPE}TotNumber \t   Rear${PARTICLE_TYPE}AveEnergy \t   Rear${PARTICLE_TYPE}TotNumber \n" >> ${OUTPUT_FILE}
#printf "# PreplasmaLength \t Density \t RampLength \t BulkLength \t ${PARTICLE_TYPE}MaxEnergy \t ${PARTICLE_TYPE}TotEnergy \t ${PARTICLE_TYPE}AveEnergy \t ${PARTICLE_TYPE}TotNumber \t Sel${PARTICLE_TYPE}AveEnergy \t Sel${PARTICLE_TYPE}TotNumber \n" >> ${OUTPUT_FILE} 

for sim in "${SIMULATION_FOLDERS[@]}"
do
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $2}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $6}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $4}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $8}'))
 CONT_LENGTH=($(echo $sim | awk -F'_' '{print $10}'))
 cd $sim
 filename=`ls -1 OUTPUT/SPECTRUM_PROTcont_*${SPEC_TIME_TO_BE_READ}.000.dat`
 if [ -f $filename ]
 then
  PROTON_MAX_ENERGY=($(tail -1 $filename | awk '{print $2}'))
  aveData=($( ${EXP_FIT_SOFTWARE} -scan -piccante $filename  ))
 else
  PROTON_MAX_ENERGY="-1"
  aveData[0]="-1"
  aveData[1]="-1"
  aveData[2]="-1"
  aveData[3]="-1"
  aveData[4]="-1"
  aveData[5]="-1"
 fi
 if [ -f "OUTPUT/diag.dat" ]
 then
  PROTON_TOT_ENERGY=($(tail -1 OUTPUT/diag.dat | awk '{print $15}'))
 else
  PROTON_TOT_ENERGY="-1"
 fi
 PROTON_AVE_ENERGY=${aveData[0]}
 PROTON_TOT_NUMBER=${aveData[1]}
 FRONT_PROTON_AVE_ENERGY=${aveData[2]}
 FRONT_PROTON_TOT_NUMBER=${aveData[3]}
 REAR_PROTON_AVE_ENERGY=${aveData[4]}
 REAR_PROTON_TOT_NUMBER=${aveData[5]}

 cd ..
 #read -p "Press [Enter] key to continue..."
 printf '        %.1f                %.1f               %.2f             %.1f            %.2f               %.4f                  %.3e                  %s                     %s                      %s                     %s                      %s                     %s' "${PREPLASMA_LENGTH}" "${DENSITY}" "${RAMP_LENGTH}" "${BULK_LENGTH}" "${CONT_LENGTH}" "${PROTON_MAX_ENERGY}" "${PROTON_TOT_ENERGY}" "${PROTON_AVE_ENERGY}" "${PROTON_TOT_NUMBER}" "${FRONT_PROTON_AVE_ENERGY}" "${FRONT_PROTON_TOT_NUMBER}" "${REAR_PROTON_AVE_ENERGY}" "${REAR_PROTON_TOT_NUMBER}" >> ${OUTPUT_FILE}
 #printf '\t      %.1f \t%.1f\t     %.1f   \t%.1f\t %.1f\t%.2f\t%.3e \t%s\t%s\t%s\t%s' "${PREPLASMA_LENGTH}" "${DENSITY}" "${RAMP_LENGTH}" "${BULK_LENGTH}" "${CONT_LENGTH}" "${PROTON_MAX_ENERGY}" "${PROTON_TOT_ENERGY}" "${PROTON_AVE_ENERGY}" "${PROTON_TOT_NUMBER}" "${SEL_PROTON_AVE_ENERGY}" "${SEL_PROTON_TOT_NUMBER}" >> ${OUTPUT_FILE}
 printf '\n'  >> ${OUTPUT_FILE}

done

if [ ${DO_SCAN} ] 
then
 for value in ${columns_values}
 do
  $SCANNER -in energy_scan.txt -out energy_scan_$value.txt -select $column $value
 done
else 
 echo "output on ${OUTPUT_FILE}"
fi


