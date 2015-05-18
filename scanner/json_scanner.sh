#! /bin/bash

#*******************************************************************************
# This file is part of tools_pic.
#
# tools_pic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# tools-pic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with tools-pic.  If not, see <http://www.gnu.org/licenses/>.
#*******************************************************************************



JOB_FILE=galileo-64.cmd
JSON_FILE=inputPiccante.json


preplasmas=$(awk 'BEGIN{for(i=1.0;i<=3.0;i+=1.0)print i}')
#preplasmas=0.0

#densities=$(awk 'BEGIN{for(i=0.5;i<=3.0;i+=0.5)print i}')
densities=1.0

#ramps=$(awk 'BEGIN{for(i=0.25;i<=0.75;i+=0.25)print i}')
ramps=0.75

centrals=$(awk 'BEGIN{for(i=2.0;i<=4.0;i+=1.0)print i}')
#centrals=$(awk 'BEGIN{for(i=2.0;i<=10.0;i+=1.0)print i}')
#centrals=2.4

#contams=$(awk 'BEGIN{for(i=0.05;i<=0.1;i+=0.01)print i}')
contams=0.08


for pre in $preplasmas
do
for dens in $densities
do
for ramp in $ramps
do
for central in $centrals
do
for contam in $contams
do

echo "pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam}" >> sim_da_fare.txt

mkdir pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam}
cd pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam}
cp ../${JOB_FILE} .

rm -f ${JSON_FILE}
touch ${JSON_FILE}

DO_DUMP="true"
DO_DUMP_EVERY=50.0

NCPU=64
DIMENSIONI=2
NUMERO_PUNTI_GRIGLIA_X=7168
NUMERO_PUNTI_GRIGLIA_Y=3584
NUMERO_PUNTI_GRIGLIA_Z=1
TMAX=50.0

XMIN=0.0
XMAX=71.68
YMIN=-35.84
YMAX=35.84
ZMIN=-1.0
ZMAX=1.0

COURANT_FRIEDRICHS_LEWY_PARAMETER=0.8

BOUNDARY_X="open"
BOUNDARY_Y="open"
BOUNDARY_Z="periodic"

LASER_POLARIZATION="P"
POSIZIONE_INIZIALE_PICCO_IMPULSO_LASER=16.51
POSIZIONE_FUOCO_LASER=33.02
LUNGHEZZA_LASER_FWHM=16.5
WAIST_LASER=6.2
PARAMETRO_ADIMENSIONALE_LASER_A0=3.0
LUNGHEZZA_ONDA_LASER=0.8

XMIN_CLOUD=$(echo ${POSIZIONE_FUOCO_LASER} + 0.01 | bc)
#XMIN_CLOUD=$(echo ${POSIZIONE_FUOCO_LASER} | bc)
XMAX_CLOUD=$(echo ${XMIN_CLOUD} + ${pre} | bc)
YMIN_CLOUD=$YMIN
YMAX_CLOUD=$YMAX
ZMIN_CLOUD=$ZMIN
ZMAX_CLOUD=$ZMAX
DENSITY_CLOUD=${dens}
NUMERO_IONIZZAZIONE_CLOUD=10
NUMERO_ATOMICO_CLOUD=27

XMIN_BULK=${XMAX_CLOUD}
XMAX_BULK=$(echo ${XMIN_BULK} + ${ramp} + ${central} | bc)
YMIN_BULK=$YMIN
YMAX_BULK=$YMAX
ZMIN_BULK=$ZMIN
ZMAX_BULK=$ZMAX
DENSITY_BULK=100.0
NUMERO_IONIZZAZIONE_BULK=10
NUMERO_ATOMICO_BULK=27

XMIN_CONT=${XMAX_BULK}
XMAX_CONT=$(echo ${XMIN_CONT} + ${contam} | bc)
YMIN_CONT=$YMIN
YMAX_CONT=$YMAX
ZMIN_CONT=$ZMIN
ZMAX_CONT=$ZMAX
DENSITY_CONT=10.0
NUMERO_IONIZZAZIONE_CONT=1
NUMERO_ATOMICO_CONT=1

NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE=4

NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE=12
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE=4
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE=9
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE=3

NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE=4
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE=4
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE=4
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE=4



### DO NOT TOUCH FROM HERE ###

printf '{\n' >> ${JSON_FILE}
printf '  \"VERSION\": 2,\n' >> ${JSON_FILE}
printf '  \"dimensions\": %s,\n' "${DIMENSIONI}" >> ${JSON_FILE}
printf '  \"masterProc\": 0,\n' >> ${JSON_FILE}
printf '  \"OutputPath\": \"OUTPUT\",\n' >> ${JSON_FILE}
printf '  \"radiationFriction\": false,\n' >> ${JSON_FILE}
printf '  \"courantFactor\": %s,\n' "${COURANT_FRIEDRICHS_LEWY_PARAMETER}" >> ${JSON_FILE}
printf '  \"nProcY\": %s,\n' "${NCPU}" >> ${JSON_FILE}
printf '  \"xRange\": [%s, %s],\n' "$XMIN" "$XMAX" >> ${JSON_FILE}
printf '  \"yRange\": [%s, %s],\n' "$YMIN" "$YMAX" >> ${JSON_FILE}
printf '  \"zRange\": [%s, %s],\n' "$ZMIN" "$ZMAX" >> ${JSON_FILE}
printf '  \"NCells\": [%s, %s, %s],\n' "${NUMERO_PUNTI_GRIGLIA_X}" "${NUMERO_PUNTI_GRIGLIA_Y}" "${NUMERO_PUNTI_GRIGLIA_Z}" >> ${JSON_FILE}
printf '  \"simulationTime\": %s,\n' "$TMAX" >> ${JSON_FILE}
printf '  \"boundaries\": [\"%s\", \"%s\", \"%s\"],\n' "${BOUNDARY_X}" "${BOUNDARY_Y}" "${BOUNDARY_Z}" >> ${JSON_FILE}
printf '  \"StretchedGrid\":{\n' >> ${JSON_FILE}
printf '    \"enabled\": false,\n' >> ${JSON_FILE}
printf '    \"x\":{\"left\": { \"limit\": -10.0, \"NCells\": 1000 },\n' >> ${JSON_FILE}
printf '      \"right\": { \"limit\":  10.0, \"NCells\": 1000 } },\n' >> ${JSON_FILE}
printf '    \"y\":{\"left\": { \"limit\": -15.0, \"NCells\": 1000 },\n' >> ${JSON_FILE}
printf '      \"right\": { \"limit\":  15.0, \"NCells\": 1000 } },\n' >> ${JSON_FILE}
printf '    \"z\":{\"left\": { \"limit\": -15.0, \"NCells\": 1000 },\n' >> ${JSON_FILE}
printf '      \"right\": { \"limit\":  15.0, \"NCells\": 1000 } }\n' >> ${JSON_FILE}
printf '  },\n' >> ${JSON_FILE}
printf '  \"MovingWindow\":{\n' >> ${JSON_FILE}
printf '    \"enabled\": false,\n' >> ${JSON_FILE}
printf '    \"start\": 0,\n' >> ${JSON_FILE}
printf '    \"NO_beta\": 1,\n' >> ${JSON_FILE}
printf '    \"NO_frequency\":20\n' >> ${JSON_FILE}
printf '  },\n' >> ${JSON_FILE}
printf '  \"restart\":{\n' >> ${JSON_FILE}
printf '    \"DumpPath\": \"DUMP\",\n' >> ${JSON_FILE}
printf '    \"doRestart\": false,\n' >> ${JSON_FILE}
printf '    \"restartFromDump\": 1,\n' >> ${JSON_FILE}
printf '    \"doDump\": %s,\n' "${DO_DUMP}" >> ${JSON_FILE}
printf '    \"dumpEvery\": %s\n' "${DO_DUMP_EVERY}" >> ${JSON_FILE}
printf '  },\n' >> ${JSON_FILE}
printf '  \"Laser\":[\n' >> ${JSON_FILE}
printf '    {\n' >> ${JSON_FILE}
printf '      \"enabled\": true,\n' >> ${JSON_FILE}
printf '      \"type\": \"GAUSSIAN\",\n' >> ${JSON_FILE}
printf '      \"polarization\": \"%s\",\n' "${LASER_POLARIZATION}" >> ${JSON_FILE}
printf '      \"durationFWHM\": %s,\n' "${LUNGHEZZA_LASER_FWHM}" >> ${JSON_FILE}
printf '      \"initialPosition\": %s,\n' "${POSIZIONE_INIZIALE_PICCO_IMPULSO_LASER}" >> ${JSON_FILE}
printf '      \"waist\": %s,\n' "${WAIST_LASER}">> ${JSON_FILE}
printf '      \"focusPosition\": %s,\n' "${POSIZIONE_FUOCO_LASER}">> ${JSON_FILE}
printf '      \"a\": %s,\n' "${PARAMETRO_ADIMENSIONALE_LASER_A0}" >> ${JSON_FILE}
printf '      \"lambda\": %s,\n' "${LUNGHEZZA_ONDA_LASER}" >> ${JSON_FILE}
printf '      \"rotation\": false,\n' >> ${JSON_FILE}
printf '      \"angle\": 10.0,\n' >> ${JSON_FILE}
printf '      \"center\": 0,\n' >> ${JSON_FILE}
printf '      \"riseTime\": 2\n' >> ${JSON_FILE}
printf '    }\n' >> ${JSON_FILE}
printf '  ],\n' >> ${JSON_FILE}
printf '  \"Plasma\":[\n' >> ${JSON_FILE}
printf '    {\n' >> ${JSON_FILE}
printf '      \"name\": \"cloud\",\n' >> ${JSON_FILE}
printf '      \"densityFunction\": \"box\",\n' >> ${JSON_FILE}
printf '      \"XRangeBox\" : [%s, %s],\n' "${XMIN_CLOUD}" "${XMAX_CLOUD}" >> ${JSON_FILE}
printf '      \"YRangeBox\" : [%s, %s],\n' "${YMIN_CLOUD}" "${YMAX_CLOUD}" >> ${JSON_FILE}
printf '      \"ZRangeBox\" : [%s, %s],\n' "${ZMIN_CLOUD}" "${ZMAX_CLOUD}" >> ${JSON_FILE}
printf '      \"DensityCoefficient\" : %s,\n' "${DENSITY_CLOUD}" >> ${JSON_FILE}
printf '      \"DensityLambda\" : %s\n' "${LUNGHEZZA_ONDA_LASER}" >> ${JSON_FILE}
printf '    },{\n' >> ${JSON_FILE}
printf '      \"name\": \"bulk\",\n' >> ${JSON_FILE}
printf '      \"densityFunction\": \"left_fixed_exp_ramp\",\n' >> ${JSON_FILE}
printf '      \"XRangeBox\" : [%s, %s],\n' "${XMIN_BULK}" "${XMAX_BULK}" >> ${JSON_FILE}
printf '      \"YRangeBox\" : [%s, %s],\n' "${YMIN_BULK}" "${YMAX_BULK}" >> ${JSON_FILE}
printf '      \"ZRangeBox\" : [%s, %s],\n' "${ZMIN_BULK}" "${ZMAX_BULK}" >> ${JSON_FILE}
printf '      \"DensityCoefficient\" : %s,\n' "${DENSITY_BULK}" >> ${JSON_FILE}
printf '      \"DensityLambda\" : %s,\n' "${LUNGHEZZA_ONDA_LASER}" >> ${JSON_FILE}
printf '      \"RampMinDensity\" : %s,\n' "${DENSITY_CLOUD}" >> ${JSON_FILE}
printf '      \"LeftRampLength\" : %s\n' "${ramp}" >> ${JSON_FILE}
printf '    },{\n' >> ${JSON_FILE}
printf '      \"name\": \"cont\",\n' >> ${JSON_FILE}
printf '      \"densityFunction\": \"box\",\n' >> ${JSON_FILE}
printf '      \"XRangeBox\" : [%s, %s],\n' "${XMIN_CONT}" "${XMAX_CONT}" >> ${JSON_FILE}
printf '      \"YRangeBox\" : [%s, %s],\n' "${YMIN_CONT}" "${YMAX_CONT}" >> ${JSON_FILE}
printf '      \"ZRangeBox\" : [%s, %s],\n' "${ZMIN_CONT}" "${ZMAX_CONT}" >> ${JSON_FILE}
printf '      \"DensityCoefficient\" : %s,\n' "${DENSITY_CONT}" >> ${JSON_FILE}
printf '      \"DensityLambda\" : %s\n' "${LUNGHEZZA_ONDA_LASER}" >> ${JSON_FILE}
printf '    }\n' >> ${JSON_FILE}
printf '  ],\n' >> ${JSON_FILE}
printf '  \"Species\":[\n' >> ${JSON_FILE}
printf '    {\n' >> ${JSON_FILE}
printf '      \"enabled\": true,\n' >> ${JSON_FILE}
printf '      \"name\": \"ELEcloud\",\n' >> ${JSON_FILE}
printf '      \"plasma\": \"cloud\",\n' >> ${JSON_FILE}
printf '      \"ParticlesPerCell\": [%s,%s,%s],\n' "${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE}" "${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE}" "${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE}" >> ${JSON_FILE}
printf '      \"type\": \"ELECTRON\",\n' >> ${JSON_FILE}
printf '      \"isMarker\": 0,\n' >> ${JSON_FILE}
printf '      \"isTest\": false,\n' >> ${JSON_FILE}
printf '      \"distribution\": \"Maxwell\",\n' >> ${JSON_FILE}
printf '      \"distributionParams\": [6.0e-4],\n' >> ${JSON_FILE}
printf '      \"distributionAddMomentum\": [0.0,0.0,0.0]\n' >> ${JSON_FILE}
printf '    },\n' >> ${JSON_FILE}
printf '    {\n' >> ${JSON_FILE}
printf '      \"enabled\": true,\n' >> ${JSON_FILE}
printf '      \"name\": \"ELEbulk\",\n' >> ${JSON_FILE}
printf '      \"plasma\": \"bulk\",\n' >> ${JSON_FILE}
printf '      \"ParticlesPerCell\": [%s,%s,%s],\n' "${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE}" "${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE}" "${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE}" >> ${JSON_FILE}
printf '      \"type\": \"ELECTRON\",\n' >> ${JSON_FILE}
printf '      \"isMarker\": 0,\n' >> ${JSON_FILE}
printf '      \"isTest\": false,\n' >> ${JSON_FILE}
printf '      \"distribution\": \"Maxwell\",\n' >> ${JSON_FILE}
printf '      \"distributionParams\": [6.0e-4],\n' >> ${JSON_FILE}
printf '      \"distributionAddMomentum\": [0.0,0.0,0.0]\n' >> ${JSON_FILE}
printf '    },\n' >> ${JSON_FILE}
printf '    {\n' >> ${JSON_FILE}
printf '      \"enabled\": true,\n' >> ${JSON_FILE}
printf '      \"name\": \"ELEcont\",\n' >> ${JSON_FILE}
printf '      \"plasma\": \"cont\",\n' >> ${JSON_FILE}
printf '      \"ParticlesPerCell\": [%s,%s,%s],\n' "${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE}" "${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE}" "${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE}" >> ${JSON_FILE}
printf '      \"type\": \"ELECTRON\",\n' >> ${JSON_FILE}
printf '      \"isMarker\": 0,\n' >> ${JSON_FILE}
printf '      \"isTest\": false,\n' >> ${JSON_FILE}
printf '      \"distribution\": \"Maxwell\",\n' >> ${JSON_FILE}
printf '      \"distributionParams\": [6.0e-4],\n' >> ${JSON_FILE}
printf '      \"distributionAddMomentum\": [0.0,0.0,0.0]\n' >> ${JSON_FILE}
printf '    },\n' >> ${JSON_FILE}
printf '    {\n' >> ${JSON_FILE}
printf '      \"enabled\": true,\n' >> ${JSON_FILE}
printf '      \"name\": \"IONcloud\",\n' >> ${JSON_FILE}
printf '      \"plasma\": \"cloud\",\n' >> ${JSON_FILE}
printf '      \"ParticlesPerCell\": [%s,%s,%s],\n' "${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE}" "${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE}" "${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE}" >> ${JSON_FILE}
printf '      \"type\": \"ION\",\n' >> ${JSON_FILE}
printf '      \"Z\": %s,\n' "${NUMERO_IONIZZAZIONE_CLOUD}" >> ${JSON_FILE}
printf '      \"A\": %s,\n' "${NUMERO_ATOMICO_CLOUD}" >> ${JSON_FILE}
printf '      \"isMarker\": 0,\n' >> ${JSON_FILE}
printf '      \"isTest\": false\n' >> ${JSON_FILE}
printf '    },\n' >> ${JSON_FILE}
printf '    {\n' >> ${JSON_FILE}
printf '      \"enabled\": true,\n' >> ${JSON_FILE}
printf '      \"name\": \"IONbulk\",\n' >> ${JSON_FILE}
printf '      \"plasma\": \"bulk\",\n' >> ${JSON_FILE}
printf '      \"ParticlesPerCell\": [%s,%s,%s],\n' "${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE}" "${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE}" "${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE}" >> ${JSON_FILE}
printf '      \"type\": \"ION\",\n' >> ${JSON_FILE}
printf '      \"Z\": %s,\n' "${NUMERO_IONIZZAZIONE_BULK}" >> ${JSON_FILE}
printf '      \"A\": %s,\n' "${NUMERO_ATOMICO_BULK}" >> ${JSON_FILE}
printf '      \"isMarker\": 0,\n' >> ${JSON_FILE}
printf '      \"isTest\": false\n' >> ${JSON_FILE}
printf '    },\n' >> ${JSON_FILE}
printf '    {\n' >> ${JSON_FILE}
printf '      \"enabled\": true,\n' >> ${JSON_FILE}
printf '      \"name\": \"PROTcont\",\n' >> ${JSON_FILE}
printf '      \"plasma\": \"cont\",\n' >> ${JSON_FILE}
printf '      \"ParticlesPerCell\": [%s,%s,%s],\n' "${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE}" "${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE}" "${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE}" >> ${JSON_FILE}
printf '      \"type\": \"ION\",\n' >> ${JSON_FILE}
printf '      \"Z\": %s,\n' "${NUMERO_IONIZZAZIONE_CONT}" >> ${JSON_FILE}
printf '      \"A\": %s,\n' "${NUMERO_ATOMICO_CONT}" >> ${JSON_FILE}
printf '      \"isMarker\": 0,\n' >> ${JSON_FILE}
printf '      \"isTest\": false\n' >> ${JSON_FILE}
printf '    }\n' >> ${JSON_FILE}
printf '  ],\n' >> ${JSON_FILE}
#printf '  \"Domains\":[\n' >> ${JSON_FILE}
#printf '    {\n' >> ${JSON_FILE}
#printf '      \"name\": \"line\",\n' >> ${JSON_FILE}
#printf '      \"freeDim\": [1,0,0],\n' >> ${JSON_FILE}
#printf '      \"pointCoord\": [10.0,0.0,0.0]\n' >> ${JSON_FILE}
#printf '    }\n' >> ${JSON_FILE}
#printf '  ],\n' >> ${JSON_FILE}
printf '  \"Output\":[\n' >> ${JSON_FILE}
#printf '    {\"type\": \"EB\", \"every\": 10.0},\n' >> ${JSON_FILE}
#printf '    {\"type\": \"Density\", \"spec\": \"ELEcloud\", \"every\":10.0},\n' >> ${JSON_FILE}
#printf '    {\"type\": \"Density\", \"spec\": \"ELEbulk\", \"every\":10.0},\n' >> ${JSON_FILE}
#printf '    {\"type\": \"Density\", \"spec\": \"ELEcloud\", \"every\":1.0, \"in\": \"line\"},\n' >> ${JSON_FILE}
#printf '    {\"type\": \"Density\", \"spec\": \"ELEbulk\", \"every\":1.0, \"in\": \"line\"},\n' >> ${JSON_FILE}
printf '    {\"type\": \"Diag\", \"every\": 1.0}\n' >> ${JSON_FILE}
printf '  ]\n' >> ${JSON_FILE}
printf '\n' >> ${JSON_FILE}
printf '}\n' >> ${JSON_FILE}
printf '\n' >> ${JSON_FILE}


qsub ${JOB_FILE}

cd ..

done
done
done
done
done


