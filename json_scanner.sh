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


JSON_FILE=input.json


rm -f ${JSON_FILE}
touch ${JSON_FILE}


### DO NOT TOUCH FROM HERE ###

printf "{\n" >> ${JSON_FILE}
printf "  \"VERSION\": 2,\n" >> ${JSON_FILE}
printf "  \"dimensions\": 2,\n" >> ${JSON_FILE}
printf "  \"masterProc\": 0,\n" >> ${JSON_FILE}
printf "  \"OutputPath\": \"OUTPUT\",\n" >> ${JSON_FILE}
printf "  \"radiationFriction\": true,\n" >> ${JSON_FILE}
printf "  \"courantFactor\": 0.8,\n" >> ${JSON_FILE}
printf "  \"nProcY\": 64,\n" >> ${JSON_FILE}
printf "  \"xRange\": [0, 43.8],\n" >> ${JSON_FILE}
printf "  \"yRange\": [-22, 22],\n" >> ${JSON_FILE}
printf "  \"NCells\": [3072, 1536, 1],\n" >> ${JSON_FILE}
printf "  \"simulationTime\": 15,\n" >> ${JSON_FILE}
printf "  \"boundaries\": [\"open\", \"open\", \"periodic\"],\n" >> ${JSON_FILE}
printf "  \"StretchedGrid\":{\n" >> ${JSON_FILE}
printf "    \"enabled\": false,\n" >> ${JSON_FILE}
printf "    \"x\":{\"left\": { \"limit\": -10.0, \"NCells\": 1000 },\n" >> ${JSON_FILE}
printf "      \"right\": { \"limit\":  10.0, \"NCells\": 1000 } },\n" >> ${JSON_FILE}
printf "    \"y\":{\"left\": { \"limit\": -15.0, \"NCells\": 1000 },\n" >> ${JSON_FILE}
printf "      \"right\": { \"limit\":  15.0, \"NCells\": 1000 } },\n" >> ${JSON_FILE}
printf "    \"z\":{\"left\": { \"limit\": -15.0, \"NCells\": 1000 },\n" >> ${JSON_FILE}
printf "      \"right\": { \"limit\":  15.0, \"NCells\": 1000 } }\n" >> ${JSON_FILE}
printf "  },\n" >> ${JSON_FILE}
printf "  \"MovingWindow\":{\n" >> ${JSON_FILE}
printf "    \"enabled\": false,\n" >> ${JSON_FILE}
printf "    \"start\": 0,\n" >> ${JSON_FILE}
printf "    \"NO_beta\": 1,\n" >> ${JSON_FILE}
printf "    \"NO_frequency\":20\n" >> ${JSON_FILE}
printf "  },\n" >> ${JSON_FILE}
printf "  \"restart\":{\n" >> ${JSON_FILE}
printf "    \"DumpPath\": \"DUMP\",\n" >> ${JSON_FILE}
printf "    \"doRestart\": false,\n" >> ${JSON_FILE}
printf "    \"restartFromDump\": 1,\n" >> ${JSON_FILE}
printf "    \"doDump\": false,\n" >> ${JSON_FILE}
printf "    \"dumpEvery\": 2.0\n" >> ${JSON_FILE}
printf "  },\n" >> ${JSON_FILE}
printf "  \"special\":{\"variabile1\": 32,\"variabile2\": 32.5},\n" >> ${JSON_FILE}
printf "  \"Laser\":[\n" >> ${JSON_FILE}
printf "    {\n" >> ${JSON_FILE}
printf "      \"enabled\": true,\n" >> ${JSON_FILE}
printf "      \"type\": \"GAUSSIAN\",\n" >> ${JSON_FILE}
printf "      \"polarization\": \"P\",\n" >> ${JSON_FILE}
printf "      \"durationFWHM\": 16.5,\n" >> ${JSON_FILE}
printf "      \"initialPosition\": 16.51,\n" >> ${JSON_FILE}
printf "      \"waist\": 6.2,\n" >> ${JSON_FILE}
printf "      \"focusPosition\": 33.01,\n" >> ${JSON_FILE}
printf "      \"a\": 3.0,\n" >> ${JSON_FILE}
printf "      \"lambda\": 0.8,\n" >> ${JSON_FILE}
printf "      \"rotation\": false,\n" >> ${JSON_FILE}
printf "      \"angle\": 10.0,\n" >> ${JSON_FILE}
printf "      \"center\": 0,\n" >> ${JSON_FILE}
printf "      \"riseTime\": 2\n" >> ${JSON_FILE}
printf "    }\n" >> ${JSON_FILE}
printf "  ],\n" >> ${JSON_FILE}
printf "  \"Plasma\":[\n" >> ${JSON_FILE}
printf "    {\n" >> ${JSON_FILE}
printf "      \"name\": \"cloud\",\n" >> ${JSON_FILE}
printf "      \"densityFunction\": \"box\",\n" >> ${JSON_FILE}
printf "      \"XRangeBox\" : [33.02, 34.02],\n" >> ${JSON_FILE}
printf "      \"YRangeBox\" : [-40, 40],\n" >> ${JSON_FILE}
printf "      \"ZRangeBox\" : [-100, 100],\n" >> ${JSON_FILE}
printf "      \"DensityCoefficient\" : 2.0,\n" >> ${JSON_FILE}
printf "      \"DensityLambda\" : 0.8\n" >> ${JSON_FILE}
printf "    },{\n" >> ${JSON_FILE}
printf "      \"name\": \"bulk\",\n" >> ${JSON_FILE}
printf "      \"densityFunction\": \"left_fixed_exp_ramp\",\n" >> ${JSON_FILE}
printf "      \"XRangeBox\" : [34.02, 37.12 ],\n" >> ${JSON_FILE}
printf "      \"YRangeBox\" : [-40, 40],\n" >> ${JSON_FILE}
printf "      \"ZRangeBox\" : [-100, 100],\n" >> ${JSON_FILE}
printf "      \"DensityCoefficient\" : 100.0,\n" >> ${JSON_FILE}
printf "      \"DensityLambda\" : 0.8,\n" >> ${JSON_FILE}
printf "      \"RampMinDensity\" : 2.0,\n" >> ${JSON_FILE}
printf "      \"LeftRampLength\" : 0.7\n" >> ${JSON_FILE}
printf "    },{\n" >> ${JSON_FILE}
printf "      \"name\": \"cont\",\n" >> ${JSON_FILE}
printf "      \"densityFunction\": \"box\",\n" >> ${JSON_FILE}
printf "      \"XRangeBox\" : [37.12, 37.2],\n" >> ${JSON_FILE}
printf "      \"YRangeBox\" : [-40, 40],\n" >> ${JSON_FILE}
printf "      \"ZRangeBox\" : [-100, 100],\n" >> ${JSON_FILE}
printf "      \"DensityCoefficient\" : 10.0,\n" >> ${JSON_FILE}
printf "      \"DensityLambda\" : 0.8\n" >> ${JSON_FILE}
printf "    }\n" >> ${JSON_FILE}
printf "  ],\n" >> ${JSON_FILE}
printf "  \"Species\":[\n" >> ${JSON_FILE}
printf "    {\n" >> ${JSON_FILE}
printf "      \"enabled\": true,\n" >> ${JSON_FILE}
printf "      \"name\": \"ELEcloud\",\n" >> ${JSON_FILE}
printf "      \"plasma\": \"cloud\",\n" >> ${JSON_FILE}
printf "      \"ParticlesPerCell\": [2,2,1],\n" >> ${JSON_FILE}
printf "      \"type\": \"ELECTRON\",\n" >> ${JSON_FILE}
printf "      \"isMarker\": 0,\n" >> ${JSON_FILE}
printf "      \"isTest\": false,\n" >> ${JSON_FILE}
printf "      \"distribution\": \"Maxwell\",\n" >> ${JSON_FILE}
printf "      \"distributionParams\": [6.0e-4],\n" >> ${JSON_FILE}
printf "      \"distributionAddMomentum\": [0.0,0.0,0.0]\n" >> ${JSON_FILE}
printf "    },\n" >> ${JSON_FILE}
printf "    {\n" >> ${JSON_FILE}
printf "      \"enabled\": true,\n" >> ${JSON_FILE}
printf "      \"name\": \"ELEbulk\",\n" >> ${JSON_FILE}
printf "      \"plasma\": \"bulk\",\n" >> ${JSON_FILE}
printf "      \"ParticlesPerCell\": [4,3,1],\n" >> ${JSON_FILE}
printf "      \"type\": \"ELECTRON\",\n" >> ${JSON_FILE}
printf "      \"isMarker\": 0,\n" >> ${JSON_FILE}
printf "      \"isTest\": false,\n" >> ${JSON_FILE}
printf "      \"distribution\": \"Maxwell\",\n" >> ${JSON_FILE}
printf "      \"distributionParams\": [6.0e-4],\n" >> ${JSON_FILE}
printf "      \"distributionAddMomentum\": [0.0,0.0,0.0]\n" >> ${JSON_FILE}
printf "    },\n" >> ${JSON_FILE}
printf "    {\n" >> ${JSON_FILE}
printf "      \"enabled\": true,\n" >> ${JSON_FILE}
printf "      \"name\": \"ELEcont\",\n" >> ${JSON_FILE}
printf "      \"plasma\": \"cont\",\n" >> ${JSON_FILE}
printf "      \"ParticlesPerCell\": [2,2,1],\n" >> ${JSON_FILE}
printf "      \"type\": \"ELECTRON\",\n" >> ${JSON_FILE}
printf "      \"isMarker\": 0,\n" >> ${JSON_FILE}
printf "      \"isTest\": false,\n" >> ${JSON_FILE}
printf "      \"distribution\": \"Maxwell\",\n" >> ${JSON_FILE}
printf "      \"distributionParams\": [6.0e-4],\n" >> ${JSON_FILE}
printf "      \"distributionAddMomentum\": [0.0,0.0,0.0]\n" >> ${JSON_FILE}
printf "    },\n" >> ${JSON_FILE}
printf "    {\n" >> ${JSON_FILE}
printf "      \"enabled\": false,\n" >> ${JSON_FILE}
printf "      \"name\": \"IONcloud\",\n" >> ${JSON_FILE}
printf "      \"plasma\": \"cloud\",\n" >> ${JSON_FILE}
printf "      \"ParticlesPerCell\": [2,2,1],\n" >> ${JSON_FILE}
printf "      \"type\": \"ION\",\n" >> ${JSON_FILE}
printf "      \"Z\": 10.0,\n" >> ${JSON_FILE}
printf "      \"A\": 27.0,\n" >> ${JSON_FILE}
printf "      \"isMarker\": 0,\n" >> ${JSON_FILE}
printf "      \"isTest\": false\n" >> ${JSON_FILE}
printf "    },\n" >> ${JSON_FILE}
printf "    {\n" >> ${JSON_FILE}
printf "      \"enabled\": false,\n" >> ${JSON_FILE}
printf "      \"name\": \"IONbulk\",\n" >> ${JSON_FILE}
printf "      \"plasma\": \"bulk\",\n" >> ${JSON_FILE}
printf "      \"ParticlesPerCell\": [1,2,1],\n" >> ${JSON_FILE}
printf "      \"type\": \"ION\",\n" >> ${JSON_FILE}
printf "      \"Z\": 10.0,\n" >> ${JSON_FILE}
printf "      \"A\": 27.0,\n" >> ${JSON_FILE}
printf "      \"isMarker\": 0,\n" >> ${JSON_FILE}
printf "      \"isTest\": false\n" >> ${JSON_FILE}
printf "    },\n" >> ${JSON_FILE}
printf "    {\n" >> ${JSON_FILE}
printf "      \"enabled\": false,\n" >> ${JSON_FILE}
printf "      \"name\": \"PROTcont\",\n" >> ${JSON_FILE}
printf "      \"plasma\": \"cont\",\n" >> ${JSON_FILE}
printf "      \"ParticlesPerCell\": [2,2,1],\n" >> ${JSON_FILE}
printf "      \"type\": \"ION\",\n" >> ${JSON_FILE}
printf "      \"Z\": 1.0,\n" >> ${JSON_FILE}
printf "      \"A\": 1.0,\n" >> ${JSON_FILE}
printf "      \"isMarker\": 0,\n" >> ${JSON_FILE}
printf "      \"isTest\": false\n" >> ${JSON_FILE}
printf "    }\n" >> ${JSON_FILE}
printf "  ],\n" >> ${JSON_FILE}
printf "  \"Domains\":[\n" >> ${JSON_FILE}
printf "    {\n" >> ${JSON_FILE}
printf "      \"name\": \"line\",\n" >> ${JSON_FILE}
printf "      \"freeDim\": [1,0,0],\n" >> ${JSON_FILE}
printf "      \"pointCoord\": [10.0,0.0,0.0]\n" >> ${JSON_FILE}
printf "    }\n" >> ${JSON_FILE}
printf "  ],\n" >> ${JSON_FILE}
printf "  \"Output\":[\n" >> ${JSON_FILE}
printf "    {\"type\": \"EB\", \"to\": 10.0, \"every\": 10.0},\n" >> ${JSON_FILE}
printf "    {\"type\": \"Density\", \"spec\": \"ELEcloud\", \"every\":10.0},\n" >> ${JSON_FILE}
printf "    {\"type\": \"Density\", \"spec\": \"ELEbulk\", \"every\":10.0},\n" >> ${JSON_FILE}
printf "    {\"type\": \"Density\", \"spec\": \"ELEcloud\", \"every\":1.0, \"in\": \"line\"},\n" >> ${JSON_FILE}
printf "    {\"type\": \"Density\", \"spec\": \"ELEbulk\", \"every\":1.0, \"in\": \"line\"},\n" >> ${JSON_FILE}
printf "    {\"type\": \"Diag\", \"every\": 1.0}\n" >> ${JSON_FILE}
printf "  ]\n" >> ${JSON_FILE}
printf "\n" >> ${JSON_FILE}
printf "}\n" >> ${JSON_FILE}
printf "\n" >> ${JSON_FILE}


