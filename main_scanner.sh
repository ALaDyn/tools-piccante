#! /bin/bash

################################################################################
# This file is part of piccante.
# 
# piccante is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# piccante is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with piccante.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

MAINFILE_NAME=main-1.cpp


rm -f ${MAINFILE_NAME}
touch ${MAINFILE_NAME}

ENABLE_OUTPUT=1
NDIM=2
NPROC_ALONG_Y=32
#NPROC_ALONG_Z=1
RESTART_FROM_DUMP=1
DO_RESTART="false"
DO_DUMP="false"
TIME_BTW_DUMP=10
DIRECTORY_OUTPUT="TEST"
DIRECTORY_DUMP="DUMP"
RANDOM_NUMBER_GENERATOR_SEED=5489
FREQUENCY_STDOUT_STATUS=5

XMIN="-50.0"
XMAX="0.0"
YMIN="-15.0"
YMAX="15.0"
#ZMIN="-1.0"
#ZMAX="1.0"
NCELLS_X=1536
NCELLS_Y=512
#NCELLS_Z=1
BND_X="xOpen"
BND_Y="yPBC"
BND_Z="zPBC"
CourantFactor="0.98"
TMAX="100.0"




### DO NOT TOUCH FROM HERE ###

printf "#define _USE_MATH_DEFINES\n" >> ${MAINFILE_NAME}
printf "#include <mpi.h>\n" >> ${MAINFILE_NAME}
printf "#include <cstdio>\n" >> ${MAINFILE_NAME}
printf "#include <iostream>\n" >> ${MAINFILE_NAME}
printf "#include <fstream>\n" >> ${MAINFILE_NAME}
printf "#include <sstream>\n" >> ${MAINFILE_NAME}
printf "#include <cmath>\n" >> ${MAINFILE_NAME}
printf "#include <iomanip>\n" >> ${MAINFILE_NAME}
printf "#include <cstring>\n" >> ${MAINFILE_NAME}
printf "#include <ctime>\n" >> ${MAINFILE_NAME}
printf "#if defined(_MSC_VER)\n" >> ${MAINFILE_NAME}
printf "#include \"gsl/gsl_rng.h\" \n" >> ${MAINFILE_NAME}
printf "#include \"gsl/gsl_randist.h\" \n" >> ${MAINFILE_NAME}
printf "#else\n" >> ${MAINFILE_NAME}
printf "#include <gsl/gsl_rng.h>\n" >> ${MAINFILE_NAME}
printf "#include <gsl/gsl_randist.h>\n" >> ${MAINFILE_NAME}
printf "#endif\n" >> ${MAINFILE_NAME}
printf "#include <cstdarg>\n" >> ${MAINFILE_NAME}
printf "#include <vector>\n" >> ${MAINFILE_NAME}
printf "#define DIMENSIONALITY $NDIM\n" >> ${MAINFILE_NAME}
printf "#include \"access.h\" \n" >> ${MAINFILE_NAME}
printf "#include \"commons.h\" \n" >> ${MAINFILE_NAME}
printf "#include \"grid.h\" \n" >> ${MAINFILE_NAME}
printf "#include \"structures.h\" \n" >> ${MAINFILE_NAME}
printf "#include \"current.h\" \n" >> ${MAINFILE_NAME}
printf "#include \"em_field.h\" \n" >> ${MAINFILE_NAME}
printf "#include \"particle_species.h\" \n" >> ${MAINFILE_NAME}
printf "#include \"output_manager.h\" \n" >> ${MAINFILE_NAME}
printf "#include \"utilities.h\" \n" >> ${MAINFILE_NAME}
printf "\n\n" >> ${MAINFILE_NAME}
printf "int main(int narg, char **args){\n" >> ${MAINFILE_NAME}
printf "  GRID grid;\n" >> ${MAINFILE_NAME}
printf "  EM_FIELD myfield;\n" >> ${MAINFILE_NAME}
printf "  CURRENT current;\n" >> ${MAINFILE_NAME}
printf "  std::vector<SPECIE*> species;\n" >> ${MAINFILE_NAME}
printf "  std::vector<SPECIE*>::const_iterator spec_iterator;\n" >> ${MAINFILE_NAME}
printf "  gsl_rng* rng = gsl_rng_alloc(gsl_rng_ranlxd1);\n" >> ${MAINFILE_NAME}
printf "  grid.setXrange(-50.0, 0.0);\n" >> ${MAINFILE_NAME}
printf "  grid.setYrange(-15.0, 15.0);\n" >> ${MAINFILE_NAME}
printf "  grid.setZrange(-1, +1);\n" >> ${MAINFILE_NAME}
printf "  grid.setNCells(1536, 512, 1);\n" >> ${MAINFILE_NAME}
printf "  grid.setNProcsAlongY(${NPROC_ALONG_Y});\n" >> ${MAINFILE_NAME}
#printf "  grid.setNProcsAlongZ(${NPROC_ALONG_Z});\n" >> ${MAINFILE_NAME}
printf "  grid.setBoundaries(xOpen | yPBC | zPBC);\n" >> ${MAINFILE_NAME}
printf "  grid.mpi_grid_initialize(&narg, args);\n" >> ${MAINFILE_NAME}
printf "  grid.setCourantFactor(0.98);\n" >> ${MAINFILE_NAME}
printf "  grid.setSimulationTime(100.0);\n" >> ${MAINFILE_NAME}
printf "  grid.with_particles = YES;//NO;\n" >> ${MAINFILE_NAME}
printf "  grid.with_current = YES;//YES;\n" >> ${MAINFILE_NAME}
printf "  grid.setStartMovingWindow(0);\n" >> ${MAINFILE_NAME}
printf "  grid.setMasterProc(0);\n" >> ${MAINFILE_NAME}
printf "  srand(time(NULL));\n" >> ${MAINFILE_NAME}
printf "  grid.initRNG(rng, ${RANDOM_NUMBER_GENERATOR_SEED});\n" >> ${MAINFILE_NAME}
printf "  grid.finalize();\n" >> ${MAINFILE_NAME}
printf "  grid.visualDiag();\n" >> ${MAINFILE_NAME}
printf "  myfield.allocate(&grid);\n" >> ${MAINFILE_NAME}
printf "  myfield.setAllValuesToZero();\n" >> ${MAINFILE_NAME}
printf "  laserPulse pulse1;\n" >> ${MAINFILE_NAME}
printf "  pulse1.setGaussianPulse();\n" >> ${MAINFILE_NAME}
printf "  pulse1.setWaist(4.0);\n" >> ${MAINFILE_NAME}
printf "  pulse1.setDurationFWHM(10.0);\n" >> ${MAINFILE_NAME}
printf "  pulse1.setNormalizedAmplitude(0.5);\n" >> ${MAINFILE_NAME}
printf "  pulse1.setCircularPolarization();\n" >> ${MAINFILE_NAME}
printf "  pulse1.setPulseInitialPosition(-10.1);\n" >> ${MAINFILE_NAME}
printf "  pulse1.setFocusPosition(0.0);\n" >> ${MAINFILE_NAME}
printf "  pulse1.setLambda(1.0);\n" >> ${MAINFILE_NAME}
printf "  pulse1.setFocusPosition(0.0);\n" >> ${MAINFILE_NAME}
printf "  myfield.addPulse(&pulse1);\n" >> ${MAINFILE_NAME}
printf "  myfield.boundary_conditions();\n" >> ${MAINFILE_NAME}
printf "  current.allocate(&grid);\n" >> ${MAINFILE_NAME}
printf "  current.setAllValuesToZero();\n" >> ${MAINFILE_NAME}
printf "  PLASMA plasma1;\n" >> ${MAINFILE_NAME}
printf "  plasma1.density_function = box;\n" >> ${MAINFILE_NAME}
printf "  plasma1.setXRangeBox(0.0, 100.0);\n" >> ${MAINFILE_NAME}
printf "  plasma1.setYRangeBox(grid.rmin[1], grid.rmax[1]);\n" >> ${MAINFILE_NAME}
printf "  plasma1.setZRangeBox(grid.rmin[2], grid.rmax[2]);\n" >> ${MAINFILE_NAME}
printf "  plasma1.setDensityCoefficient(0.0025);\n" >> ${MAINFILE_NAME}
printf "  SPECIE  electrons1(&grid);\n" >> ${MAINFILE_NAME}
printf "  electrons1.plasma = plasma1;\n" >> ${MAINFILE_NAME}
printf "  electrons1.setParticlesPerCellXYZ(2, 2, 2);\n" >> ${MAINFILE_NAME}
printf "  electrons1.setName(\"ELE1\");\n" >> ${MAINFILE_NAME}
printf "  electrons1.type = ELECTRON;\n" >> ${MAINFILE_NAME}
printf "  electrons1.creation();\n" >> ${MAINFILE_NAME}
printf "  species.push_back(&electrons1);\n" >> ${MAINFILE_NAME}
printf "  SPECIE ions1(&grid);\n" >> ${MAINFILE_NAME}
printf "  ions1.plasma = plasma1;\n" >> ${MAINFILE_NAME}
printf "  ions1.setParticlesPerCellXYZ(2, 2, 2);\n" >> ${MAINFILE_NAME}
printf "  ions1.setName(\"ION1\");\n" >> ${MAINFILE_NAME}
printf "  ions1.type = ION;\n" >> ${MAINFILE_NAME}
printf "  ions1.Z = 6.0;\n" >> ${MAINFILE_NAME}
printf "  ions1.A = 12.0;\n" >> ${MAINFILE_NAME}
printf "  for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) (*spec_iterator)->printParticleNumber();\n" >> ${MAINFILE_NAME}
printf "  OUTPUT_MANAGER manager(&grid, &myfield, &current, species);\n" >> ${MAINFILE_NAME}
printf "  manager.addEBFieldFrom(0.0, 10.0);\n" >> ${MAINFILE_NAME}
printf "  manager.addSpeciesDensityFrom(electrons1.name, 0.0, 10.0);\n" >> ${MAINFILE_NAME}
printf "  manager.addSpeciesPhaseSpaceFrom(electrons1.name, 0.0, 10.0);\n" >> ${MAINFILE_NAME}
printf "  manager.addDiagFrom(0.0, 1.0);\n" >> ${MAINFILE_NAME}
printf "  manager.initialize(${DIRECTORY_OUTPUT});\n" >> ${MAINFILE_NAME}
printf "  grid.setDumpPath(${DIRECTORY_DUMP});\n" >> ${MAINFILE_NAME}
printf "  if (grid.myid == grid.master_proc){\n" >> ${MAINFILE_NAME}
printf "    printf(\"----- START temporal cicle -----\\n\");\n" >> ${MAINFILE_NAME}
printf "    fflush(stdout);\n" >> ${MAINFILE_NAME}
printf "  }\n" >> ${MAINFILE_NAME}
printf "  int Nstep = grid.getTotalNumberOfTimesteps();\n" >> ${MAINFILE_NAME}
printf "  int dumpID = 1, dumpEvery;\n" >> ${MAINFILE_NAME}
printf "  if (${DO_DUMP}){\n" >> ${MAINFILE_NAME}
printf "    dumpEvery = (int)${TIME_BTW_DUMP} / grid.dt;\n" >> ${MAINFILE_NAME}
printf "  }\n" >> ${MAINFILE_NAME}
printf "  grid.istep = 0;\n" >> ${MAINFILE_NAME}
printf "  if (${DO_RESTART}){\n" >> ${MAINFILE_NAME}
printf "    dumpID = ${RESTART_FROM_DUMP};\n" >> ${MAINFILE_NAME}
printf "    restartFromDump(&dumpID, &grid, &myfield, species);\n" >> ${MAINFILE_NAME}
printf "  }\n" >> ${MAINFILE_NAME}
printf "  while (grid.istep <= Nstep) {\n" >> ${MAINFILE_NAME}
printf "    grid.printTStepEvery(${FREQUENCY_STDOUT_STATUS});\n" >> ${MAINFILE_NAME}
printf "    manager.callDiags(grid.istep);\n" >> ${MAINFILE_NAME}
printf "    myfield.openBoundariesE_1();\n" >> ${MAINFILE_NAME}
printf "    myfield.new_halfadvance_B();\n" >> ${MAINFILE_NAME}
printf "    myfield.boundary_conditions();\n" >> ${MAINFILE_NAME}
printf "    current.setAllValuesToZero();\n" >> ${MAINFILE_NAME}
printf "    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) (*spec_iterator)->current_deposition_standard(&current);\n" >> ${MAINFILE_NAME}
printf "    current.pbc();\n" >> ${MAINFILE_NAME}
printf "    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) (*spec_iterator)->position_parallel_pbc();\n" >> ${MAINFILE_NAME}
printf "    myfield.openBoundariesB();\n" >> ${MAINFILE_NAME}
printf "    myfield.new_advance_E(&current);\n" >> ${MAINFILE_NAME}
printf "    myfield.boundary_conditions();\n" >> ${MAINFILE_NAME}
printf "    myfield.openBoundariesE_2();\n" >> ${MAINFILE_NAME}
printf "    myfield.new_halfadvance_B();\n" >> ${MAINFILE_NAME}
printf "    myfield.boundary_conditions();\n" >> ${MAINFILE_NAME}
printf "    for (spec_iterator = species.begin(); spec_iterator != species.end(); spec_iterator++) (*spec_iterator)->momenta_advance(&myfield);\n" >> ${MAINFILE_NAME}
printf "    grid.time += grid.dt;\n" >> ${MAINFILE_NAME}
printf "    moveWindow(&grid, &myfield, species);\n" >> ${MAINFILE_NAME}
printf "    grid.istep++;\n" >> ${MAINFILE_NAME}
printf "    if (${DO_DUMP}){\n" >> ${MAINFILE_NAME}
printf "      if (grid.istep != 0 && !(grid.istep %% (dumpEvery))) {\n" >> ${MAINFILE_NAME}
printf "        dumpFilesForRestart(&dumpID, &grid, &myfield, species);\n" >> ${MAINFILE_NAME}
printf "      }\n" >> ${MAINFILE_NAME}
printf "    }\n" >> ${MAINFILE_NAME}
printf "  }\n" >> ${MAINFILE_NAME}
printf "  manager.close();\n" >> ${MAINFILE_NAME}
printf "  MPI_Finalize();\n" >> ${MAINFILE_NAME}
printf "  exit(1);\n" >> ${MAINFILE_NAME}
printf "}\n" >> ${MAINFILE_NAME}
printf "\n\n" >> ${MAINFILE_NAME}


