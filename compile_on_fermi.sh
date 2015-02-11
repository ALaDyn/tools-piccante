module purge
module load bgq-xl

make -f makefile.fermi.xl all
cp newReader lightReader titan frogReader MPItitan ~/bin
