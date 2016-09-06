# 2016 / 09 / 06
# Javier Hernando-Ayuso
# see https://gcc.gnu.org/wiki/GFortranGettingStarted

clear
clear


# the -c option is used to compile code split in different files
# the -ffree-line-length-none is necessary as long lines are contained in the code
gfortran -c  -ffree-line-length-none factorial.f90  CollisionProbability_GPHA.f90 CollisionProbability_Chan.f90 CollisionProbability_Alfano.f90 CollisionProbability_Serra.f90 CollisionProbability_Foster.f90 test.f90
gfortran -o  test factorial.o CollisionProbability_GPHA.o CollisionProbability_Chan.o CollisionProbability_Alfano.o CollisionProbability_Serra.o CollisionProbability_Foster.o test.o

# delete the object files
rm factorial.o  CollisionProbability_GPHA.o CollisionProbability_Chan.o CollisionProbability_Alfano.o CollisionProbability_Serra.o CollisionProbability_Foster.o test.o
