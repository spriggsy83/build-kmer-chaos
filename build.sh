#Requires Boost C++ Libraries and OpenMPI
#module load boost
#module load openmpi

g++ -std=c++0x -o buildKmerChaos buildKmerChaos.cpp KmerChaosGen.cpp SeqReader.cpp -fopenmp -lboost_iostreams -lz
