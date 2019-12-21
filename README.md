## README
  
Build chaos game representation of kmers of sequence, as per:
https://towardsdatascience.com/chaos-game-representation-of-a-genetic-sequence-4681f1a67e14
  
C++ program prerequisites: Boost C++ Libraries (www.boost.org) and OpenMPI (www.open-mpi.org/)  
Ubuntu preparation:
```
sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev
sudo apt-get install libboost-all-dev libboost-dev libboost-doc libboost-container-dev
sudo apt-get install zlib1g-dev zlib1g
```
  
Compile on Unix systems with:  
`bash build.sh`
  
Run as:
`buildKmerChaos <infile> <outfile> <kmer size>`
Input file may be fasta or fastq and may be .gz compressed.
  

Code by Andrew Spriggs, CSIRO Ag&Food (www.csiro.au)  
andrew.spriggs@csiro.au  
github.com/spriggsy83
