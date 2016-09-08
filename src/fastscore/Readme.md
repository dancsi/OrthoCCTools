# A C implementation of Amy Keating CC scoring function # 

## Build the test ##

g++ main.cpp Interaction.cpp -o test

g++ score-fasta.cpp Interaction.cpp -o score-fasta

## Swig wrapper ##

Just run ./make.sh