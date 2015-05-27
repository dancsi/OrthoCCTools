CC=g++
FLAGS=-O3 -std=c++11

solver.exe: src/solver.cpp
	$(CC) $(FLAGS) src/solver.cpp -o solver.exe

.PHONY: clean

clean:
	rm -f solver solver.exe