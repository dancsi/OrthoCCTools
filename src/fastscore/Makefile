CC = g++
ALL_CFLAGS = -std=c++11 -O3 -fopenmp -g $(CFLAGS)
SOURCES=main.cpp Interaction.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=fastscore

ifeq ($(OS),Windows_NT)
	SOURCES += mman.c
endif

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(ALL_CFLAGS)  $(OBJECTS) -o $@

.cpp.o:
	$(CC) -c $(ALL_CFLAGS) $< -o $@