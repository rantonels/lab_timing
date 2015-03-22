MKDIR_P = mkdir -p
DIRS = bin build out tmp
FLAGS = -std=c++0x -O3
GCC = g++ $(FLAGS)

.PHONY: directories

all: bin/comptonfit bin/root2csv bin/timegaussians

#datafile

build/datafile.o: src/datafile.cxx
	$(GCC) -c src/datafile.cxx -o build/datafile.o

#comptonfit

build/comptonfit.o:	src/comptonfit.cxx src/datafile.h
	$(GCC) -c src/comptonfit.cxx -o build/comptonfit.o

bin/comptonfit: build/comptonfit.o build/datafile.o
	$(GCC) build/datafile.o build/comptonfit.o -o bin/comptonfit

#root2csv

R_CFLAGS = $(shell root-config --cflags)
R_LIBS   = $(shell root-config --libs)
R_GLIBS  = $(shell root-config --glibs)

bin/root2csv: src/root2csv.cxx
	$(GCC) $(R_CFLAGS) src/root2csv.cxx -o bin/root2csv $(R_LIBS)

#timegaussians

bin/timegaussians: src/timegaussians.cxx
	$(GCC) $(R_CFLAGS) src/timegaussians.cxx -o bin/timegaussians $(R_LIBS)

#cleanup

clean:
	rm -f bin/*
	rm -f build/*

#creazione cartelle

$(shell mkdir -p $(DIRS))
