MKDIR_P = mkdir -p
DIRS = bin build out tmp
.PHONY: directories

all: bin/comptonfit

#datafile

build/datafile.o: src/datafile.cxx
	g++ -c src/datafile.cxx -o build/datafile.o

#comptonfit

build/comptonfit.o:	src/comptonfit.cxx src/datafile.h
	g++ -c src/comptonfit.cxx -o build/comptonfit.o

bin/comptonfit: build/comptonfit.o build/datafile.o
	g++ build/datafile.o build/comptonfit.o -o bin/comptonfit

#cleanup

clean:
	rm -f bin/*
	rm -f build/*

#creazione cartelle

$(shell mkdir -p $(DIRS))
