build/datafile.o: src/datafile.cxx
	g++ -c src/datafile.cxx -o build/datafile.o

bin/comptonfit.o:	src/comptonfit.cxx src/datafile.h
	g++ -c src/comptonfit.cxx -o bin/comptonfit

bin/comptonfit: build/comptonfit.o build/datafile.o
	g++ build/datafile.o build/comptonfit.o -o bin/comptonfit

clean:
	rm -f bin/*
	rm -f build/*
