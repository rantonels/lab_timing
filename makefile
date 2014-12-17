bin/comptonfit:	src/comptonfit.cxx
	g++ src/comptonfit.cxx -o bin/comptonfit

clean:
	rm -f bin/*
