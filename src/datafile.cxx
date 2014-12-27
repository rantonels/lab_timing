#include "datafile.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <stdexcept>

//nuova funzione per caricare i file in un array
int loadfile(double * a, string fname, int buffersize)
{
	double buffer[buffersize];
	ifstream f(fname.c_str());
	if(!f.is_open())
	{
		cerr << "datafile ERRORE: errore nella lettura di " << fname << ": file non aperto" << endl;
		exit(1);
	}
	int n = 0;
	while(n < buffersize)
	{
		string line;
		f >> line;
		
		if (line == "")
			continue;
	
		if (line[0] == '#')
			continue;
	
		try {
			buffer[n] = stod(line);
			n++;
		}
		catch (const invalid_argument & ia) {
			cerr << "datafile ERRORE: (linea " << n << " del file " << fname << ") conversione stringa > double fallita" << endl;
			cerr << "Stringa di errore: " << ia.what() << endl;
			exit(1);
		}
	
		
		if(f.eof())
			break;
		n++;
	}
	f.close();

	a = new double[n];
	for(int i=0; i<n; i++)
		a[i] = buffer[i];

	return n;
}
