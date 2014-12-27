#include "datafile.h"
#include <fstream>
#include <iostream>
#include <cstdlib>

//nuova funzione per caricare i file in un array
int loadfile(double * a, string fname, int buffersize)
{
	double buffer[buffersize];
	ifstream f(fname.c_str());
	if(!f.is_open())
	{
		cout << "comptonfit ERRORE: errore nella lettura di " << fname << ": file non aperto" << endl;
		exit(1);
	}
	int n = 0;
	while(n < buffersize)
	{
		string line;
		f >> string;
		
		if (f == "")
			continue;
		
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
