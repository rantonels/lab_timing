#ifndef __DATAFILE_H__
#define __DATAFILE_H__

#include <string>

using namespace std;

const string datafname = "data/?.dat";

struct FileArray{
	public:
		int N;
		double * array;
};

FileArray loadfile(string fname=datafname, int buffersize=100000);






#endif
