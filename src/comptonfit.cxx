#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <climits>
#include <string>
#include "datafile.h"

using namespace std;

//Compton edge rispettivamente dei fotoni di 511 keV e 1275 keV
const int comp_edge1 = 340;
const int comp_edge2 = 1062;

/*Contiene i fattori del fit: energia del canale 0, rapporto eventi fotone 511 e 1275 keV, 
  fattore di scala x e y, sigma gaussiana che tiene conto della risoluzione energetica. 
  S è il chi2*/
class fattori {
	public:
		double E;
		double k;
		double x;
		double y;
		double sigma;
		double S;

		const fattori operator+(const fattori &other) const;
		const fattori operator-(const fattori &other) const;
		const fattori operator*(const double &other) const;
};

const fattori fattori::operator+(const fattori &other) const {
	fattori result = *this;   
	result.E += other.E;
	result.k += other.k;
	result.x += other.x;
	result.y += other.y;
	result.sigma += other.sigma;

	return result;  
}

const fattori fattori::operator-(const fattori &other) const {
	fattori result = *this;   
	result.E -= other.E;
	result.k -= other.k;
	result.x -= other.x;
	result.y -= other.y;
	result.sigma -= other.sigma;

	return result;  
}

const fattori fattori::operator*(const double &other) const {
	fattori result = *this;   
	result.E *= other;
	result.k *= other;
	result.x *= other;
	result.y *= other;
	result.sigma *= other;

	return result;  
}



//!!! Questa funzione e' stata sostituita con load_file da datafile.h, con allocazione dinamica
//ficca in Dati[] i dati registrati
int leggi(double Dati[], int max, string fname=datafname)
{
	long int N = 0;
	ifstream f(datafname.c_str());
	while(N < max)
	{
		f >> Dati[N];
		if(f.eof()) 
			break;
		N++;
	}
	f.close();
	return N;
}


//calcola la posizione dell'elemento corrispondente al compton edge del 511 keV 
int edge1(fattori fit, int N) {
	double a;
	if (int((comp_edge1-fit.E)/fit.x)<N )      a= int((comp_edge1-fit.E)/fit.x);
	else a= N;
	return a;
}

//analogo per 1275
int edge2(fattori fit, int N) {
	double a;
	if (int((comp_edge2-fit.E)/fit.x)<N)       a= int((comp_edge2-fit.E)/fit.x);
	else a= N;
	return a;
}

/*Riempie l'array A con il profilo compton generato dai due fotoni, 
  k definisce il rapporto tra le probabilità che un evento registrato 
  sia generato da un fotone a 511 keV invece che a 1275 keV */
void profilo(int N, double A[], fattori fit)
{
	int i;
	double E;
	for (i=0; i<N; i++)
		A[i] = 0;
	for (i=0; i<edge1(fit,N); i++){
		E = fit.E + fit.x*i;
		A[i] = 2-2*E/(511.0-E)+E*E/((511.0-E)*(511.0-E))+E*E/(511*(511.0-E));
		A[i]= A[i]/511.0;}
	double t = 1275/511.0;
	for (i=0; i<edge2(fit,N); i++)
		A[i] = A[i]+fit.k*(2-2*E/(t*(1275.0-E))+E*E/(t*t*(1275.0-E)*(1275.0-E))+E*E/(1275*(1275.0-E)))/(t*1275.0);	
	for (i=0; i<N; i++)
		A[i] = A[i]*fit.y;
}

/*Ogni elemento A[i] genera una gaussiana di centro i, area A[i], e sigma data da parametro.
  B[i] è dato dalla somma di tutti i valori che assumono le gaussiane in i; per motivi pratici
  tronchiamo le gaussiane a 5*fit.sigma*/
void convoluzione_alt(int N, double A[], double B[], fattori * fitp)
{
	fattori fit = (*fitp);
	int i,j;
	for (i=0; i<N; i++)
		B[i] = 0;
	for (i=1; i<edge2(fit,N); i++){
		if ( (i>=5*fit.sigma) and ((i+5*fit.sigma)<N) ) 
			for (j=(i-5*fit.sigma); j<(i+5*fit.sigma); j++)
				B[j] = B[j] + A[i]/(sqrt(2*M_PI)*fit.sigma)*exp(-1.0*(j-i)*(j-i)/(fit.sigma*fit.sigma*2));
		else{
			if ( (i<5*fit.sigma) and ((i+5*fit.sigma)<N) ) 
				for (j=0; j<(i+5*fit.sigma); j++)
					B[j] = B[j] + A[i]/(sqrt(2*M_PI)*fit.sigma)*exp(-1.0*(j-i)*(j-i)/(fit.sigma*fit.sigma*2));
			else{
				for (j=(i-5*fit.sigma); j<N; j++)
					B[j] = B[j] + A[i]/(sqrt(2*M_PI)*fit.sigma)*exp(-1.0*(j-i)*(j-i)/(fit.sigma*fit.sigma*2));}
		}
	}

}	

void convoluzione(int N, double A[], double B[], fattori * fitp)
{
	int maxin = (int)(5*(fitp->sigma))+1;
	//	cout << "MAXIN = " << maxin << endl;

	double * gauss = new double [maxin];
	for(int i=0; i<maxin; i++)
	{
		gauss[i] = 1/(sqrt(2*M_PI)*(fitp->sigma)) * exp(-1.0*(i)*(i)/((fitp->sigma)*(fitp->sigma)*2));
	};



	for(int i = 0; i < N; i++)
	{
		B[i] = 0;
		for(int j = - min(i,maxin); j < min(N-i,maxin); j++)
		{
			B[i] += A[i+j] * gauss[abs(j)];
		}
	};

}


void minimize_chi2(int N, double Dati[], fattori fit1, fattori fit2, fattori & best, double passo)
{
	double A[1000];
	double B[1000];
	fattori bo, delta;
	best.S = HUGE_VAL; //Nan??
	//delta.E = (fit2.E-fit1.E)/passo;
	//delta.x = (fit2.x-fit1.x)/passo;
	//delta.y = (fit2.y-fit1.y)/passo;
	//delta.k = (fit2.k-fit1.k)/passo;

	delta = (fit2-fit1)*(1/passo);

	delta.sigma = (fit2.sigma-fit1.sigma)/passo;
	for (bo.E=fit1.E; bo.E<=fit2.E; bo.E= bo.E + delta.E)
		for (bo.x=fit1.x; bo.x<=fit2.x; bo.x=bo.x +delta.x)
			for (bo.y=fit1.y; bo.y<=fit2.y; bo.y=bo.y + delta.y)
				for (bo.k=fit1.k; bo.k<=fit2.k; bo.k=bo.k + delta.k)
					for (bo.sigma=fit1.sigma; bo.sigma<=fit2.sigma; bo.sigma=bo.sigma + delta.sigma){
						//cout << bo.sigma << endl;
						profilo(N,A,bo);
						convoluzione(N,A,B,&bo);
						bo.S = 0;
						for (int i=0; i<1000; i++)
							bo.S = bo.S + (B[i]-Dati[i])*(B[i]-Dati[i])/(Dati[i]+1); //S e' il chi^2. Somma delle differenze al quadrato diviso la varianza al quadrato. La sigma e' sqrt(N) e N=Dati[i]
						//cout << bo.S << endl;
						if (bo.S<best.S){
							best.E = bo.E;
							best.x = bo.x;
							best.y = bo.y;
							best.k = bo.k;
							best.sigma = bo.sigma;
							best.S = bo.S;
						}
					}
	//cout << best.S << endl;
}	


void fit(int N, double Dati[], fattori fitbeg, fattori fitend, fattori * best, double npassi, int iterazioni)
{
	fattori lbounds, ubounds, tmpbest, semidelta;

	cout << "* fit iterativo..." << endl;
	cout << "num passi: " << npassi << "; iterazioni: "<< iterazioni << endl;


	lbounds = fitbeg;
	ubounds = fitend;

	for(int it=0; it<iterazioni; it++)
	{
		cout << "* Iterazione di fit " << it << "/" << iterazioni << endl;
		minimize_chi2(N,Dati, lbounds, ubounds, tmpbest, npassi);

		semidelta = (ubounds-lbounds)* (0.5/npassi);

		lbounds = tmpbest - semidelta;
		ubounds = tmpbest + semidelta;

	}
	
	(*best) = tmpbest;
}

void salva_array(int N, double arr[], string filename, string comment = "")
{
	ofstream file(filename.c_str());

	file << "#" << comment << endl;
	for(int i=0; i<N; i++)
		file << i << "\t" << arr[i] << endl;

	file.close();

}

void salva_fattori(fattori tos, string filename, string comment)
{
		ofstream f(filename.c_str());
		f << "#" << comment;
		f << tos.E << endl;
		f << tos.x << endl;
		f << tos.y << endl;
		f << tos.k << endl;
		f << tos.sigma << endl;
		f << tos.S << endl;
		
		f.close();
}

void analisi(string fname)
{
	cout << "Fit compton di " << fname << "..." << endl;
	cout << "* caricamento file...";
	double * Dati;
	int N = loadfile(Dati,fname);
	cout << " " << N << " dati caricati." << endl;
	
	
	fattori bestfat;
	
	//fattori di bound a caso, CAMBIARE!!!
	fattori fit1,fit2;
		fit1.E=50;
		fit1.x = 0.5;
		fit1.y = 100000000;
		fit1.k = 0.0005;
		fit1.sigma = 5;
		fit1.S = 1000;
		fit2.E=300;
		fit2.x = 1.5;
		fit2.y = 400000000;
		fit2.k = 0.01;
		fit2.sigma = 50;
		fit2.S = 1000;
	
	
	fit( N, Dati, fit1, fit2, &bestfat, 3, 5);
	
	string outfname = "tmp/" + fname + ".cfit";
	
	cout << "fit terminato. Salvo in " << outfname << endl;
	
	salva_fattori(bestfat, outfname, "fattori risultanti dal fit di "+ fname);
	
}

int test_random()
{
	int i,N = 1000;
	double A[1000];
	double B[1000];
	fattori fit,best,fit1,fit2;
	fit.E=100;
	fit.x = 1;
	fit.y = 200000000;
	fit.k = 0.001;
	fit.sigma = 25;
	fit.S = 1000;
	profilo(1000,A,fit);
	convoluzione(1000,A,B,&fit);

	for(i=0;i<N;i++){
		B[i] = (rand() % 11 + 95)/100.0*(B[i]);}
	//cout << i << "     " << B[i] << endl;}

	salva_array(N,B, "tmp/randomcompton", "curva Compton simulata da bin/comptonfit. Generata con bin/comptonfit test");



		fit1.E=50;
		fit1.x = 0.5;
		fit1.y = 100000000;
		fit1.k = 0.0005;
		fit1.sigma = 5;
		fit1.S = 1000;
		fit2.E=300;
		fit2.x = 1.5;
		fit2.y = 400000000;
		fit2.k = 0.01;
		fit2.sigma = 50;
		fit2.S = 1000;
		minimize_chi2(N,B,fit1,fit2,best,3.0);
		cout << best.S << "     " << best.E << "    " << best.x << "     " << best.y << "    " << best.k <<"      " << best.sigma << endl;
		profilo(1000,A,best);
		convoluzione(1000,A,B,&best);

		for(i=0;i<N;i++){
			cout << i << "     " << B[i] << endl;
		}


return 0;
}

int test_convoluzione()
{
	fattori fit;
	fit.E=100;
	fit.x = 1;
	fit.y = 200000000;
	fit.k = 0.001;
	fit.sigma = 25;
	fit.S = 1000;

	double cost [1000];
	for (int i = 0; i<1000; i++)
		cost[i] = 0;
	cost[500] = 1000;

	double convcost [1000];

	convoluzione(1000, cost, convcost, &fit);

	for (int j=0; j<1000; j++)
		cout << j << "\t" << convcost[j] << endl;

	double count = 0;
	for (int j=0; j<1000; j++)
		count += convcost[j];

	cout << "Somma: " << count << endl;

	return 0;
}


const string helpstring = "programma di fit del profilo compton\n\nUtilizzo:\n\ncomptonfit XXX \t\t\tcarica i dati da XXX, esegui il fit e salva in tmp/XXX.cfit\n\ncomptonfit test\t\t\tesegui un test con dati generati (salvati in tmp/randomcompton) e stampa a terminale i risultati\n\ncomptonfit h\t\t\tstampa questo messaggio.";

int main(int argcount, char* argv[])
{
	if (argcount<=1)
	{
		cout << "comptonfit ERRORE: sono un po' troppo pochi gli argomenti... prova bin/comptonfit analisi nomefiledati. E non demoralizzarti. Sei bravissimo." << endl;
		exit(0);
	}
	else
	{
		string command = string(argv[1]);		

		if(command == "analisi")
		{
			if (argcount <=2)
			{
				cout << "comptonfit ERRORE: file dati non specificato." << endl;
				exit(1);
			}
			else
			{
				analisi(argv[2]);
				exit(0);
			}
		}

		if(argcount >= 3)
		{
			cout << "comptonfit ERRORE: troppi argomenti in linea di comando." << endl;
			exit(1);
		};


		if(command == "test")
		{
			test_random();
			exit(0);
		}
		if((command == "h") or (command == "-h") or (command == "--help") or (command == "?") or (command == "-?"))
		{
			cout << helpstring << endl;
			exit(0);
		}
		//default
		cout << "comptonfit ERRORE: comando non riconosciuto: \"" << command << "\"" << endl;
		exit(1);

	}

}



