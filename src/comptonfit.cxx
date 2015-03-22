#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <climits>
#include <string>
#include "datafile.h"
#include <sstream>

using namespace std;

//Compton edge rispettivamente dei fotoni di 511 keV e 1275 keV
const int comp_edge1 = 340;
const int comp_edge2 = 1062;

float multiplier;

#define LORENTZIAN 1 //commentare per tornare alla gaussiana, decommentare per la lorentziana

double sqr(double r)
{
	return r*r;
}


//tutto cio' che riguarda le gaussiane

const int GAUSS_STEPS 	= 500;	//numero di entrate per stdev nella lookup table
const int GAUSS_STDEVS 	= 5;	//numero di deviazioni standard nella tabella
const int GAUSS_TSIZE 	= GAUSS_STEPS * GAUSS_STDEVS;

float gaussiana(float t, float sigma)
{
	return exp(-0.5 * sqr(t/sigma))/(sigma*sqrt(2*M_PI));
}


class Gaussianatore {
	private:
		double lutb[GAUSS_TSIZE];
	public:
		Gaussianatore();
		double gaussiamelo(float t, float stdev);

};

Gaussianatore::Gaussianatore(){
	cout << "Inizializzo il gaussianatore..." << endl;
	for (int i = 0; i<GAUSS_TSIZE; i++)
		{
			lutb[i] = gaussiana(float(i)/GAUSS_STEPS,1);
		};
	cout << "Inizializzazione lookup table completata." << endl;
}

double Gaussianatore::gaussiamelo(float t,float stdev)
{

//	cout << "richiesti servigi del gaussianatore: " << t << " " << stdev << endl << flush;

	float cstdev = max(stdev,1.F);

	double result =  lutb[
		min(GAUSS_TSIZE-1,
				(int)(abs(t)/cstdev*GAUSS_STEPS)
				)
	]/cstdev;

//	cout << "computato resultato del gauss...ntore" << endl << flush;
	return result;
};

Gaussianatore gaussianator;

/*Contiene i fattori del fit: energia del canale 0, rapporto eventi fotone 511 e 1275 keV, 
  fattore di scala x e y, sigma gaussiana che tiene conto della risoluzione energetica. 
  S è il chi2*/
class fattori {
	public:
		//double E;
		double e1;
		double e2;
		double k;

		double y;
		double sigma;
		double S;


		const fattori operator+(const fattori &other) const;
		const fattori operator-(const fattori &other) const;
		const fattori operator*(const double &other) const;

		double E();
		double x();
};

double fattori::E()
{
	return (comp_edge1 - e1/e2 * comp_edge2) / (1 - e1/e2);
}

double fattori::x()
{
	return ( comp_edge1 - E())/e1;
}

const fattori fattori::operator+(const fattori &other) const {
	fattori result = *this;   
	result.e1 += other.e1;
	result.k += other.k;
	result.e2 += other.e2;
	result.y += other.y;
	result.sigma += other.sigma;

	return result;  
}

const fattori fattori::operator-(const fattori &other) const {
	fattori result = *this;   
	result.e1 -= other.e1;
	result.k -= other.k;
	result.e2 -= other.e2;
	result.y -= other.y;
	result.sigma -= other.sigma;

	return result;  
}

const fattori fattori::operator*(const double &other) const {
	fattori result = *this;   
	result.e1 *= other;
	result.k *= other;
	result.e2 *= other;
	result.y *= other;
	result.sigma *= other;

	return result;  
}

string fattorirep(fattori f)
{
	ostringstream strs;
	
	strs << "e1=" << f.e1 << ",\te2=" << f.e2 << ",\tk=" << f.k << ",\ty=" << f.y << ",\tsigma/gamma=" << f.sigma;
	
	return strs.str();
	
}

fattori trim(fattori f)
{
	fattori o = {
				max(0.,f.e1), 
				max(0.,f.e2), 
				max(0.,f.k), 
				max(0.,f.y), 
				max(0.,f.sigma),
				f.S
			};
	return o;
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


////calcola la posizione dell'elemento corrispondente al compton edge del 511 keV 
//int edge1(fattori fit, int N) {
	//double a;
	//if (int((comp_edge1-fit.E)/fit.x)<N )      a= int((comp_edge1-fit.E)/fit.x);
	//else a= N;
	//return a;
//}

////analogo per 1275
//int edge2(fattori fit, int N) {
	//double a;
	//if (int((comp_edge2-fit.E)/fit.x)<N)       a= int((comp_edge2-fit.E)/fit.x);
	//else a= N;
	//return a;
//}


/*Riempie l'array A con il profilo compton generato dai due fotoni, 
  k definisce il rapporto tra le probabilità che un evento registrato 
  sia generato da un fotone a 511 keV invece che a 1275 keV */
double * profilo(int N, fattori fit)
{
	double * A = new double[N];
	
	double fitE = fit.E();
	double fitx = fit.x();
	
	for (int i=0; i<N; i++)
		A[i] = 0;
		
	for (int i=0; i<min((int)fit.e1,N); i++)
		{
			double E = fitE + fitx*i;
			double tmp1 = 2-2*E/(511.0-E)+E*E/((511.0-E)*(511.0-E))+E*E/(511*(511.0-E));
			A[i] += tmp1;
		}
		
		
	//double t = 1275/511.0;
	//for (i=0; i<min(edge2(fit,N),N); i++)
	//	A[i] = A[i]+fit.k*(2-2*E/(t*(1275.0-E))+E*E/(sqr(t)*(1275.0-E)*(1275.0-E))+sqr(E)/(1275*(1275.0-E)))/(t*1275.0);	

	for (int i=0; i<min((int)fit.e2,N); i++)
		{
			double E = fitE + fitx*i;
			double tmp2 = 2-2*E/(1275.0-E)+E*E/((1275.0-E)*(1275.0-E))+E*E/(1275.0*(1275.0-E));
			A[i] += fit.k*tmp2;
		}


	for (int i=0; i<N; i++)
		A[i] *= fit.y;
		

//	if ((A[i]<0)or std::isnan(std::abs(A[i])))
//		A[i] = 0;
		
	return A;
}



/*Ogni elemento A[i] genera una gaussiana di centro i, area A[i], e sigma data da parametro.
  B[i] è dato dalla somma di tutti i valori che assumono le gaussiane in i; per motivi pratici
  tronchiamo le gaussiane a 5*fit.sigma*/
//void convoluzione_alt(int N, double A[], double B[], fattori * fitp)
//{
	//fattori fit = (*fitp);
	//int i,j;
	//for (i=0; i<N; i++)
		//B[i] = 0;
	//for (i=1; i<edge2(fit,N); i++){
		//if ( (i>=5*fit.sigma) and ((i+5*fit.sigma)<N) ) 
			//for (j=(i-5*fit.sigma); j<(i+5*fit.sigma); j++)
				//B[j] = B[j] + A[i]/(sqrt(2*M_PI)*fit.sigma)*exp(-1.0*(j-i)*(j-i)/(fit.sigma*fit.sigma*2));
		//else{
			//if ( (i<5*fit.sigma) and ((i+5*fit.sigma)<N) ) 
				//for (j=0; j<(i+5*fit.sigma); j++)
					//B[j] = B[j] + A[i]/(sqrt(2*M_PI)*fit.sigma)*exp(-1.0*(j-i)*(j-i)/(fit.sigma*fit.sigma*2));
			//else{
				//for (j=(i-5*fit.sigma); j<N; j++)
					//B[j] = B[j] + A[i]/(sqrt(2*M_PI)*fit.sigma)*exp(-1.0*(j-i)*(j-i)/(fit.sigma*fit.sigma*2));}
		//}
	//}

//}	

double * convoluzione(int N, double * A, fattori * fitp)
{
	double * B = new double[N];
//	#ifdef LORENTZIAN
//		int maxin = (int)(10*(fitp->sigma))+1;
//	#else
//		int maxin = (int)(3*(fitp->sigma))+1;
//	#endif
//	
//	//	cout << "MAXIN = " << maxin << endl;
//
//	double * gauss = new double [maxin];
//	for(int i=0; i<maxin; i++)
//	{
//		#ifdef LORENTZIAN
//			gauss[i] = ((fitp->sigma)/M_PI )/( sqr(i) + sqr(fitp->sigma)  );
//		#else
//			gauss[i] =	0.8*		 1/(sqrt(2*M_PI)*(fitp->sigma)) * exp(-1.0*(i)*(i)/((fitp->sigma)*(fitp->sigma)*2));
//			gauss[i] += 0.2*		((fitp->sigma)/M_PI )/( sqr(i) + sqr(fitp->sigma)  );
//	//		cout << "gauss " << i << "\t" << gauss[i] << endl;
//		#endif
//	};
//


	for(int i = 0; i < N; i++)
	{
		B[i] = 0;
		float stdev = (fitp->sigma)*i;
		int maxin = (int)(GAUSS_STDEVS*stdev)+1;
		for(int j = - min(i,maxin-1); j < min(N-i,maxin); j++)
		{
			if((i+j<0) or (i+j >= N) or (abs(j) >= maxin))
			{	cerr << "comptonfit ERRORE: fuori range in convoluzione\n"; exit(1);}
			B[i] += A[i+j] * gaussianator.gaussiamelo(j,stdev);
		}
	};

	return B;

}


void minimize_chi2(int N, double Dati[], fattori fit1, fattori fit2, fattori & best, double passo)
{
	double * A;
	double * B;
	fattori bo, delta, limit;
	best.S = HUGE_VAL; //Nan??
	
	delta = (fit2-fit1)*(1/passo);

	limit = fit2 + delta*0.1; // questo permette di evitare confronti float pelo pelo
	
	
	int bcounter = 0;
	int bMAX = pow(passo+1,4) * (2*passo + 1);

	
	for (bo.e1=fit1.e1; bo.e1<=limit.e1; bo.e1= bo.e1 + delta.e1)
	for (bo.e2=fit1.e2; bo.e2<=limit.e2; bo.e2=bo.e2 +delta.e2)
	for (bo.k=fit1.k; bo.k<=limit.k; bo.k=bo.k + delta.k)
	for (bo.sigma=fit1.sigma; bo.sigma<=limit.sigma; bo.sigma=bo.sigma + delta.sigma)
	{

		bo.y = 1.;
		A = profilo(N,bo);
		B = convoluzione(N,A,&bo);
		delete [] A;

		for (bo.y=fit1.y; bo.y<=limit.y; bo.y=bo.y + 0.5*delta.y)
		{	
			//da qui in poi il profilo � bo.y*B;	
			
			
			bo.S = 0;
			for (int i=int(100*multiplier); i<min(N,800); i++)
				bo.S = bo.S + sqr(bo.y*B[i]-Dati[i]);///max( Dati[i] , 1.0); 
				//S e' il chi^2. Somma delle differenze al quadrato diviso la varianza al quadrato. La sigma e' sqrt(N) e N=Dati[i]
		

			if (bo.S<best.S){							
				
				best.e1 = bo.e1;
				best.e2 = bo.e2;
				best.y = bo.y;
				best.k = bo.k;
				best.sigma = bo.sigma;
				best.S = bo.S;
			}
			
			if (bcounter%50 == 0)
			{
				cout << "\rchi^2 = " << best.S << " " << (100*bcounter)/bMAX << "%       ";
				cout.flush();
			}
			
			
			bcounter ++;
		}


		delete [] B;
	
	}

	cout << endl;
	if (best.S == HUGE_VAL)
	{
		cerr << "comptonfit ERRORE: chi^2 del fit rimane infinito!\n";
		exit(1);
	}	
	
	cout << "parametri migliori parziali: " << fattorirep(best) << endl;
	cout << "chi^2 migliore parziale: " << best.S << endl;
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
		
		cout << endl << "* Iterazione di fit " << it+1 << "/" << iterazioni << endl;
		
		cout << "limiti per il fit:" << endl;
		cout << fattorirep(lbounds) << endl;
		cout << fattorirep(ubounds) << endl;
		
		minimize_chi2(N,Dati, lbounds, ubounds, tmpbest, npassi);

		semidelta = (ubounds-lbounds)* (0.5/npassi);

		lbounds = trim(tmpbest - semidelta);
		ubounds = tmpbest + semidelta;

	}
	
	(*best) = tmpbest;
}

void salva_array(int N, double arr[], string filename, string comment = "")
{
	ofstream file(filename.c_str());

	if (not file.is_open())
	{
		cerr << "comptonfit ERRORE: impossibile aprire " << filename << " per la scrittura.\n";
		exit(1);
	}

	file << "#" << comment << endl;
	for(int i=0; i<N; i++)
		file << i << "\t" << arr[i] << endl;

	file.close();

}

void salva_fattori(fattori tos, string filename, string comment)
{
		ofstream f(filename.c_str());
		
		if (not f.is_open())
		{
		cerr << "comptonfit ERRORE: impossibile aprire " << filename << " per la scrittura.\n";
			exit(1);
	}
		
		f << "#" << comment << endl;
		f << tos.e1 << endl;
		f << tos.e2 << endl;
		f << tos.k << endl;
		f << tos.y << endl;
		f << tos.sigma << endl;
		f << tos.S << endl;
		
		f.close();
}

void analisi(string fname, float gain = 1)
{
	cout << "Fit compton di " << fname << "..." << endl;
	cout << "Gain: " << gain << endl;
	cout << "* caricamento file...";
	cout.flush();
	double * Dati;
	FileArray farr = loadfile(fname);
	int N = farr.N;
	Dati = farr.array;
	cout << " " << N << " dati caricati." << endl;
	
	//for(int i=0;i<N;i++)
		//cout << i << "\t\t" << Dati[i] << endl;
	//exit(1);	
	
	//stima parametri
	
	//massimo
	double maxval = -HUGE_VAL;
	int maxpos = -1;
	
	for (int i=50; i<N; i++)
	{
		if(Dati[i] > maxval)
		{
			maxval = Dati[i];
			maxpos = i;
		}
	}
	
	if(maxval < 0)
	{
		cerr << "comptonfit ERRORE: il massimo dei dati e' negativo, o non e' stato trovato.\n";
		cerr << "massimo in posizione " << maxpos << " e valore " << maxval << "\n";
		exit(1);
	}
	
	cout << "Massimo dei dati: " << maxval << " al canale " << maxpos << endl;
	
	multiplier = maxpos/207.0;

	cout << "Moltiplicatore: " << multiplier << endl;

	fattori bestfat;
	
	fattori fit1,fit2;
	
	//INSERIRE RANGE MANUALI QUI

	//rapporto altezza picchi
	fit1.k = 0.003;
	fit2.k = 0.004;

	//coefficiente di sigma = coefficiente * x
	fit1.sigma = 0.12 * multiplier;
	fit2.sigma = .18 * multiplier;

	//scala y
	fit1.y = maxval* 0.5 * 0.8 * gain;
	fit2.y = maxval* 0.5 * 2.0 * gain;
	
	//posizione canale primo picco
	fit1.e1 = 225 * multiplier;
	fit2.e1 = 260 * multiplier;

	//posizione canale secondo picco
	fit1.e2 = 680 * multiplier;
	fit2.e2 = 730 * multiplier;
	
	#ifdef LORENTZIAN
	
	//fit2.E = 40;
	//fit1.k = 0.02;
	//fit1.sigma = 20;
	//fit1.x = 200.0/maxpos * 1.1;
	//fit2.x = 200.0/maxpos * 1.4;
	
	#endif
	
	
	fit( N, Dati, fit1, fit2, &bestfat, 
			5,	//STEP 
			3	//ITERAZIONI
			);
	
	string outfname =  fname + ".cfit";
	
	cout << "fit terminato. Salvo fattori ottimali in " << outfname << endl;
	
	salva_fattori(bestfat, outfname, "fattori risultanti dal fit di "+ fname);
	
	cout << "salvo curva di fit in " << fname + ".ccurve" << endl;
	
	double * outprof = profilo(N,bestfat);
	outprof = convoluzione(N,outprof,&bestfat);
	
	salva_array(N, outprof, fname + ".ccurve", "curva compton di fit per il file "+fname);
	
	delete [] outprof;

	cout << "salvo profilo non convoluto in " << fname + ".craw" << endl;
	
	double * outprofraw = profilo(N,bestfat);
	
	salva_array(N, outprof, fname + ".craw", "curva compton grezza (non convoluta) per il file "+fname);
	
	delete [] outprofraw;

}

int test_random()
{
	int i,N = 1000;
	double * A;
	double * B;
	fattori fit,best,fit1,fit2;
	//fit.E=-50;
	//fit.x = 1.68;
	fit.e1 = 225;
	fit.e2 = 675;
	fit.y = 19000.2;// 8700000;// 200000000;
	fit.k = 0.005;
	fit.sigma = 0.1;
	fit.S = 1000;
	A = profilo(1000,fit);
	
	salva_array(N,A, "tmp/randomcompton_raw", "curva Compton grezza simulata da bin/comptonfit. Generata con bin/comptonfit test");
	
	B = convoluzione(1000,A,&fit);

	salva_array(N,B, "tmp/randomcompton_con", "curva Compton convoluta simulata da bin/comptonfit. Generata con bin/comptonfit test");

	for(i=0;i<N;i++){
		B[i] = (rand() % 11 + 95)/100.0*(B[i]);}
	//cout << i << "     " << B[i] << endl;}

	salva_array(N,B, "tmp/randomcompton_ran", "curva Compton con rumore simulata da bin/comptonfit. Generata con bin/comptonfit test");



		//fit1.E=50;
		//fit1.x = 0.5;
		fit1.y = 100000000;
		fit1.k = 0.0005;
		fit1.sigma = 5;
		fit1.S = 1000;
		//fit2.E=300;
		//fit2.x = 1.5;
		fit2.y = 400000000;
		fit2.k = 0.01;
		fit2.sigma = 50;
		fit2.S = 1000;
		minimize_chi2(N,B,fit1,fit2,best,3.0);
		cout << fattorirep(best) << endl;
		A = profilo(1000,best);
		B = convoluzione(1000,A,&best);

		for(i=0;i<N;i++){
			cout << i << "     " << B[i] << endl;
		}


return 0;
}

int test_convoluzione()
{
	fattori fit;
	//fit.E=100;
	//fit.x = 1;
	fit.y = 200000000;
	fit.k = 0.001;
	fit.sigma = 25;
	fit.S = 1000;

	double cost [1000];
	for (int i = 0; i<1000; i++)
		cost[i] = 0;
	cost[500] = 1000;

	double * convcost;

	convcost = convoluzione(1000, cost, &fit);

	for (int j=0; j<1000; j++)
		cout << j << "\t" << convcost[j] << endl;

	double count = 0;
	for (int j=0; j<1000; j++)
		count += convcost[j];

	cout << "Somma: " << count << endl;

	return 0;
}


const string helpstring = "programma di fit del profilo compton\n\nUtilizzo:\n\ncomptonfit analisi XXX \t\t\tcarica i dati da XXX, esegui il fit e salva in tmp/XXX.cfit\n\ncomptonfit test\t\t\tesegui un test con dati generati (salvati in tmp/randomcompton) e stampa a terminale i risultati\n\ncomptonfit h\t\t\tstampa questo messaggio.";

int main(int argcount, char* argv[])
{
	if (argcount<=1)
	{
		cerr << "comptonfit ERRORE: sono un po' troppo pochi gli argomenti... prova bin/comptonfit analisi nomefiledati. E non demoralizzarti. Sei bravissimo." << endl;
		exit(0);
	}
	else
	{
		string command = string(argv[1]);		

		if(command == "analisi")
		{
			if (argcount <=2)
			{
				cerr << "comptonfit ERRORE: file dati non specificato." << endl;
				exit(1);
			}
			else
			{
				if (argcount == 3)
					analisi(argv[2]);
				else
					analisi(argv[2],atof(argv[3]));
				exit(0);
			}
		}

		if(argcount >= 3)
		{
			cerr << "comptonfit ERRORE: troppi argomenti in linea di comando." << endl;
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
		cerr << "comptonfit ERRORE: comando non riconosciuto: \"" << command << "\"" << endl;
		exit(1);

	}

}



