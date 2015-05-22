#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1.h"

#include <iostream>

#include <string>

#include <fstream>


int main(int argc, char** argv){
	
	int numBins = atoi(argv[2]);
	int binSize = atoi(argv[3]);
	int chNum = atoi(argv[4]);

	std::ofstream outf;
	outf.open(argv[5],std::ofstream::out);

	TH1F *h1; 
	
	TFile * fp = new TFile(argv[1]);

	TTree *pjmca = (TTree*)fp->Get("pjmca");

	float ch1,ch2,ch3,ch0,temp,ene;
	pjmca->SetBranchAddress("ch0",&ch0);
	pjmca->SetBranchAddress("ch1",&ch1);
	//pjmca->SetBranchAddress("ch2",&ch2);
	//pjmca->SetBranchAddress("ch3",&ch3);

	float chan;
	if (chNum == 1)
		pjmca->SetBranchAddress("ch2",&chan);
	else if (chNum == 2)
		pjmca->SetBranchAddress("ch3",&chan);

	for (int b = 0; b < numBins; b++)
	{

	h1 = new TH1F("h1","h1 ",2000,0,2000);

	int numEntries = pjmca->GetEntries();

	int counter = 0;

	for (int i=0; i<numEntries; i++)
	{
		pjmca->GetEntry(i);

		if (chan>(b*binSize) && chan<((b+1)*binSize))
		{
			h1->Fill(ch1);
			counter++;
		}
	};


	printf("%d entries loaded.\t",counter);

	printf("fitting...\t");

	h1->Fit("gaus","Q");
	TF1 * gfun = h1->GetFunction("gaus");

	Double_t sigma = gfun->GetParameter(2);
	Double_t sigma_err = gfun->GetParError(2);

	printf("sigma = %f\n",sigma);
	printf("sigma_err = %f\n",sigma);

	outf << sigma << "," << sigma_err << std::endl;

	delete h1;

	};


	outf.close();
}
