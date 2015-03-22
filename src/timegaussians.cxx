#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1.h"

#include <iostream>

#include <string>


const int numBins = 10;
const int binSize = 50;

int main(int argc, char** argv){
	
	int numBins = atoi(argv[2]);
	int binSize = atoi(argv[3]);


	TH1F *h1; 
	
	TFile * fp = new TFile(argv[1]);

	TTree *pjmca = (TTree*)fp->Get("pjmca");

	float ch1,ch2,ch3,ch0,temp,ene;
	pjmca->SetBranchAddress("ch0",&ch0);
	pjmca->SetBranchAddress("ch1",&ch1);
	pjmca->SetBranchAddress("ch2",&ch2);
	pjmca->SetBranchAddress("ch3",&ch3);

	for (int b = 0; b < numBins; b++)
	{

	h1 = new TH1F("h1","h1 ",2000,0,2000);

	int numEntries = pjmca->GetEntries();

	int counter = 0;

	for (int i=0; i<numEntries; i++)
	{
		pjmca->GetEntry(i);

		if (ch3>(b*binSize) && ch3<((b+1)*binSize))
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

	printf("sigma = %f\n",sigma);

	delete h1;

	}
}
