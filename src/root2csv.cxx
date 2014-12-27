#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"

#include <cstdlib>
#include <iostream>

#include <fstream>
#include <string>


const char * HNAMES [] = { "h0", "h1", "h2", "h3"};
const char * CNAMES [] = { "ch0", "ch1", "ch2", "ch3"};

int main(int argc, char** argv) {

	if(argc<2)
	{
		std::cerr << "root2csv ERRORE: argomenti insufficienti.\n";
		exit(1);
	}
	
	if(argc>4)
	{
		std::cerr << "root2csv ERRORE: troppi argomenti in linea di comando.\n";
		exit(1);
	}
	
	bool CLIMODE;
	
	CLIMODE = argc < 3;
	

    TFile file(argv[1]);
    TNtuple* tuple = reinterpret_cast<TNtuple*>(file.Get("pjmca"));
    Float_t ch0;
    Float_t ch1;
    Float_t ch2;
    Float_t ch3;
    tuple->SetBranchAddress("ch0", &ch0);
    tuple->SetBranchAddress("ch1", &ch1);
    tuple->SetBranchAddress("ch2", &ch2);
    tuple->SetBranchAddress("ch3", &ch3);
    int n = tuple->GetEntries();

	// preparazione file di output
	std::ofstream outf;
	std::string outfnames [4];
	if (not CLIMODE)
	{
		for (int i =0; i<4; i++)
			outfnames[i] = std::string(argv[2])+"_"+std::string(HNAMES[i]);
	}
	

    for (int faggot = 0; faggot < 4; ++faggot)
    {
		
		if (CLIMODE)
			std::cout << "###" << CNAMES[faggot] << "###" << std::endl;
        
        TH1D histogram(HNAMES[faggot], HNAMES[faggot], 2000, 0, 2000);
        tuple->Project(HNAMES[faggot], CNAMES[faggot]);
        
        
        if (CLIMODE)
        {
			for (int i = 0; i < histogram.GetNbinsX(); ++i) {
				std::cout << i << " " << histogram.GetBinContent(i) << std::endl;
			}
        }    
        else
        {
			outf.open(outfnames[faggot].c_str());
		
			for(int i= 0; i<histogram.GetNbinsX(); ++i)
				outf << histogram.GetBinContent(i) << std::endl;
				
			outf.close();
		}
		
	
    }
    //for (int i = 0; i < n; ++i) {
    //int read = tuple->GetEntry(i);
    //std::cout << i << " " << ch0 << " " << ch1 << " " << ch2 << " " << ch3 << std::endl;
    //}
    file.Close();
    return 0;
}
