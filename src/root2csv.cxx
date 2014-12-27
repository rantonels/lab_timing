#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"

#include <cstdlib>
#include <iostream>

const char * HNAMES [] = { "h0", "h1", "h2", "h3"};
const char * CNAMES [] = { "ch0", "ch1", "ch2", "ch3"};

int main(int argc, char** argv) {
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

    for (int faggot = 0; faggot < 4; ++faggot)
    {



        std::cout << "###" << CNAMES[faggot] << "###" << std::endl;
        TH1D histogram(HNAMES[faggot], HNAMES[faggot], 2000, 0, 2000);
        tuple->Project(HNAMES[faggot], CNAMES[faggot]);
        for (int i = 0; i < histogram.GetNbinsX(); ++i) {
            std::cout << i << " " << histogram.GetBinContent(i) << std::endl;
        }            
    }
    //for (int i = 0; i < n; ++i) {
    //int read = tuple->GetEntry(i);
    //std::cout << i << " " << ch0 << " " << ch1 << " " << ch2 << " " << ch3 << std::endl;
    //}
    file.Close();
    return 0;
}
