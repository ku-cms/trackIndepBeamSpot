// phi_fit.cpp

#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

void fit(std::string input_file, std::string input_dir, std::string plot_dir)
{
    std::cout << "Fitting for input file: " << input_file << std::endl;
    TFile *a = new TFile(Form("%s/%s.root", input_dir.c_str(), input_file.c_str()), "READ");
    TF1   *f = new TF1("f", "[0]*sin(x - [1]) + [2]", -3.14159265, 3.14159265);
    f->SetParNames("amp", "shift", "avg");
    f->SetParameters(1000, 0, 1.0e4);
    f->SetParLimits(0, 0, 1e8);
    f->SetParLimits(1, -3.14159265, 3.14159265);
    
    double amp[64], s[64], avg[64];
    TGraph  *g[64];
    TCanvas *c[64];
  
    for(int i = 0; i < 64; ++i)
    {
        g[i] = (TGraph*) a->Get(Form("gr_phi_occ_ring_%d",i));
        c[i] = new TCanvas(Form("%s gr_phi_occ_ring_%d", input_file.c_str(), i), Form("%s gr_phi_occ_ring_%d", input_file.c_str(), i));
        
        g[i]->Fit(f, "R");
        g[i]->SetTitle(Form("%s gr_phi_occ_ring_%d", input_file.c_str(), i));
        g[i]->SetMarkerStyle(20);
        //g[i]->GetYaxis()->SetRangeUser(0.0,  20000.0);
        g[i]->GetYaxis()->SetRangeUser(0.0, 100000.0);
   
        c[i]->cd();
        g[i]->Draw("AP");

        amp[i] = f->GetParameter(0);
        s[i]   = f->GetParameter(1);
        avg[i] = f->GetParameter(2);
    }
    
    // create pdf
    c[0]->Print(Form("%s/%s.pdf(", plot_dir.c_str(), input_file.c_str()));
    for(int i = 1; i < 63; ++i)
    {
        c[i]->Print(Form("%s/%s.pdf", plot_dir.c_str(), input_file.c_str()));
    }
    c[63]->Print(Form("%s/%s.pdf)", plot_dir.c_str(), input_file.c_str()));
    
    // delete canvases
    for(int i = 0; i < 64; ++i)
    {
        delete c[i];
    }
    
    delete a;
    delete f;
}

void loop()
{
    std::string input_dir = "output";
    std::string plot_dir  = "phi_fit_plots";
    std::vector<std::string> input_files;
    input_files.push_back("TTBar_AllClusters_zsmear");
    input_files.push_back("SingleMuon_AllClusters");
    //input_files.push_back("output/TTBar_AllClusters_zsmear.root");
    //input_files.push_back("output/SingleMuon_AllClusters.root");
    for (int i = 0; i < input_files.size(); ++i)
    {
        fit(input_files[i], input_dir, plot_dir);
    }
}

int main()
{
    std::cout << "Run phi fit." << std::endl;
    loop();
    return 0;
}

