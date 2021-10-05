// phi_fit.cpp

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"

namespace fs = std::filesystem;

void fit(std::string input_file, std::string input_dir, std::string plot_dir, double y_min, double y_max)
{
    std::cout << "Fitting for input file: " << input_file << std::endl;

    // create directory for plots
    fs::create_directory(plot_dir);
    
    gStyle->SetOptFit(111);
    gStyle->SetStatW(0.1);                
    gStyle->SetStatH(0.1);
    gStyle->SetStatY(0.99);
    gStyle->SetTitleFontSize(0.04);

    double pi = 3.14159265;
    
    TFile *a = new TFile(Form("%s/%s.root", input_dir.c_str(), input_file.c_str()), "READ");
    TF1   *f = new TF1("f", "[0]*sin(x - [1]) + [2]", -1 * pi, pi);
    f->SetParNames("amp", "shift", "avg");
    
    f->SetParLimits(0, 0, 1e8);
    f->SetParLimits(1, -1 * pi, pi);
    
    double amp[64], s[64], avg[64];
    TGraph  *g[64];
    TCanvas *c[64];

    std::string tag = "";
    int option = 0;

    if (option == 1)
    {
        tag = "subtracted";
    }
    else if (option == 2)
    {
        tag = "postcut";
    }
    else
    {
        tag = "standard";
    }
  
    for(int i = 0; i < 64; ++i)
    {
        // set parameter starting values for each fit
        f->SetParameters(1000, 0, 1.0e4);
        
        // standard
        //g[i] = (TGraph*) a->Get(Form("gr_phi_occ_ring_%d",i));
        //c[i] = new TCanvas(Form("%s gr_phi_occ_ring_%d", input_file.c_str(), i), Form("%s gr_phi_occ_ring_%d", input_file.c_str(), i));
        //g[i]->Fit(f, "R");
        //g[i]->SetTitle(Form("%s gr_phi_occ_ring_%d", input_file.c_str(), i));
        
        // subtracted
        //g[i] = (TGraph*) a->Get(Form("gr_phi_occ_ring_subtracted_%d",i));
        //c[i] = new TCanvas(Form("%s gr_phi_occ_ring_subtracted_%d", input_file.c_str(), i), Form("%s gr_phi_occ_ring_subtracted_%d", input_file.c_str(), i));
        //g[i]->Fit(f, "R");
        //g[i]->SetTitle(Form("%s gr_phi_occ_ring_subtracted_%d", input_file.c_str(), i));
        
        // postcut
        //g[i] = (TGraph*) a->Get(Form("gr_phi_occ_ring_postcut_%d",i));
        //c[i] = new TCanvas(Form("%s gr_phi_occ_ring_postcut_%d", input_file.c_str(), i), Form("%s gr_phi_occ_ring_postcut_%d", input_file.c_str(), i));
        //g[i]->Fit(f, "R");
        //g[i]->SetTitle(Form("%s gr_phi_occ_ring_postcut_%d", input_file.c_str(), i));
        
        //std::string ring_name = "gr_phi_occ_ring_" + tag + "_" + std::to_string(i);
        const char* ring_name = Form("gr_phi_occ_ring_%s_%d", tag.c_str(), i);
        const char* file_ring_name = Form("%s gr_phi_occ_ring_%s_%d", input_file.c_str(), tag.c_str(), i);
        std::cout << ring_name << std::endl;
        std::cout << file_ring_name << std::endl;
        
        //g[i] = (TGraph*) a->Get(Form("gr_phi_occ_ring_%s_%d", tag.c_str(), i));
        //c[i] = new TCanvas(Form("%s gr_phi_occ_ring_%s_%d", input_file.c_str(), tag.c_str(), i), Form("%s gr_phi_occ_ring_%s_%d", input_file.c_str(), tag.c_str(), i));
        //g[i]->Fit(f, "R");
        //g[i]->SetTitle(Form("%s gr_phi_occ_ring_%s_%d", input_file.c_str(), tag.c_str(), i));
        
        g[i] = (TGraph*) a->Get(ring_name);
        c[i] = new TCanvas(file_ring_name, file_ring_name);
        g[i]->Fit(f, "R");
        g[i]->SetTitle(file_ring_name);
        
        g[i]->SetMarkerStyle(20);
        g[i]->GetYaxis()->SetRangeUser(y_min, y_max);
   
        c[i]->cd();
        g[i]->Draw("AP");

        amp[i] = f->GetParameter(0);
        s[i]   = f->GetParameter(1);
        avg[i] = f->GetParameter(2);
    }
    
    // create pdf
    
    c[0]->Print(Form("%s/%s_%s.pdf(", plot_dir.c_str(), input_file.c_str(), tag.c_str()));
    for(int i = 1; i < 63; ++i)
    {
        c[i]->Print(Form("%s/%s_%s.pdf", plot_dir.c_str(), input_file.c_str(), tag.c_str()));
    }
    c[63]->Print(Form("%s/%s_%s.pdf)", plot_dir.c_str(), input_file.c_str(), tag.c_str()));
    
    //const char* file_name_first  = Form("%s/%s.pdf(", plot_dir.c_str(), input_file.c_str());
    //const char* file_name_middle = Form("%s/%s.pdf", plot_dir.c_str(), input_file.c_str());
    //const char* file_name_last   = Form("%s/%s.pdf)", plot_dir.c_str(), input_file.c_str());
    //
    //c[0]->Print(file_name_first);
    //for(int i = 1; i < 63; ++i)
    //{
    //    c[i]->Print(file_name_middle);
    //}
    //c[63]->Print(file_name_last);
    
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
    //input_files.push_back("TTBar_AllClusters_zsmear");
    //input_files.push_back("SingleMuon_AllClusters");
    
    //input_files.push_back("ZeroBias_2017B_AllClusters");
    //input_files.push_back("ZeroBias_2017B_ClusterSize2_AllClusters");
    //input_files.push_back("ZeroBias_2017B_ClusterSize2_NumberClusters2000_AllClusters");
    
    //input_files.push_back("ZeroBias_2017B_MoreEvents_AllClusters");
    //input_files.push_back("ZeroBias_2017B_MoreEvents_ClusterSize2_AllClusters");
    //input_files.push_back("SingleMuon_2017B_MoreEvents_AllClusters");
    //input_files.push_back("SingleMuon_2017B_MoreEvents_ClusterSize2_AllClusters");

    // Legacy 2017
    input_files.push_back("ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters");
    input_files.push_back("SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters");
    
    std::vector<double> y_min_vals;
    std::vector<double> y_max_vals;
    
    y_min_vals.push_back(0.0); 
    y_min_vals.push_back(0.0); 
    y_max_vals.push_back(300000.0);
    y_max_vals.push_back(300000.0);
    
    //y_min_vals.push_back(-100000.0); 
    //y_min_vals.push_back(-100000.0); 
    //y_max_vals.push_back(100000.0); 
    //y_max_vals.push_back(100000.0); 

    for (int i = 0; i < input_files.size(); ++i)
    {
        fit(input_files[i], input_dir, plot_dir, y_min_vals[i], y_max_vals[i]);
    }
}

int main()
{
    std::cout << "Run phi fits." << std::endl;
    loop();
    return 0;
}

