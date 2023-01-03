// phi_fit.cpp

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <math.h>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"

namespace fs = std::filesystem;

// draw and save histogram
void draw(TH1F &hist, std::string output_name, std::string x_label, std::string y_label, double x_min, double x_max, double y_min, double y_max)
{
    TCanvas c;
    c.cd();
    hist.GetXaxis()->SetTitle(x_label.c_str());
    hist.GetYaxis()->SetTitle(y_label.c_str());
    hist.GetXaxis()->SetRangeUser(x_min, x_max);
    hist.GetYaxis()->SetRangeUser(y_min, y_max);
    hist.Draw();
    c.SaveAs(Form("%s.pdf", output_name.c_str()));
}

// run fit
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
    
    double chisq[64];
    double amp[64];
    double shift[64];
    double offset[64];
    double amp_err[64];
    double shift_err[64];
    double offset_err[64];

    //printf("M_PI = %f\n", M_PI);
    
    TGraph  *g[64];
    TCanvas *c[64];
    TH1F    *h_dummy[64];
    
    TH1F *h_num_phi     = new TH1F("h_num_phi",     "Number of phi points", 64, 0, 64);
    TH1F *h_chisq       = new TH1F("h_chisq",       "Fit chi squares", 64, 0, 64);
    TH1F *h_amp         = new TH1F("h_amp",         "Fit amplitudes", 64, 0, 64);
    TH1F *h_shift       = new TH1F("h_shift",       "Fit phi shifts", 64, 0, 64);
    TH1F *h_offset      = new TH1F("h_offset",      "Fit offsets", 64, 0, 64);
    TH1F *h_num_phi_cut = new TH1F("h_num_phi_cut", "Number of phi points after cut", 64, 0, 64);
    TH1F *h_chisq_cut   = new TH1F("h_chisq_cut",   "Fit chi squares after cut", 64, 0, 64);
    TH1F *h_amp_cut     = new TH1F("h_amp_cut",     "Fit amplitudes after cut", 64, 0, 64);
    TH1F *h_shift_cut   = new TH1F("h_shift_cut",   "Fit phi shifts after cut", 64, 0, 64);
    TH1F *h_offset_cut  = new TH1F("h_offset_cut",  "Fit offsets after cut", 64, 0, 64);
    
    TFile *a = new TFile(Form("%s/%s.root", input_dir.c_str(), input_file.c_str()), "READ");
    TF1   *f = new TF1("f", "[0]*sin(x - [1]) + [2]", -M_PI, M_PI);
    
    // set parameter names and limits
    f->SetParNames("amp", "shift", "offset");
    f->SetParLimits(0, 0, 1e10);
    f->SetParLimits(1, -M_PI, M_PI);
    f->SetParLimits(2, 0, 1e10);

    std::string tag = "";
    int option = 2;

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

    std::string base_name = plot_dir + "/" + input_file + "_" + tag;
    
    // number of phi after cut
    TGraph *num_phi = (TGraph*) a->Get("gr_num_phi_per_ring");
    if (num_phi) printf("PASS: loaded num_phi\n");
    else         printf("FAIL: did not load num_phi\n");

    for(int i = 0; i < 64; ++i)
    {
        // set parameter starting values for each fit
        // starting values for 2017B data
        //f->SetParameters(10000, 1.0, 1.0e5);
        // starting values for 2022F data
        f->SetParameters(10000, 1.0, 1.0e5);
        
        const char* ring_name = Form("gr_phi_occ_ring_%s_%d", tag.c_str(), i);
        const char* file_ring_name = Form("%s gr_phi_occ_ring_%s_%d", input_file.c_str(), tag.c_str(), i);
        
        c[i] = new TCanvas(file_ring_name, file_ring_name);
        g[i] = (TGraph*) a->Get(ring_name);
        h_dummy[i] = new TH1F("h_dummy", "h_dummy", 2, -M_PI, M_PI);
        g[i]->Fit(f, "R");
        g[i]->SetTitle(ring_name);
        
        // test dummy
        //h_dummy[i]->SetBinContent(1, 80000.0);
        //h_dummy[i]->SetBinContent(2, 80000.0);

        c[i]->cd();
        g[i]->Draw("AP");
        
        // TODO: TGraph auto zooms in on x-axis
        // may need to use dummy hist to set x-axis range to [-pi, pi]
        g[i]->SetMarkerStyle(20);
        g[i]->GetXaxis()->SetRangeUser(-M_PI, M_PI);
        g[i]->GetYaxis()->SetRangeUser(y_min, y_max);
        //g[i]->GetXaxis()->SetRangeUser(-10, 10);
        //g[i]->GetYaxis()->SetRangeUser(0, 1e6);
        //h_dummy[i]->GetXaxis()->SetRangeUser(-M_PI, M_PI);
        //h_dummy[i]->GetYaxis()->SetRangeUser(y_min, y_max);
   
        //c[i]->cd();
        //g[i]->Draw("AP");
        //h_dummy[i]->Draw("hist same");

        g[i]->Draw("AP");
        c[i]->Update();

        double num_phi_x = num_phi->GetPointX(i);
        double num_phi_y = num_phi->GetPointY(i);
        printf("ring %d: number of phi (x, y) = (%f, %f)\n", i, num_phi_x, num_phi_y);
        h_num_phi->SetBinContent(i + 1, num_phi_y);

        chisq[i]        = f->GetChisquare();
        amp[i]          = f->GetParameter(0);
        shift[i]        = f->GetParameter(1);
        offset[i]       = f->GetParameter(2);
        amp_err[i]      = f->GetParError(0);
        shift_err[i]    = f->GetParError(1);
        offset_err[i]   = f->GetParError(2);

        h_chisq->SetBinContent(i + 1,  chisq[i]);
        h_amp->SetBinContent(i + 1,    amp[i]);
        h_shift->SetBinContent(i + 1,  shift[i]);
        h_offset->SetBinContent(i + 1, offset[i]);
        h_amp->SetBinError(i + 1,      amp_err[i]);
        h_shift->SetBinError(i + 1,    shift_err[i]);
        h_offset->SetBinError(i + 1,   offset_err[i]);
        
        if(num_phi_y >= 9)
        {
            h_num_phi_cut->SetBinContent(i + 1, num_phi_y);
            h_chisq_cut->SetBinContent(i + 1,   chisq[i]);
            h_amp_cut->SetBinContent(i + 1,     amp[i]);
            h_shift_cut->SetBinContent(i + 1,   shift[i]);
            h_offset_cut->SetBinContent(i + 1,  offset[i]);
            h_amp_cut->SetBinError(i + 1,       amp_err[i]);
            h_shift_cut->SetBinError(i + 1,     shift_err[i]);
            h_offset_cut->SetBinError(i + 1,    offset_err[i]);
        }
    }

    
    // create pdf
    c[0]->Print(Form("%s.pdf(", base_name.c_str()));
    for(int i = 1; i < 63; ++i)
    {
        c[i]->Print(Form("%s.pdf", base_name.c_str()));
    }
    c[63]->Print(Form("%s.pdf)", base_name.c_str()));
    
    // draw histograms

    // limits for ttbar
    //draw(*h_chisq,  base_name + "_chisq",  "ring", "chi sq.",   0, 64, 0, 1e6);
    //draw(*h_amp,    base_name + "_amp",    "ring", "amplitude", 0, 64, 0, 2e3);
    //draw(*h_shift,  base_name + "_shift",  "ring", "shift",     0, 64, -M_PI, M_PI);
    //draw(*h_offset, base_name + "_offset", "ring", "offset",    0, 64, 0, 2e4);
    
    // limits for legacy 2017 data
    //draw(*h_num_phi,        base_name + "_num_phi",     "ring", "num_phi",   0, 64, 0, 20);
    //draw(*h_chisq,          base_name + "_chisq",       "ring", "chi sq.",   0, 64, 0, 2e10);
    //draw(*h_amp,            base_name + "_amp",         "ring", "amplitude", 0, 64, 0, 5e4);
    //draw(*h_shift,          base_name + "_shift",       "ring", "shift",     0, 64, -M_PI, M_PI);
    //draw(*h_offset,         base_name + "_offset",      "ring", "offset",    0, 64, 0, 3e5);
    
    //draw(*h_num_phi_cut,    base_name + "_num_phi_cut", "ring", "num_phi",   0, 64, 0, 20);
    //draw(*h_chisq_cut,      base_name + "_chisq_cut",   "ring", "chi sq.",   0, 64, 0, 2e10);
    //draw(*h_amp_cut,        base_name + "_amp_cut",     "ring", "amplitude", 0, 64, 0, 5e4);
    //draw(*h_shift_cut,      base_name + "_shift_cut",   "ring", "shift",     0, 64, -M_PI, M_PI);
    //draw(*h_offset_cut,     base_name + "_offset_cut",  "ring", "offset",    0, 64, 0, 3e5);
    
    // limits for legacy 2022 data (v1)
    draw(*h_num_phi,        base_name + "_num_phi",     "ring", "num_phi",   0, 64, 0, 20);
    draw(*h_chisq,          base_name + "_chisq",       "ring", "chi sq.",   0, 64, 0, 1e5);
    draw(*h_amp,            base_name + "_amp",         "ring", "amplitude", 0, 64, 0, 1e3);
    draw(*h_shift,          base_name + "_shift",       "ring", "shift",     0, 64, -M_PI, M_PI);
    draw(*h_offset,         base_name + "_offset",      "ring", "offset",    0, 64, 0, 1e4);
    
    draw(*h_num_phi_cut,    base_name + "_num_phi_cut", "ring", "num_phi",   0, 64, 0, 20);
    draw(*h_chisq_cut,      base_name + "_chisq_cut",   "ring", "chi sq.",   0, 64, 0, 1e5);
    draw(*h_amp_cut,        base_name + "_amp_cut",     "ring", "amplitude", 0, 64, 0, 1e3);
    draw(*h_shift_cut,      base_name + "_shift_cut",   "ring", "shift",     0, 64, -M_PI, M_PI);
    draw(*h_offset_cut,     base_name + "_offset_cut",  "ring", "offset",    0, 64, 0, 1e4);
    
    // delete canvases
    for(int i = 0; i < 64; ++i)
    {
        delete c[i];
        delete h_dummy[i];
    }
    
    delete a;
    delete f;
    delete h_num_phi;
    delete h_chisq;
    delete h_amp;
    delete h_shift;
    delete h_offset;
    delete h_num_phi_cut;
    delete h_chisq_cut;
    delete h_amp_cut;
    delete h_shift_cut;
    delete h_offset_cut;
}

void loop()
{
    std::string input_dir = "output";
    std::string plot_dir  = "phi_fit_plots";
    
    std::vector<std::string> input_files;
    // input files: do not include ".root", as this is added later
    
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
    //input_files.push_back("ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters");
    //input_files.push_back("SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters");
    // 2022 (v1)
    input_files.push_back("ZeroBias_2022F_nFiles1_NoCuts_weighted_avg");
    input_files.push_back("ZeroBias_2022F_nFiles1_ClustSize2_nClust2000_weighted_avg");
    
    std::vector<double> y_min_vals;
    std::vector<double> y_max_vals;
    
    // limits for ttbar
    //y_min_vals.push_back(0.0); 
    //y_max_vals.push_back(30000.0);
    
    // limits for legacy 2017 data
    //y_min_vals.push_back(0.0); 
    //y_min_vals.push_back(0.0); 
    //y_max_vals.push_back(300000.0);
    //y_max_vals.push_back(300000.0);
    // limits for 2022 data (v1)
    y_min_vals.push_back(0.0); 
    y_min_vals.push_back(0.0); 
    y_max_vals.push_back(5000.0);
    y_max_vals.push_back(5000.0);
    
    // limits for subtracted version
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

