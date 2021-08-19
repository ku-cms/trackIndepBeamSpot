// fit 64 phi ring TGraphs

void multiple_phifit_v2() {
  
  gStyle->SetOptFit(111);
  gStyle->SetStatW(0.1);                
  gStyle->SetStatH(0.1);
  gStyle->SetStatY(0.99);
  gStyle->SetStatY(0.99);
  gStyle->SetTitleFontSize(0.04);

  //string h  = "TTBar_pileup_0p2_AllClusters";
  //string h  = "TTBar_AllClusters_zsmear";
  string h  = "SingleMuon_AllClusters";

  TFile *a = new TFile(Form("%s.root", h.c_str()), "READ");

  TF1 *f = new TF1("f", "[0]*sin(x - [1]) + [2]", -3.1415926, 3.1415926);
  f->SetParNames("amp", "shift", "avg");
  f->SetParameters(1000, 0, 1.0e4);
  f->SetParLimits(0, 0, 1e8);
  f->SetParLimits(1, -3.1415926, 3.1415926);
  
  double amp[64], s[64], avg[64];
  TGraph *g[64];
  TCanvas *c[64];
  for(int i=0; i<64; i++) {
    g[i] = (TGraph*) a->Get(Form("gr_phi_occ_ring_%d",i));
    c[i] = new TCanvas(Form("%s gr_phi_occ_ring_%d", h.c_str(), i), Form("%s gr_phi_occ_ring_%d", h.c_str(), i));
    
    g[i]->Fit(f, "R");
    g[i]->SetTitle(Form("%s gr_phi_occ_ring_%d", h.c_str(), i));
    g[i]->SetMarkerStyle(20);
    //g[i]->GetYaxis()->SetRangeUser(0.0,  20000.0);
    g[i]->GetYaxis()->SetRangeUser(0.0, 100000.0);
   
    c[i]->cd();
    g[i]->Draw("AP");

    amp[i] = f->GetParameter(0);
    s[i]   = f->GetParameter(1);
    avg[i] = f->GetParameter(2);
  }

   c[0]->Print(Form("%s.pdf(", h.c_str()));
   for(int i=1; i<63; i++) {
     c[i]->Print(Form("%s.pdf", h.c_str()));
   }
   c[63]->Print(Form("%s.pdf)", h.c_str()));
   for(int i=0; i<64; i++) {
     delete c[i];
   }

   for(int i=0; i<64; i++) {
     cout << amp[i] << endl;
   }
   cout << endl;
   for(int i=0; i<64; i++) {
     cout << s[i] << endl;
   }
   cout << endl;
   for(int i=0; i<64; i++) {
     cout << avg[i] << endl;
   }
}


