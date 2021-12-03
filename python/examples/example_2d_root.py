# example_2d_root.py

import ROOT

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)

plot_name = "plots/example_2d_hist"
c1 = ROOT.TCanvas("c1", "c1", 800, 800)
pad = c1.cd(1)
h2d = ROOT.TH2F("h2d", "h2d", 10, 0.0, 10.0, 10, 0.0, 10.0)

# fill
for i in range(0, 10, 2):
    h2d.Fill(i, i)
    h2d.Fill(i, 2 * i)
h2d.SetBinContent(5, 1, 3.0)
h2d.SetStats(0)
h2d.Draw("colz")
c1.Update()
c1.SaveAs(plot_name + ".pdf")
c1.SaveAs(plot_name + ".png")

