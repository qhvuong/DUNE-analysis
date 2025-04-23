#include <DUNEStyle.h>
#include <TLegend.h>


// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}



double GetCriticalChi2(TH1D* h, double confidenceLevel) {
    if (!h) {
        std::cerr << "Histogram is null!" << std::endl;
        return -1;
    }

    // Total number of entries in the histogram
    double totalEntries = h->Integral();

    // Desired cumulative number of entries corresponding to the confidence level
    double targetEntries = confidenceLevel * totalEntries;

    // Initialize the cumulative sum
    double cumulativeEntries = 0.0;

    // Loop over the bins to find the bin where the cumulative sum exceeds targetEntries
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        cumulativeEntries += h->GetBinContent(bin);
        
        if (cumulativeEntries >= targetEntries) {
            // Found the critical bin
            double criticalChi2 = h->GetBinLowEdge(bin) + h->GetBinWidth(bin) / 2;
            return criticalChi2;
        }
    }

    // If we didn't find it, return the highest bin center as a fallback
    return h->GetBinLowEdge(h->GetNbinsX()) + h->GetBinWidth(h->GetNbinsX()) / 2;
}



void plot_from_root()
{
  const char* fnames[4] = {"all", "stat", "flux", "zeroSeeding"};
  const char* onames[4] = {"stat + flux + sys", "stat only", "stat + flux", "zeroSeeding"};

  TH1D *h[4];

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetLogy();
  c->SetGrid();
  TLegend *leg = MakeLegend(0.45, 0.55, 0.8, 0.75);

  for(int run=0; run<4; run++){
    TFile *f = new TFile(Form("FCtot_FINAL_%s.root",fnames[run]));
    h[run] = (TH1D*)f->Get("h");
    h[run]->Scale(1.0 / h[run]->Integral());
    auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, run==0 ? 0 : -1);
    h[run]->SetLineColor(color);
    h[run]->SetLineWidth(2);
    dunestyle::CenterTitles(h[run]);
    if(run==3) continue;
    h[run]->Draw(run==0 ? "hist" : "hist same");
    leg->AddEntry(h[run], Form("%s", onames[run]));
  }
/*
  c->RedrawAxis();  // otherwise the last histogram drawn overlaps with the frame
  leg->Draw();
  h[0]->SetMaximum(5.); // make some space for the watermark
  dunestyle::WIP();
  c->SaveAs("FC_uncertainties.png");

  c->Clear();
  c->SetGrid();
  h[3]->Draw("hist");
  h[0]->Draw("hist same");
  c->RedrawAxis();  // otherwise the last histogram drawn overlaps with the frame
  TLegend *leg1 = MakeLegend(0.45, 0.55, 0.8, 0.75);
  leg1->AddEntry(h[0], "seeding algorithm");
  leg1->AddEntry(h[3], "zero seeding");
  leg1->Draw();
  h[3]->SetMaximum(5.);
  dunestyle::WIP();
  c->SaveAs("FC_seedings.png");
*/

  auto color1 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 1);
  auto color2 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 2);
  auto color3 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 3);


  double level1 = GetCriticalChi2(h[0], 0.682);
  double level2 = GetCriticalChi2(h[0], 0.954);
  double level3 = GetCriticalChi2(h[0], 0.997);

  std::cout << 0.682 << "\t" << GetCriticalChi2(h[0], 0.682) << "\n";
  std::cout << 0.90 << "\t" << GetCriticalChi2(h[0], 0.90) << "\n";
  std::cout << 0.95 << "\t" << GetCriticalChi2(h[0], 0.95) << "\n";
  std::cout << 0.954 << "\t" << GetCriticalChi2(h[0], 0.954) << "\n";
  std::cout << 0.997 << "\t" << GetCriticalChi2(h[0], 0.997) << "\n";
  

  // Draw Confidence Level Lines
  TLine *line1 = new TLine(level1, 0, level1, h[0]->GetMaximum());
  TLine *line2 = new TLine(level2, 0, level2, h[0]->GetMaximum());
  TLine *line3 = new TLine(level3, 0, level3, h[0]->GetMaximum());

  line1->SetLineColor(color1);
  line2->SetLineColor(color2);
  line3->SetLineColor(color3);

  line1->SetLineStyle(kSolid);
  line2->SetLineStyle(kSolid);
  line3->SetLineStyle(kSolid);

  line1->SetLineWidth(2);
  line2->SetLineWidth(2);
  line3->SetLineWidth(2);

/*
  c->Clear();
  c->SetGrid();
  c->SetLogy();
  h[0]->Draw("HIST");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  TLegend *legend = MakeLegend(0.45, 0.55, 0.8, 0.75);
  legend->AddEntry(line1, Form("1#sigma (68.2%%) #Delta#chi^{2} = %.2f",level1), "l");
  legend->AddEntry(line2, Form("2#sigma (95.4%%) #Delta#chi^{2} = %.2f",level2), "l");
  legend->AddEntry(line3, Form("3#sigma (99.7%%) #Delta#chi^{2} = %.2f",level3), "l");
  legend->Draw();
  dunestyle::CenterTitles(h[0]);
  dunestyle::WIP();
  //c->SaveAs("FCdchi2.png");
  c->SaveAs("FC_FINAL.png");
*/
}
