#include <DUNEStyle.h>
#include <TLegend.h>


// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}


void plot_3pars()
{
  
  double ue42[4], um42[4], ut42[4], dm2[4], chi2[4], nochi2[4], dchi2[4], ratio[4];
  int N = 100000;
  int count = 0;


  double xMin = 1e-5;
  double xMax = 0.7;
  double mMin = 1E-3;
  double mMax = 1E3;

  int nbins = 100;
  double logXMin = TMath::Log10(xMin);
  double logXMax = TMath::Log10(xMax);
  double binWidth = (logXMax - logXMin) / nbins;
  double logMMin = TMath::Log10(mMin);
  double logMMax = TMath::Log10(mMax);
  double binWidthM = (logMMax - logMMin) / nbins;

  double binEdges[nbins + 1], binEdgesM[nbins+1];
  for (int i = 0; i <= nbins; ++i) {
    binEdges[i] = TMath::Power(10, logXMin + i * binWidth);
    binEdgesM[i] = TMath::Power(10, logMMin + i * binWidthM);
  }


  TH1D *h[4], *hd[4];
  TH2D *hme[4], *hmt[4], *hdme[4], *hdmm[4], *hdmt[4];

  for(int j=0; j<4; j++){
  h[j]  = new TH1D(Form("h%d",j), "", nbins, 0, 500);
  hd[j] = new TH1D(Form("hd%d",j), "", nbins, 0, 50);
  hme[j] = new TH2D(Form("hme%d",j), "", nbins, binEdges, nbins, binEdges);
  hmt[j] = new TH2D(Form("hmt%d",j), "", nbins, binEdges, nbins, binEdges);
  hdme[j] = new TH2D(Form("hdme%d",j), "", nbins, binEdges, nbins, binEdgesM);
  hdmm[j] = new TH2D(Form("hdmm%d",j), "", nbins, binEdges, nbins, binEdgesM);
  hdmt[j] = new TH2D(Form("hdmt%d",j), "", nbins, binEdges, nbins, binEdgesM);
  }

  const char data_path[] = "/pnfs/dune/scratch/users/qvuong/output";
  const char out_path[] = "/exp/dune/app/users/qvuong/data/lownu/Feldman_Cousins";

  for(int i=0; i<N; i++) {
    if(i%100==0) std::cout << i*100./N << " percent\n";

    ifstream f0(Form("%s/FC_3pars/ut42_0/output_%d.txt",data_path,i));
    ifstream f1(Form("%s/FC_3pars/ut42_1/output_%d.txt",data_path,i));
    ifstream f2(Form("%s/FC_3pars/ut42_2/output_%d.txt",data_path,i));
    //ifstream f3(Form("%s/FC_3pars/test/output_%d.txt",data_path,i));
    ifstream f3(Form("%s/FC_4pars/0827/output_%d.txt",data_path,i));
    if(f0 && f1 && f2 && f3) {
      count+=1;
      f0 >> ue42[0] >> um42[0] >> ut42[0] >> dm2[0] >> chi2[0] >> nochi2[0] >> ratio[0];
      f1 >> ue42[1] >> um42[1] >> ut42[1] >> dm2[1] >> chi2[1] >> nochi2[1] >> ratio[1];
      f2 >> ue42[2] >> um42[2] >> ut42[2] >> dm2[2] >> chi2[2] >> nochi2[2] >> ratio[2];
      f3 >> ue42[3] >> um42[3] >> ut42[3] >> dm2[3] >> chi2[3] >> nochi2[3] >> ratio[3];
      //double dchi2 = nochi2[0]-chi2[0];
      //if(dchi2 > 20)
      //std::cout << i << "\t" << nochi2[0] << "\t" << chi2[0] << "\t" << ue42[0] << "\t" << um42[0] << "\t" << ut42[0] << "\t" << dm2[0] << "\n";
      
      for(int j=0; j<4; j++){
      dchi2[j] = nochi2[j] - chi2[j];
      //if(dchi2[j] > 40) std::cout << i << "\t" << j << "\t" << dchi2[j] << "\t" << ue42[j] << "\t" << um42[j] << "\t" << ut42[j] << "\t" << dm2[j] << "\n";
      h[j]->Fill(chi2[j]);
      hd[j]->Fill(dchi2[j]);
      hme[j]->Fill(um42[j], ue42[j]);
      hmt[j]->Fill(um42[j], ut42[j]);
      hdme[j]->Fill(ue42[j], dm2[j]);
      hdmm[j]->Fill(um42[j], dm2[j]);
      hdmt[j]->Fill(ut42[j], dm2[j]);
      }
      
      f0.close();
      f1.close();
      f2.close();
      f3.close();
    }

  }

  std::cout << "count = " << count << "\n";


  hd[0]->SetTitle("Feldman-Cousins #Delta#chi^{2} distribution");
  hd[0]->GetXaxis()->SetTitle("#Delta#chi^{2}");
  hd[0]->GetYaxis()->SetTitle("Number of universes");

  TCanvas c;
  c.Clear();
  c.cd();
  c.SetGrid();
  gPad->SetLogy();
  TLegend * leg = MakeLegend(0.50, 0.55, 0.90, 0.80);
  TH1 * hFirst = nullptr;
  for (std::size_t histIdx = 0; histIdx < 4; histIdx++)
  {
    //if(histIdx == 1 || histIdx == 2) continue;
    TH1D *h = hd[histIdx];
    auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, histIdx==0 ? 0 : -1);
    h->SetLineColor(color);
    h->SetFillStyle(0);
    dunestyle::CenterTitles(h);
    auto newh = h->DrawCopy(histIdx == 0 ? "" : "same");  // need to leak it so it doesn't disappear
    if (!hFirst)
      hFirst = newh;

    // we do this the hard way so the legend has the top-most histogram in the stack first
    //leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(newh, Form("true #nu < %.1f GeV", nu_cuts[histIdx]), "l"));
  }
  leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(hd[0], Form("fixed U_{#tau4}^{2} = 0.0"), "l"));
  leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(hd[1], Form("fixed U_{#tau4}^{2} = 0.3"), "l"));
  leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(hd[2], Form("fixed U_{#tau4}^{2} = 0.7"), "l"));
  leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(hd[3], Form("free U_{#tau4}^{2}"), "l"));
  c.RedrawAxis();  // otherwise the last histogram drawn overlaps with the frame
  leg->Draw();
  hFirst->SetMaximum(hFirst->GetMaximum()*3.); // make some space for the watermark
  dunestyle::WIP();
  c.SaveAs("dchi2_100bins_0909.png");






/*
  h[0]->SetTitle("Best-fit #chi^{2} distribution");
  h[0]->GetXaxis()->SetTitle("#chi^{2}");
  h[0]->GetYaxis()->SetTitle("Entries");
  


  int colors[8] = {kRed, kBlue, kGreen, kBlack, kMagenta, kCyan, kOrange, kYellow};

  for (int i = 0; i < 2; ++i) {
    h[i]->SetLineColor(colors[i]);
    hd[i]->SetLineColor(colors[i]);
    
    h[i]->SetLineWidth(2.);
    hd[i]->SetLineWidth(2.);

    h[i]->SetStats(0);
    hd[i]->SetStats(0);
    hme[i]->SetStats(0);
    hmt[i]->SetStats(0);
    hdme[i]->SetStats(0);
    hdmm[i]->SetStats(0);
    hdmt[i]->SetStats(0);
  }

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogy();
  TLegend *lg = new TLegend(0.65,0.75,0.9,0.9);
  h[0]->Draw();
  h[1]->Draw("same");
  //h[2]->Draw("same");
  //lg->AddEntry(h[1],Form("stats+flux"), "f");
  lg->AddEntry(h[0],Form("fixed Ut42 = 0.1"), "f");
  lg->AddEntry(h[1],Form("free Ut42"), "f");
  lg->Draw();
  c->SaveAs(Form("%s/0815_chi2.png",out_path));
  
  TCanvas *cd = new TCanvas("cd","",800,600);
  cd->SetGrid();
  cd->SetLogy();
  TLegend *lgd = new TLegend(0.65,0.75,0.9,0.9);
  hd[0]->Draw();
  //hd[2]->Draw("same");
  hd[1]->Draw("same");
  //lgd->AddEntry(hd[0],Form("stats only"), "f");
  //lgd->AddEntry(hd[1],Form("stats+flux"), "f");  
  lgd->AddEntry(hd[0],Form("fixed Ut42 = 0.1"), "f");
  lgd->AddEntry(hd[1],Form("free Ut42"), "f");
  lgd->Draw();
  cd->SaveAs(Form("%s/0815_dchi2.png",out_path));
*/

}
