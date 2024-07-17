void plotContour()
{
  double Uee2[3], Umm2[3], dm2[3], chi2[3], nochi2[3], bf_Uee2[3], bf_Umm2[3], bf_dm2[3], dchi2[3];
  int N=10000;


  double xMin = 1e-4;
  double xMax = 0.3;
  double mMin = 1e-4;
  double mMax = 1e+3;

  int nBins = 100;
  double logXMin = TMath::Log10(xMin);
  double logXMax = TMath::Log10(xMax);
  double binWidth = (logXMax - logXMin) / nBins;

  double logmMin = TMath::Log10(mMin);
  double logmMax = TMath::Log10(mMax);
  double mbinWidth = (logmMax - logmMin) / nBins;

  std::vector<double> binEdges(nBins + 1, 0);
  std::vector<double> mbinEdges(nBins + 1, 0);
  for (int i = 0; i <= nBins; ++i) {
    binEdges[i] = TMath::Power(10, logXMin + i * binWidth);
    mbinEdges[i] = TMath::Power(10, logmMin + i * mbinWidth);
  }

  TH1D *h0 = new TH1D("h0", "", 100,0,1000);
  TH1D *h1 = new TH1D("h1", "", 100,0,1000);
  TH1D *h2 = new TH1D("h2", "", 100,0,1000);
  TH1D *h0NO = new TH1D("h0NO", "", 100,0,1000);
  TH1D *h1NO = new TH1D("h1NO", "", 100,0,1000);
  TH1D *h2NO = new TH1D("h2NO", "", 100,0,1000);

  TH2D *hdm1   = new TH2D("hdm1",  "", nBins, &binEdges[0], nBins, &binEdges[0]);
  TH2D *hdm10  = new TH2D("hdm10", "", nBins, &binEdges[0], nBins, &binEdges[0]);
  TH2D *hdm100 = new TH2D("hdm100","", nBins, &binEdges[0], nBins, &binEdges[0]);

  const char data_path[] = "/pnfs/dune/scratch/users/qvuong/output/contour_grid";
  const char out_path[] = "/exp/dune/app/users/qvuong/data/lownu/contours";

  for(int i=0; i<N; i++) {
    ifstream f0(Form("%s/dm1/output_%d.txt",data_path,i));
    ifstream f1(Form("%s/dm10/output_%d.txt",data_path,i));
    ifstream f2(Form("%s/dm100/output_%d.txt",data_path,i));

    if(i%100==0) std::cout << i*100./N << " percent\n";

    if(f0) {
      f0 >> Uee2[0] >> Umm2[0] >> dm2[0] >> bf_Uee2[0] >> bf_Umm2[0] >> bf_dm2[0] >> chi2[0] >> nochi2[0];
      dchi2[0] = nochi2[0] - chi2[0];
      h0->Fill(chi2[0]);
      h0NO->Fill(nochi2[0]);
      hdm1->Fill(Umm2[0], Uee2[0], dchi2[0]);
    }

    if(f1) {
      f1 >> Uee2[1] >> Umm2[1] >> dm2[1] >> bf_Uee2[1] >> bf_Umm2[1] >> bf_dm2[1] >> chi2[1] >> nochi2[1];
      dchi2[1] = nochi2[1] - chi2[1];
      h1->Fill(chi2[1]);
      h1NO->Fill(nochi2[1]);
      hdm10->Fill(Umm2[1], Uee2[1], dchi2[1]);
    }

    if(f2) {
      f2 >> Uee2[2] >> Umm2[2] >> dm2[2] >> bf_Uee2[2] >> bf_Umm2[2] >> bf_dm2[2] >> chi2[2] >> nochi2[2];
      dchi2[2] = nochi2[2] - chi2[2];
      h2->Fill(chi2[2]);
      h2NO->Fill(nochi2[2]);
      hdm100->Fill(Umm2[2], Uee2[2], dchi2[2]);
    }

/*
    if(!f0 || !f1 || !f2) continue;

    f0 >> Uee2[0] >> Umm2[0] >> dm2[0] >> chi2[0] >> nochi2[0];
    f1 >> Uee2[1] >> Umm2[1] >> dm2[1] >> chi2[1] >> nochi2[1];
    f2 >> Uee2[2] >> Umm2[2] >> dm2[2] >> chi2[2] >> nochi2[2];
*/    

    f0.close();
    f1.close();
    f2.close();
  }

  h0->SetStats(0);
  h0->SetLineColor(kBlack);
  h0->SetLineWidth(2.);
  h1->SetStats(0);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(2.);
  h2->SetStats(0);
  h2->SetLineColor(kRed);
  h2->SetLineWidth(2.);
  h0->SetTitle("Best-fit #chi^{2} distribution");
  h0->GetXaxis()->SetTitle("#chi^{2}");
  h0->GetYaxis()->SetTitle("Entries");

  h0NO->SetStats(0);
  h0NO->SetLineColor(kBlack);
  h0NO->SetLineWidth(2.);
  h1NO->SetStats(0);
  h1NO->SetLineColor(kBlue);
  h1NO->SetLineWidth(2.);
  h2NO->SetStats(0);
  h2NO->SetLineColor(kRed);
  h2NO->SetLineWidth(2.);
  h0NO->SetTitle("#chi^{2} at no oscillation distribution");
  h0NO->GetXaxis()->SetTitle("#chi^{2} NO");
  h0NO->GetYaxis()->SetTitle("Entries");

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogy();
  h0->Draw();
  h1->Draw("same");
  h2->Draw("same");
  TLegend *lg = new TLegend(0.55,0.70,0.9,0.9);
  lg->AddEntry(h0,"#Deltam^{2}=1.0 eV^{2}");
  lg->AddEntry(h1,"#Deltam^{2}=10.0 eV^{2}");
  lg->AddEntry(h2,"#Deltam^{2}=100.0 eV^{2}");
  lg->Draw();
  c->SaveAs(Form("%s/chi2.png",out_path));

  TCanvas *c0 = new TCanvas("c0","",800,600);
  c0->SetGrid();
  c0->SetLogy();
  h0NO->Draw();
  h1NO->Draw("same");
  h2NO->Draw("same");
  TLegend *lg0 = new TLegend(0.55,0.70,0.9,0.9);
  lg0->AddEntry(h0NO,"#Deltam^{2}=1.0 eV^{2}");
  lg0->AddEntry(h1NO,"#Deltam^{2}=10.0 eV^{2}");
  lg0->AddEntry(h2NO,"#Deltam^{2}=100.0 eV^{2}");
  lg0->Draw();
  c0->SaveAs(Form("%s/nochi2.png",out_path));

  TFile *fout = new TFile(Form("%s/contour.root",out_path), "RECREATE");
  hdm1->Write();
  hdm10->Write();
  hdm100->Write();
  fout->Close();
}
