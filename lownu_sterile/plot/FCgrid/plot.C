void plot()
{
  char seed[10]="11";
  
  double Uee2[3], Umm2[3], dm2[3], chi2[3], nochi2[3], dchi2[3];
  int N=10000;


  double xMin = 1e-5;
  double xMax = 1.0;
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


  for(int i=0; i<N; i++) {
    ifstream f0(Form("/pnfs/dune/scratch/users/qvuong/output/stat_new/s%s/output_%d.txt",seed,i));
    ifstream f1(Form("/pnfs/dune/scratch/users/qvuong/output/flux_new/s%s/output_%d.txt",seed,i));
    ifstream f2(Form("/pnfs/dune/scratch/users/qvuong/output/sys_new/s%s/output_%d.txt",seed,i));

    if(i%100==0) std::cout << i*100./N << " percent\n";

    if(!f0 || !f1 || !f2) continue;

    f0 >> Uee2[0] >> Umm2[0] >> dm2[0] >> chi2[0] >> nochi2[0];
    f1 >> Uee2[1] >> Umm2[1] >> dm2[1] >> chi2[1] >> nochi2[1];
    f2 >> Uee2[2] >> Umm2[2] >> dm2[2] >> chi2[2] >> nochi2[2];
    
    for(int ii=0; ii<3; ii++){
      dchi2[ii] = nochi2[ii] - chi2[ii];
    }

    h0->Fill(chi2[0]);
    h1->Fill(chi2[1]);
    h2->Fill(chi2[2]);
  
    h0NO->Fill(nochi2[0]);
    h1NO->Fill(nochi2[1]);
    h2NO->Fill(nochi2[2]);

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
  lg->AddEntry(h0,"stat only");
  lg->AddEntry(h1,"stat+flux");
  lg->AddEntry(h2,"stat+flux+sig");
  lg->Draw();
  c->SaveAs(Form("new/chi2_%s.png",seed));

  TCanvas *c0 = new TCanvas("c0","",800,600);
  c0->SetGrid();
  c0->SetLogy();
  h0NO->Draw();
  h1NO->Draw("same");
  h2NO->Draw("same");
  TLegend *lg0 = new TLegend(0.55,0.70,0.9,0.9);
  lg0->AddEntry(h0NO,"stat only");
  lg0->AddEntry(h1NO,"stat+flux");
  lg0->AddEntry(h2NO,"stat+flux+sig");
  lg0->Draw();
  c0->SaveAs(Form("new/nochi2_%s.png",seed));
}
