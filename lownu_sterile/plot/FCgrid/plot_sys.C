void plot_sys()
{
  double Uee2[7], Umm2[7], dm2[7], chi2[7], nochi2[7], dchi2[7];
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

  TH1D *h00 = new TH1D("h00", "", 100,0,1000);
  TH1D *h11 = new TH1D("h11", "", 100,0,1000);
  TH1D *h12 = new TH1D("h12", "", 100,0,1000);
  TH1D *h13 = new TH1D("h13", "", 100,0,1000);
  TH1D *h21 = new TH1D("h21", "", 100,0,1000);
  TH1D *h22 = new TH1D("h22", "", 100,0,1000);
  TH1D *h23 = new TH1D("h23", "", 100,0,1000);
  
  TH1D *h00d = new TH1D("h00d", "", 100,0,100);
  TH1D *h11d = new TH1D("h11d", "", 100,0,100);
  TH1D *h12d = new TH1D("h12d", "", 100,0,100);
  TH1D *h13d = new TH1D("h13d", "", 100,0,100);
  TH1D *h21d = new TH1D("h21d", "", 100,0,100);
  TH1D *h22d = new TH1D("h22d", "", 100,0,100);
  TH1D *h23d = new TH1D("h23d", "", 100,0,100);
/*
  TH1D *h00NO = new TH1D("h00NO", "", 100,0,1000);
  TH1D *h11NO = new TH1D("h11NO", "", 100,0,1000);
  TH1D *h12NO = new TH1D("h12NO", "", 100,0,1000);
  TH1D *h13NO = new TH1D("h13NO", "", 100,0,1000);
  TH1D *h21NO = new TH1D("h21NO", "", 100,0,1000);
  TH1D *h22NO = new TH1D("h22NO", "", 100,0,1000);
  TH1D *h23NO = new TH1D("h23NO", "", 100,0,1000);
*/

  for(int i=0; i<N; i++) {
    ifstream f00(Form("/pnfs/dune/scratch/users/qvuong/output/sys_new/s00/output_%d.txt",i));

    ifstream f11(Form("/pnfs/dune/scratch/users/qvuong/output/sys_new/s11/output_%d.txt",i));
    ifstream f12(Form("/pnfs/dune/scratch/users/qvuong/output/sys_new/s12/output_%d.txt",i));
    ifstream f13(Form("/pnfs/dune/scratch/users/qvuong/output/sys_new/s13/output_%d.txt",i));

    ifstream f21(Form("/pnfs/dune/scratch/users/qvuong/output/sys_new/s21/output_%d.txt",i));
    ifstream f22(Form("/pnfs/dune/scratch/users/qvuong/output/sys_new/s22/output_%d.txt",i));
    ifstream f23(Form("/pnfs/dune/scratch/users/qvuong/output/sys_new/s23/output_%d.txt",i));

    if(i%100==0) std::cout << i*100./N << " percent\n";

    //if(!f00 || !f11 || !f12 || !f13 || !f21 || !f22 || !f23) continue;
    //if(!f00 || !f11 || !f12 || !f13) continue;

    if(f00) f00 >> Uee2[0] >> Umm2[0] >> dm2[0] >> chi2[0] >> nochi2[0];
    if(f11) f11 >> Uee2[1] >> Umm2[1] >> dm2[1] >> chi2[1] >> nochi2[1];
    if(f12) f12 >> Uee2[2] >> Umm2[2] >> dm2[2] >> chi2[2] >> nochi2[2];
    if(f13) f13 >> Uee2[3] >> Umm2[3] >> dm2[3] >> chi2[3] >> nochi2[3];

    if(f21) f21 >> Uee2[4] >> Umm2[4] >> dm2[4] >> chi2[4] >> nochi2[4];
    if(f22) f22 >> Uee2[5] >> Umm2[5] >> dm2[5] >> chi2[5] >> nochi2[5];
    if(f23) f23 >> Uee2[6] >> Umm2[6] >> dm2[6] >> chi2[6] >> nochi2[6];
    
    for(int ii=0; ii<7; ii++){
      dchi2[ii] = nochi2[ii] - chi2[ii];
    }

    h00->Fill(chi2[0]);
    h11->Fill(chi2[1]);
    h12->Fill(chi2[2]);
    h13->Fill(chi2[3]);
    h21->Fill(chi2[4]);
    h22->Fill(chi2[5]);
    h23->Fill(chi2[6]);

    h00d->Fill(dchi2[0]);
    h11d->Fill(dchi2[1]);
    h12d->Fill(dchi2[2]);
    h13d->Fill(dchi2[3]);
    h21d->Fill(dchi2[4]);
    h22d->Fill(dchi2[5]);
    h23d->Fill(dchi2[6]);
/*  
    h00NO->Fill(nochi2[0]);
    h11NO->Fill(nochi2[1]);
    h12NO->Fill(nochi2[2]);
    h13NO->Fill(nochi2[3]);
    h21NO->Fill(nochi2[4]);
    h22NO->Fill(nochi2[5]);
    h23NO->Fill(nochi2[6]);
*/
    f00.close();
    f11.close();
    f12.close();
    f13.close();
    f21.close();
    f22.close();
    f23.close();
  }

  h00->SetStats(0);
  h00->SetLineColor(kBlack);
  h00->SetLineWidth(2.);
  h00->SetTitle("Best-fit #chi^{2} distribution");
  h00->GetXaxis()->SetTitle("#chi^{2}");
  h00->GetYaxis()->SetTitle("Entries");
  
  h00d->SetStats(0);
  h00d->SetLineColor(kBlack);
  h00d->SetLineWidth(2.);
  h00d->SetTitle("Best-fit #Delta#chi^{2} distribution");
  h00d->GetXaxis()->SetTitle("#Delta#chi^{2}");
  h00d->GetYaxis()->SetTitle("Entries");

  h11->SetStats(0);
  h11->SetLineColor(kBlue);
  h11->SetLineWidth(2.);
  h12->SetStats(0);
  h12->SetLineColor(kRed);
  h12->SetLineWidth(2.);
  h13->SetStats(0);
  h13->SetLineColor(kGreen);
  h13->SetLineWidth(2.);
  
  h21->SetStats(0);
  h21->SetLineColor(kBlue);
  h21->SetLineWidth(2.);
  h22->SetStats(0);
  h22->SetLineColor(kRed);
  h22->SetLineWidth(2.);
  h23->SetStats(0);
  h23->SetLineColor(kGreen);
  h23->SetLineWidth(2.);
  
  
  h11d->SetStats(0);
  h11d->SetLineColor(kBlue);
  h11d->SetLineWidth(2.);
  h12d->SetStats(0);
  h12d->SetLineColor(kRed);
  h12d->SetLineWidth(2.);
  h13d->SetStats(0);
  h13d->SetLineColor(kGreen);
  h13d->SetLineWidth(2.);

  h21d->SetStats(0);
  h21d->SetLineColor(kBlue);
  h21d->SetLineWidth(2.);
  h22d->SetStats(0);
  h22d->SetLineColor(kRed);
  h22d->SetLineWidth(2.);
  h23d->SetStats(0);
  h23d->SetLineColor(kGreen);
  h23d->SetLineWidth(2.);

/*
  h00NO->SetStats(0);
  h00NO->SetLineColor(kBlack);
  h00NO->SetLineWidth(2.);
  h00NO->SetTitle("#chi^{2} at no oscillation distribution");
  h00NO->GetXaxis()->SetTitle("#chi^{2} NO");
  h00NO->GetYaxis()->SetTitle("Entries");

  h11NO->SetStats(0);
  h11NO->SetLineColor(kBlue);
  h11NO->SetLineWidth(2.);
  h12NO->SetStats(0);
  h12NO->SetLineColor(kRed);
  h12NO->SetLineWidth(2.);
  h13NO->SetStats(0);
  h13NO->SetLineColor(kGreen);
  h13NO->SetLineWidth(2.);
  
  h21NO->SetStats(0);
  h21NO->SetLineColor(kBlue);
  h21NO->SetLineWidth(2.);
  h22NO->SetStats(0);
  h22NO->SetLineColor(kRed);
  h22NO->SetLineWidth(2.);
  h23NO->SetStats(0);
  h23NO->SetLineColor(kGreen);
  h23NO->SetLineWidth(2.);
*/
/*
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogy();
  h00->Draw();
  h11->Draw("same");
  h12->Draw("same");
  h13->Draw("same");
  TLegend *lg = new TLegend(0.55,0.70,0.9,0.9);
  lg->AddEntry(h00,"seeds = (0.0, 0.0, 0.0)");
  lg->AddEntry(h11,"seeds = (0.01, 0.01, 1.0)");
  lg->AddEntry(h12,"seeds = (0.01, 0.01, 10.0)");
  lg->AddEntry(h13,"seeds = (0.01, 0.01, 100.0)");
  lg->Draw();
  c->SaveAs(Form("new/chi2_1.png"));
*/
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->SetGrid();
  c1->SetLogy();
  h00->Draw();
  h11->Draw("same");
  h12->Draw("same");
  h13->Draw("same");
  TLegend *lg1 = new TLegend(0.55,0.70,0.9,0.9);
  lg1->AddEntry(h00,"seeds = (0.0, 0.0, 0.0)");
  lg1->AddEntry(h11,"seeds = (0.01, 0.01, 1.0)");
  lg1->AddEntry(h12,"seeds = (0.01, 0.01, 10.0)");
  lg1->AddEntry(h13,"seeds = (0.01, 0.01, 100.0)");
  lg1->Draw();
  c1->SaveAs(Form("new/chi2_1.png"));

  TCanvas *c1d = new TCanvas("c1d","",800,600);
  c1d->SetGrid();
  c1d->SetLogy();
  h00d->Draw();
  h11d->Draw("same");
  h12d->Draw("same");
  h13d->Draw("same");
  TLegend *lg1d = new TLegend(0.55,0.70,0.9,0.9);
  lg1d->AddEntry(h00d,"seeds = (0.0, 0.0, 0.0)");
  lg1d->AddEntry(h11d,"seeds = (0.01, 0.01, 1.0)");
  lg1d->AddEntry(h12d,"seeds = (0.01, 0.01, 10.0)");
  lg1d->AddEntry(h13d,"seeds = (0.01, 0.01, 100.0)");
  lg1d->Draw();
  c1d->SaveAs(Form("new/dchi2_1.png"));
  

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetGrid();
  c2->SetLogy();
  h00->Draw();
  h21->Draw("same");
  h22->Draw("same");
  h23->Draw("same");
  TLegend *lg2 = new TLegend(0.55,0.70,0.9,0.9);
  lg2->AddEntry(h00,"seeds = (0.0, 0.0, 0.0)");
  lg2->AddEntry(h21,"seeds = (0.05, 0.05, 1.0)");
  lg2->AddEntry(h22,"seeds = (0.05, 0.05, 10.0)");
  lg2->AddEntry(h23,"seeds = (0.05, 0.05, 100.0)");
  lg2->Draw();
  c2->SaveAs(Form("new/chi2_2.png"));

  TCanvas *c2d = new TCanvas("c2d","",800,600);
  c2d->SetGrid();
  c2d->SetLogy();
  h00d->Draw();
  h21d->Draw("same");
  h22d->Draw("same");
  h23d->Draw("same");
  TLegend *lg2d = new TLegend(0.55,0.70,0.9,0.9);
  lg2d->AddEntry(h00d,"seeds = (0.0, 0.0, 0.0)");
  lg2d->AddEntry(h21d,"seeds = (0.05, 0.05, 1.0)");
  lg2d->AddEntry(h22d,"seeds = (0.05, 0.05, 10.0)");
  lg2d->AddEntry(h23d,"seeds = (0.05, 0.05, 100.0)");
  lg2d->Draw();
  c2d->SaveAs(Form("new/dchi2_2.png"));
/*
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetGrid();
  c2->SetLogy();
  h00NO->Draw();
  h11NO->Draw("same");
  h12NO->Draw("same");
  h13NO->Draw("same");
  TLegend *lg2 = new TLegend(0.55,0.70,0.9,0.9);
  lg2->AddEntry(h00NO,"seeds = (0.0, 0.0, 0.0)");
  lg2->AddEntry(h11NO,"seeds = (0.01, 0.01, 1.0)");
  lg2->AddEntry(h12NO,"seeds = (0.01, 0.01, 10.0)");
  lg2->AddEntry(h13NO,"seeds = (0.01, 0.01, 100.0)");
  lg2->Draw();
  c2->SaveAs(Form("new/nochi2_1.png"));
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetGrid();
  c3->SetLogy();
  h00NO->Draw();
  h21NO->Draw("same");
  h22NO->Draw("same");
  h23NO->Draw("same");
  TLegend *lg3 = new TLegend(0.55,0.70,0.9,0.9);
  lg3->AddEntry(h00NO,"seeds = (0.0, 0.0, 0.0)");
  lg3->AddEntry(h21NO,"seeds = (0.05, 0.05, 1.0)");
  lg3->AddEntry(h22NO,"seeds = (0.05, 0.05, 10.0)");
  lg3->AddEntry(h23NO,"seeds = (0.05, 0.05, 100.0)");
  lg3->Draw();
  c3->SaveAs(Form("new/nochi2_2.png"));
*/
}
