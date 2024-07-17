void plot()
{
  
  double Uee2[8], Umm2[8], dm2[8], chi2[8], nochi2[8], dchi2[8];
  double Uee2_t[8], Umm2_t[8], dm2_t[8], chi2_t[8], nochi2_t[8], dchi2_t[8];
  int N=10000;

  const char *dm[] = {"1.0", "5.0", "10.0", "100.0"};

  double xMin = 1e-3;
  double xMax = 50;

  int nBins = 100;
  double logXMin = TMath::Log10(xMin);
  double logXMax = TMath::Log10(xMax);
  double binWidth = (logXMax - logXMin) / nBins;

  double binEdges[nBins + 1];
  for (int i = 0; i <= nBins; ++i) {
    binEdges[i] = TMath::Power(10, logXMin + i * binWidth);
  }

  TH1D *h[5];
  TH1D *hd[5];

  int colors[8] = {kRed, kBlue, kGreen, kBlack, kMagenta, kCyan, kOrange, kYellow};


  for (int i = 0; i < 3; ++i) {
    h[i]    = new TH1D(Form("h%d", i), "", 100, 0, 500);
    hd[i]   = new TH1D(Form("hd%d", i), "", 100, 0, 40);
  }

  const char data_path[] = "/pnfs/dune/scratch/users/qvuong/output";
  const char out_path[] = "/exp/dune/app/users/qvuong/data/lownu/Feldman_Cousins";

  for(int i=0; i<N; i++) {
    if(i%100==0) std::cout << i*100./N << " percent\n";

    ifstream f(Form("%s/FCflux/seeding_algorithm/output_%d.txt",data_path,i));
    ifstream f1(Form("%s/FC_2pars/FCflux_dm1/output_%d.txt",data_path,i));
    ifstream f2(Form("%s/FCstat_new/seeding_algorithm/output_%d.txt",data_path,i));
/*
    ifstream f5(Form("%s/FC_2pars/FCflux_dm5/output_%d.txt",data_path,i));
    ifstream f10(Form("%s/FC_2pars/FCflux_dm10/output_%d.txt",data_path,i));
    ifstream f100(Form("%s/FC_2pars/FCflux_dm100/output_%d.txt",data_path,i));
*/
    if(!f) continue;
    //if(!f1 || !f5 || !f10 || !f100) continue;
    if(!f1 || !f2) continue;

    f >> Uee2[0] >> Umm2[0] >> dm2[0] >> chi2[0] >> nochi2[0];
/*    
    f5 >> Uee2[1] >> Umm2[1] >> dm2[1] >> chi2[1] >> nochi2[1];
    
    f10 >> Uee2[2] >> Umm2[2] >> dm2[2] >> chi2[2] >> nochi2[2];
    
    f100 >> Uee2[3] >> Umm2[3] >> dm2[3] >> chi2[3] >> nochi2[3];
*/    
    f1 >> Uee2[1] >> Umm2[1] >> dm2[1] >> chi2[1] >> nochi2[1];
    f2 >> Uee2[2] >> Umm2[2] >> dm2[2] >> chi2[2] >> nochi2[2];
    
    
    for(int ii=0; ii<3; ii++){
      dchi2[ii] = nochi2[ii] - chi2[ii];

      h[ii]->Fill(chi2[ii]);
      
      hd[ii]->Fill(dchi2[ii]);

    }
    f.close();
    f1.close();
    f2.close();
/*
    f5.close();
    f10.close();
    f100.close();
*/
  }

  
  for (int i = 0; i < 3; ++i) {
    h[i]->SetLineColor(colors[i]);
    hd[i]->SetLineColor(colors[i]);
    
    h[i]->SetLineWidth(2.);
    hd[i]->SetLineWidth(2.);

    h[i]->SetStats(0);
    hd[i]->SetStats(0);
  }

  h[0]->SetTitle("Best-fit #chi^{2} distribution");
  h[0]->GetXaxis()->SetTitle("#chi^{2}");
  h[0]->GetYaxis()->SetTitle("Entries");
  
  hd[0]->SetTitle("Best-fit #Delta#chi^{2} distribution");
  hd[0]->GetXaxis()->SetTitle("#Delta#chi^{2}");
  hd[0]->GetYaxis()->SetTitle("Entries");

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogy();
  TLegend *lg = new TLegend(0.55,0.60,0.9,0.9);
  h[0]->Draw();
  h[2]->Draw("same");
  h[1]->Draw("same");
  lg->AddEntry(h[0],Form("stats+flux"), "f");
  lg->AddEntry(h[2],Form("stats only"), "f");
  lg->AddEntry(h[1],Form("#Deltam^{2} = 1.0"), "f");
  lg->Draw();
  c->SaveAs(Form("%s/chi2_01.png",out_path));
  
  TCanvas *cd = new TCanvas("cd","",800,600);
  cd->SetGrid();
  cd->SetLogy();
  TLegend *lgd = new TLegend(0.55,0.60,0.9,0.9);
  hd[0]->Draw();
  hd[2]->Draw("same");
  hd[1]->Draw("same");
  lgd->AddEntry(hd[0],Form("stats+flux"), "f");
  lgd->AddEntry(hd[2],Form("stats only"), "f");
  lgd->AddEntry(hd[1],Form("#Deltam^{2} = 1.0"), "f");
  lgd->Draw();
  cd->SaveAs(Form("%s/dchi2_01.png",out_path));

}
