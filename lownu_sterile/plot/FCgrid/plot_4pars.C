void plot_4pars()
{
  
  double ue42[3], um42[3], ut42[3], dm2[3], chi2[3], nochi2[3], dchi2[3], ratio[3];
  int N=500;


  double xMin = 1e-5;
  double xMax = 0.7;
  double mMin = 1E-3;
  double mMax = 1E3;

  int nbins = 10;
  double logXMin = TMath::Log10(xMin);
  double logXMax = TMath::Log10(xMax);
  double binWidth = (logXMax - logXMin) / nbins;
  double logMMin = TMath::Log10(mMin);
  double logMMax = TMath::Log10(mMax);
  double binWidthM = (logMMax - logMMin) / nbins;

  double binEdges[nbins + 1], binEdgesM[nbins+1];
  for (int i = 0; i <= nbins; ++i) {
    binEdges[i] = TMath::Power(10, logXMin + i * binWidth);
    binEdgesM[i] = TMath::Power(10, logMMin + i * binWidth);
  }


  TH1D *h[3], *hd[3];
  TH2D *hme[3], *hmt[3], *hdme[3], *hdmm[3], *hdmt[3];

  for(int j=0; j<3; j++){
  h[j]  = new TH1D(Form("h%d",j), "", nbins, 0, 500);
  hd[j] = new TH1D(Form("hd%d",j), "", nbins, 0, 40);
  /*
  hme[j] = new TH2D(Form("hme%s",j), "", nbins, binEdges, nbins, binEdges);
  hmt[j] = new TH2D(Form("hmt%s",j), "", nbins, binEdges, nbins, binEdges);
  hdme[j] = new TH2D(Form("hdme%s",j), "", nbins, binEdges, nbins, binEdgesM);
  hdmm[j] = new TH2D(Form("hdmm%s",j), "", nbins, binEdges, nbins, binEdgesM);
  hdmt[j] = new TH2D(Form("hdmt%s",j), "", nbins, binEdges, nbins, binEdgesM);
  */
  }

  const char data_path[] = "/pnfs/dune/scratch/users/qvuong/output/FC_4pars_new";
  const char out_path[] = "/exp/dune/app/users/qvuong/data/lownu/Feldman_Cousins";

  for(int i=0; i<N; i++) {
    if(i%100==0) std::cout << i*100./N << " percent\n";

    ifstream f0(Form("%s/alg/output_%d.txt",data_path,i));
    //ifstream f1(Form("%s/flux_alg/output_%d.txt",data_path,i));
    //ifstream f2(Form("%s/alg/output_%d.txt",data_path,i));
    if(f0) {
      f0 >> ue42[0] >> um42[0] >> ut42[0] >> dm2[0] >> chi2[0] >> nochi2[0] >> ratio[0];
      //f1 >> ue42[1] >> um42[1] >> ut42[1] >> dm2[1] >> chi2[1] >> nochi2[1] >> ratio[1];
      //f2 >> ue42[2] >> um42[2] >> ut42[2] >> dm2[2] >> chi2[2] >> nochi2[2] >> ratio[2];
      for(int j=0; j<1; j++){
      dchi2[j] = nochi2[j] - chi2[j];
      h[j]->Fill(chi2[j]);
      hd[j]->Fill(dchi2[j]);
      /*
      hme[j]->Fill(um42[j], ue42[j]);
      hmt[j]->Fill(um42[j], ut42[j]);
      hdme[j]->Fill(ue42[j], dm2[j]);
      hdmm[j]->Fill(um42[j], dm2[j]);
      hdmt[j]->Fill(ut42[j], dm2[j]);
      */
      }
      f0.close();
      //f1.close();
      //f2.close();
      
    }

  }


  //h[0]->Draw();

  h[0]->SetTitle("Best-fit #chi^{2} distribution");
  h[0]->GetXaxis()->SetTitle("#chi^{2}");
  h[0]->GetYaxis()->SetTitle("Entries");
  
  hd[0]->SetTitle("Best-fit #Delta#chi^{2} distribution");
  hd[0]->GetXaxis()->SetTitle("#Delta#chi^{2}");
  hd[0]->GetYaxis()->SetTitle("Entries");

  int colors[8] = {kRed, kBlue, kGreen, kBlack, kMagenta, kCyan, kOrange, kYellow};

  for (int i = 0; i < 1; ++i) {
    h[i]->SetLineColor(colors[i]);
    hd[i]->SetLineColor(colors[i]);
    
    h[i]->SetLineWidth(2.);
    hd[i]->SetLineWidth(2.);

    h[i]->SetStats(0);
    hd[i]->SetStats(0);
    /*
    hme[i]->SetStats(0);
    hmt[i]->SetStats(0);
    hdme[i]->SetStats(0);
    hdmm[i]->SetStats(0);
    hdmt[i]->SetStats(0);
    */
  }

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogy();
  TLegend *lg = new TLegend(0.65,0.60,0.9,0.9);
  h[0]->Draw();
  h[1]->Draw("same");
  h[2]->Draw("same");
  lg->AddEntry(h[0],Form("stats only"), "f");
  lg->AddEntry(h[1],Form("stats+flux"), "f");
  lg->AddEntry(h[2],Form("stats+flux+sys"), "f");
  lg->Draw();
  c->SaveAs(Form("%s/tot_chi2_4pars.png",out_path));
  
  TCanvas *cd = new TCanvas("cd","",800,600);
  cd->SetGrid();
  cd->SetLogy();
  TLegend *lgd = new TLegend(0.65,0.60,0.9,0.9);
  hd[0]->Draw();
  hd[2]->Draw("same");
  hd[1]->Draw("same");
  lgd->AddEntry(hd[0],Form("stats only"), "f");
  lgd->AddEntry(hd[1],Form("stats+flux"), "f");
  lgd->AddEntry(hd[2],Form("stats+flux+sys"), "f");
  lgd->Draw();
  cd->SaveAs(Form("%s/tot_dchi2_4pars.png",out_path));

}
