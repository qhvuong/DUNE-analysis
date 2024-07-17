void plot()
{
  char seed_alg[]="seeding_algorithm";
  char seed_true[]="seeding_true";
  
  double Uee2[8], Umm2[8], dm2[8], chi2[8], nochi2[8], dchi2[8];
  double Uee2_t[8], Umm2_t[8], dm2_t[8], chi2_t[8], nochi2_t[8], dchi2_t[8];
  int N=100000;


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

  TH1D *h[8];
  TH1D *h_t[8];
  TH1D *hd[8];
  TH1D *hd_t[8];
  TH2D *hchi2_cor[8];

  int colors[8] = {kRed, kBlue, kGreen, kBlack, kMagenta, kCyan, kOrange, kYellow};

  std::vector<std::string> names = {"stats only", "stats+flux", "stats+flux+sig", "stats+flux+CCQE", 
                                    "stats+flux+BeRPA_A", "stats+flux+BeRPA_B", "stats+flux+BeRPA_D", "stats+flux+MnGaus"};


  for (int i = 0; i < 8; ++i) {
    h[i]    = new TH1D(Form("h%d", i), "", 100, 0, 500);
    h_t[i]  = new TH1D(Form("h_t%d", i), "", 100, 0, 500);
    hd[i]   = new TH1D(Form("hd%d", i), "", 100, 0, 40);
    hd_t[i] = new TH1D(Form("hd_t%d", i), "", 100, 0, 40);
    hchi2_cor[i] = new TH2D(Form("hchi2_cor%d",i), "", nBins, binEdges, nBins, binEdges);}

  const char data_path[] = "/pnfs/dune/scratch/users/qvuong/output";
  const char out_path[] = "/exp/dune/app/users/qvuong/data/lownu/Feldman_Cousins";

  int j=0;

  for(int i=0; i<N; i++) {
    if(i%100==0) std::cout << i*100./N << " percent\n";

    ifstream fstat_true(Form("%s/FCstat_new/%s/output_%d.txt",data_path,seed_true,i));
    ifstream fstat_alg(Form("%s/FCstat_new/%s/output_%d.txt",data_path,seed_alg,i));
    ifstream fflux_true(Form("%s/FCflux/%s/output_%d.txt",data_path,seed_true,i));
    ifstream fflux_alg(Form("%s/FCflux/%s/output_%d.txt",data_path,seed_alg,i));
/*
    ifstream fgrid_true(Form("%s/FCgrid_new/%s/output_%d.txt",data_path,seed_true,i));
    ifstream fgrid_alg(Form("%s/FCgrid_new/%s/output_%d.txt",data_path,seed_alg,i));
    ifstream fgrid_MaCCQE_true(Form("%s/FCgrid_MaCCQE/%s/output_%d.txt",data_path,seed_true,i));
    ifstream fgrid_MaCCQE_alg(Form("%s/FCgrid_MaCCQE/%s/output_%d.txt",data_path,seed_alg,i));
    ifstream fgrid_BeRPA_A_true(Form("%s/FCgrid_BeRPA_A/%s/output_%d.txt",data_path,seed_true,i));
    ifstream fgrid_BeRPA_A_alg(Form("%s/FCgrid_BeRPA_A/%s/output_%d.txt",data_path,seed_alg,i));
    ifstream fgrid_BeRPA_B_true(Form("%s/FCgrid_BeRPA_B/%s/output_%d.txt",data_path,seed_true,i));
    ifstream fgrid_BeRPA_B_alg(Form("%s/FCgrid_BeRPA_B/%s/output_%d.txt",data_path,seed_alg,i));
    ifstream fgrid_BeRPA_D_true(Form("%s/FCgrid_BeRPA_D/%s/output_%d.txt",data_path,seed_true,i));
    ifstream fgrid_BeRPA_D_alg(Form("%s/FCgrid_BeRPA_D/%s/output_%d.txt",data_path,seed_alg,i));
    ifstream fgrid_MnGaus_true(Form("%s/FCgrid_MnGaus/%s/output_%d.txt",data_path,seed_true,i));
    ifstream fgrid_MnGaus_alg(Form("%s/FCgrid_MnGaus/%s/output_%d.txt",data_path,seed_alg,i));
*/
    if(!fstat_true || !fstat_alg) continue;
    if(!fflux_true || !fflux_alg) continue;
    j++;
/*
    if(!fgrid_true || !fgrid_alg) continue;
    if(!fgrid_MaCCQE_true || !fgrid_MaCCQE_alg) continue;
    if(!fgrid_BeRPA_A_true || !fgrid_BeRPA_A_alg) continue;
    if(!fgrid_BeRPA_B_true || !fgrid_BeRPA_B_alg) continue;
    if(!fgrid_BeRPA_D_true || !fgrid_BeRPA_D_alg) continue;
    if(!fgrid_MnGaus_true || !fgrid_MnGaus_alg) continue;
*/
    

    fstat_true >> Uee2_t[0] >> Umm2_t[0] >> dm2_t[0] >> chi2_t[0] >> nochi2_t[0];
    fstat_alg >> Uee2[0] >> Umm2[0] >> dm2[0] >> chi2[0] >> nochi2[0];
    
    fflux_true >> Uee2_t[1] >> Umm2_t[1] >> dm2_t[1] >> chi2_t[1] >> nochi2_t[1];
    fflux_alg >> Uee2[1] >> Umm2[1] >> dm2[1] >> chi2[1] >> nochi2[1];
/*    
    fgrid_true >> Uee2_t[2] >> Umm2_t[2] >> dm2_t[2] >> chi2_t[2] >> nochi2_t[2];
    fgrid_alg >> Uee2[2] >> Umm2[2] >> dm2[2] >> chi2[2] >> nochi2[2];
    
    fgrid_MaCCQE_true >> Uee2_t[3] >> Umm2_t[3] >> dm2_t[3] >> chi2_t[3] >> nochi2_t[3];
    fgrid_MaCCQE_alg >> Uee2[3] >> Umm2[3] >> dm2[3] >> chi2[3] >> nochi2[3];
    
    fgrid_BeRPA_A_true >> Uee2_t[4] >> Umm2_t[4] >> dm2_t[4] >> chi2_t[4] >> nochi2_t[4];
    fgrid_BeRPA_A_alg >> Uee2[4] >> Umm2[4] >> dm2[4] >> chi2[4] >> nochi2[4];
    
    fgrid_BeRPA_B_true >> Uee2_t[5] >> Umm2_t[5] >> dm2_t[5] >> chi2_t[5] >> nochi2_t[5];
    fgrid_BeRPA_B_alg >> Uee2[5] >> Umm2[5] >> dm2[5] >> chi2[5] >> nochi2[5];
    
    fgrid_BeRPA_D_true >> Uee2_t[6] >> Umm2_t[6] >> dm2_t[6] >> chi2_t[6] >> nochi2_t[6];
    fgrid_BeRPA_D_alg >> Uee2[6] >> Umm2[6] >> dm2[6] >> chi2[6] >> nochi2[6];
    
    fgrid_MnGaus_true >> Uee2_t[7] >> Umm2_t[7] >> dm2_t[7] >> chi2_t[7] >> nochi2_t[7];
    fgrid_MnGaus_alg >> Uee2[7] >> Umm2[7] >> dm2[7] >> chi2[7] >> nochi2[7];
*/    
    for(int ii=0; ii<1; ii++){
      dchi2[ii] = nochi2[ii] - chi2[ii];
      dchi2_t[ii] = nochi2_t[ii] - chi2_t[ii];

      h[ii]->Fill(chi2[ii]);
      h_t[ii]->Fill(chi2_t[ii]);
      
      hd[ii]->Fill(dchi2[ii]);
      hd_t[ii]->Fill(dchi2_t[ii]);

      hchi2_cor[ii]->Fill(dchi2[ii], dchi2_t[ii]);

    }
    fstat_true.close();
    fstat_alg.close();
    fflux_true.close();
    fflux_alg.close();
/*
    fgrid_true.close();
    fgrid_alg.close();
    fgrid_MaCCQE_true.close();
    fgrid_MaCCQE_alg.close();
    fgrid_BeRPA_A_true.close();
    fgrid_BeRPA_A_alg.close();
    fgrid_BeRPA_B_true.close();
    fgrid_BeRPA_B_alg.close();
    fgrid_BeRPA_D_true.close();
    fgrid_BeRPA_D_alg.close();
    fgrid_MnGaus_true.close();
    fgrid_MnGaus_alg.close();
*/
  }

  
  for (int i = 0; i < 1; ++i) {
    h[i]->SetLineColor(colors[i]);
    h_t[i]->SetLineColor(colors[i]);
    hd[i]->SetLineColor(colors[i]);
    hd_t[i]->SetLineColor(colors[i]);
    
    h[i]->SetLineWidth(2.);
    h_t[i]->SetLineWidth(2.);
    hd[i]->SetLineWidth(2.);
    hd_t[i]->SetLineWidth(2.);

    h[i]->SetStats(0);
    h_t[i]->SetStats(0);
    hd[i]->SetStats(0);
    hd_t[i]->SetStats(0);

    hchi2_cor[i]->SetStats(0);
    hchi2_cor[i]->GetXaxis()->SetTitle("seeding_algorithm #Delta#chi^{2}");
    hchi2_cor[i]->GetYaxis()->SetTitle("seeding_true #Delta#chi^{2}");
  }

  hchi2_cor[0]->SetTitle("seeding_true vs seeding_algorithm #Delta#chi^{2} (stats only)");
  hchi2_cor[1]->SetTitle("seeding_true vs seeding_algorithm #Delta#chi^{2} (stats + flux)");
/*
  hchi2_cor[2]->SetTitle("seeding_true vs seeding_algorithm #Delta#chi^{2} (stats + flux + sig)");
  hchi2_cor[3]->SetTitle("seeding_true vs seeding_algorithm #Delta#chi^{2} (stats + flux + MaCCQE)");
  hchi2_cor[4]->SetTitle("seeding_true vs seeding_algorithm #Delta#chi^{2} (stats + flux + BeRPA_A)");
  hchi2_cor[5]->SetTitle("seeding_true vs seeding_algorithm #Delta#chi^{2} (stats + flux + BeRPA_B)");
  hchi2_cor[6]->SetTitle("seeding_true vs seeding_algorithm #Delta#chi^{2} (stats + flux + BeRPA_D)");
  hchi2_cor[7]->SetTitle("seeding_true vs seeding_algorithm #Delta#chi^{2} (stats + flux + MnGaus)");
*/
  h[0]->SetTitle("Best-fit #chi^{2} distribution (seeding_alg)");
  h[0]->GetXaxis()->SetTitle("#chi^{2}");
  h[0]->GetYaxis()->SetTitle("Entries");
  
  h_t[0]->SetTitle("Best-fit #chi^{2} distribution (seeding_true)");
  h_t[0]->GetXaxis()->SetTitle("#chi^{2}");
  h_t[0]->GetYaxis()->SetTitle("Entries");
  
  hd[0]->SetTitle("Best-fit #Delta#chi^{2} distribution (seeding_alg)");
  hd[0]->GetXaxis()->SetTitle("#Delta#chi^{2}");
  hd[0]->GetYaxis()->SetTitle("Entries");
  
  hd_t[0]->SetTitle("Best-fit #Delta#chi^{2} distribution (seeding_true)");
  hd_t[0]->GetXaxis()->SetTitle("#Delta#chi^{2}");
  hd_t[0]->GetYaxis()->SetTitle("Entries");

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogy();
  TLegend *lg = new TLegend(0.55,0.60,0.9,0.9);
  h[0]->Draw();
  lg->AddEntry(h[7],names[7].c_str(), "f");
  for (int i = 6; i >= 0; i--) {
  h[i]->Draw("same");
  lg->AddEntry(h[i],names[i].c_str(), "f");}
  lg->Draw();
  c->SaveAs(Form("%s/chi2.png",out_path));
  
  TCanvas *c_t = new TCanvas("c_t","",800,600);
  c_t->SetGrid();
  c_t->SetLogy();
  TLegend *lg_t = new TLegend(0.55,0.60,0.9,0.9);
  h_t[0]->Draw();
  lg_t->AddEntry(h_t[7],names[7].c_str(), "f");
  for (int i = 6; i >= 0; i--) {
  h_t[i]->Draw("same");
  lg_t->AddEntry(h_t[i],names[i].c_str(), "f");}
  lg_t->Draw();
  c_t->SaveAs(Form("%s/chi2_t.png",out_path));
  
  TCanvas *cd = new TCanvas("cd","",800,600);
  cd->SetGrid();
  cd->SetLogy();
  TLegend *lgd = new TLegend(0.55,0.60,0.9,0.9);
  hd[0]->Draw();
  lgd->AddEntry(hd[7],names[7].c_str(), "f");
  for (int i = 6; i >= 0; i--) {
  hd[i]->Draw("same");
  lgd->AddEntry(hd[i],names[i].c_str(), "f");}
  lgd->Draw();
  cd->SaveAs(Form("%s/dchi2.png",out_path));
  
  TCanvas *cd_t = new TCanvas("cd_t","",800,600);
  cd_t->SetGrid();
  cd_t->SetLogy();
  TLegend *lgd_t = new TLegend(0.55,0.60,0.9,0.9);
  hd_t[0]->Draw();
  lgd_t->AddEntry(hd_t[7],names[7].c_str(), "f");
  for (int i = 6; i >= 0; i--) {
  hd_t[i]->Draw("same");
  lgd_t->AddEntry(hd_t[i],names[i].c_str(), "f");}
  lgd_t->Draw();
  cd_t->SaveAs(Form("%s/dchi2_t.png",out_path));


  gStyle->SetNumberContours(999);
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  for (int i = 0; i < 8; ++i) {
    TString canvasName = Form("canvas_%d", i);
    TCanvas *canvas = new TCanvas(canvasName, canvasName, 800, 800);
    canvas->SetGrid();
    canvas->SetLogx();
    canvas->SetLogy();
    canvas->SetLogz();
    hchi2_cor[i]->Draw("colz");
    canvas->SaveAs(Form("%s/hchi2_cor%d.png",out_path,i));
  }

}
