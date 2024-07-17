void plot()
{
/*
  double ue42, um42, ut42, dm2, bf_ue42, bf_um42, bf_ut42, bf_dm2, chi2, nochi2, dchi2;
  const char out_path[] = "/pnfs/dune/scratch/users/qvuong/output/contour_4pars";
  int dm[4] = {1, 5, 10, 100};
  int N_tot = 10000;

  double U_max = 0.7;
  double U_min = 0.0001;

  int N = 100;
  double log_min = log10(U_min);
  double log_max = log10(U_max);
  double binWidth = (log_max - log_min) / N;

  double bin_edges[N+1], bin_center[N];
  for(int i=0; i<N+1; i++){
    bin_edges[i] = pow(10, (log_min + i * binWidth));
    //bin_edges[i+1] = U_min + pow(10, (log_min + (i+1) * binWidth));
    std::cout << bin_edges[i] << "\t";
  }

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  for(int j=1; j<4; j++){
  std::cout << "file = " << dm[j] << "\n";
  TH2D *h = new TH2D("h","",N,bin_edges,N,bin_edges);
  h->SetStats(0);
  h->SetTitle("Sensitivity Contour");

  for(int i=0; i<N_tot; i++){
    ifstream f(Form("%s/dm%d_true/output_%d.txt",out_path,dm[j],i));
    if(i%50==0) std::cout << i*100./N_tot << " percent...\n";
    if(!f) {
      std::cout << "failed: " << i << "\n";
      continue;
    }
    else{
      f >> ue42 >> um42 >> ut42 >> dm2 >> bf_ue42 >> bf_um42 >> bf_ut42 >> bf_dm2 >> chi2 >> nochi2;
      dchi2 = sqrt(nochi2 - chi2);
      h->Fill(um42, ue42, dchi2);
      f.close();
    }
  }

  TFile *out = new TFile(Form("contour_dm%d.root", dm[j]), "RECREATE");
  h->Write();
  out->Close();

  }
*/

  
  TFile *f0 = new TFile(Form("contour_dm1.root"),"READ");
  TFile *f1 = new TFile(Form("contour_dm5.root"),"READ");
  TFile *f2 = new TFile(Form("contour_dm10.root"),"READ");
  TFile *f3 = new TFile(Form("contour_dm100.root"),"READ");

  TH2D *h0 = (TH2D*)f0->Get("h");
  TH2D *h1 = (TH2D*)f1->Get("h");
  TH2D *h2 = (TH2D*)f2->Get("h");
  TH2D *h3 = (TH2D*)f3->Get("h");

  double contours[1];
  contours[0] = 3; 

  h0->SetContour(1, contours);
  h1->SetContour(1, contours);
  h2->SetContour(1, contours);
  h3->SetContour(1, contours);

  h0->SetStats(0);
  h1->SetStats(0);
  h2->SetStats(0);
  h3->SetStats(0);

  h0->GetXaxis()->SetRangeUser(1e-4,0.7);
  h0->GetYaxis()->SetRangeUser(1e-4,0.7);
  h1->GetXaxis()->SetRangeUser(1e-4,0.7);
  h1->GetYaxis()->SetRangeUser(1e-4,0.7);
  h2->GetXaxis()->SetRangeUser(1e-4,0.7);
  h2->GetYaxis()->SetRangeUser(1e-4,0.7);
  h3->GetXaxis()->SetRangeUser(1e-4,0.7);
  h3->GetYaxis()->SetRangeUser(1e-4,0.7);
 
  h0->SetLineColor(kBlack);
  h0->SetLineStyle(2);
  h0->SetLineWidth(2);
  h1->SetLineColor(kRed);
  h1->SetLineStyle(1);
  h1->SetLineWidth(2);
  h2->SetLineColor(kBlue);
  h2->SetLineStyle(9);
  h2->SetLineWidth(2);
  h3->SetLineColor(kGreen);
  //h3->SetLineStyle(9);
  h3->SetLineWidth(2);

  h0->SetTitle("3#sigma Contours");
  h0->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  h0->GetYaxis()->SetTitle("U_{e4}^{2}");

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogx();
  c->SetLogy();
  c->SetLogz();
  h0->Draw("cont3 same");
  h1->Draw("cont3 same");
  h2->Draw("cont3 same");
  h3->Draw("cont3 same");
  TLegend *lg = new TLegend(0.65,0.75,0.9,0.9);
  lg->AddEntry(h0,"#Deltam^{2} = 1.0");
  lg->AddEntry(h1,"#Deltam^{2} = 5.0");
  lg->AddEntry(h2,"#Deltam^{2} = 10.0");
  lg->AddEntry(h3,"#Deltam^{2} = 100.0");
  lg->Draw();
  c->SaveAs(Form("3sigma.png"));

}
