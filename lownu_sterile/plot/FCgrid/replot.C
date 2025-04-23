void replot()
{
  char seed_alg[]="seeding_algorithm";

  
  double Uee2, Umm2, dm2, chi2, nochi2, dchi2;
  int N=1000;


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

  TH1D *hd[8];

  for (int i = 0; i < 1; ++i) {
    hd[i]   = new TH1D(Form("hd%d", i), "", 100, 0, 40);}

  const char data_path[] = "/pnfs/dune/scratch/users/qvuong/output";
  const char out_path[] = "/exp/dune/app/users/qvuong/data/lownu/Feldman_Cousins";

  int j=0;

  for(int i=0; i<N; i++) {
    if(i%100==0) std::cout << i*100./N << " percent\n";
        
    ifstream fgrid_alg(Form("%s/FCgrid_new/%s/output_%d.txt",data_path,seed_alg,i));

    fgrid_alg >> Uee2 >> Umm2 >> dm2 >> chi2 >> nochi2;
    dchi2 = nochi2 - chi2;

      
    hd[0]->Fill(dchi2);
    fgrid_alg.close();

  }

  
  for (int i = 0; i < 1; ++i) {

    hd[i]->SetLineColor(kBlack);

    hd[i]->SetLineWidth(2.);

    hd[i]->SetStats(0);
  }

  
  hd[0]->SetTitle("Best-fit #Delta#chi^{2} distribution (seeding_alg)");
  hd[0]->GetXaxis()->SetTitle("#Delta#chi^{2}");
  hd[0]->GetYaxis()->SetTitle("Entries");

  
  TCanvas *cd = new TCanvas("cd","",800,600);
  cd->SetGrid();
  cd->SetLogy();
  hd[0]->Draw();
  cd->SaveAs(Form("%s/dchi2.png",out_path));
  


}
