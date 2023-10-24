void plot()
{
  int j = 1;

  double Uee2, Umm2, dm2, sen;
/*
  double min=log(0.0001);
  double maxx=log(0.2);
  double maxy=log(0.05);
*/
  //double min=-4;
  //double maxx=-0.7;
  //double maxy=-1.3;
  int N=100;
/*
  double widthx = (maxx-min)/N;
  double widthy = (maxy-min)/N;
  double binx[N],biny[N];
  for(int i=0; i<N; i++){
    binx[i] = TMath::Power(10, min + i * widthx); 
    biny[i] = TMath::Power(10, min + i * widthy); 

    std::cout << binx[i] << "\t" << biny[i] << "\n";
  }
*/
  TH2D *h = new TH2D("h","",N,0.0001,0.1,N,0.0001,0.1);

  for(int i=0; i<N*N; i++) {
    ifstream f(Form("/pnfs/dune/scratch/users/qvuong/output/noNuCut/100dm%d/output_%d.txt",j,i));

    if(i%100==0) std::cout << i*1./100 << " percent" << "\n";

    //if (f.good())
    if (f) {
      f >> Uee2 >> Umm2 >> dm2 >> sen;
      //if(sen<0) std::cout << i << "\n";

      h->Fill(Uee2,Umm2,sen);    
    }
    
    else std::cout << i << "\n"; 

    f.close();
  }

  //h->Draw();

  h->SetStats(0);
  h->SetMaximum(20);
  //h->SetMinimum(0.1);
  h->SetTitle(Form("Sensitivity Contour (dm2=%.2f)",dm2));
  h->GetXaxis()->SetTitle("U_{e4}^{2}");
  h->GetYaxis()->SetTitle("U_{m4}^{2}");

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  TCanvas *cL = new TCanvas("cL","",800,800);
  cL->SetLogx();
  cL->SetLogy();
  //cL->SetLogz();
  h->Draw("colz");
  //cL->SaveAs(Form("dm2_%d_contour_Log.png",j));
/*
  TCanvas *c = new TCanvas("c","",800,800);
  h->Draw("colz");
  c->SaveAs(Form("dm2_%d_contour.png",j));
*/    
  
  TFile *out = new TFile(Form("noNuCut_dm%d.root",j),"RECREATE");
  h->Write();
  out->Close();
}
