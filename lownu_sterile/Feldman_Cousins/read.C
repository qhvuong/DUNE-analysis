void read()
{
  TFile *f = new TFile("FC3_500000_2.root","READ");

  TMatrixD* mtr = dynamic_cast<TMatrixD*>(f->Get("scales"));

  int N = 500000;
  int nbins = 250;
  TMatrixD scales(N, nbins);

  TH2D* hscales = new TH2D("hscales","",N,0,N,nbins,0,nbins);

  for(int i=0; i<N; i++){
    for(int j=0; j<nbins; j++){
      hscales->SetBinContent(i+1,j+1, (*mtr)(i,j));
    }
  }

  TFile *fout = new TFile("FC3_2.root","RECREATE");
  hscales->Write();
  fout->Close(); 

}
