void test()
{
  TFile *f = new TFile("LE.root","READ");
  TH2D *h = (TH2D*)f->Get("h");
  TH1D *hL[480];

  for(int i = 0; i < 10; i++) {
    hL[i] = (TH1D*)h->ProjectionY("",i+1,i+1);
    if(i<40) hL[i]->Fit("expo");
  }

}
