void nue()
{
  TFile *f_nue = new TFile("/exp/dune/app/users/qvuong/data/lownu/input_dfiles/nue_output_test.root");

  TH1D *h = (TH1D*)f_nue->Get("hElep0");

  TCanvas *c = new TCanvas("c","",800,600);
  h->Draw();
  c->SaveAs("hElep_8bins.png");


}
