void test()
{
  TH2D *h = new TH2D("h","",10,0,10,10,0,10);
  h->SetTitle("#Delta#chi^{2}");
  h->GetXaxis()->SetTitle("U_{#mu4}^{2}"); 
  h->GetYaxis()->SetTitle("#Deltam^{2}"); 
  h->Draw();
}
