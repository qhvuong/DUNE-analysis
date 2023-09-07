void test()
{

   //TAxis *axis = h->GetXaxis();
   int bins = 100;

   Double_t from = -6;
   Double_t to = -0.7;
   Double_t width = (to - from) / bins;
   //Double_t *new_bins = new Axis_t[bins + 1];
   Double_t new_bins[100];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
   }
   //axis->Set(bins, new_bins);
   //delete new_bins;
   
   TH1D *h = new TH1D("h","",100,new_bins);

   
   for (int i = 0; i <= bins; i++) {
     //new_bins[i] = TMath::Power(10, from + i * width);
     h->Fill(new_bins[i]);
   }
   h->Draw();

}
