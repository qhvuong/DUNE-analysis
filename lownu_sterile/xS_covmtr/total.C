static const int nbins_CC = 58;
static const int nbins_nue = 8;
static const int nbins = 2*nbins_CC + nbins_nue;

void total()
{
  int cutNu = 3;

  TH2D *hcv  = new TH2D("hcv","",nbins,0,nbins,nbins,0,nbins);
  TH2D *frhcv  = new TH2D("frhcv","",nbins,0,nbins,nbins,0,nbins);

  std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"};
  const char *name[] = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"};

  int N_wgt = namelist.size();
  std::cout << N_wgt << "\n";

  char name1[100] = "wgt_Mnv2p2hGaussEnhancement";


  for(int iw=0; iw<N_wgt; iw++){
    if(iw%10==0) std::cout << iw*100.0/N_wgt << "\n";

    TFile *f = new TFile(Form("%s_covmtr%d.root",name[iw],cutNu),"READ");

    TH2D *cv2 = (TH2D*)f->Get("hcv2");    
    TH2D *cv  = (TH2D*)f->Get("hcv");    
    TH2D *frcv2 = (TH2D*)f->Get("hfrcv2");    
    TH2D *frcv  = (TH2D*)f->Get("hfrcv");    

    if(name[iw]==name1){
      hcv->Add(cv2);
      frhcv->Add(frcv2);}
    else{
      hcv->Add(cv);
      frhcv->Add(frcv);}
  }

  hcv->SetTitle("Total Covariance Matrix");
  frhcv->SetTitle("Fractional Covariance Matrix");

  hcv->SetStats(0);
  frhcv->SetStats(0);


  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  TCanvas *c = new TCanvas("c","",1600,600);
  c->Divide(2,1);
  c->cd(1);
  hcv->Draw("colz");
  c->cd(2);
  frhcv->Draw("colz");
  c->SaveAs(Form("sigma%d_5sig.png",cutNu));

  TCanvas *c1 = new TCanvas("c1","",1600,600);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogz();
  hcv->Draw("colz");
  c1->cd(2);
  frhcv->Draw("colz");
  c1->SaveAs(Form("sigma%d_5sig_logz.png",cutNu));


  TFile *out = new TFile(Form("total_sigmtr%d_5sig.root",cutNu),"RECREATE");
  hcv->Write();
  frhcv->Write();
  out->Close();

}
