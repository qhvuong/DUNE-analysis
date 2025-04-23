void xS()
{
  TFile *fl = new TFile("flux_50bins.root","READ");
  TH1D *mfl = (TH1D*)fl->Get("h");
  TH1D *efl = (TH1D*)fl->Get("he");
  double totalPOT = dynamic_cast<TParameter<double>*>(fl->Get("total_pot"))->GetVal();
  double scalePOT = 1.1E21/totalPOT;
  std::cout << totalPOT << "\n";

  TFile *f = new TFile("../CC/eff.root","READ");

  TH1D *xSm[5], *xSe[5];

  TH1D *mEff[5], *eEff[5];

  double nA = 6.0221408e+23;
  double nn = 5E7*nA;
  double factor = 10*1.E42/nn / scalePOT;

  for(int i=0; i<5; i++){
    xSm[i] = (TH1D*)f->Get(Form("hm%d",i));
    xSe[i] = (TH1D*)f->Get(Form("he%d",i));

    mEff[i] = (TH1D*)f->Get(Form("hmEff%d",i));
    eEff[i] = (TH1D*)f->Get(Form("heEff%d",i));

    xSm[i]->Divide(mEff[i]);
    xSe[i]->Divide(eEff[i]);

    xSm[i]->Divide(mfl);
    xSe[i]->Divide(efl);

    xSm[i]->Scale(factor);
    xSe[i]->Scale(factor);
  }

  //xSm[4]->Draw("hist");

  TFile *out = new TFile("xS.root","RECREATE");
  for(int i=0; i<5; i++){
    xSm[i]->Write();
    xSe[i]->Write();
  }
  out->Close();

  /*
  TH1D *hEv0_nc = (TH1D*)f->Get("hEv0_nc");

  TH1D *hEv0 = (TH1D*)f->Get("hEv0");

  TH1D *hEv_nc_t0 = (TH1D*)f->Get("hEv_nc_t0");
  TH1D *hEv_nc_t1 = (TH1D*)f->Get("hEv_nc_t1");
  TH1D *hEv_nc_t2 = (TH1D*)f->Get("hEv_nc_t2");
  TH1D *hEv_nc_n0 = (TH1D*)f->Get("hEv_nc_n0");
  TH1D *hEv_nc_n1 = (TH1D*)f->Get("hEv_nc_n1");
  TH1D *hEv_nc_n2 = (TH1D*)f->Get("hEv_nc_n2");
  
  TH1D *hEv_t0 = (TH1D*)f->Get("hEv_t0");
  TH1D *hEv_t1 = (TH1D*)f->Get("hEv_t1");
  TH1D *hEv_t2 = (TH1D*)f->Get("hEv_t2");
  TH1D *hEv_n0 = (TH1D*)f->Get("hEv_n0");
  TH1D *hEv_n1 = (TH1D*)f->Get("hEv_n1");
  TH1D *hEv_n2 = (TH1D*)f->Get("hEv_n2");

  TH1D *hEff    = (TH1D*)hEv0->Clone();
  TH1D *hEff_t0 = (TH1D*)hEv_t0->Clone();
  TH1D *hEff_t1 = (TH1D*)hEv_t1->Clone();
  TH1D *hEff_t2 = (TH1D*)hEv_t2->Clone();
  TH1D *hEff_n0 = (TH1D*)hEv_n0->Clone();
  TH1D *hEff_n1 = (TH1D*)hEv_n1->Clone();
  TH1D *hEff_n2 = (TH1D*)hEv_n2->Clone();

  hEff->Divide(hEv0_nc);
  hEff_t0->Divide(hEv_nc_t0);
  hEff_t1->Divide(hEv_nc_t1);
  hEff_t2->Divide(hEv_nc_t2);
  hEff_n0->Divide(hEv_nc_n0);
  hEff_n1->Divide(hEv_nc_n1);
  hEff_n2->Divide(hEv_nc_n2);

  TH1D *xs = (TH1D*)hEv0->Clone();
  TH1D *xs_t0 = (TH1D*)hEv_t0->Clone();
  TH1D *xs_t1 = (TH1D*)hEv_t1->Clone();
  TH1D *xs_t2 = (TH1D*)hEv_t2->Clone();
  TH1D *xs_n0 = (TH1D*)hEv_n0->Clone();
  TH1D *xs_n1 = (TH1D*)hEv_n1->Clone();
  TH1D *xs_n2 = (TH1D*)hEv_n2->Clone();

  xs->Divide(hEff);
  xs_t0->Divide(hEff_t0);
  xs_t1->Divide(hEff_t1);
  xs_t2->Divide(hEff_t2);
  xs_n0->Divide(hEff_n0);
  xs_n1->Divide(hEff_n1);
  xs_n2->Divide(hEff_n2);

  xs->Divide(hfl);
  xs_t0->Divide(hfl);
  xs_t1->Divide(hfl);
  xs_t2->Divide(hfl);
  xs_n0->Divide(hfl);
  xs_n1->Divide(hfl);
  xs_n2->Divide(hfl);

  
  xs->Scale(factor);
  xs_t0->Scale(factor);
  xs_t1->Scale(factor);
  xs_t2->Scale(factor);
  xs_n0->Scale(factor);
  xs_n1->Scale(factor);
  xs_n2->Scale(factor);

  xs->SetStats(0);
  xs_t0->SetStats(0);
  xs_t1->SetStats(0);
  xs_t2->SetStats(0);
  xs_n0->SetStats(0);
  xs_n1->SetStats(0);
  xs_n2->SetStats(0);

  xs->SetMaximum(3.5);

  xs->SetLineColor(1);
  xs_t0->SetLineColor(2);
  xs_t1->SetLineColor(3);
  xs_t2->SetLineColor(4);

  xs->SetLineWidth(2);
  xs_t0->SetLineWidth(2);
  xs_t1->SetLineWidth(2);
  xs_t2->SetLineWidth(2);

  xs_n0->SetMarkerStyle(21);
  xs_n1->SetMarkerStyle(21);
  xs_n2->SetMarkerStyle(21);
  xs_n0->SetMarkerColor(7);
  xs_n1->SetMarkerColor(6);
  xs_n2->SetMarkerColor(5);

  xs->GetXaxis()->SetTitle("Ev (GeV)");
  xs->GetYaxis()->SetTitle("#sigma (#times 10^{-38} cm^2/nucleon)");


  TCanvas *cxs = new TCanvas("cxs","",800,600);
  gPad->SetGrid();
  xs->Draw("hist");
  xs_t0->Draw("hist same");
  xs_t1->Draw("hist same");
  xs_t2->Draw("hist same");
  xs_n0->Draw("hist same p");
  xs_n1->Draw("hist same p");
  xs_n2->Draw("hist same p");
  TLegend *lxs = new TLegend(0.10,0.65,0.4,0.9);
  lxs->AddEntry(xs,"no nu cut");
  lxs->AddEntry(xs_t0,"true_nu < 2.0GeV");
  lxs->AddEntry(xs_t1,"true_nu < 0.8GeV");
  lxs->AddEntry(xs_t2,"true_nu < 0.3GeV");
  lxs->AddEntry(xs_n0,"reco_nu < 2.0GeV");
  lxs->AddEntry(xs_n1,"reco_nu < 0.8GeV");
  lxs->AddEntry(xs_n2,"reco_nu < 0.3GeV");
  lxs->Draw();
  cxs->SaveAs("cutSigma_50bins.png");
  */
  
}
