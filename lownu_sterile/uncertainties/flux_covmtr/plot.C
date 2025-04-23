#include "DUNEStyle.h"

static const int nbins_CC = 56;
static const int nbins_nue = 8;
static const int nbins = 2*nbins_CC + nbins_nue;

const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu/uncertainties/flux_covmtr";



void plotMtr(TH2D* cov, const char *name, TCanvas *c, TLine *l1, TLine *l2, TLine *l3, TLine *l4){
  cov->GetXaxis()->SetTitle("Bin number");
  cov->GetYaxis()->SetTitle("Bin number");
  cov->SetStats(0);

  c->Clear();
  cov->Draw("colz");
  l1->Draw("same");
  l2->Draw("same");
  l3->Draw("same");
  l4->Draw("same");
  c->SaveAs(Form("%s.png",name));
}


void plot()
{
  TFile *f = new TFile(Form("%s/flux_covmtr_TEST.root",data_path), "READ");

  const int N_nucut = 5;
  double nucut[N_nucut] = {10., 3., 1., 0.5, 0.3};

  TH2D *hcv[N_nucut];
  TH2D *hfrcv[N_nucut];
  TH2D *hcr[N_nucut]; TH2D *hcr_log[N_nucut];

  //gStyle->SetPalette(kColorPrintableOnGrey); 
  TColor::InvertPalette();
  TCanvas *c = new TCanvas("c","",800,600);

  for(int j=0; j<N_nucut; j++){
    hcv[j] = (TH2D*)f->Get(Form("hcv%d",j));

    hfrcv[j] = (TH2D*)f->Get(Form("hfrcv%d",j));

    hcr[j] = (TH2D*)f->Get(Form("hcr%d",j));
    hcr_log[j] = (TH2D*)f->Get(Form("hcr%d",j));

    hcv[j]->SetMinimum(1E-3);
    hcv[j]->SetMaximum(1E11);

    hfrcv[j]->SetMinimum(1E-5);
    hfrcv[j]->SetMaximum(1.);

    hcr[j]->SetMinimum(0.);
    hcr[j]->SetMaximum(1.);
    hcr_log[j]->SetMinimum(0.05);
    hcr_log[j]->SetMaximum(1.);
  }

  //makes lines to divide matrix to show the three different samples
  double x = hcv[0]->GetXaxis()->GetBinUpEdge(nbins_CC);
  double y1 = hcv[0]->GetYaxis()->GetBinLowEdge(1);
  double y2 = hcv[0]->GetYaxis()->GetBinUpEdge(nbins);

  double y = hcv[0]->GetYaxis()->GetBinUpEdge(nbins_CC);
  double x1 = hcv[0]->GetXaxis()->GetBinLowEdge(1);
  double x2 = hcv[0]->GetXaxis()->GetBinUpEdge(nbins);

  double w = hcv[0]->GetXaxis()->GetBinUpEdge(2*nbins_CC);
  double v1 = hcv[0]->GetYaxis()->GetBinLowEdge(1);
  double v2 = hcv[0]->GetYaxis()->GetBinUpEdge(nbins);

  double v = hcv[0]->GetYaxis()->GetBinUpEdge(2*nbins_CC);
  double w1 = hcv[0]->GetXaxis()->GetBinLowEdge(1);
  double w2 = hcv[0]->GetXaxis()->GetBinUpEdge(nbins);

  TLine *l1 = new TLine(x, y1, x, y2);  //splits muCC from the rest in x axis 
  TLine *l2 = new TLine(x1, y, x2, y);  //splits muCC from the rest in y axis
  TLine *l3 = new TLine(w, v1, w, v2);  //splits nu+e from the rest in x axis
  TLine *l4 = new TLine(w1, v, w2, v);  //splits nu+e from the rest in y axis

  l1->SetLineWidth(2);
  l1->SetLineColor(kBlack);
  l2->SetLineWidth(2);
  l2->SetLineColor(kBlack);
  l3->SetLineWidth(2);
  l3->SetLineColor(kBlack);
  l4->SetLineWidth(2);
  l4->SetLineColor(kBlack);


  // PLOT   
  
  gPad->SetLogz(0);
  gPad->Update();
  for(int j=0; j<N_nucut; j++){
    plotMtr(hcr[j], Form("%s/flux_cor%d",data_path,j), c, l1, l2, l3, l4);
  }
  

  gPad->SetLogz();
  gPad->Update();
  for(int j=0; j<N_nucut; j++){
    //plotMtr(hcv[j], Form("%s/flux_cov%d_log",data_path,j), c, l1, l2, l3, l4);
    plotMtr(hfrcv[j], Form("%s/flux_frcov%d_log",data_path,j), c, l1, l2, l3, l4);
    //plotMtr(hcr_log[j], Form("%s/flux_cor%d_log",data_path,j), c, l1, l2, l3, l4);
  }

  c->Close();
}