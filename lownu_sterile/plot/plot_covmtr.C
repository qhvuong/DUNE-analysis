#include "DUNEStyle.h"


static const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu";

static const int nbins_CC = 56;
static const int nbins_nue = 8;
static const int nbins = nbins_CC + nbins_CC + nbins_nue;

void plot_covmtr()
{
    //TFile *CC_f  = new TFile(Form("%s/input_dfiles/CC_output_56bins.root",data_path),"READ");
    //TFile *nue_f = new TFile(Form("%s/input_dfiles/nue_output_8bins.root",data_path),"READ");
    TFile *f_sig = new TFile(Form("%s/uncertainties/xS_covmtr/xS_unc_5.root",data_path), "READ");
    TFile *f_fl  = new TFile(Form("%s/uncertainties/flux_covmtr/flux_covmtr_1104.root",data_path),"READ");
    TFile *f_det = new TFile(Form("%s/uncertainties/det_covmtr/det_covmtr.root",data_path), "READ");

    //TH1D *

    TH2D *fl_cov  = (TH2D*)f_fl->Get("hcv");
    TH2D *sig_cov = (TH2D*)f_sig->Get("hcv_tot");
    TH2D *det_cov = (TH2D*)f_det->Get("hcov4");

    fl_cov->SetStats(0);
    sig_cov->SetStats(0);
    det_cov->SetStats(0);

    fl_cov->SetTitle("Absolute Flux Uncertainty");
    sig_cov->SetTitle("Absolute Cross Section Uncertainty");
    det_cov->SetTitle("Absolute Detector Uncertainty");

    fl_cov->SetMaximum(1E11);
    sig_cov->SetMaximum(1E11);
    det_cov->SetMaximum(1E11);
    
    fl_cov->SetMinimum(1.);
    sig_cov->SetMinimum(1.);
    det_cov->SetMinimum(1.);

    //makes lines to divide matrix to show the three different samples
    double x = fl_cov->GetXaxis()->GetBinUpEdge(nbins_CC);
    double y1 = fl_cov->GetYaxis()->GetBinUpEdge(0);
    double y2 = fl_cov->GetYaxis()->GetBinUpEdge(nbins);

    double y = fl_cov->GetYaxis()->GetBinUpEdge(nbins_CC);
    double x1 = fl_cov->GetXaxis()->GetBinUpEdge(0);
    double x2 = fl_cov->GetXaxis()->GetBinUpEdge(nbins);

    double w = fl_cov->GetXaxis()->GetBinUpEdge(2*nbins_CC);
    double v1 = fl_cov->GetYaxis()->GetBinUpEdge(0);
    double v2 = fl_cov->GetYaxis()->GetBinUpEdge(nbins);

    double v = fl_cov->GetYaxis()->GetBinUpEdge(2*nbins_CC);
    double w1 = fl_cov->GetXaxis()->GetBinUpEdge(0);
    double w2 = fl_cov->GetXaxis()->GetBinUpEdge(nbins);

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

    //gStyle->SetPalette(kColorPrintableOnGrey); 
    TColor::InvertPalette(); 
    //gStyle->SetNumberContours(999);

    TCanvas *c = new TCanvas("c", "", 2100, 500);
    c->Divide(3,1);
    c->cd(1);
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    fl_cov->Draw("colz");
    l1->Draw("same");
    l2->Draw("same");
    l3->Draw("same");
    l4->Draw("same");
    c->cd(2);
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    sig_cov->Draw("colz");
    l1->Draw("same");
    l2->Draw("same");
    l3->Draw("same");
    l4->Draw("same");
    c->cd(3);
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    det_cov->Draw("colz");
    l1->Draw("same");
    l2->Draw("same");
    l3->Draw("same");
    l4->Draw("same");
    c->SaveAs(Form("cov_tot.png"));
}
