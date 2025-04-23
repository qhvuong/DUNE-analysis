#include <iostream>
#include <string>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TTree.h>
#include <TRandom3.h>
#include <THStack.h>
#include <TChain.h>
#include <TLegend.h>
#include <TParameter.h>
#include <cmath>

int mu_eCC_nue_covars()
{
    TFile *file = TFile::Open("./CCnue_output_1.root", "update");  //this is the file made in sample_prep.cpp

//get the histograms (nominal, pos and neg shifts for mu and e) -------------------------------------------------------
    TH1D *hnom = nullptr;
    TH1D *hpmu = nullptr;
    TH1D *hnmu = nullptr;
    TH1D *hpel = nullptr;
    TH1D *hnel = nullptr;

    file->GetObject("hnom", hnom);
    file->GetObject("hpmu", hpmu);
    file->GetObject("hnmu", hnmu);
    file->GetObject("hpel", hpel);
    file->GetObject("hnel", hnel);

//useful constants ----------------------------------------------------------------------------------------------------
    const int nbins = hnom->GetNbinsX();  //or 2*nbinsCC + nbinsY or 120
    const int nbinsCC = 56;
    const int nbinsY = 8;

//average covariance histograms ---------------------------------------------------------------------------------------
    TH2D* hmu_cov = new TH2D("hmu_cov", "2% Mu Reco Energy Shift", nbins, 0, nbins, nbins, 0, nbins);
    TH2D* hel_cov = new TH2D("hel_cov", "2.5% electron Energy Shift", nbins, 0, nbins, nbins, 0, nbins);

//function computes average fractional covariance ---------------------------------------------------------------------
    auto covariance = [hnom, nbins](TH1D* hp, TH1D* hn, TH2D* covhist){
        for( int ii = 1; ii < nbins + 1; ++ii){
            for( int jj = 1; jj < nbins + 1; ++jj){
                //cov_ij = delta_i * delta_j
                //delta_i = s_i - n_i (difference in number of entries in bin i for nominal and shift)

                double n_i = hnom->GetBinContent(ii);
                double ps_i = hp->GetBinContent(ii);  //positive shift histogram (hpmu or hpel)
                double ns_i = hn->GetBinContent(ii);  //negative shift histogram (hnmu or hnel)

                double n_j = hnom->GetBinContent(jj);
                double ps_j = hp->GetBinContent(jj);
                double ns_j = hn->GetBinContent(jj);

                double prod = n_i*n_j;

                //C++ doesn't flag dividing by zero so need to check for it
                if (prod == 0){
                    covhist->SetBinContent(ii, jj, 0.);         
                }
                else{ 
                    double pdelta_i = ps_i - n_i;
                    double pdelta_j = ps_j - n_j;
                    double ndelta_i = ns_i - n_i;
                    double ndelta_j = ns_j - n_j;
              
                    double pfcov_ij = pdelta_i*pdelta_j / prod;  //positive shift covariance
                    double nfcov_ij = ndelta_i*ndelta_j / prod;  //negative shift covariance
                    double afcov_ij = (pfcov_ij+nfcov_ij) / 2;   //average covariance which ends up in the histograms

                    covhist->SetBinContent(ii, jj, afcov_ij);  //hmu_cov or hel_cov
                }
            }
        }
    };
//apply the function---------------------------------------------------------------------------------------------------
    covariance(hpel, hnel, hel_cov);
    covariance(hpmu, hnmu, hmu_cov);

//---------------------------------------------------------------------------------------------------------------------

//running this part of the script will cause the script to never finish 
//If you run this but put the output into a text file and not a histogram then it computes instantly
/*
//total fractional covariance -----------------------------------------------------------------------------------------

    TH2D* htot_cov = new TH2D("htot_cov", "muCC, eCC, nu+e e Covariances", nbins, 0, nbins, nbins, 0, nbins);

    //this SHOULD add the covariances for mu and electron in quadrature
    for(int ii = 1; ii < nbins + 1; ++ii){
        for(int jj = 1; jj < nbins + 1; ++ii){
            double mu = hmu_cov->GetBinContent(ii, jj);
            double el = hel_cov->GetBinContent(ii, jj);

            double quad = sqrt(mu*mu + el*el);

            htot_cov->SetBinContent(ii, jj, quad);
        }
    } 
//--------------------------------------------------------------------------------------------------------------------- 
*/

//plot ----------------------------------------------------------------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "Covars Shift", 800, 800);

    //makes the red to white to blue gradient 
    Double_t Red[3] = { 0.0, 1.0, 1.0 }; 
    Double_t Green[3] = { 0.0, 1.0, 0.0 }; 
    Double_t Blue[3]  = { 1.0, 1.0, 0.0 }; 
    Double_t Stops[3] = { 0.0, 0.5, 1.0 }; 
    TColor::CreateGradientColorTable(3, Stops, Red, Green, Blue, 999);
    gStyle->SetNumberContours(999);
    gPad->SetRightMargin(0.15);

    //makes lines to divide matrix to show the three different samples
    double x = hmu_cov->GetXaxis()->GetBinCenter(nbinsCC);
    double y1 = hmu_cov->GetYaxis()->GetBinLowEdge(0);
    double y2 = hmu_cov->GetYaxis()->GetBinUpEdge(nbins);

    double y = hmu_cov->GetYaxis()->GetBinCenter(nbinsCC);
    double x1 = hmu_cov->GetXaxis()->GetBinLowEdge(0);
    double x2 = hmu_cov->GetXaxis()->GetBinUpEdge(nbins);

    double w = hmu_cov->GetXaxis()->GetBinCenter(2*nbinsCC);
    double v1 = hmu_cov->GetYaxis()->GetBinLowEdge(0);
    double v2 = hmu_cov->GetYaxis()->GetBinUpEdge(nbins);

    double v = hmu_cov->GetYaxis()->GetBinCenter(2*nbinsCC);
    double w1 = hmu_cov->GetXaxis()->GetBinLowEdge(0);
    double w2 = hmu_cov->GetXaxis()->GetBinUpEdge(nbins);

    TLine *l1 = new TLine(x, y1, x, y2);  //splits muCC from the rest in x axis 
    TLine *l2 = new TLine(x1, y, x2, y);  //splits muCC from the rest in y axis
    TLine *l3 = new TLine(w, v1, w, v2);  //splits nu+e from the rest in x axis
    TLine *l4 = new TLine(w1, v, w2, v);  //splits nu+e from the rest in y axis

    hmu_cov->SetStats(0);
    hmu_cov->SetTitle("2% Mu Energy Shift");
    hmu_cov->Draw("colz");
    hmu_cov->GetZaxis()->SetRangeUser(-0.01, 0.01);  //must have this to make zero look white in the plot
    l1->SetLineWidth(2);
    l1->SetLineColor(kBlack);
    l1->Draw("same");
    l2->SetLineWidth(2);
    l2->SetLineColor(kBlack);
    l2->Draw("same");
    l3->SetLineWidth(2);
    l3->SetLineColor(kBlack);
    l3->Draw("same");
    l4->SetLineWidth(2);
    l4->SetLineColor(kBlack);
    l4->Draw("same");
    c1->SaveAs("./hmu_cov.png");

    c1->Clear();

    hel_cov->SetStats(0);
    hel_cov->SetTitle("2.5% electron Energy Shift");
    hel_cov->Draw("colz");
    hel_cov->GetZaxis()->SetRangeUser(-0.01, 0.01);
    l1->SetLineWidth(2);
    l1->SetLineColor(kBlack);
    l1->Draw("same");
    l2->SetLineWidth(2);
    l2->SetLineColor(kBlack);
    l2->Draw("same");
    l3->SetLineWidth(2);
    l3->SetLineColor(kBlack);
    l3->Draw("same");
    l4->SetLineWidth(2);
    l4->SetLineColor(kBlack);
    l4->Draw("same");
    c1->SaveAs("./hel_cov.png"); 

/* Doesn't work right now since total covariance doesn't work    
    htot_cov->SetStats(0);
    htot_cov->SetTitle("2% muCC and 2.5% e shifts");
    htot_cov->Draw("colz");
    htot_cov->GetZaxis()->SetRangeUser(-0.01, 0.01);  //this step makes sure that 0 is white
    l1->SetLineWidth(2);
    l1->SetLineColor(kBlack);
    l1->Draw("same");
    l2->SetLineWidth(2);
    l2->SetLineColor(kBlack);
    l2->Draw("same");
    l3->SetLineWidth(2);
    l3->SetLineColor(kBlack);
    l3->Draw("same");
    l4->SetLineWidth(2);
    l4->SetLineColor(kBlack);
    l4->Draw("same");
    c1->SaveAs("./htot_cov.png");
*/

//---------------------------------------------------------------------------------------------------------------------

//save to ROOT file ---------------------------------------------------------------------------------------------------
    hmu_cov->Write();
    hel_cov->Write();
    //htot_cov->Write();
//---------------------------------------------------------------------------------------------------------------------
    file->Close();

    return(0);
}