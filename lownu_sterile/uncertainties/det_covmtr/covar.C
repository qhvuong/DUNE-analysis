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
#include "DUNEStyle.h"

static const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu";

static const int nbins_CC = 56;
static const int nbins_nue = 8;
static const int nbins = nbins_CC + nbins_CC + nbins_nue;


void plot1D_pn(TH1D* hp, TH1D *hn, TH1D *ho, const char* name){
    ho->SetLineColor(kBlack);
    hp->SetLineColor(kBlue);
    hn->SetLineColor(kRed);
    ho->SetLineWidth(2.);
    hp->SetLineWidth(2.);
    hn->SetLineWidth(2.);
    ho->SetStats(0);
    hp->SetStats(0);
    hn->SetStats(0);
    TCanvas *c = new TCanvas("c", "", 800, 600);
    c->SetLogy();
    c->SetGrid();
    hp->Draw();
    hn->Draw("same");
    ho->Draw("same");
    TLegend *lg = new TLegend(0.65,0.75,0.9,0.9);
    lg->AddEntry(ho, "nominal");
    lg->AddEntry(hp, "positive");
    lg->AddEntry(hn, "negative");
    lg->Draw("same");
    c->SaveAs(Form("%s.png",name));
}

void plot1D(TH1D* ho, TH1D *hp, const char* name){
    ho->SetLineColor(kBlack);
    hp->SetLineColor(kRed);
    ho->SetLineWidth(2.);
    hp->SetLineWidth(2.);
    ho->SetStats(0);
    hp->SetStats(0);
    TCanvas *c = new TCanvas("c", "", 800, 600);
    c->SetLogy();
    c->SetGrid();
    hp->Draw();
    ho->Draw("same");
    TLegend *lg = new TLegend(0.65,0.75,0.9,0.9);
    lg->AddEntry(ho, "nominal");
    lg->AddEntry(hp, "shifted");
    lg->Draw("same");
    c->SaveAs(Form("%s.png",name));
}



void covariance_pn(TH1D* hnom, TH1D* hp, TH1D *hn, TH2D* cov, TH2D* frCov, TCanvas *c, const char *name){
    int nbins = hnom->GetNbinsX();
    //std::cout << "nbins = " << nbins << "\n";
    frCov->GetXaxis()->SetTitle("Bin number");
    frCov->GetYaxis()->SetTitle("Bin number");

    for(int ii=0; ii<nbins; ii++){
        double n_i = hnom->GetBinContent(ii+1);
        double ps_i = hp->GetBinContent(ii+1);  //positive shift histogram (hpmu or hpel)
        double ns_i = hn->GetBinContent(ii+1);  //negative shift histogram (hnmu or hnel)

        double pdelta_i = ps_i - n_i;
        double ndelta_i = ns_i - n_i;

        for(int jj=0; jj<nbins; jj++){
            double n_j = hnom->GetBinContent(jj+1);
            double ps_j = hp->GetBinContent(jj+1);
            double ns_j = hn->GetBinContent(jj+1);

            double pdelta_j = ps_j - n_j;
            double ndelta_j = ns_j - n_j;

            double pfcov_ij = abs(pdelta_i) * abs(pdelta_j);  //positive shift covariance
            double nfcov_ij = abs(ndelta_i) * abs(ndelta_j);  //negative shift covariance
            double afcov_ij = (pfcov_ij+nfcov_ij) / 2;   //average covariance which ends up in the histograms
            cov->SetBinContent(ii+1, jj+1, afcov_ij);
            
            afcov_ij = afcov_ij/(n_i*n_j);
            frCov->SetBinContent(ii+1, jj+1, afcov_ij);
        }
    }
    //std::cout << cov->GetBinContent(1,1);

    c->Clear();
    frCov->SetMinimum(1E-5);
    frCov->SetMaximum(1.);
    frCov->SetStats(0);
    c->SetLogz();
    frCov->Draw("colz");
    c->SaveAs(Form("%s/uncertainties/det_covmtr/plots/fractional_%s_logz.png",data_path,name));
/*
    c->Clear();
    cov->SetMinimum(1E-3);
    cov->SetMaximum(1E11);
    cov->SetStats(0);
    c->SetLogz();
    cov->Draw("colz");
    //c->SaveAs(Form("%s/uncertainties/det_covmtr/plots/%s_logz.png",data_path,name));
*/
}



void covariance(TH1D* hnom, TH1D* h, TH2D* cov, TH2D* frCov, TCanvas *c, const char *name){
    int nbins = hnom->GetNbinsX();
    //std::cout << "nbins = " << nbins << "\n";
    frCov->GetXaxis()->SetTitle("Bin number");
    frCov->GetYaxis()->SetTitle("Bin number");

    for(int ii=0; ii<nbins; ii++){
        double n_i = hnom->GetBinContent(ii+1);
        double s_i = h->GetBinContent(ii+1);  // shift histogram (hpmu or hpel)
        double delta_i = s_i - n_i;

        for(int jj=0; jj<nbins; jj++){
            double n_j = hnom->GetBinContent(jj+1);
            double s_j = h->GetBinContent(jj+1);
            double delta_j = s_j - n_j;
            double cov_ij = abs(delta_i) * abs(delta_j);  //positive shift covariance
            //cov->SetBinContent(ii+1, jj+1, cov_ij);  //hmu_cov or hel_cov
            cov_ij = cov_ij/(n_i*n_j);
            frCov->SetBinContent(ii+1, jj+1, cov_ij);
        }
    }
    std::cout << frCov->GetBinContent(1,1);


    c->Clear();
    frCov->SetMinimum(1E-5);
    frCov->SetMaximum(1.);
    frCov->SetStats(0);
    c->SetLogz();
    frCov->Draw("colz");
    c->SaveAs(Form("%s/uncertainties/det_covmtr/plots/fractional_%s_logz.png",data_path,name));
/*
    c->Clear();
    cov->SetMinimum(1E-3);
    cov->SetMaximum(1E11);
    cov->SetStats(0);
    c->SetLogz();
    cov->Draw("colz");
    c->SaveAs(Form("%s/uncertainties/det_covmtr/plots/%s_logz.png",data_path,name));
*/
}


int covar()
{
    TFile *fCC = new TFile(Form("%s/uncertainties/det_covmtr/outputCC_detCov_0102.root",data_path), "READ");  //this is the file made in sample_prep.cpp
    //TFile *fCC = new TFile("output_CC.root","READ");
    TFile *fnue = new TFile(Form("%s/input_dfiles/nue_output.root",data_path), "READ");

    const int N_nucut = 5;
    double nucut[N_nucut] = {100., 10., 1., 0.5, 0.3};

    TH1D *hn = (TH1D*)fnue->Get("hElep");
    TH1D *hnp = (TH1D*)fnue->Get("hElep_p");
    TH1D *hnn = (TH1D*)fnue->Get("hElep_n");
    TH1D *hn_resp = (TH1D*)fnue->Get("hEres_p");
    TH1D *hn_resn = (TH1D*)fnue->Get("hEres_n");

    //gStyle->SetPalette(kColorPrintableOnGrey); 
    //TColor::InvertPalette();
    TCanvas *c = new TCanvas("c","",800,600);
    TFile *out = new TFile(Form("%s/uncertainties/det_covmtr/det_covmtr.root",data_path),"RECREATE");

    // Get the histograms (nominal, pos and neg shifts for mu and e) -------------------------------------------------------
    TH1D* hm[N_nucut]; 
    TH1D* hmp[N_nucut]; TH1D* hmn[N_nucut];
    TH1D* hm_resp[N_nucut]; TH1D* hm_resn[N_nucut];

    TH1D* hm_nup[N_nucut]; TH1D* hm_nun[N_nucut];
    TH1D* hm_eRecoNp[N_nucut]; TH1D* hm_eRecoNn[N_nucut];
    TH1D* hm_eRecoPp[N_nucut]; TH1D* hm_eRecoPn[N_nucut];
    TH1D* hm_eRecoPip[N_nucut]; TH1D* hm_eRecoPin[N_nucut];
    TH1D* hm_LeNp[N_nucut]; TH1D* hm_LeNn[N_nucut];
    TH1D* hm_LeN[N_nucut];


    TH1D* he[N_nucut]; 
    TH1D* hep[N_nucut]; TH1D* hen[N_nucut];
    TH1D* he_resp[N_nucut]; TH1D* he_resn[N_nucut];

    TH1D* he_nup[N_nucut]; TH1D* he_nun[N_nucut];
    TH1D* he_eRecoNp[N_nucut]; TH1D* he_eRecoNn[N_nucut];
    TH1D* he_eRecoPp[N_nucut]; TH1D* he_eRecoPn[N_nucut];
    TH1D* he_eRecoPip[N_nucut]; TH1D* he_eRecoPin[N_nucut];
    TH1D* he_LeNp[N_nucut]; TH1D* he_LeNn[N_nucut];
    TH1D* he_LeN[N_nucut];


    // Forming 1d histograms with all 3 samples -----------------------------------------------------------------------------------
    TH1D* hnom[N_nucut];
    TH1D* hmp_tot[N_nucut]; TH1D* hmn_tot[N_nucut];
    TH1D* hep_tot[N_nucut]; TH1D* hen_tot[N_nucut];
    TH1D* hresp[N_nucut]; TH1D* hresn[N_nucut];
    TH1D* hnup[N_nucut]; TH1D* hnun[N_nucut];
    TH1D* hnup_eRecoN[N_nucut]; TH1D* hnun_eRecoN[N_nucut];
    TH1D* hnup_eRecoP[N_nucut]; TH1D* hnun_eRecoP[N_nucut];
    TH1D* hnup_eRecoPi[N_nucut]; TH1D* hnun_eRecoPi[N_nucut];
    TH1D* hLeNp[N_nucut]; TH1D* hLeNn[N_nucut];
    TH1D* hLeN[N_nucut];

    // Covariance matrices -------------------------------------------------------------------------------------------------------------
    TH2D* hcov[N_nucut];
    TH2D* hcov_had[N_nucut]; TH2D* hcov_lep[N_nucut];

    TH2D* hfrcov[N_nucut];
    TH2D* hfrcov_had[N_nucut]; TH2D* hfrcov_lep[N_nucut];


    TH2D* hmu_cov[N_nucut]; TH2D* hmu_frcov[N_nucut];
    TH2D* hel_cov[N_nucut]; TH2D* hel_frcov[N_nucut];
    TH2D* hEres_cov[N_nucut]; TH2D* hEres_frcov[N_nucut];
    TH2D* hnu_cov[N_nucut]; TH2D* hnu_frcov[N_nucut];
    TH2D* heRecoN_cov[N_nucut]; TH2D* heRecoN_frcov[N_nucut];
    TH2D* heRecoP_cov[N_nucut]; TH2D* heRecoP_frcov[N_nucut];
    TH2D* heRecoPi_cov[N_nucut]; TH2D* heRecoPi_frcov[N_nucut];
    TH2D* hLeNp_cov[N_nucut]; TH2D* hLeNp_frcov[N_nucut];
    TH2D* hLeNn_cov[N_nucut]; TH2D* hLeNn_frcov[N_nucut];
    TH2D* hLeN_cov[N_nucut]; TH2D* hLeN_frcov[N_nucut];



    for(int j=0; j<N_nucut; j++){
        hm[j] = (TH1D*)fCC->Get(Form("hm%d",j));
        hmp[j] = (TH1D*)fCC->Get(Form("hmp%d",j));
        hmn[j] = (TH1D*)fCC->Get(Form("hmn%d",j));
        hm_resp[j] = (TH1D*)fCC->Get(Form("hm_resp%d",j));
        hm_resn[j] = (TH1D*)fCC->Get(Form("hm_resn%d",j));

        hm_nup[j] = (TH1D*)fCC->Get(Form("hm_nup%d",j));
        hm_nun[j] = (TH1D*)fCC->Get(Form("hm_nun%d",j));
        hm_eRecoNp[j] = (TH1D*)fCC->Get(Form("hm_eRecoNp%d",j));
        hm_eRecoNn[j] = (TH1D*)fCC->Get(Form("hm_eRecoNn%d",j));
        hm_eRecoPp[j] = (TH1D*)fCC->Get(Form("hm_eRecoPp%d",j));
        hm_eRecoPn[j] = (TH1D*)fCC->Get(Form("hm_eRecoPn%d",j));
        hm_eRecoPip[j] = (TH1D*)fCC->Get(Form("hm_eRecoPip%d",j));
        hm_eRecoPin[j] = (TH1D*)fCC->Get(Form("hm_eRecoPin%d",j));
        hm_LeNp[j] = (TH1D*)fCC->Get(Form("hm_LeNp%d",j));
        hm_LeNn[j] = (TH1D*)fCC->Get(Form("hm_LeNn%d",j));
        hm_LeN[j] = (TH1D*)fCC->Get(Form("hm_LeN%d",j));


        he[j] = (TH1D*)fCC->Get(Form("he%d",j));
        hep[j] = (TH1D*)fCC->Get(Form("hep%d",j));
        hen[j] = (TH1D*)fCC->Get(Form("hen%d",j));
        he_resp[j] = (TH1D*)fCC->Get(Form("he_resp%d",j));
        he_resn[j] = (TH1D*)fCC->Get(Form("he_resn%d",j));

        he_nup[j] = (TH1D*)fCC->Get(Form("he_nup%d",j));
        he_nun[j] = (TH1D*)fCC->Get(Form("he_nun%d",j));
        he_eRecoNp[j] = (TH1D*)fCC->Get(Form("he_eRecoNp%d",j));
        he_eRecoNn[j] = (TH1D*)fCC->Get(Form("he_eRecoNn%d",j));
        he_eRecoPp[j] = (TH1D*)fCC->Get(Form("he_eRecoPp%d",j));
        he_eRecoPn[j] = (TH1D*)fCC->Get(Form("he_eRecoPn%d",j));
        he_eRecoPip[j] = (TH1D*)fCC->Get(Form("he_eRecoPip%d",j));
        he_eRecoPin[j] = (TH1D*)fCC->Get(Form("he_eRecoPin%d",j));
        he_LeNp[j] = (TH1D*)fCC->Get(Form("he_LeNp%d",j));
        he_LeNn[j] = (TH1D*)fCC->Get(Form("he_LeNn%d",j));
        he_LeN[j] = (TH1D*)fCC->Get(Form("he_LeN%d",j));


        
        hnom[j] = new TH1D(Form("hnom%d",j), "", nbins, 0, nbins);

        hmp_tot[j] = new TH1D(Form("hmp_tot%d",j), "", nbins, 0, nbins);
        hmn_tot[j] = new TH1D(Form("hmn_tot%d",j), "", nbins, 0, nbins);

        hep_tot[j] = new TH1D(Form("hep_tot%d",j), "", nbins, 0, nbins);
        hen_tot[j] = new TH1D(Form("hen_tot%d",j), "", nbins, 0, nbins);

        hresp[j] = new TH1D(Form("hresp%d",j), "", nbins, 0, nbins);
        hresn[j] = new TH1D(Form("hresn%d",j), "", nbins, 0, nbins);

        hnup[j] = new TH1D(Form("hnup%d",j), "", nbins, 0, nbins);
        hnun[j] = new TH1D(Form("hnun%d",j), "", nbins, 0, nbins);

        hnup_eRecoN[j] = new TH1D(Form("hnup_eRecoN%d",j), "", nbins, 0, nbins);
        hnun_eRecoN[j] = new TH1D(Form("hnun_eRecoN%d",j), "", nbins, 0, nbins);

        hnup_eRecoP[j] = new TH1D(Form("hnup_eRecoP%d",j), "", nbins, 0, nbins);
        hnun_eRecoP[j] = new TH1D(Form("hnun_eRecoP%d",j), "", nbins, 0, nbins);

        hnup_eRecoPi[j] = new TH1D(Form("hnup_eRecoPi%d",j), "", nbins, 0, nbins);
        hnun_eRecoPi[j] = new TH1D(Form("hnun_eRecoPi%d",j), "", nbins, 0, nbins);

        hLeN[j] = new TH1D(Form("hLeN%d",j), "", nbins, 0, nbins);
        hLeNp[j] = new TH1D(Form("hLeNp%d",j), "", nbins, 0, nbins);
        hLeNn[j] = new TH1D(Form("hLeNn%d",j), "", nbins, 0, nbins);


        hcov[j] = new TH2D(Form("hcov%d",j), Form("Total Detector Uncertainty (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        hfrcov[j] = new TH2D(Form("hfrcov%d",j), Form("Fractional Detector Uncertainty (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        hcov_had[j] = new TH2D(Form("hcov_had%d",j), Form("Detector Uncertainty from Hadrons (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        hfrcov_had[j] = new TH2D(Form("hfrcov_had%d",j), Form("Fractional Detector Uncertainty from Hadrons (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        hcov_lep[j] = new TH2D(Form("hcov_lep%d",j), Form("Detector Uncertainty from Leptons (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        hfrcov_lep[j] = new TH2D(Form("hfrcov_lep%d",j), Form("Fractional Detector Uncertainty from Leptons (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        hmu_cov[j] = new TH2D(Form("hmu_cov%d",j), Form("2%% Mu Reco Energy Scaling (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        hmu_frcov[j] = new TH2D(Form("hmu_frcov%d",j), Form("2%% Mu Reco Energy Scaling (fractional) (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        hel_cov[j] = new TH2D(Form("hel_cov%d",j), Form("2.5%% electron Energy Scaling (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        hel_frcov[j] = new TH2D(Form("hel_frcov%d",j), Form("2.5%% electron Energy Scaling (fractional) (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        hEres_cov[j] = new TH2D(Form("hEres_cov%d",j), Form("10%% Energy Resolution Scaling (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        hEres_frcov[j] = new TH2D(Form("hEres_frcov%d",j), Form("10%% Energy Resolution Scaling (fractional) (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        hnu_cov[j] = new TH2D(Form("hnu_cov%d",j), Form("2%% #nu cut Scaling (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        hnu_frcov[j] = new TH2D(Form("hnu_frcov%d",j), Form("2%% #nu cut Scaling (fractional) (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        heRecoN_cov[j] = new TH2D(Form("heRecoN_cov%d",j), Form("30%% Neutron Energy Scaling (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        heRecoN_frcov[j] = new TH2D(Form("heRecoN_frcov%d",j), Form("30%% Neutron Energy Scaling (fractional) (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        heRecoP_cov[j] = new TH2D(Form("heRecoP_cov%d",j), Form("5%% Proton Energy Scaling (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        heRecoP_frcov[j] = new TH2D(Form("heRecoP_frcov%d",j), Form("5%% Proton Energy Scaling (fractional) (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        heRecoPi_cov[j] = new TH2D(Form("heRecoPi_cov%d",j), Form("5%% Pion Energy Scaling (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        heRecoPi_frcov[j] = new TH2D(Form("heRecoPi_frcov%d",j), Form("5%% Pion Energy Scaling (fractional) (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        hLeNp_cov[j] = new TH2D(Form("hLeNp_cov%d",j), Form("Low-energy Neutron (#nu<%.2f GeV) (eRecoN < eN)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        hLeNp_frcov[j] = new TH2D(Form("hLeNp_frcov%d",j), Form("Low-energy Neutron (fractional) (#nu<%.2f GeV) (eRecoN < eN)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        hLeNn_cov[j] = new TH2D(Form("hLeNn_cov%d",j), Form("Low-energy Neutron (#nu<%.2f GeV) (eRecoN >= eN)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        hLeNn_frcov[j] = new TH2D(Form("hLeNn_frcov%d",j), Form("Low-energy Neutron (fractional) (#nu<%.2f GeV) (eRecoN >= eN)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        hLeN_cov[j] = new TH2D(Form("hLeN_cov%d",j), Form("Low-energy Neutron (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);
        hLeN_frcov[j] = new TH2D(Form("hLeN_frcov%d",j), Form("Low-energy Neutron (fractional) (#nu<%.2f GeV)",nucut[j]), nbins, 0, nbins, nbins, 0, nbins);

        for(int i=0; i<nbins; i++){
            if(i<nbins_CC){
                //std::cout << hmo->GetBinContent(i+1) << "\t" << hm->GetBinContent(i+1) << "\n";
                hnom[j]->SetBinContent(i+1, hm[j]->GetBinContent(i+1));
                
                hmp_tot[j]->SetBinContent(i+1, hmp[j]->GetBinContent(i+1));
                hmn_tot[j]->SetBinContent(i+1, hmn[j]->GetBinContent(i+1));

                hep_tot[j]->SetBinContent(i+1, hm[j]->GetBinContent(i+1));
                hen_tot[j]->SetBinContent(i+1, hm[j]->GetBinContent(i+1));

                hresp[j]->SetBinContent(i+1, hm_resp[j]->GetBinContent(i+1));
                hresn[j]->SetBinContent(i+1, hm_resn[j]->GetBinContent(i+1));

                hnup[j]->SetBinContent(i+1, hm_nup[j]->GetBinContent(i+1));
                hnun[j]->SetBinContent(i+1, hm_nun[j]->GetBinContent(i+1));

                hnup_eRecoN[j]->SetBinContent(i+1, hm_eRecoNp[j]->GetBinContent(i+1));
                hnun_eRecoN[j]->SetBinContent(i+1, hm_eRecoNn[j]->GetBinContent(i+1));

                hnup_eRecoP[j]->SetBinContent(i+1, hm_eRecoPp[j]->GetBinContent(i+1));
                hnun_eRecoP[j]->SetBinContent(i+1, hm_eRecoPn[j]->GetBinContent(i+1));

                hnup_eRecoPi[j]->SetBinContent(i+1, hm_eRecoPip[j]->GetBinContent(i+1));
                hnun_eRecoPi[j]->SetBinContent(i+1, hm_eRecoPin[j]->GetBinContent(i+1));

                //hLeNp[j]->SetBinContent(i+1, hm_LeNp[j]->GetBinContent(i+1));
                //hLeNn[j]->SetBinContent(i+1, hm_LeNn[j]->GetBinContent(i+1));
                hLeN[j]->SetBinContent(i+1, hm_LeN[j]->GetBinContent(i+1));

                //std::cout << i << "\t" << hnom->GetBinContent(i+1) << "\t" << hresp->GetBinContent(i+1) << "\n";
                
            }

            else if(i>=nbins_CC && i<2*nbins_CC){
                hnom[j]->SetBinContent(i+1, he[j]->GetBinContent(i-nbins_CC+1));

                hmp_tot[j]->SetBinContent(i+1, he[j]->GetBinContent(i-nbins_CC+1));
                hmn_tot[j]->SetBinContent(i+1, he[j]->GetBinContent(i-nbins_CC+1));

                hep_tot[j]->SetBinContent(i+1, hep[j]->GetBinContent(i-nbins_CC+1));
                hen_tot[j]->SetBinContent(i+1, hen[j]->GetBinContent(i-nbins_CC+1));

                hresp[j]->SetBinContent(i+1, he_resp[j]->GetBinContent(i-nbins_CC+1));
                hresn[j]->SetBinContent(i+1, he_resn[j]->GetBinContent(i-nbins_CC+1));

                hnup[j]->SetBinContent(i+1, he_nup[j]->GetBinContent(i-nbins_CC+1));
                hnun[j]->SetBinContent(i+1, he_nun[j]->GetBinContent(i-nbins_CC+1));

                hnup_eRecoN[j]->SetBinContent(i+1, he_eRecoNp[j]->GetBinContent(i-nbins_CC+1));
                hnun_eRecoN[j]->SetBinContent(i+1, he_eRecoNn[j]->GetBinContent(i-nbins_CC+1));

                hnup_eRecoP[j]->SetBinContent(i+1, he_eRecoPp[j]->GetBinContent(i-nbins_CC+1));
                hnun_eRecoP[j]->SetBinContent(i+1, he_eRecoPn[j]->GetBinContent(i-nbins_CC+1));

                hnup_eRecoPi[j]->SetBinContent(i+1, he_eRecoPip[j]->GetBinContent(i-nbins_CC+1));
                hnun_eRecoPi[j]->SetBinContent(i+1, he_eRecoPin[j]->GetBinContent(i-nbins_CC+1));

                //hLeNp[j]->SetBinContent(i+1, he_LeNp[j]->GetBinContent(i-nbins_CC+1));
                //hLeNn[j]->SetBinContent(i+1, he_LeNn[j]->GetBinContent(i-nbins_CC+1));
                hLeN[j]->SetBinContent(i+1, he_LeN[j]->GetBinContent(i-nbins_CC+1));

                //std::cout << i << "\t" << hnom->GetBinContent(i+1) << "\t" << hresp->GetBinContent(i+1) << "\n";
                
            }

            else{
                hnom[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));

                hmp_tot[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));
                hmn_tot[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));

                hep_tot[j]->SetBinContent(i+1, hnp->GetBinContent(i-2*nbins_CC+1));
                hen_tot[j]->SetBinContent(i+1, hnn->GetBinContent(i-2*nbins_CC+1));

                hresp[j]->SetBinContent(i+1, hn_resp->GetBinContent(i-2*nbins_CC+1));
                hresn[j]->SetBinContent(i+1, hn_resn->GetBinContent(i-2*nbins_CC+1));

                hnup[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));
                hnun[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));

                hnup_eRecoN[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));
                hnun_eRecoN[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));

                hnup_eRecoP[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));
                hnun_eRecoP[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));

                hnup_eRecoPi[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));
                hnun_eRecoPi[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));

                //hLeNp[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));
                //hLeNn[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));
                hLeN[j]->SetBinContent(i+1, hn->GetBinContent(i-2*nbins_CC+1));

                //std::cout << i << "\t" << hnom->GetBinContent(i+1) << "\t" << hresp->GetBinContent(i+1) << "\n";
                
            }
        }

        


        covariance_pn(hnom[j], hmp_tot[j], hmn_tot[j], hmu_cov[j], hmu_frcov[j], c, Form("hmu%d",j));
        covariance_pn(hnom[j], hep_tot[j], hen_tot[j], hel_cov[j], hel_frcov[j], c, Form("hel%d",j));
        covariance_pn(hnom[j], hresp[j], hresn[j], hEres_cov[j], hEres_frcov[j], c, Form("hEres%d",j));
        covariance_pn(hnom[j], hnup[j], hnun[j], hnu_cov[j], hnu_frcov[j], c, Form("hnu%d",j));
        covariance_pn(hnom[j], hnup_eRecoN[j], hnun_eRecoN[j], heRecoN_cov[j], heRecoN_frcov[j], c, Form("heRecoN%d",j));
        covariance_pn(hnom[j], hnup_eRecoP[j], hnun_eRecoP[j], heRecoP_cov[j], heRecoP_frcov[j], c, Form("heRecoP%d",j));
        covariance_pn(hnom[j], hnup_eRecoPi[j], hnun_eRecoPi[j], heRecoPi_cov[j], heRecoPi_frcov[j], c, Form("heRecoPi%d",j));

        
        //covariance(hnom[j], hLeNp[j], hLeNp_cov[j], hLeNp_frcov[j], c, Form("hLeNp_%d",j));
        //covariance(hnom[j], hLeNn[j], hLeNn_cov[j], hLeNn_frcov[j], c, Form("hLeNn_%d",j));
          
        covariance(hnom[j], hLeN[j], hLeN_cov[j], hLeN_frcov[j], c, Form("hLeN_%d",j));
        
        
        for(int ii = 0; ii < nbins; ++ii){
            for(int jj = 0; jj < nbins; ++jj){
                double mu = hmu_frcov[j]->GetBinContent(ii+1, jj+1);
                double el = hel_frcov[j]->GetBinContent(ii+1, jj+1);
                double res = hEres_frcov[j]->GetBinContent(ii+1, jj+1);
                double nu = hnu_frcov[j]->GetBinContent(ii+1, jj+1);
                double eRecoN = heRecoN_frcov[j]->GetBinContent(ii+1, jj+1);
                double eRecoP = heRecoP_frcov[j]->GetBinContent(ii+1, jj+1);
                double eRecoPi = heRecoPi_frcov[j]->GetBinContent(ii+1, jj+1);
                double LeN = hLeN_frcov[j]->GetBinContent(ii+1, jj+1);
                //double LeNn = hLeNn_frcov[j]->GetBinContent(ii+1, jj+1);

                hfrcov[j]->SetBinContent(ii+1, jj+1, mu + el + res + nu + eRecoN + eRecoP + eRecoPi + LeN);
                hfrcov_had[j]->SetBinContent(ii+1, jj+1, nu + eRecoN + eRecoP + eRecoPi + LeN);
                hfrcov_lep[j]->SetBinContent(ii+1, jj+1, mu + el + res);
            }
        } 

        
        for(int ii = 0; ii < nbins; ++ii){
            for(int jj = 0; jj < nbins; ++jj){
                double mu = hmu_cov[j]->GetBinContent(ii+1, jj+1);
                double el = hel_cov[j]->GetBinContent(ii+1, jj+1);
                double res = hEres_cov[j]->GetBinContent(ii+1, jj+1);
                double nu = hnu_cov[j]->GetBinContent(ii+1, jj+1);
                double eRecoN = heRecoN_cov[j]->GetBinContent(ii+1, jj+1);
                double eRecoP = heRecoP_cov[j]->GetBinContent(ii+1, jj+1);
                double eRecoPi = heRecoPi_cov[j]->GetBinContent(ii+1, jj+1);
                double LeN = hLeN_cov[j]->GetBinContent(ii+1, jj+1);

                hcov[j]->SetBinContent(ii+1, jj+1, mu + el + res + nu + eRecoN + eRecoP + eRecoPi + LeN);
                hcov_had[j]->SetBinContent(ii+1, jj+1, nu + eRecoN + eRecoP + eRecoPi + LeN);
                hcov_lep[j]->SetBinContent(ii+1, jj+1, mu + el + res);
            }
        } 
        
        
    }
    //covariance(hnom[4], hLeN[4], hLeN_cov[4], hLeN_frcov[4], c, Form("hLeN_4"));
    //plot1D(hnom[4], hLeN[4], "LeN4");
    //covariance_pn(hnom[4], hnup_eRecoP[4], hnun_eRecoP[4], heRecoP_cov[4], heRecoP_frcov[4], c, Form("eRecoP4_cov"));
    //covariance(hnom[1], hLeNp[1], hLeNp_cov[1], hLeNp_frcov[1], c, Form("hLeNp_1"));
   

    //makes lines to divide matrix to show the three different samples
    double x = hcov[0]->GetXaxis()->GetBinUpEdge(nbins_CC);
    double y1 = hcov[0]->GetYaxis()->GetBinLowEdge(1);
    double y2 = hcov[0]->GetYaxis()->GetBinUpEdge(nbins);

    double y = hcov[0]->GetYaxis()->GetBinUpEdge(nbins_CC);
    double x1 = hcov[0]->GetXaxis()->GetBinLowEdge(1);
    double x2 = hcov[0]->GetXaxis()->GetBinUpEdge(nbins);

    double w = hcov[0]->GetXaxis()->GetBinUpEdge(2*nbins_CC);
    double v1 = hcov[0]->GetYaxis()->GetBinLowEdge(1);
    double v2 = hcov[0]->GetYaxis()->GetBinUpEdge(nbins);

    double v = hcov[0]->GetYaxis()->GetBinUpEdge(2*nbins_CC);
    double w1 = hcov[0]->GetXaxis()->GetBinLowEdge(1);
    double w2 = hcov[0]->GetXaxis()->GetBinUpEdge(nbins);

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

    for(int j=0; j<N_nucut; j++){
        c->Clear();
        gPad->SetLogz();
        gPad->SetRightMargin(0.15);
        hfrcov[j]->SetStats(0);
        hfrcov[j]->SetMaximum(1.);
        hfrcov[j]->SetMinimum(1E-4);
        hfrcov[j]->Draw("colz");
        l1->Draw("same");
        l2->Draw("same");
        l3->Draw("same");
        l4->Draw("same");
        c->SaveAs(Form("%s/uncertainties/det_covmtr/plots/fractional_detcov%d.png",data_path,j));

        c->Clear();
        //gPad->SetLogz();
        gPad->SetRightMargin(0.15);
        hfrcov_had[j]->SetStats(0);
        hfrcov_had[j]->SetMaximum(1.);
        hfrcov_had[j]->SetMinimum(1E-5);
        hfrcov_had[j]->Draw("colz");
        l1->Draw("same");
        l2->Draw("same");
        l3->Draw("same");
        l4->Draw("same");
        c->SaveAs(Form("%s/uncertainties/det_covmtr/plots/fractional_detcov_had%d.png",data_path,j));

        c->Clear();
        //gPad->SetLogz();
        gPad->SetRightMargin(0.15);
        hfrcov_lep[j]->SetStats(0);
        hfrcov_lep[j]->SetMaximum(1.);
        hfrcov_lep[j]->SetMinimum(1E-5);
        hfrcov_lep[j]->Draw("colz");
        l1->Draw("same");
        l2->Draw("same");
        l3->Draw("same");
        l4->Draw("same");
        c->SaveAs(Form("%s/uncertainties/det_covmtr/plots/fractional_detcov_lep%d.png",data_path,j));
    }

    c->Close();



    
    for(int j=0; j<N_nucut; j++){
        hcov[j]->Write();

        hfrcov[j]->Write();
        hfrcov_had[j]->Write();
        hfrcov_lep[j]->Write();

        hmu_frcov[j]->Write();
        hel_frcov[j]->Write();
        hEres_frcov[j]->Write();
        hnu_frcov[j]->Write();
        heRecoN_frcov[j]->Write();
        heRecoP_frcov[j]->Write();
        heRecoPi_frcov[j]->Write();
        hLeN_frcov[j]->Write();
        //hLeNn_frcov[j]->Write();
    }
    out->Close();



    return(0);
}