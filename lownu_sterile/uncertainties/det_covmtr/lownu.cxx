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
//#include "DUNEStyle.h"

static const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu/uncertainties/det_covmtr";

void plot1D_pn(TH1D* ho, TH1D *hp, TH1D *hn, const char* name){
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
    hp->SetLineColor(kBlue);
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
    lg->AddEntry(hp, "shift");
    lg->Draw("same");
    c->SaveAs(Form("%s.png",name));
}

void plot1_1D(TH1D *h, const char *name){
    h->SetStats(0);
    TCanvas *c = new TCanvas("c", "", 800, 600);
    c->SetGrid();
    c->SetLogy();
    h->Draw();
    c->SaveAs(Form("%s.png",name));
}
void plot1_2D(TH2D *h, const char *name){
    h->SetStats(0);
    TCanvas *c = new TCanvas("c", "", 800, 600);
    c->SetGrid();
    c->SetLogx();
    c->SetLogy();
    h->Draw("colz");
    c->SaveAs(Form("%s.png",name));
}

int main()
{
// CC and nu+e ROOT files ---------------------------------------------------------------------------------------------
    TChain * tree = new TChain( "cafTree", "cafTree" );
    TChain * meta = new TChain( "meta", "meta" );

    for(int i = 0; i<400;i++){
        tree->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%02d.root", i) );
        meta->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%02d.root", i) );
    }

    //tree->Print();

//protons on target for CC and nu+e------------------------------------------------------------------------------------
    double total_pot = 0.;
    double pot;
    meta->SetBranchAddress( "pot", &pot );
    const int Nfiles = meta->GetEntries();
    for( int ii = 0; ii < Nfiles; ++ii ) {
        meta->GetEntry(ii);
        total_pot += pot;
    }
    double yrPOT = 1.1E21;
    double scalePOT = yrPOT/total_pot;

    const int N_nucut = 5;
    double nucut[N_nucut] = {100., 10., 1., 0.5, 0.3};

//CC branches ---------------------------------------------------------------------------------------------------------
    double vtx_x, vtx_y, vtx_z;
    double Ev_reco, Elep_reco;
    int reco_numu, reco_nue;
    int muon_contained, muon_tracker, muon_ecal;
    double LepE, Ev;
    double LepNuAngle;
    double eN, eP, eDepN, ePip, ePim, ePi0, eOther;
    double eRecoN, eRecoP, eRecoPip, eRecoPim, eRecoPi0, eRecoOther;

    tree->SetBranchAddress( "vtx_x", &vtx_x );
    tree->SetBranchAddress( "vtx_y", &vtx_y );
    tree->SetBranchAddress( "vtx_z", &vtx_z );
    tree->SetBranchAddress( "Ev_reco", &Ev_reco );
    tree->SetBranchAddress( "Elep_reco", &Elep_reco );
    tree->SetBranchAddress( "reco_numu", &reco_numu );
    tree->SetBranchAddress( "reco_nue", &reco_nue );
    tree->SetBranchAddress( "muon_contained", &muon_contained );
    tree->SetBranchAddress( "muon_tracker", &muon_tracker );
    tree->SetBranchAddress( "muon_ecal", &muon_ecal );
    tree->SetBranchAddress( "LepE", &LepE );
    tree->SetBranchAddress( "LepNuAngle", &LepNuAngle );


    tree->SetBranchAddress( "Ev", &Ev );
    tree->SetBranchAddress( "eP", &eP );
    tree->SetBranchAddress( "eN", &eN );
    //tree->SetBranchAddress( "eDepN", &eDepN );
    tree->SetBranchAddress( "ePip", &ePip );
    tree->SetBranchAddress( "ePim", &ePim );
    tree->SetBranchAddress( "ePi0", &ePi0 );
    tree->SetBranchAddress( "eOther", &eOther );
    tree->SetBranchAddress( "eRecoP", &eRecoP );
    tree->SetBranchAddress( "eRecoN", &eRecoN );
    //tree->SetBranchAddress( "eDepN", &eDepN );
    tree->SetBranchAddress( "eRecoPip", &eRecoPip );
    tree->SetBranchAddress( "eRecoPim", &eRecoPim );
    tree->SetBranchAddress( "eRecoPi0", &eRecoPi0 );
    tree->SetBranchAddress( "eRecoOther", &eRecoOther );

//---------------------------------------------------------------------------------------------------------------------


//CC binning and smear setup. small hist ------------------------------------------------------------------------------
    const Int_t nbinsCC = 56;
    double CCEdges[nbinsCC+1];
    CCEdges[0] = 0.;
    CCEdges[1] = 0.3;

    for(int i=1; i<nbinsCC+1; i++)
    {
      if(i<38)               CCEdges[i+1] = CCEdges[i] + 0.1;
      else if(i>=38 && i<43) CCEdges[i+1] = CCEdges[i] + 0.2;
      else if(i>=43 && i<48) CCEdges[i+1] = CCEdges[i] + 0.4;
      else if(i>=48 && i<53) CCEdges[i+1] = CCEdges[i] + 0.8;
      else if(i>=53 && i<55) CCEdges[i+1] = CCEdges[i] + 1.5;
      else                   CCEdges[i+1] = CCEdges[i] + 2.0;
    }

    TRandom3 *r = new TRandom3(12345);

    TH1D* hm[N_nucut]; 
    TH1D* hmp[N_nucut]; TH1D* hmn[N_nucut];
    TH1D* hm_resp[N_nucut]; TH1D* hm_resn[N_nucut];

    TH1D* hm_nup[N_nucut]; TH1D* hm_nun[N_nucut];
    TH1D* hm_eRecoNp[N_nucut]; TH1D* hm_eRecoNn[N_nucut];
    TH1D* hm_eRecoPp[N_nucut]; TH1D* hm_eRecoPn[N_nucut];
    TH1D* hm_eRecoPip[N_nucut]; TH1D* hm_eRecoPin[N_nucut];
    TH1D* hm_LeNp[N_nucut]; TH1D* hm_LeNn[N_nucut]; 
    TH1D* hm_LeN[N_nucut]; TH1D* hm_HeN[N_nucut];

    TH1D* he[N_nucut]; 
    TH1D* hep[N_nucut]; TH1D* hen[N_nucut];
    TH1D* he_resp[N_nucut]; TH1D* he_resn[N_nucut];

    TH1D* he_nup[N_nucut]; TH1D* he_nun[N_nucut];
    TH1D* he_eRecoNp[N_nucut]; TH1D* he_eRecoNn[N_nucut];
    TH1D* he_eRecoPp[N_nucut]; TH1D* he_eRecoPn[N_nucut];
    TH1D* he_eRecoPip[N_nucut]; TH1D* he_eRecoPin[N_nucut];
    TH1D* he_LeNp[N_nucut]; TH1D* he_LeNn[N_nucut]; 
    TH1D* he_LeN[N_nucut]; TH1D* he_HeN[N_nucut];


    for(int inu=0; inu<N_nucut; inu++){
        hm[inu] = new TH1D(Form("hm%d",inu), Form("nominal mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hmp[inu] = new TH1D(Form("hmp%d",inu), Form("+2%% lepton energy mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hmn[inu] = new TH1D(Form("hmn%d",inu), Form("-2%% lepton energy mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_resp[inu] = new TH1D(Form("hm_resp%d",inu), Form("+10.0%% energy resolution mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_resn[inu] = new TH1D(Form("hm_resn%d",inu), Form("-10.0%% energy resolution mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);

        hm_nup[inu] = new TH1D(Form("hm_nup%d",inu), Form("+2.0%% #nu cut mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_nun[inu] = new TH1D(Form("hm_nun%d",inu), Form("-2.0%% #nu cut mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_eRecoNp[inu] = new TH1D(Form("hm_eRecoNp%d",inu), Form("+30.0%% eRecoN mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_eRecoNn[inu] = new TH1D(Form("hm_eRecoNn%d",inu), Form("-30.0%% eRecoN mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_eRecoPp[inu] = new TH1D(Form("hm_eRecoPp%d",inu), Form("+5.0%% eRecoP mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_eRecoPn[inu] = new TH1D(Form("hm_eRecoPn%d",inu), Form("-5.0%% eRecoP mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_eRecoPip[inu] = new TH1D(Form("hm_eRecoPip%d",inu), Form("+5.0%% eRecoPi mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_eRecoPin[inu] = new TH1D(Form("hm_eRecoPin%d",inu), Form("-5.0%% eRecoPi mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_LeNp[inu] = new TH1D(Form("hm_LeNp%d",inu), Form("Low energy neutron effect mCC (#nu < %.2f GeV) (+1#sigma)",nucut[inu]), nbinsCC, CCEdges);
        hm_LeNn[inu] = new TH1D(Form("hm_LeNn%d",inu), Form("Low energy neutron effect mCC (#nu < %.2f GeV) (-1#sigma)",nucut[inu]), nbinsCC, CCEdges);
        hm_LeN[inu] = new TH1D(Form("hm_LeN%d",inu), Form("Low energy neutron effect mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hm_HeN[inu] = new TH1D(Form("hm_HeN%d",inu), Form("High energy neutron effect mCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);



        he[inu] = new TH1D(Form("he%d",inu), Form("nominal eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hep[inu] = new TH1D(Form("hep%d",inu), Form("+2%% lepton energy eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        hen[inu] = new TH1D(Form("hen%d",inu), Form("-2%% lepton energy eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_resp[inu] = new TH1D(Form("he_resp%d",inu), Form("+10.0%% energy resolution eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_resn[inu] = new TH1D(Form("he_resn%d",inu), Form("-10.0%% energy resolution eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);

        he_nup[inu] = new TH1D(Form("he_nup%d",inu), Form("+2.0%% #nu cut eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_nun[inu] = new TH1D(Form("he_nun%d",inu), Form("-2.0%% #nu cut eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_eRecoNp[inu] = new TH1D(Form("he_eRecoNp%d",inu), Form("+30.0%% eRecoN eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_eRecoNn[inu] = new TH1D(Form("he_eRecoNn%d",inu), Form("-30.0%% eRecoN eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_eRecoPp[inu] = new TH1D(Form("he_eRecoPp%d",inu), Form("+5.0%% eRecoP eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_eRecoPn[inu] = new TH1D(Form("he_eRecoPn%d",inu), Form("-5.0%% eRecoP eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_eRecoPip[inu] = new TH1D(Form("he_eRecoPip%d",inu), Form("+5.0%% eRecoPi eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_eRecoPin[inu] = new TH1D(Form("he_eRecoPin%d",inu), Form("-5.0%% eRecoPi eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_LeNp[inu] = new TH1D(Form("he_LeNp%d",inu), Form("Low energy neutron effect eCC (#nu < %.2f GeV) (+1#sigma)",nucut[inu]), nbinsCC, CCEdges);
        he_LeNn[inu] = new TH1D(Form("he_LeNn%d",inu), Form("Low energy neutron effect eCC (#nu < %.2f GeV) (-1#sigma)",nucut[inu]), nbinsCC, CCEdges);
        he_LeN[inu] = new TH1D(Form("he_LeN%d",inu), Form("Low energy neutron effect eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        he_HeN[inu] = new TH1D(Form("he_HeN%d",inu), Form("High energy neutron effect eCC (#nu < %.2f GeV)",nucut[inu]), nbinsCC, CCEdges);
        
    }    
    


    TH1D *hnu = new TH1D("hnu","",100,0,20);
    TH1D *hnu_sum = new TH1D("hnu_sum","",100,0,20);
    TH1D *hnu_sum_reco = new TH1D("hnu_sum_reco","",100,0,20);
    TH1D *hN = new TH1D("hN","",100,0,50);
    TH1D *hRecoN = new TH1D("hRecoN","",100,0,50);
    TH2D *hCorN = new TH2D("hCorN","",100,0.,2.,100,0.,1.5);



//fill CC samples (mu and eCC) ----------------------------------------------------------------------------------------
    const int N = tree->GetEntries();
    //const int N = 20;
    for(int ii = 0; ii < N; ++ii){
        if( ii % 100000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
        tree->GetEntry(ii);
        if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;

        double nu = Ev_reco - Elep_reco;
        double nu_true = eP + eN + ePip + ePim + ePi0 + eOther;
        double nu_reco = eRecoP + eRecoN + eRecoPip + eRecoPim + eRecoPi0 + eRecoOther;

        //Smearing true energy to get reco energy (5%) because of existing problem in CAF
        double LepE_sm = r->Gaus(LepE, 0.05); 
        double Ev_sm = r->Gaus(Ev, 0.05);


        // Energy resolution scaling 10%
        double Eresp = LepE_sm + (1.+0.1)*(LepE_sm-LepE);
        double Eresn = LepE_sm + (1.-0.1)*(LepE_sm-LepE);


        // Nu cut scaling 2%
        double nup = nu_reco * (1.+0.02);
        double nun = nu_reco * (1.-0.02);


        // Nu cut affected by eRecoN scaling 30%
        double nup_eN = nu_reco + eRecoN * 0.3;
        double nun_eN = nu_reco - eRecoN * 0.3;

        // Nu cut affected by eRecoP scaling 5%
        double nup_eP = nu_reco + eRecoP * 0.05;
        double nun_eP = nu_reco - eRecoP * 0.05;

        // Nu cut affected by eRecoPi scaling 5%
        double nup_ePi = nu_reco + (eRecoPip + eRecoPim) * 0.05;
        double nun_ePi = nu_reco - (eRecoPip + eRecoPim) * 0.05;


        hN->Fill(eN);
        hRecoN->Fill(eRecoN);
        hCorN->Fill(eN, eRecoN/eN);
        
        //check for numu
        if( reco_numu && (muon_contained || muon_tracker || muon_ecal )){
            
            // Low nu cut for numuCC 
            for(int j=0; j<N_nucut; j++){
                if(nu_reco < nucut[j]){
                    // Nominal
                    hm[j]->Fill(LepE_sm);

                    // Energy resolution scaling 10%
                    hm_resp[j]->Fill(Eresp);
                    hm_resn[j]->Fill(Eresn);

                    // Energy scaling 2% for muon
                    double m_Ep = LepE_sm * (1.+0.02);
                    double m_En = LepE_sm * (1.-0.02);
                    if(muon_contained || muon_tracker){
                        hmp[j]->Fill(m_Ep);
                        hmn[j]->Fill(m_En);
                    }
                }

                // Nu cut scaling 2%
                if(nup < nucut[j]) hm_nup[j]->Fill(LepE_sm);
                if(nun < nucut[j]) hm_nun[j]->Fill(LepE_sm);

                // eRecoN scaling 30%
                if(nup_eN < nucut[j]) hm_eRecoNp[j]->Fill(LepE_sm);
                if(nun_eN < nucut[j]) hm_eRecoNn[j]->Fill(LepE_sm);

                // eRecoP scaling 5%
                if(nup_eP < nucut[j]) hm_eRecoPp[j]->Fill(LepE_sm);
                if(nun_eP < nucut[j]) hm_eRecoPn[j]->Fill(LepE_sm);

                // eRecoPi scaling 5%
                if(nup_ePi < nucut[j]) hm_eRecoPip[j]->Fill(LepE_sm);
                if(nun_ePi < nucut[j]) hm_eRecoPin[j]->Fill(LepE_sm);


                
                //printf("eRecoP: nucut = %.2f, nu_reco = %f, nup_eP = %f, nun_eP = %f\n", nucut[j], nu_reco, nup_eP, nun_eP);

                // Low eN
                if(eN >= 0.2){
                    if(nu_reco < nucut[j]){
                        hm_HeN[j]->Fill(LepE_sm);
                        //hm_LeNp[j]->Fill(LepE_sm);
                        //hm_LeNn[j]->Fill(LepE_sm);
                    }
                }
                else if(eN < 0.2){
                    //if(nu_reco < nucut[j]) hm_LeN[j]->Fill(LepE_sm);
                    double nu_LeN = nu_reco - eRecoN + r->Uniform(0, eN);
                    if(nu_LeN < nucut[j]) hm_LeN[j]->Fill(LepE_sm);
                    /*
                    // Dial selection low eN
                    if(eRecoN < eN){
                        double nu_LeN = nu_reco - eRecoN + r->Uniform(eRecoN, eN);
                        if(nu_LeN < nucut[j]) hm_LeNp[j]->Fill(LepE_sm);
                    }
                    else{
                        double nu_LeN = nu_reco - eRecoN + r->Uniform(0., eN);
                        if(nu_LeN < nucut[j]) hm_LeNn[j]->Fill(LepE_sm);
                    }
                    */
                }
            }
        }


        //now do nue 
        else if(reco_nue){
            
            // Smearing true energy (electron energy and nue energy) to get reco energy (5%) because of some existing problem with this energy in CAF
            // Low nu cut for nueCC
            for(int j=0; j<N_nucut; j++){
                if(nu_reco < nucut[j]){
                    // Nominal
                    if(LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he[j]->Fill(LepE_sm);

                    // Energy resolution scaling 10%
                    if(Eresp*LepNuAngle*LepNuAngle*1E3 > 3.) he_resp[j]->Fill(Eresp);
                    if(Eresn*LepNuAngle*LepNuAngle*1E3 > 3.) he_resn[j]->Fill(Eresn);

                    // Energy scaling 2.5% for electron
                    double e_Ep = LepE_sm * (1+0.025);
                    double e_En = LepE_sm * (1-0.025);
                    if(e_Ep*LepNuAngle*LepNuAngle*1E3 > 3.) hep[j]->Fill(e_Ep);
                    if(e_En*LepNuAngle*LepNuAngle*1E3 > 3.) hen[j]->Fill(e_En);
                }

                // Nu cut scaling 2%
                if(nup < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_nup[j]->Fill(LepE_sm);
                if(nun < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_nun[j]->Fill(LepE_sm);

                // eRecoN scaling 30%
                if(nup_eN < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_eRecoNp[j]->Fill(LepE_sm);
                if(nun_eN < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_eRecoNn[j]->Fill(LepE_sm);

                // eRecoP scaling 5%
                if(nup_eP < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_eRecoPp[j]->Fill(LepE_sm);
                if(nun_eP < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_eRecoPn[j]->Fill(LepE_sm);

                // eRecoPi scaling 5%
                if(nup_ePi < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_eRecoPip[j]->Fill(LepE_sm);
                if(nun_ePi < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_eRecoPin[j]->Fill(LepE_sm);


                // Low eN
                if(eN >= 0.2){
                    if(nu_reco < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.){
                        he_HeN[j]->Fill(LepE_sm);
                        //he_LeNp[j]->Fill(LepE_sm);
                        //he_LeNn[j]->Fill(LepE_sm);
                    }
                }
                else if(eN < 0.2){
                    //if(nu_reco < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_LeN[j]->Fill(LepE_sm);
                    double nu_LeN = nu_reco - eRecoN + r->Uniform(0, eN);
                    if(nu_LeN < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_LeN[j]->Fill(LepE_sm);
                    /*
                    if(eRecoN < eN){
                        double nu_LeN = nu_reco - eRecoN + r->Uniform(eRecoN, eN);
                        if(nu_LeN < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_LeNp[j]->Fill(LepE_sm);
                    }
                    else{
                        double nu_LeN = nu_reco - eRecoN + r->Uniform(0., eN);
                        if(nu_LeN < nucut[j] && LepE_sm*LepNuAngle*LepNuAngle*1E3 > 3.) he_LeNn[j]->Fill(LepE_sm);
                    }
                    */
                }


            }
        }
    }


    /*
    plot1_1D(hN, "N");
    plot1_1D(hRecoN, "RecoN");
    plot1_2D(hCorN, "CorN");
    
    
    hnu->SetLineColor(kBlue);
    hnu->SetLineWidth(2.);
    hnu_sum->SetLineColor(kRed);
    hnu_sum->SetLineWidth(2.);
    hnu_sum_reco->SetLineColor(kBlack);
    hnu_sum_reco->SetLineWidth(2.);
    hnu->SetStats(0);
    hnu_sum->SetStats(0);
    hnu_sum_reco->SetStats(0);
    hnu->SetTitle("Energy Transfer nu");
    hnu->GetXaxis()->SetTitle("nu (GeV)");
    TCanvas *c = new TCanvas("c","",800,600);
    c->SetGrid();
    c->SetLogy();
    hnu->Draw();
    hnu_sum->Draw("same");
    hnu_sum_reco->Draw("same");
    TLegend *lg = new TLegend(0.65,0.75,0.9,0.9);
    lg->AddEntry(hnu, "Enu_reco - Elep_reco");
    lg->AddEntry(hnu_sum, "energy sum");
    lg->AddEntry(hnu_sum_reco, "reco energy sum");
    lg->Draw("same");
    c->SaveAs("nu.png");

    plot_test(hN, "eN");
    plot_test(hP, "eP");
    //plot_test(hDepN, "eDepN");
    */

//---------------------------------------------------------------------------------------------------------------------



// Scaling
for(int j=0; j<N_nucut; j++){
    hm_LeN[j]->Add(hm_HeN[j]);
    he_LeN[j]->Add(he_HeN[j]);

    hm[j]->Scale(scalePOT);
    hmp[j]->Scale(scalePOT);  //positive mu shift
    hmn[j]->Scale(scalePOT);  //negative mu shift
    hm_resp[j]->Scale(scalePOT);
    hm_resn[j]->Scale(scalePOT);

    hm_nup[j]->Scale(scalePOT);
    hm_nun[j]->Scale(scalePOT);
    hm_eRecoNp[j]->Scale(scalePOT);
    hm_eRecoNn[j]->Scale(scalePOT);
    hm_eRecoPp[j]->Scale(scalePOT);
    hm_eRecoPn[j]->Scale(scalePOT);
    hm_eRecoPip[j]->Scale(scalePOT);
    hm_eRecoPin[j]->Scale(scalePOT);
    //hm_LeNp[j]->Scale(scalePOT);
    //hm_LeNn[j]->Scale(scalePOT);
    hm_LeN[j]->Scale(scalePOT);

    

    he[j]->Scale(scalePOT);  //unshifted e energies
    hep[j]->Scale(scalePOT);  //positive electron shift
    hen[j]->Scale(scalePOT);  //negative electron shift
    he_resp[j]->Scale(scalePOT);
    he_resn[j]->Scale(scalePOT);

    he_nup[j]->Scale(scalePOT);
    he_nun[j]->Scale(scalePOT);
    he_eRecoNp[j]->Scale(scalePOT);
    he_eRecoNn[j]->Scale(scalePOT);
    he_eRecoPp[j]->Scale(scalePOT);
    he_eRecoPn[j]->Scale(scalePOT);
    he_eRecoPip[j]->Scale(scalePOT);
    he_eRecoPin[j]->Scale(scalePOT);
    //he_LeNp[j]->Scale(scalePOT);
    //he_LeNn[j]->Scale(scalePOT);
    he_LeN[j]->Scale(scalePOT);
}

//plot1D(hm[1], hm_eRecoPp[1], hm_eRecoPn[1], "eRecoP");
//plot1D(hm[4], hm_LeN[4], "test");


//save to new ROOT file -----------------------------------------------------------------------------------------------
    TFile *out = new TFile(Form("%s/outputCC_detCov_0102.root",data_path),"RECREATE");
    //hmo->Write();  //unshifted mu energies
    
    for(int j=0; j<N_nucut; j++){
        hm[j]->Write();
        hmp[j]->Write();  //positive mu shift
        hmn[j]->Write();  //negative mu shift
        hm_resp[j]->Write();
        hm_resn[j]->Write();

        hm_nup[j]->Write();
        hm_nun[j]->Write();
        hm_eRecoNp[j]->Write();
        hm_eRecoNn[j]->Write();
        hm_eRecoPp[j]->Write();
        hm_eRecoPn[j]->Write();
        hm_eRecoPip[j]->Write();
        hm_eRecoPin[j]->Write();
        //hm_LeNp[j]->Write();
        //hm_LeNn[j]->Write();
        hm_LeN[j]->Write();

        

        he[j]->Write();  //unshifted e energies
        hep[j]->Write();  //positive electron shift
        hen[j]->Write();  //negative electron shift
        he_resp[j]->Write();
        he_resn[j]->Write();

        he_nup[j]->Write();
        he_nun[j]->Write();
        he_eRecoNp[j]->Write();
        he_eRecoNn[j]->Write();
        he_eRecoPp[j]->Write();
        he_eRecoPn[j]->Write();
        he_eRecoPip[j]->Write();
        he_eRecoPin[j]->Write();
        //he_LeNp[j]->Write();
        //he_LeNn[j]->Write();
        he_LeN[j]->Write();
    }

    //total POT is saved for CC and nu+e samples to the file
    TParameter<double> totalPOT("total_pot", total_pot);
    totalPOT.Write();
    out->Close();

    //plot1D(hmo, hmp, hmn, "muon");
    //plot1D(heo, hep, hen, "electron");
    

//---------------------------------------------------------------------------------------------------------------------


    return(0);
}