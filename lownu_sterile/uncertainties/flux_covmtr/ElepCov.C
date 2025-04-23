#include "DUNEStyle.h"



static const int N = 10000; // number of universes

static const int nbins_CC = 56;
static const int nbins_nue = 8;
static const int nbins_tot = 2*nbins_CC + nbins_nue;

static const int n_mu = 19; // number of muon bins in covariance
static const int n_e = 7; // number of electron bins in covariance
static const int nbins = 2*n_mu + 2*n_e;

// bin edges -- fix these to be whatever they actually are
const double mubins[20] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,12.,16.,20.,40.,100.};
const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};
  

TMatrixD covmx( nbins, nbins );
TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, nbins );

TMatrixD E_m(N, nbins_CC);
TMatrixD E_e(N, nbins_CC);
TMatrixD E_nue(N, nbins_nue);
  
TMatrixD ECovars_mm   ( nbins_CC, nbins_CC );
TMatrixD ECovars_ee   ( nbins_CC, nbins_CC );
TMatrixD ECovars_nuenue  ( nbins_nue, nbins_nue );
TMatrixD ECovars_me  ( nbins_CC, nbins_CC );
TMatrixD ECovars_mnue( nbins_CC, nbins_nue );
TMatrixD ECovars_em  ( nbins_CC, nbins_CC );
TMatrixD ECovars_enue( nbins_CC, nbins_nue );
TMatrixD ECovars_nuem( nbins_nue, nbins_CC );
TMatrixD ECovars_nuee( nbins_nue, nbins_CC );

TMatrixD frECovars_mm   ( nbins_CC, nbins_CC );
TMatrixD frECovars_ee   ( nbins_CC, nbins_CC );
TMatrixD frECovars_nuenue  ( nbins_nue, nbins_nue );
TMatrixD frECovars_me  ( nbins_CC, nbins_CC );
TMatrixD frECovars_mnue( nbins_CC, nbins_nue );
TMatrixD frECovars_em  ( nbins_CC, nbins_CC );
TMatrixD frECovars_enue( nbins_CC, nbins_nue );
TMatrixD frECovars_nuem( nbins_nue, nbins_CC );
TMatrixD frECovars_nuee( nbins_nue, nbins_CC );

TMatrixD ECorrel_mm   ( nbins_CC, nbins_CC );
TMatrixD ECorrel_ee   ( nbins_CC, nbins_CC );
TMatrixD ECorrel_nuenue  ( nbins_nue, nbins_nue );
TMatrixD ECorrel_me  ( nbins_CC, nbins_CC );
TMatrixD ECorrel_mnue( nbins_CC, nbins_nue );
TMatrixD ECorrel_em  ( nbins_CC, nbins_CC );
TMatrixD ECorrel_enue( nbins_CC, nbins_nue );
TMatrixD ECorrel_nuem( nbins_nue, nbins_CC );
TMatrixD ECorrel_nuee( nbins_nue, nbins_CC );


TMatrixD ECovars( nbins_tot, nbins_tot );
TMatrixD frECovars( nbins_tot, nbins_tot );
TMatrixD ECorrel( nbins_tot, nbins_tot );


/*
void plotCov(TH2D* cov, const char *name, TCanvas *c){
  //int nbins = hnom->GetNbinsX();
  //std::cout << "nbins = " << nbins << "\n";

  cov->GetXaxis()->SetTitle("Bin number");
  cov->GetYaxis()->SetTitle("Bin number");

  c->Clear();
  
  //makes lines to divide matrix to show the three different samples
    double x = cov->GetXaxis()->GetBinUpEdge(nbins_CC);
    double y1 = cov->GetYaxis()->GetBinLowEdge(1);
    double y2 = cov->GetYaxis()->GetBinUpEdge(nbins_tot);

    double y = cov->GetYaxis()->GetBinUpEdge(nbins_CC);
    double x1 = cov->GetXaxis()->GetBinLowEdge(1);
    double x2 = cov->GetXaxis()->GetBinUpEdge(nbins_tot);

    double w = cov->GetXaxis()->GetBinUpEdge(2*nbins_CC);
    double v1 = cov->GetYaxis()->GetBinLowEdge(1);
    double v2 = cov->GetYaxis()->GetBinUpEdge(nbins_tot);

    double v = cov->GetYaxis()->GetBinUpEdge(2*nbins_CC);
    double w1 = cov->GetXaxis()->GetBinLowEdge(1);
    double w2 = cov->GetXaxis()->GetBinUpEdge(nbins_tot);

    TLine *l1 = new TLine(x, y1, x, y2);  //splits muCC from the rest in x axis 
    TLine *l2 = new TLine(x1, y, x2, y);  //splits muCC from the rest in y axis
    TLine *l3 = new TLine(w, v1, w, v2);  //splits nu+e from the rest in x axis
    TLine *l4 = new TLine(w1, v, w2, v);  //splits nu+e from the rest in y axis

    //TCanvas *c = new TCanvas("c", "", 800, 600);
    cov->SetMinimum(1E-3);
    cov->SetMaximum(1E11);
    cov->SetStats(0);
    cov->Draw("colz");
    c->SaveAs(Form("%s.png",name));
    c->SetLogz();
    cov->Draw("colz");
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
    c->SaveAs(Form("%s_logz.png",name));
}





void plotfrCov(TH2D* cov, const char *name, TCanvas *c){
  //int nbins = hnom->GetNbinsX();
  //std::cout << "nbins = " << nbins << "\n";

  cov->GetXaxis()->SetTitle("Bin number");
  cov->GetYaxis()->SetTitle("Bin number");

  c->Clear();
  
  //makes lines to divide matrix to show the three different samples
    double x = cov->GetXaxis()->GetBinUpEdge(nbins_CC);
    double y1 = cov->GetYaxis()->GetBinLowEdge(1);
    double y2 = cov->GetYaxis()->GetBinUpEdge(nbins_tot);

    double y = cov->GetYaxis()->GetBinUpEdge(nbins_CC);
    double x1 = cov->GetXaxis()->GetBinLowEdge(1);
    double x2 = cov->GetXaxis()->GetBinUpEdge(nbins_tot);

    double w = cov->GetXaxis()->GetBinUpEdge(2*nbins_CC);
    double v1 = cov->GetYaxis()->GetBinLowEdge(1);
    double v2 = cov->GetYaxis()->GetBinUpEdge(nbins_tot);

    double v = cov->GetYaxis()->GetBinUpEdge(2*nbins_CC);
    double w1 = cov->GetXaxis()->GetBinLowEdge(1);
    double w2 = cov->GetXaxis()->GetBinUpEdge(nbins_tot);

    TLine *l1 = new TLine(x, y1, x, y2);  //splits muCC from the rest in x axis 
    TLine *l2 = new TLine(x1, y, x2, y);  //splits muCC from the rest in y axis
    TLine *l3 = new TLine(w, v1, w, v2);  //splits nu+e from the rest in x axis
    TLine *l4 = new TLine(w1, v, w2, v);  //splits nu+e from the rest in y axis

    //TCanvas *c = new TCanvas("c", "", 800, 600);
    cov->SetMinimum(1E-3);
    cov->SetMaximum(1.);
    cov->SetStats(0);
    cov->Draw("colz");
    c->SaveAs(Form("%s.png",name));
    c->SetLogz();
    cov->Draw("colz");
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
    c->SaveAs(Form("%s_logz.png",name));
}



void plotCor(TH2D* cov, const char *name, TCanvas *c){
  //int nbins = hnom->GetNbinsX();
  //std::cout << "nbins = " << nbins << "\n";

  cov->GetXaxis()->SetTitle("Bin number");
  cov->GetYaxis()->SetTitle("Bin number");

  c->Clear();
  
  //makes lines to divide matrix to show the three different samples
    double x = cov->GetXaxis()->GetBinUpEdge(nbins_CC);
    double y1 = cov->GetYaxis()->GetBinLowEdge(1);
    double y2 = cov->GetYaxis()->GetBinUpEdge(nbins_tot);

    double y = cov->GetYaxis()->GetBinUpEdge(nbins_CC);
    double x1 = cov->GetXaxis()->GetBinLowEdge(1);
    double x2 = cov->GetXaxis()->GetBinUpEdge(nbins_tot);

    double w = cov->GetXaxis()->GetBinUpEdge(2*nbins_CC);
    double v1 = cov->GetYaxis()->GetBinLowEdge(1);
    double v2 = cov->GetYaxis()->GetBinUpEdge(nbins_tot);

    double v = cov->GetYaxis()->GetBinUpEdge(2*nbins_CC);
    double w1 = cov->GetXaxis()->GetBinLowEdge(1);
    double w2 = cov->GetXaxis()->GetBinUpEdge(nbins_tot);

    TLine *l1 = new TLine(x, y1, x, y2);  //splits muCC from the rest in x axis 
    TLine *l2 = new TLine(x1, y, x2, y);  //splits muCC from the rest in y axis
    TLine *l3 = new TLine(w, v1, w, v2);  //splits nu+e from the rest in x axis
    TLine *l4 = new TLine(w1, v, w2, v);  //splits nu+e from the rest in y axis

    //TCanvas *c = new TCanvas("c", "", 800, 600);
    cov->SetMinimum(-1.);
    cov->SetMaximum(1.);
    cov->SetStats(0);
    cov->Draw("colz");
    c->SaveAs(Form("%s.png",name));
    c->SetLogz();
    cov->Draw("colz");
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
    c->SaveAs(Form("%s_logz.png",name));
}
*/



const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu";

void ElepCov()
{

  // covariance matrix
  TFile * covfile = new TFile( Form("%s/uncertainties/flux_covmtr/total_covariance_DUNE_opt.root",data_path), "OLD" );
  TH2D * hcovmx = (TH2D*) covfile->Get( "total_covariance" );

  // only need the ND FHC part, which is the first 52 bins probably
  for( int x = 0; x < nbins; ++x ) {
    for( int y = 0; y < nbins; ++y ) {
      covmx[x][y] = hcovmx->GetBinContent( x+1, y+1 );
    }
  }

    // cholesky decomposition
  TDecompChol decomp( covmx );
  if( !decomp.Decompose() ) {
    printf( "Main covariance matrix failed Cholesky decomposition\n" );
    return;
  }
  const TMatrixD chol = decomp.GetU();

  TRandom3 * rand = new TRandom3(12345);

  // make random number matrix
  for( int i = 0; i < N; ++i ) {
    double mean = 0.;
    for( int j = 0; j < nbins; ++j ) {
      double val = rand->Gaus( 0., 1. );
      scales[i][j] = val;
      mean += val;
    }

    mean /= N;
    // force each column mean to be exactly 0, this just eliminates tiny statistical fluctuations in the mean weight
    for( int j = 0; j < nbins; ++j ) {
      double old = scales[i][j];
      scales[i][j] = old - mean;
    }
  }

  for( int i = 0; i < nbins; ++i ) { // columns
    for( int j = 0; j < nbins; ++j ) { // columns
      // compute column covariance
      double covar = 0.;
      for( int k = 0; k < N; ++k ) { // rows
        covar += (scales[k][i] * scales[k][j]); // means are 0 already by construction
      }
      scaleCovars[i][j] = covar/N;
    }
  }

  TDecompChol scaleDecomp( scaleCovars );
  if( !scaleDecomp.Decompose() ) printf( "Scale matrix didn't decompolse\n" );
  TMatrixD toInvert = scaleDecomp.GetU();
  TMatrixD inverse = toInvert.Invert();
  scales *= inverse;

  scales *= chol;


  double CCEdges[nbins_CC+1];
  CCEdges[0]=0., CCEdges[1]=0.3;
  for(int i=1; i<nbins_CC+1; i++)
  {
    if(i<38)               CCEdges[i+1] = CCEdges[i] + 0.1;
    else if(i>=38 && i<43) CCEdges[i+1] = CCEdges[i] + 0.2;
    else if(i>=43 && i<48) CCEdges[i+1] = CCEdges[i] + 0.4;
    else if(i>=48 && i<53) CCEdges[i+1] = CCEdges[i] + 0.8;
    else if(i>=53 && i<55) CCEdges[i+1] = CCEdges[i] + 1.5;
    else                   CCEdges[i+1] = CCEdges[i] + 2.0;
  }

  const double nueEdges[9] = {0., 0.3, 0.6, 0.92, 1.3, 1.75, 2.45, 3.9, 16.0};

  TH1D *m     = new TH1D("m","",nbins_CC,CCEdges);
  TH1D *e     = new TH1D("e","",nbins_CC,CCEdges);
  TH1D *m_nue = new TH1D("m_nue","",nbins_nue,nueEdges);
  TH1D *e_nue = new TH1D("e_nue","",nbins_nue,nueEdges);
  TH1D *nue   = new TH1D("nue","",nbins_nue,nueEdges);
  TH1D *CC_m_nom  = new TH1D("CC_m_nom","",nbins_CC,CCEdges);
  TH1D *CC_e_nom  = new TH1D("CC_e_nom","",nbins_CC,CCEdges);
  TH1D *nue_m_nom = new TH1D("nue_m_nom","",nbins_nue,nueEdges);
  TH1D *nue_e_nom = new TH1D("nue_e_nom","",nbins_nue,nueEdges);
  TH1D *nue_nom   = new TH1D("nue_nom","",nbins_nue,nueEdges);

  TH1D *tp_m[n_mu];
  TH1D *tp_e[n_e];
  TH1D *tp_m_nue[n_mu];
  TH1D *tp_e_nue[n_e];


  TFile *f     = new TFile(Form("%s/input_dfiles/CCtest_output.root",data_path), "READ");
  TFile *f_nue = new TFile(Form("%s/input_dfiles/nue_output.root",data_path), "READ");

  TH2D *CC_m  = (TH2D*)f->Get(Form("mElepRecoVsEv4_cov"));
  TH2D *CC_e  = (TH2D*)f->Get(Form("eElepRecoVsEv4_cov"));

  TH2D *nue_m = (TH2D*)f_nue->Get(Form("mElepRecoVsEv_cov"));
  TH2D *nue_e = (TH2D*)f_nue->Get(Form("eElepRecoVsEv_cov"));

  const int N_nucut = 5;
  double nucut[N_nucut] = {50., 5., 1., 0.5, 0.3};

  TH2D *hcv[N_nucut];
  TH2D *hfrcv[N_nucut];
  TH2D *hcr[N_nucut];


  //gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  //TCanvas *c = new TCanvas("c","",800,600);


  for(int inu=0; inu<N_nucut; inu++){

    CC_m = (TH2D*)f->Get(Form("mElepRecoVsEv%d_cov",inu));
    CC_e = (TH2D*)f->Get(Form("eElepRecoVsEv%d_cov",inu));


    for(int mb=0; mb<n_mu; mb++){
      tp_m[mb]     = (TH1D*)CC_m->ProjectionY(Form("m_bin%d",mb+1),mb+1,mb+1);
      tp_m_nue[mb] = (TH1D*)nue_m->ProjectionY(Form("m_bin_nue%d",mb+1),mb+1,mb+1);
    }
    for(int eb=0; eb<n_e; eb++){
      tp_e[eb]     = (TH1D*)CC_e->ProjectionY(Form("e_bin%d",eb+1),eb+1,eb+1);
      tp_e_nue[eb] = (TH1D*)nue_e->ProjectionY(Form("e_bin_nue%d",eb+1),eb+1,eb+1);
    }


    for( int u = 0; u < N; ++u ) {
      m->Reset();
      e->Reset();
      m_nue->Reset();
      e_nue->Reset();
      nue->Reset();
      CC_m_nom->Reset();
      CC_e_nom->Reset();
      nue_m_nom->Reset();
      nue_e_nom->Reset();
      nue_nom->Reset();
      
      for(int mb=0; mb<n_mu; mb++) {
        int fluxbin = mb+1;
        double evtwgt = scales[u][fluxbin];

        m->Add(tp_m[mb],1.+evtwgt);
        m_nue->Add(tp_m_nue[mb],1.+evtwgt);

        CC_m_nom->Add(tp_m[mb], 1.);
        nue_m_nom->Add(tp_m_nue[mb],1.);
      }

      for(int eb=0; eb<n_e; eb++) {
        int fluxbin = 38+eb+1;
        double evtwgt = scales[u][fluxbin];

        e->Add(tp_e[eb],1.+evtwgt);
        e_nue->Add(tp_e_nue[eb],1.+evtwgt);

        CC_e_nom->Add(tp_e[eb], 1.);
        nue_e_nom->Add(tp_e_nue[eb],1.);
      }
      
      nue->Add(m_nue);         nue->Add(e_nue);
      nue_nom->Add(nue_m_nom); nue_nom->Add(nue_e_nom);

      for(int i=0; i<nbins_CC; i++){
        E_m[u][i]     = m->GetBinContent(i+1);
        E_e[u][i]     = e->GetBinContent(i+1);
      }
      for(int i=0; i<nbins_nue; i++){
        E_nue[u][i]   = nue->GetBinContent(i+1);
      }
    }


    for( int i = 0; i < nbins_tot; ++i ) { // columns
      for( int j = 0; j < nbins_tot; ++j ) { // columns
        // compute column covariance
        double covar_mm = 0.;
        double covar_ee = 0.;
        double covar_nuenue = 0.;
        double covar_me = 0.;
        double covar_mnue = 0.;
        double covar_em = 0.;
        double covar_enue = 0.;
        double covar_nuem = 0.;
        double covar_nuee = 0.;

        double var_m_i = 0.;
        double var_m_j = 0.;
        double var_e_i = 0.;
        double var_e_j = 0.;
        double var_nue_i = 0.;
        double var_nue_j = 0.;

        for( int k = 0; k < N; ++k ) { // rows

          if(i<nbins_CC){
            if(j<nbins_CC)                  covar_mm   += (E_m[k][i] - CC_m_nom->GetBinContent(i+1)) * (E_m[k][j]              - CC_m_nom->GetBinContent(j+1));
            if(j>=nbins_CC && j<2*nbins_CC) covar_me   += (E_m[k][i] - CC_m_nom->GetBinContent(i+1)) * (E_e[k][j-nbins_CC]     - CC_e_nom->GetBinContent(j-nbins_CC+1));
            if(j>=2*nbins_CC)               covar_mnue += (E_m[k][i] - CC_m_nom->GetBinContent(i+1)) * (E_nue[k][j-2*nbins_CC] - nue_nom->GetBinContent(j-2*nbins_CC+1));}

          if(i>=nbins_CC && i<2*nbins_CC){
            if(j<nbins_CC)                  covar_em     += (E_e[k][i-nbins_CC] - CC_e_nom->GetBinContent(i-nbins_CC+1)) * (E_m[k][j]              - CC_m_nom->GetBinContent(j+1));
            if(j>=nbins_CC && j<2*nbins_CC) covar_ee     += (E_e[k][i-nbins_CC] - CC_e_nom->GetBinContent(i-nbins_CC+1)) * (E_e[k][j-nbins_CC]     - CC_e_nom->GetBinContent(j-nbins_CC+1));
            if(j>=2*nbins_CC)               covar_enue   += (E_e[k][i-nbins_CC] - CC_e_nom->GetBinContent(i-nbins_CC+1)) * (E_nue[k][j-2*nbins_CC] - nue_nom->GetBinContent(j-2*nbins_CC+1));}

          if(i>=2*nbins_CC){
            if(j<nbins_CC)                  covar_nuem   += (E_nue[k][i-2*nbins_CC] - nue_nom->GetBinContent(i-2*nbins_CC+1))  * (E_m[k][j]              - CC_m_nom->GetBinContent(j+1));
            if(j>=nbins_CC && j<2*nbins_CC) covar_nuee   += (E_nue[k][i-2*nbins_CC] - nue_nom->GetBinContent(i-2*nbins_CC+1))  * (E_e[k][j-nbins_CC]     - CC_e_nom->GetBinContent(j-nbins_CC+1));
            if(j>=2*nbins_CC)               covar_nuenue += (E_nue[k][i-2*nbins_CC] - nue_nom->GetBinContent(i-2*nbins_CC+1))  * (E_nue[k][j-2*nbins_CC] - nue_nom->GetBinContent(j-2*nbins_CC+1));}

          if(i<nbins_CC)                  var_m_i   += (E_m[k][i]              - CC_m_nom->GetBinContent(i+1))           * (E_m[k][i]              - CC_m_nom->GetBinContent(i+1));
          if(i>=nbins_CC && i<2*nbins_CC) var_e_i   += (E_e[k][i-nbins_CC]     - CC_e_nom->GetBinContent(i-nbins_CC+1))  * (E_e[k][i-nbins_CC]     - CC_e_nom->GetBinContent(i-nbins_CC+1));
          if(i>=2*nbins_CC)               var_nue_i += (E_nue[k][i-2*nbins_CC] - nue_nom->GetBinContent(i-2*nbins_CC+1)) * (E_nue[k][i-2*nbins_CC] - nue_nom->GetBinContent(i-2*nbins_CC+1));

          if(j<nbins_CC)                  var_m_j   += (E_m[k][j]              - CC_m_nom->GetBinContent(j+1))           * (E_m[k][j]              - CC_m_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) var_e_j   += (E_e[k][j-nbins_CC]     - CC_e_nom->GetBinContent(j-nbins_CC+1))  * (E_e[k][j-nbins_CC]     - CC_e_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               var_nue_j += (E_nue[k][j-2*nbins_CC] - nue_nom->GetBinContent(j-2*nbins_CC+1)) * (E_nue[k][j-2*nbins_CC] - nue_nom->GetBinContent(j-2*nbins_CC+1));
        }


        //This is covariance matrix
        if(i<nbins_CC){
          if(j<nbins_CC)                  ECovars_mm[i][j]              = covar_mm  /N;
          if(j>=nbins_CC && j<2*nbins_CC) ECovars_me[i][j-nbins_CC]     = covar_me  /N;
          if(j>=2*nbins_CC)               ECovars_mnue[i][j-2*nbins_CC] = covar_mnue/N;}

        if(i>=nbins_CC && i<2*nbins_CC){
          if(j<nbins_CC)                  ECovars_em[i-nbins_CC][j]              = covar_em  /N;
          if(j>=nbins_CC && j<2*nbins_CC) ECovars_ee[i-nbins_CC][j-nbins_CC]     = covar_ee  /N;
          if(j>=2*nbins_CC)               ECovars_enue[i-nbins_CC][j-2*nbins_CC] = covar_enue/N;}

        if(i>=2*nbins_CC){
          if(j<nbins_CC)                  ECovars_nuem[i-2*nbins_CC][j]              = covar_nuem  /N;
          if(j>=nbins_CC && j<2*nbins_CC) ECovars_nuee[i-2*nbins_CC][j-nbins_CC]     = covar_nuee  /N;
          if(j>=2*nbins_CC)               ECovars_nuenue[i-2*nbins_CC][j-2*nbins_CC] = covar_nuenue/N;}

        
        //This is fractional covariance matrix
        if(i<nbins_CC){
          if(j<nbins_CC)                  frECovars_mm[i][j]              = covar_mm  /(N * CC_m_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) frECovars_me[i][j-nbins_CC]     = covar_me  /(N * CC_m_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               frECovars_mnue[i][j-2*nbins_CC] = covar_mnue/(N * CC_m_nom->GetBinContent(i+1)  * nue_nom->GetBinContent(j-2*nbins_CC+1));}

        if(i>=nbins_CC && i<2*nbins_CC){
          if(j<nbins_CC)                  frECovars_em[i-nbins_CC][j]              = covar_em  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1)  * CC_m_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) frECovars_ee[i-nbins_CC][j-nbins_CC]     = covar_ee  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1)  * CC_e_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               frECovars_enue[i-nbins_CC][j-2*nbins_CC] = covar_enue/(N * CC_e_nom->GetBinContent(i-nbins_CC+1)  * nue_nom->GetBinContent(j-2*nbins_CC+1));}

        if(i>=2*nbins_CC){
          if(j<nbins_CC)                  frECovars_nuem[i-2*nbins_CC][j]              = covar_nuem  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * CC_m_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) frECovars_nuee[i-2*nbins_CC][j-nbins_CC]     = covar_nuee  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * CC_e_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               frECovars_nuenue[i-2*nbins_CC][j-2*nbins_CC] = covar_nuenue/(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * nue_nom->GetBinContent(j-2*nbins_CC+1));}



        //This is correlation matrix
        if(i<nbins_CC){
          if(j<nbins_CC)                  ECorrel_mm[i][j]              = covar_mm  /sqrt(var_m_i * var_m_j);
          if(j>=nbins_CC && j<2*nbins_CC) ECorrel_me[i][j-nbins_CC]     = covar_me  /sqrt(var_m_i * var_e_j);
          if(j>=2*nbins_CC)               ECorrel_mnue[i][j-2*nbins_CC] = covar_mnue/sqrt(var_m_i * var_nue_j);}

        if(i>=nbins_CC && i<2*nbins_CC){
          if(j<nbins_CC)                  ECorrel_em[i-nbins_CC][j]              = covar_em  /sqrt(var_e_i * var_m_j);
          if(j>=nbins_CC && j<2*nbins_CC) ECorrel_ee[i-nbins_CC][j-nbins_CC]     = covar_ee  /sqrt(var_e_i * var_e_j);
          if(j>=2*nbins_CC)               ECorrel_enue[i-nbins_CC][j-2*nbins_CC] = covar_enue/sqrt(var_e_i * var_nue_j);}

        if(i>=2*nbins_CC){
          if(j<nbins_CC)                  ECorrel_nuem[i-2*nbins_CC][j]              = covar_nuem  /sqrt(var_nue_i * var_m_j);
          if(j>=nbins_CC && j<2*nbins_CC) ECorrel_nuee[i-2*nbins_CC][j-nbins_CC]     = covar_nuee  /sqrt(var_nue_i * var_e_j);
          if(j>=2*nbins_CC)               ECorrel_nuenue[i-2*nbins_CC][j-2*nbins_CC] = covar_nuenue/sqrt(var_nue_i * var_nue_j);}
        
      }
    }

    for(int i = 0; i < nbins_tot; i++) {
      for(int j = 0; j < nbins_tot; j++) {

        if(i<nbins_CC){
          if(j<nbins_CC)                  ECovars[i][j] = ECovars_mm[i][j];
          if(j>=nbins_CC && j<2*nbins_CC) ECovars[i][j] = ECovars_me[i][j-nbins_CC];
          if(j>=2*nbins_CC)               ECovars[i][j] = ECovars_mnue[i][j-2*nbins_CC];}

        if(i>=nbins_CC && i<2*nbins_CC){
          if(j<nbins_CC)                  ECovars[i][j] = ECovars_em[i-nbins_CC][j];
          if(j>=nbins_CC && j<2*nbins_CC) ECovars[i][j] = ECovars_ee[i-nbins_CC][j-nbins_CC];
          if(j>=2*nbins_CC)               ECovars[i][j] = ECovars_enue[i-nbins_CC][j-2*nbins_CC];}

        if(i>=2*nbins_CC){
          if(j<nbins_CC)                  ECovars[i][j] = ECovars_nuem[i-2*nbins_CC][j];
          if(j>=nbins_CC && j<2*nbins_CC) ECovars[i][j] = ECovars_nuee[i-2*nbins_CC][j-nbins_CC];
          if(j>=2*nbins_CC)               ECovars[i][j] = ECovars_nuenue[i-2*nbins_CC][j-2*nbins_CC];}

        


        if(i<nbins_CC){
          if(j<nbins_CC)                  frECovars[i][j] = frECovars_mm[i][j];
          if(j>=nbins_CC && j<2*nbins_CC) frECovars[i][j] = frECovars_me[i][j-nbins_CC];
          if(j>=2*nbins_CC)               frECovars[i][j] = frECovars_mnue[i][j-2*nbins_CC];}

        if(i>=nbins_CC && i<2*nbins_CC){
          if(j<nbins_CC)                  frECovars[i][j] = frECovars_em[i-nbins_CC][j];
          if(j>=nbins_CC && j<2*nbins_CC) frECovars[i][j] = frECovars_ee[i-nbins_CC][j-nbins_CC];
          if(j>=2*nbins_CC)               frECovars[i][j] = frECovars_enue[i-nbins_CC][j-2*nbins_CC];}

        if(i>=2*nbins_CC){
          if(j<nbins_CC)                  frECovars[i][j] = frECovars_nuem[i-2*nbins_CC][j];
          if(j>=nbins_CC && j<2*nbins_CC) frECovars[i][j] = frECovars_nuee[i-2*nbins_CC][j-nbins_CC];
          if(j>=2*nbins_CC)               frECovars[i][j] = frECovars_nuenue[i-2*nbins_CC][j-2*nbins_CC];}




        if(i<nbins_CC){
          if(j<nbins_CC)                  ECorrel[i][j] = ECorrel_mm[i][j];
          if(j>=nbins_CC && j<2*nbins_CC) ECorrel[i][j] = ECorrel_me[i][j-nbins_CC];
          if(j>=2*nbins_CC)               ECorrel[i][j] = ECorrel_mnue[i][j-2*nbins_CC];}

        if(i>=nbins_CC && i<2*nbins_CC){
          if(j<nbins_CC)                  ECorrel[i][j] = ECorrel_em[i-nbins_CC][j];
          if(j>=nbins_CC && j<2*nbins_CC) ECorrel[i][j] = ECorrel_ee[i-nbins_CC][j-nbins_CC];
          if(j>=2*nbins_CC)               ECorrel[i][j] = ECorrel_enue[i-nbins_CC][j-2*nbins_CC];}

        if(i>=2*nbins_CC){
          if(j<nbins_CC)                  ECorrel[i][j] = ECorrel_nuem[i-2*nbins_CC][j];
          if(j>=nbins_CC && j<2*nbins_CC) ECorrel[i][j] = ECorrel_nuee[i-2*nbins_CC][j-nbins_CC];
          if(j>=2*nbins_CC)               ECorrel[i][j] = ECorrel_nuenue[i-2*nbins_CC][j-2*nbins_CC];}
        
      }
    }

    hcv[inu] = new TH2D(Form("hcv%d",inu), Form("Total Flux Covariance Matrix (#nu < %.2f GeV)", nucut[inu]),nbins_tot,0,nbins_tot,nbins_tot,0,nbins_tot);
    hfrcv[inu] = new TH2D(Form("hfrcv%d",inu), Form("Fractional Flux Covariance Matrix (#nu < %.2f GeV)", nucut[inu]),nbins_tot,0,nbins_tot,nbins_tot,0,nbins_tot);
    hcr[inu] = new TH2D(Form("hcr%d",inu), Form("Total Flux Correlation Matrix (#nu < %.2f GeV)", nucut[inu]),nbins_tot,0,nbins_tot,nbins_tot,0,nbins_tot);

    for(int i=0; i<nbins_tot; i++) {
      for(int j=0; j<nbins_tot; j++) {
        hcv[inu]->SetBinContent(i+1, j+1, ECovars[i][j]);
        hfrcv[inu]->SetBinContent(i+1, j+1, frECovars[i][j]);
        hcr[inu]->SetBinContent(i+1, j+1, ECorrel[i][j]);
      }
    }
    /*
    plotCov(hcv[inu], Form("%s/uncertainties/flux_covmtr/flux_totCovmtr%d",data_path,inu), c);
    plotfrCov(hfrcv[inu], Form("%s/uncertainties/flux_covmtr/flux_frCovmtr%d",data_path,inu), c);
    plotCor(hcr[inu], Form("%s/uncertainties/flux_covmtr/flux_totCormtr%d",data_path,inu), c);
    */
  }

  //c->Close();


  TFile *out = new TFile(Form("%s/uncertainties/flux_covmtr/flux_covmtr_TEST.root", data_path),"RECREATE");
  for(int inu=0; inu<N_nucut; inu++){
    hcv[inu]->Write();
    hfrcv[inu]->Write();
    hcr[inu]->Write();
  }
  out->Close();


}

