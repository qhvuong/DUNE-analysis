static const int N = 10000; // number of universes
static const int nbins = 52;
static const int nbins_CC = 100;
static const int nbins_nue = 50;
static const int nbins_tot = 2*nbins_CC + nbins_nue;

const int n_mu = 19; // number of muon bins in covariance
const int n_e = 7; // number of electron bins in covariance

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
TMatrixD ECovars_m_nue( nbins_CC, nbins_nue );
TMatrixD ECovars_e_nue( nbins_CC, nbins_nue );
TMatrixD ECovars_nuenue  ( nbins_nue, nbins_nue );

TMatrixD ECovars_me  ( nbins_CC, nbins_CC );
TMatrixD ECovars_mnue( nbins_CC, nbins_nue );
TMatrixD ECovars_em  ( nbins_CC, nbins_CC );
TMatrixD ECovars_enue( nbins_CC, nbins_nue );
TMatrixD ECovars_nuem( nbins_nue, nbins_CC );
TMatrixD ECovars_nuee( nbins_nue, nbins_CC );

TMatrixD ECorrel_mm ( nbins_CC, nbins_CC );
TMatrixD ECorrel_ee ( nbins_CC, nbins_CC );
TMatrixD ECorrel_nuenue( nbins_nue, nbins_nue );

TMatrixD ECorrel_me  ( nbins_CC, nbins_CC );
TMatrixD ECorrel_mnue( nbins_CC, nbins_nue );
TMatrixD ECorrel_em  ( nbins_CC, nbins_CC );
TMatrixD ECorrel_enue( nbins_CC, nbins_nue );
TMatrixD ECorrel_nuem( nbins_nue, nbins_CC );
TMatrixD ECorrel_nuee( nbins_nue, nbins_CC );

TMatrixD ECovars( nbins_tot, nbins_tot );
TMatrixD ECorrel( nbins_tot, nbins_tot );

void test()
{

  // covariance matrix
  TFile * covfile = new TFile( "total_covariance_DUNE_opt.root", "OLD" );
  TH2D * hcovmx = (TH2D*) covfile->Get( "total_covariance" );

  // only need the ND FHC part, which is the first 52 bins probably
  for( int x = 0; x < nbins; ++x ) {
    for( int y = 0; y < nbins; ++y ) {
      covmx[x][y] = hcovmx->GetBinContent( x+1, y+1 );

      //std::cout << x << "\t" << y << "\t" << covmx[x][y] << "\n";
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
    //if( i>0 && i%50==0 ) std::cout << "Progress: " << i/N*100. << "..." << "\n";
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

      //std::cout << i << "\t" << j << "\t" << scaleCovars[i][j] << "\n";
    }
  }

  TDecompChol scaleDecomp( scaleCovars );
  if( !scaleDecomp.Decompose() ) printf( "Scale matrix didn't decompolse\n" );
  TMatrixD toInvert = scaleDecomp.GetU();
  TMatrixD inverse = toInvert.Invert();
  scales *= inverse;

  scales *= chol;

  char name[20] = "ElepReco";
  int cutNu = 3;

  //for(para = 1; para <3; para++) {
  TFile *f     = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/data/lownu/nue_output.root");
  TH2D *CC_m  = (TH2D*)f->Get(Form("m_h%sVsEv%d_cov",name,cutNu));
  TH2D *CC_e  = (TH2D*)f->Get(Form("e_h%sVsEv%d_cov",name,cutNu));
  TH2D *nue_m = (TH2D*)f_nue->Get(Form("m_h%sVsEv0_cov",name));
  TH2D *nue_e = (TH2D*)f_nue->Get(Form("e_h%sVsEv0_cov",name));

  TH1D *tp_m[n_mu];
  TH1D *tp_e[n_e];
  TH1D *tp_m_nue[n_mu];
  TH1D *tp_e_nue[n_e];

  for(int mb=0; mb<n_mu; mb++){
    tp_m[mb]     = (TH1D*)CC_m->ProjectionY(Form("m_bin%d",mb+1),mb+1,mb+1);
    tp_m_nue[mb] = (TH1D*)nue_m->ProjectionY(Form("m_bin_nue%d",mb+1),mb+1,mb+1);
  }
  for(int eb=0; eb<n_e; eb++){
    tp_e[eb]     = (TH1D*)CC_e->ProjectionY(Form("e_bin%d",eb+1),eb+1,eb+1);
    tp_e_nue[eb] = (TH1D*)nue_e->ProjectionY(Form("e_bin_nue%d",eb+1),eb+1,eb+1);
  }

  TH1D *m     = new TH1D("m","",nbins_CC,0,16);
  TH1D *e     = new TH1D("e","",nbins_CC,0,16);
  TH1D *m_nue = new TH1D("m_nue","",nbins_nue,0,16);
  TH1D *e_nue = new TH1D("e_nue","",nbins_nue,0,16);
  TH1D *nue   = new TH1D("nue","",nbins_nue,0,16);
  TH1D *CC_m_nom  = new TH1D("CC_m_nom","",nbins_CC,0,16);
  TH1D *CC_e_nom  = new TH1D("CC_e_nom","",nbins_CC,0,16);
  TH1D *nue_m_nom = new TH1D("nue_m_nom","",nbins_nue,0,16);
  TH1D *nue_e_nom = new TH1D("nue_e_nom","",nbins_nue,0,16);
  TH1D *nue_nom   = new TH1D("nue_nom","",nbins_nue,0,16);

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

/*
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
*/

//This is fractional covariance matrix

      if(i<nbins_CC){
        if(j<nbins_CC)                  ECovars_mm[i][j]              = covar_mm  /(N * CC_m_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_me[i][j-nbins_CC]     = covar_me  /(N * CC_m_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               ECovars_mnue[i][j-2*nbins_CC] = covar_mnue/(N * CC_m_nom->GetBinContent(i+1)  * nue_nom->GetBinContent(j-2*nbins_CC+1));}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  ECovars_em[i-nbins_CC][j]              = covar_em  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1)  * CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_ee[i-nbins_CC][j-nbins_CC]     = covar_ee  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1)  * CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               ECovars_enue[i-nbins_CC][j-2*nbins_CC] = covar_enue/(N * CC_e_nom->GetBinContent(i-nbins_CC+1)  * nue_nom->GetBinContent(j-2*nbins_CC+1));}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  ECovars_nuem[i-2*nbins_CC][j]              = covar_nuem  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_nuee[i-2*nbins_CC][j-nbins_CC]     = covar_nuee  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               ECovars_nuenue[i-2*nbins_CC][j-2*nbins_CC] = covar_nuenue/(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * nue_nom->GetBinContent(j-2*nbins_CC+1));}


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

  TH2D *hcv = new TH2D("hcv","",nbins_tot,0,nbins_tot,nbins_tot,0,nbins_tot);
  TH2D *hcr = new TH2D("hcr","",nbins_tot,0,nbins_tot,nbins_tot,0,nbins_tot);
  for(int i=0; i<nbins_tot; i++) {
    for(int j=0; j<nbins_tot; j++) {
      hcv->SetBinContent(i+1, j+1, ECovars[i][j]);
      hcr->SetBinContent(i+1, j+1, ECorrel[i][j]);
      //if(i>nbins_CC || j>nbins_CC) std::cout << ECorrel[i][j] << "\t";  
    }
  }
 
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();


  hcv->SetStats(0);

  hcr->SetStats(0);

  TCanvas *c = new TCanvas("c","",700,700);
  hcv->Draw("colz");
  c->SaveAs(Form("flux_frCovmtr%d_%d.png",cutNu,N));
  

  const Int_t Number = 3;
  Double_t Red[Number]    = { 0.00, 1.00, 1.00};
  Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  Double_t Blue[Number]   = { 1.00, 1.00, 0.00};
  Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  Int_t nb=50;
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

  TCanvas *c1 = new TCanvas("c1","",700,700);
  hcr->GetZaxis()->SetRangeUser(-1., 1.);
  hcr->Draw("colz");
  c1->SaveAs(Form("flux_Cormtr%d_%d.png",cutNu,N));

  //gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  TFile *out = new TFile(Form("flux_frCovmtr%d_%d.root",cutNu,N),"RECREATE");
  hcv->Write();
  hcr->Write();
  out->Close();

}

