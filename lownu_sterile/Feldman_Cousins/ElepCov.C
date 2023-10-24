static const int N = 500; 
static const int nbins_CC = 100;
static const int nbins_nue = 50;
static const int nbins = 2*nbins_CC + nbins_nue;



TMatrixD covmx( nbins, nbins );
TMatrixD sysmx( nbins, nbins );
TMatrixD statmx( nbins, nbins );
TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, nbins );
TMatrixD scalesChol( N, nbins );
TMatrixD scalesSVD( N, nbins );

TMatrixD E_m(N, nbins_CC);
TMatrixD E_e(N, nbins_CC);
TMatrixD E_nue(N, nbins_nue);
  
TMatrixD ECovars_mm    ( nbins_CC, nbins_CC );
TMatrixD ECovars_ee    ( nbins_CC, nbins_CC );
TMatrixD ECovars_nuenue  ( nbins_nue, nbins_nue );

TMatrixD ECovars_me  ( nbins_CC, nbins_CC );
TMatrixD ECovars_mnue( nbins_CC, nbins_nue );
TMatrixD ECovars_em  ( nbins_CC, nbins_CC );
TMatrixD ECovars_enue( nbins_CC, nbins_nue );
TMatrixD ECovars_nuem( nbins_nue, nbins_CC );
TMatrixD ECovars_nuee( nbins_nue, nbins_CC );

TMatrixD ECovars( nbins, nbins );

void ElepCov_stats()
{
  int cutNu = 3;
	
  const char name[20] = "wgt_MaCCQE";

  TFile *f     = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root", "READ");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/data/lownu/nue_output.root", "READ");
  TFile *f_sys = new TFile(Form("../xS_covmtr/%s_sigmtr%d.root",name,cutNu), "READ");
  TFile *f_fl  = new TFile(Form("../flux_covmtr/flux_covmtr%d_10000.root",cutNu),"READ");

  TH2D *hsys   = (TH2D*)f_sys->Get( "hcv" );
  TH2D *hfl    = (TH2D*)f_fl->Get("hcv");

  TH1D *CC_m_nom = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma3",name,cutNu));
  TH1D *CC_e_nom = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma3",name,cutNu));
  TH1D *nue_nom  = (TH1D*)f_nue->Get("hElep0");
  

  // only need the ND FHC part, which is the first 52 bins probably
  for( int x = 0; x < nbins; ++x ) {
    for( int y = 0; y < nbins; ++y ) {
      //sysmx[x][y]  = hfl->GetBinContent(x+1, y+1) + hsys->GetBinContent(x+1, y+1);
      sysmx[x][y]  = hfl->GetBinContent(x+1, y+1);
      statmx[x][y] = 0.;

      if(x==y){
        if(x<nbins_CC)                  statmx[x][y] = CC_m_nom->GetBinContent( x+1 );
        if(x>=nbins_CC && x<2*nbins_CC) statmx[x][y] = CC_e_nom->GetBinContent( x-nbins_CC+1 );
        if(x>=2*nbins_CC)               statmx[x][y] = nue_nom->GetBinContent( x-2*nbins_CC+1 );
      }

    }
  }

/*
  for( int i = 0; i < nbins; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      sysmx[i][j]                 = (hfl->GetBinContent(i+1, j+1)         + hsys->GetBinContent(i+1, j+1))                 * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1);
      sysmx[i][j+nbins]           = (hfl->GetBinContent(i+1, j+nbins+1)   + hsys->GetBinContent(i+1, j+nbins+1))           * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1);
      sysmx[i][j+2*nbins]         = (hfl->GetBinContent(i+1, j+2*nbins+1) + hsys->GetBinContent(i+1, j+2+nbins+1))         * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1);

      sysmx[i+nbins][j]           = (hfl->GetBinContent(i+nbins+1, j+1)         + hsys->GetBinContent(i+nbins+1, j+1))           * CC_e_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1);
      sysmx[i+nbins][j+nbins]     = (hfl->GetBinContent(i+nbins+1, j+nbins+1)   + hsys->GetBinContent(i+nbins+1, j+nbins+1))     * CC_e_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1);
      sysmx[i+nbins][j+2*nbins]   = (hfl->GetBinContent(i+nbins+1, j+2*nbins+1) + hsys->GetBinContent(i+nbins+1, j+2*nbins+1))   * CC_e_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1);

      sysmx[i+2*nbins][j]         = (hfl->GetBinContent(i+2*nbins+1, j+1)         + hsys->GetBinContent(i+2*nbins+1, j+1))         * nue_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1);
      sysmx[i+2*nbins][j+nbins]   = (hfl->GetBinContent(i+2*nbins+1, j+nbins+1)   + hsys->GetBinContent(i+2*nbins+1, j+nbins+1))   * nue_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j+1);
      sysmx[i+2*nbins][j+2*nbins] = (hfl->GetBinContent(i+2*nbins+1, j+2*nbins+1) + hsys->GetBinContent(i+2*nbins+1, j+2*nbins+1)) * nue_nom->GetBinContent(i+1)  * nue_nom->GetBinContent(j+1);
  }
  } 
*/
  

  TH2D *hcovmx = new TH2D("hcovmx","",nbins,0,nbins,nbins,0,nbins);
  TH2D *hsysmx = new TH2D("hsysmx","",nbins,0,nbins,nbins,0,nbins);
  TH2D *hstatmx = new TH2D("hstatmx","",nbins,0,nbins,nbins,0,nbins);
  for( int i = 0; i < nbins; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      //covmx[i][j] = statmx[i][j] + sysmx[i][j];
      covmx[i][j] = statmx[i][j];
      hsysmx->SetBinContent(i+1, j+1, sysmx[i][j]);
      hstatmx->SetBinContent(i+1, j+1, statmx[i][j]);
      hcovmx->SetBinContent(i+1, j+1, covmx[i][j]);
  }
  } 


/*
  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  int MaxAttempts = 1000;
  int iAttempt = 0;
  bool CanDecomp = false;
  TDecompChol chdcmp;

  for (iAttempt=0;iAttempt<MaxAttempts;iAttempt++) {
    chdcmp = TDecompChol(covmx);
    if (chdcmp.Decompose()) {
      CanDecomp = true;
      break;
    } else {
      for (int iVar=0;iVar<3*nbins;iVar++) {
        (covmx)(iVar,iVar) += pow(10,-9);
      }
    }
  }

  if (!CanDecomp) {
    std::cerr << "Tried " << MaxAttempts << " times to shift diagonal but still can not decompose the matrix" << std::endl;
    std::cerr << "This indicates that something is wrong with the input matrix" << std::endl;
    throw;
  }

  std::cout << "Had to shift diagonal " << iAttempt << " time(s) to allow the covariance matrix to be decomposed" << std::endl;
*/


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
  //std::cout << scales[1][1] << "\n";


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
  scales = scales*inverse;
  scales *= chol;



  TH1D *m     = new TH1D("m","",nbins_CC,0,16);
  TH1D *e     = new TH1D("e","",nbins_CC,0,16);
  TH1D *nue   = new TH1D("nue","",nbins_nue,0,16);

  
  for( int u = 0; u < N; ++u ) {
    //if(u%100==0) std::cout << u << "\n";

    m->Reset();
    e->Reset();
    nue->Reset();

    for(int b=0; b<nbins_CC; b++) {
      int fluxbin = b;
      double evtwgt = scales[u][fluxbin];
      double bincontent = CC_m_nom->GetBinContent(b+1);

      m->SetBinContent(b+1, bincontent*(1.+evtwgt));
    }
    for(int b=0; b<nbins_CC; b++) {
      int fluxbin = b+nbins_CC;
      double evtwgt = scales[u][fluxbin];
      double bincontent = CC_e_nom->GetBinContent(b+1);

      e->SetBinContent(b+1, bincontent*(1.+evtwgt));
    }
    for(int b=0; b<nbins_nue; b++) {
      int fluxbin = b+2*nbins_CC;
      double evtwgt = scales[u][fluxbin];
      double bincontent = nue_nom->GetBinContent(b+1);

      nue->SetBinContent(b+1, bincontent*(1.+evtwgt));
    }


    for(int i=0; i<nbins_CC; i++){
      E_m[u][i]     = m->GetBinContent(i+1);
      E_e[u][i]     = e->GetBinContent(i+1);
    }
    for(int i=0; i<nbins_nue; i++){
      E_nue[u][i]   = nue->GetBinContent(i+1);
    }
   }


  for( int i = 0; i < nbins; ++i ) { // columns
    for( int j = 0; j < nbins; ++j ) { // columns
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

      }

      if(i<nbins_CC){
        if(j<nbins_CC)                  ECovars_mm[i][j]              = covar_mm  /(N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_me[i][j-nbins_CC]     = covar_me  /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               ECovars_mnue[i][j-2*nbins_CC] = covar_mnue/(N * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}

      if(i>=nbins_CC && i<2*nbins_CC){ 
        if(j<nbins_CC)                  ECovars_em[i-nbins_CC][j]              = covar_em  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1) * CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_ee[i-nbins_CC][j-nbins_CC]     = covar_ee  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1) * CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               ECovars_enue[i-nbins_CC][j-2*nbins_CC] = covar_enue/(N * CC_e_nom->GetBinContent(i-nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  ECovars_nuem[i-2*nbins_CC][j]              = covar_nuem  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_nuee[i-2*nbins_CC][j-nbins_CC]     = covar_nuee  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               ECovars_nuenue[i-2*nbins_CC][j-2*nbins_CC] = covar_nuenue/(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * nue_nom->GetBinContent(j-2*nbins_CC+1));}

    }
  }

  for(int i = 0; i < nbins; i++) {
    for(int j = 0; j < nbins; j++) {
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

    }
  }

  TH2D *hcv = new TH2D("cv","",nbins,0,nbins,nbins,0,nbins);
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      hcv->SetBinContent(i+1, j+1, ECovars[i][j]);
    }
  }
 
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

/*
  hcovmx->SetStats(0);
  hsysmx->SetStats(0);
  hstatmx->SetStats(0);
  hcovmx->SetTitle("total covmtr");
  hsysmx->SetTitle("systematic covmtr");
  hstatmx->SetTitle("statistical covmtr");

  TCanvas *c = new TCanvas("c","",1200,600);
  c->Divide(3,2);
  c->cd(1);
  gPad->SetLogz();
  hsysmx->Draw("colz");
  c->cd(2);
  gPad->SetLogz();
  hstatmx->Draw("colz");
  c->cd(3);
  gPad->SetLogz();
  hcovmx->Draw("colz");
  c->cd(4);
  hsysmx->Draw("colz");
  c->cd(5);
  hstatmx->Draw("colz");
  c->cd(6);
  hcovmx->Draw("colz");
  c->SaveAs(Form("covmtr%d_3sig.png",cutNu));
*/

  hcovmx->SetStats(0);
  hcovmx->SetTitle("stat mtr");
  hcv->SetStats(0);
  hcv->SetTitle(Form("FC stat %d universes",N));

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetLogz();
  hcovmx->Draw("colz");
  c1->cd(2);
  gPad->SetLogz();
  hcv->Draw("colz");
  c1->cd(3);
  hcovmx->Draw("colz");
  c1->cd(4);
  hcv->Draw("colz");
  c1->SaveAs(Form("FC_stat_%d_%d.png",cutNu,N));

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();


  TFile *out = new TFile(Form("FC_stat_%d_%d.root",cutNu,N),"RECREATE");
  hcv->Write();
  out->Close();

  //} 
}

