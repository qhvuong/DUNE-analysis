static const int N = 500; // number of universes
static const int nbins = 100;


TMatrixD covmx( 3*nbins, 3*nbins );
TMatrixD scaleCovars( 3*nbins, 3*nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, 3*nbins );

TMatrixD E_m(N, nbins);
TMatrixD E_e(N, nbins);
TMatrixD E_nue(N, nbins);
  
TMatrixD ECovars_mm    ( nbins, nbins );
TMatrixD ECovars_ee    ( nbins, nbins );
TMatrixD ECovars_nuenue  ( nbins, nbins );

TMatrixD ECovars_me  ( nbins, nbins );
TMatrixD ECovars_mnue( nbins, nbins );
TMatrixD ECovars_em  ( nbins, nbins );
TMatrixD ECovars_enue( nbins, nbins );
TMatrixD ECovars_nuem( nbins, nbins );
TMatrixD ECovars_nuee( nbins, nbins );

TMatrixD ECovars( 3*nbins, 3*nbins );

void ElepCov_stats()
{
  int cutNu, cutEv;
  const char name[20] = "wgt_MaCCQE";

  //for(para = 1; para <3; para++) {
  TFile *f     = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/data/lownu/nue_output.root");
  //for(cutNu = 0; cutNu < 4; cutNu++) {
  //if(cutNu != 0 && cutNu != 3 ) continue;
  cutNu = 3;		// cut nu<0.3GeV
  //std::cout << cutNu << "\n";
  //for(cutEv = 0; cutEv < 3; cutEv++) {
  cutEv = 2;		// cut Etheta<0.5MeV

  TH1D *CC_m_nom = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma3",name,cutNu));
  TH1D *CC_e_nom = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma3",name,cutNu));
  TH1D *nue_nom  = (TH1D*)f_nue->Get(Form("hElep%d",cutEv));


  // only need the ND FHC part, which is the first 52 bins probably
  for( int x = 0; x < nbins; ++x ) {
    for( int y = 0; y < nbins; ++y ) {
      if(x==y){
        covmx[x][y]                 = CC_m_nom->GetBinContent( x+1, y+1 );
        covmx[x+nbins][y+nbins]     = CC_e_nom->GetBinContent( x+1, y+1 );
        covmx[x+2*nbins][y+2*nbins] = nue_nom->GetBinContent( x+1, y+1 );
        std::cout << covmx[x][y] << "\t" << covmx[x+nbins][y+nbins] << "\t" << covmx[x+2*nbins][y+2*nbins] << "\n";
      }

      else{
        covmx[x][y]                     = 0.;
        covmx[x+nbins][y+nbins]         = 0.;
        covmx[x+2*nbins][y+2*nbins]     = 0.;}
    }
  }

  TH2D *hcovmx = new TH2D("hcovmx","",300,0,300,300,0,300);
  for( int i = 0; i < nbins; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
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
    for( int j = 0; j < 3*nbins; ++j ) {
      double val = rand->Gaus( 0., 1. );
      scales[i][j] = val;
      mean += val;
    }

    mean /= N;
    // force each column mean to be exactly 0, this just eliminates tiny statistical fluctuations in the mean weight
    for( int j = 0; j < 3*nbins; ++j ) {
      double old = scales[i][j];
      scales[i][j] = old - mean;
    }
  }
  std::cout << scales[1][1] << "\n";


  for( int i = 0; i < 3*nbins; ++i ) { // columns
    for( int j = 0; j < 3*nbins; ++j ) { // columns

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
  std::cout << scales[1][1] << "\n";

/*
  TH1D *m     = new TH1D("m","",100,0,16);
  TH1D *e     = new TH1D("e","",100,0,16);
  TH1D *nue   = new TH1D("nue","",100,0,16);

  
  for( int u = 0; u < N; ++u ) {
    //if(u%100==0) std::cout << u << "\n";

    m->Reset();
    e->Reset();
    nue->Reset();

    for(int CCmb=0; CCmb<nbins; CCmb++) {
      int fluxbin = CCmb+1;
      double evtwgt = scales[u][fluxbin];

      m->Fill(CC_m_nom->GetBinContent(CCmb+1),1.+evtwgt);
    }
    for(int CCeb=0; CCeb<nbins; CCeb++) {
      int fluxbin = CCeb+101;
      double evtwgt = scales[u][fluxbin];

      e->Fill(CC_e_nom->GetBinContent(CCeb+1),1.+evtwgt);
    }
    for(int nueb=0; nueb<nbins; nueb++) {
      int fluxbin = nueb+201;
      double evtwgt = scales[u][fluxbin];

      nue->Fill(nue_nom->GetBinContent(nueb+1),1.+evtwgt);
    }


    for(int i=0; i<nbins; i++){
      E_m[u][i]     = m->GetBinContent(i+1);
      E_e[u][i]     = e->GetBinContent(i+1);
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
        covar_mm     += (E_m[k][i]   - CC_m_nom->GetBinContent(i+1)) * (E_m[k][j]   - CC_m_nom->GetBinContent(j+1));
        covar_ee     += (E_e[k][i]   - CC_e_nom->GetBinContent(i+1)) * (E_e[k][j]   - CC_e_nom->GetBinContent(j+1)); 
        covar_nuenue += (E_nue[k][i] - nue_nom->GetBinContent(i+1))  * (E_nue[k][j] - nue_nom->GetBinContent(j+1)); 

        covar_me     += (E_m[k][i]   - CC_m_nom->GetBinContent(i+1)) * (E_e[k][j]   - CC_e_nom->GetBinContent(j+1));
        covar_mnue   += (E_m[k][i]   - CC_m_nom->GetBinContent(i+1)) * (E_nue[k][j] - nue_nom->GetBinContent(j+1));
        covar_em     += (E_e[k][i]   - CC_e_nom->GetBinContent(i+1)) * (E_m[k][j]   - CC_m_nom->GetBinContent(j+1)); 
        covar_enue   += (E_e[k][i]   - CC_e_nom->GetBinContent(i+1)) * (E_nue[k][j] - nue_nom->GetBinContent(j+1)); 
        covar_nuem   += (E_nue[k][i] - nue_nom->GetBinContent(i+1))  * (E_m[k][j]   - CC_m_nom->GetBinContent(j+1)); 
        covar_nuee   += (E_nue[k][i] - nue_nom->GetBinContent(i+1))  * (E_e[k][j]   - CC_e_nom->GetBinContent(j+1)); 

      }

      ECovars_mm[i][j]     = covar_mm    /(N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      ECovars_ee[i][j]     = covar_ee    /(N * CC_e_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      ECovars_nuenue[i][j] = covar_nuenue/(N * nue_nom->GetBinContent(i+1)  * nue_nom->GetBinContent(j+1));

      ECovars_me[i][j]   = covar_me  /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      ECovars_mnue[i][j] = covar_mnue/(N * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      ECovars_em[i][j]   = covar_em  /(N * CC_e_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      ECovars_enue[i][j] = covar_enue/(N * CC_e_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      ECovars_nuem[i][j] = covar_nuem/(N * nue_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
      ECovars_nuee[i][j] = covar_nuee/(N * nue_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j+1));
    }
  }

  for(int i = 0; i < nbins; i++) {
    for(int j = 0; j < nbins; j++) {
      ECovars[i][j]           = ECovars_mm[i][j];
      ECovars[i][j+nbins]   = ECovars_me[i][j];
      ECovars[i][j+2*nbins] = ECovars_mnue[i][j];

      ECovars[i+nbins][j]           = ECovars_em[i][j];
      ECovars[i+nbins][j+nbins]   = ECovars_ee[i][j];
      ECovars[i+nbins][j+2*nbins] = ECovars_enue[i][j];

      ECovars[i+2*nbins][j]           = ECovars_nuem[i][j];
      ECovars[i+2*nbins][j+nbins]   = ECovars_nuee[i][j];
      ECovars[i+2*nbins][j+2*nbins] = ECovars_nuenue[i][j];

    }
  }

  TH2D *hcv = new TH2D("hcv","",300,0,300,300,0,300);
  for(int i=0; i<300; i++) {
    for(int j=0; j<300; j++) {
      hcv->SetBinContent(i+1, j+1, ECovars[i][j]);
    }
  }
 
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  double nu, Ev;
  if(cutNu == 0)      nu = 10.0;
  else if(cutNu == 3) nu = 0.3;
  if(cutEv == 0)      Ev = 3.0;
  else if(cutEv == 1) Ev = 0.8;
  else if(cutEv == 2) Ev = 0.5;

  hcv->SetStats(0);
  hcv->SetTitle(Form(" Covariance (nu<%.1fGeV && Etheta2<%.1fMeV)",nu,Ev));

  //hcr->SetStats(0);
  //hcr->SetTitle(Form("%s Correlation (nu<%.1fGeV && Etheta2<%.1fMeV)",name,nu,Ev));

  TCanvas *c0 = new TCanvas("c0","",900,800);
  hcovmx->Draw("colz");
  c0->SaveAs("hcovmx.png");

  TCanvas *ccv = new TCanvas("ccv","",900,800);
  hcv->Draw("colz");
  ccv->SaveAs(Form("FC_stats_%d.png",N));

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
*/


/*
  TFile *out = new TFile(Form("/dune/app/users/qvuong/lownu/Feldman_Cousins/FC_tot3_%d.root",N),"RECREATE");
  hcv->Write();
  out->Close();
*/
  
}

