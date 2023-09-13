static const int N = 10; // number of universes
static const int nb = 2;
static const int nbins = nb*3;


TMatrixD covmx( nbins, nbins );
TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, nbins );

TMatrixD E_m(N, nb);
TMatrixD E_e(N, nb);
TMatrixD E_nue(N, nb);

 
TMatrixD ECovars_m    ( nb, nb );
TMatrixD ECovars_e    ( nb, nb );
TMatrixD ECovars_nue  ( nb, nb );

TMatrixD ECovars_me  ( nb, nb );
TMatrixD ECovars_mnue( nb, nb );
TMatrixD ECovars_em  ( nb, nb );
TMatrixD ECovars_enue( nb, nb );
TMatrixD ECovars_nuem( nb, nb );
TMatrixD ECovars_nuee( nb, nb );

TMatrixD ECovars( nbins, nbins );


void test()
{

  TH1D *CC_m_nom = new TH1D("CC_m_nom", "", nb,0,nb);
  TH1D *CC_e_nom = new TH1D("CC_e_nom", "", nb,0,nb);;
  TH1D *nue_nom  = new TH1D("nue_nom", "", nb,0,nb);;


  for( int x = 0; x < nb; ++x ) {
    CC_m_nom->SetBinContent(x+1, x+10);
    CC_e_nom->SetBinContent(x+1, x+1);
    nue_nom->SetBinContent(x+1, x+2);


    for( int y = 0; y < nb; ++y ) {
      if(x==y){
        covmx[x][y]           = CC_m_nom->GetBinContent( x+1 );
        covmx[x+nb][y+nb]     = CC_e_nom->GetBinContent( x+1 );
        covmx[x+2*nb][y+2*nb] = nue_nom->GetBinContent( x+1 );
        //std::cout << covmx[x][y] << "\t" << covmx[x+nbins][y+nbins] << "\t" << covmx[x+2*nbins][y+2*nbins] << "\n";
      }

      else{
        covmx[x][y]               = 0.;
        covmx[x+nb][y+nb]         = 0.;
        covmx[x+2*nb][y+2*nb]     = 0.;}
    }
  }
  
  //TCanvas *c = new TCanvas("c","",800,800);
  //CC_m_nom->Draw();

  TH2D *hcovmx = new TH2D("hcovmx","",nbins,0,nbins,nbins,0,nbins);
  TH2D *hchol = new TH2D("hchol","",nbins,0,nbins,nbins,0,nbins);
  for( int i = 0; i < nbins; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      hcovmx->SetBinContent(i+1, j+1, covmx[i][j]);
    }
  }
  for( int i = 0; i < nbins; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      std::cout << covmx[i][j] << "\t";
    }
      std::cout << "\n";
  }

    
  // cholesky decomposition
  TDecompChol decomp( covmx );
  if( !decomp.Decompose() ) {
    printf( "Main covariance matrix failed Cholesky decomposition\n" );
    return;
  }
  const TMatrixD chol = decomp.GetU();
  for( int i = 0; i < nbins; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      std::cout << chol[i][j] << "\t";
    }
      std::cout << "\n";
  }
  
  

  TRandom3 * rand = new TRandom3(12345);
  TH2D *hscales = new TH2D("hscales", "", nbins,0,nbins, nbins,0,nbins);

  // make random number matrix
  for( int i = 0; i < N; ++i ) {
    double mean = 0.;
    for( int j = 0; j < nbins; ++j ) {
      double val = rand->Gaus( 0., 1. );
      scales[i][j] = val;
      mean += val;
    }
/*
    mean /= N;
    for( int j = 0; j < nbins; ++j ) {
      double old = scales[i][j];
      scales[i][j] = old - mean;
    }
*/
  }


  printf("scales\n");
  for( int i = 0; i < N; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      if(abs(scales[i][j])>1)
      std::cout << scales[i][j] << "\t";
    }
      std::cout << "\n";
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

  printf("scaleCovars\n");
  for( int i = 0; i < nbins; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      std::cout << scaleCovars[i][j] << "\t";
    }
      std::cout << "\n";
  }

  TDecompChol scaleDecomp( scaleCovars );
  if( !scaleDecomp.Decompose() ) printf( "Scale matrix didn't decompolse\n" );
  TMatrixD toInvert = scaleDecomp.GetU();
  TMatrixD inverse = toInvert.Invert();
  scales *= inverse;

  scales *= chol;

  printf("scales_final\n");
  for( int i = 0; i < N; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      if(abs(scales[i][j])>1)
      std::cout << scales[i][j] << "\t";
    }
      std::cout << "\n";
  }


  TH1D *m     = new TH1D("m","",nb,0,nb);
  TH1D *e     = new TH1D("e","",nb,0,nb);
  TH1D *nue   = new TH1D("nue","",nb,0,nb);

  
  for( int u = 0; u < N; ++u ) {
    m->Reset();
    e->Reset();
    nue->Reset();
    
    for(int b=0; b<nb; b++) {
      int bin = b;
      double evtwgt = scales[u][bin];
      double bincontent = CC_m_nom->GetBinContent(b+1);
      m->SetBinContent(b+1, bincontent*(1.+evtwgt));
      //std::cout << bincontent << "\t" << evtwgt << "\t" << m->GetBinContent(b+1) << "\n";
    }
    for(int b=0; b<nb; b++) {
      int bin = b+nb;
      double evtwgt = scales[u][bin];
      double bincontent = CC_e_nom->GetBinContent(b+1);
      e->SetBinContent(b+1, bincontent*(1.+evtwgt));
    }
    for(int b=0; b<nb; b++) {
      int bin = b+2*nb;
      double evtwgt = scales[u][bin];
      double bincontent = nue_nom->GetBinContent(b+1);
      nue->SetBinContent(b+1, bincontent*(1.+evtwgt));
    }

    
    for(int i=0; i<nb; i++){
      E_m[u][i]     = m->GetBinContent(i+1);
      E_e[u][i]     = e->GetBinContent(i+1);
      E_nue[u][i]   = nue->GetBinContent(i+1);
    }

   }


  for( int i = 0; i < nb; ++i ) { // columns
    for( int j = 0; j < nb; ++j ) { // columns
      // compute column covariance
      double covar_m = 0.;
      double covar_e = 0.;
      double covar_nue = 0.;
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
        covar_m     += (E_m[k][i]     - CC_m_nom->GetBinContent(i+1))  * (E_m[k][j]     - CC_m_nom->GetBinContent(j+1));
        covar_e     += (E_e[k][i]     - CC_e_nom->GetBinContent(i+1))  * (E_e[k][j]     - CC_e_nom->GetBinContent(j+1)); 
        covar_nue   += (E_nue[k][i]   - nue_nom->GetBinContent(i+1))   * (E_nue[k][j]   - nue_nom->GetBinContent(j+1)); 

        covar_me   += (E_m[k][i]   - CC_m_nom->GetBinContent(i+1)) * (E_e[k][j]   - CC_e_nom->GetBinContent(j+1));
        covar_mnue += (E_m[k][i]   - CC_m_nom->GetBinContent(i+1)) * (E_nue[k][j] - nue_nom->GetBinContent(j+1));
        covar_em   += (E_e[k][i]   - CC_e_nom->GetBinContent(i+1)) * (E_m[k][j]   - CC_m_nom->GetBinContent(j+1)); 
        covar_enue += (E_e[k][i]   - CC_e_nom->GetBinContent(i+1)) * (E_nue[k][j] - nue_nom->GetBinContent(j+1)); 
        covar_nuem += (E_nue[k][i] - nue_nom->GetBinContent(i+1))  * (E_m[k][j]   - CC_m_nom->GetBinContent(j+1)); 
        covar_nuee += (E_nue[k][i] - nue_nom->GetBinContent(i+1))  * (E_e[k][j]   - CC_e_nom->GetBinContent(j+1)); 

      }

      ECovars_m[i][j]     = covar_m    /(N * CC_m_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
      ECovars_e[i][j]     = covar_e    /(N * CC_e_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j+1));
      ECovars_nue[i][j]   = covar_nue  /(N * nue_nom->GetBinContent(i+1)   * nue_nom->GetBinContent(j+1));

      ECovars_me[i][j]   = covar_me  /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      ECovars_mnue[i][j] = covar_mnue/(N * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      ECovars_em[i][j]   = covar_em  /(N * CC_e_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      ECovars_enue[i][j] = covar_enue/(N * CC_e_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      ECovars_nuem[i][j] = covar_nuem/(N * nue_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
      ECovars_nuee[i][j] = covar_nuee/(N * nue_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j+1));

    }
  }



  for(int i = 0; i < nb; i++) {
    for(int j = 0; j < nb; j++) {
      ECovars[i][j]         = ECovars_m[i][j];
      ECovars[i][j+nb]   = ECovars_me[i][j];
      ECovars[i][j+2*nb] = ECovars_mnue[i][j];

      ECovars[i+nb][j]         = ECovars_em[i][j];
      ECovars[i+nb][j+nb]   = ECovars_e[i][j];
      ECovars[i+nb][j+2*nb] = ECovars_enue[i][j];

      ECovars[i+2*nb][j]         = ECovars_nuem[i][j];
      ECovars[i+2*nb][j+nb]   = ECovars_nuee[i][j];
      ECovars[i+2*nb][j+2*nb] = ECovars_nue[i][j];
    }
  }



  TH2D *hcv = new TH2D("hcv","",nbins,0,nbins,nbins,0,nbins);
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      hcv->SetBinContent(i+1, j+1, ECovars[i][j]);
    }
  }
 
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  hcv->SetStats(0);

  TCanvas *ccv = new TCanvas("ccv","",900,800);
  hcv->SetMaximum(1.0);
  hcv->Draw("colz");
  //ccv->SaveAs(Form("Cov%d%d_%d.png",cutNu,cutEv,N));

/*
  TFile *out = new TFile(Form("/dune/app/users/qvuong/lownu/cov_matrix/%s_covmtr%d%d_%d.root",name,cutNu,cutEv,N),"RECREATE");
  hcv->Write();
  hcr->Write();
  out->Close();
  }
  //}
  //}
*/
}

