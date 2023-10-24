static const int N = 1000000; 
static const int nbins_CC = 100;
static const int nbins_nue = 50;
static const int nbins = 2*nbins_CC + nbins_nue;


TMatrixD covmx( nbins, nbins );
TMatrixD sysmx( nbins, nbins );
TMatrixD statmx( nbins, nbins );
TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, nbins );
  

void test()
{

  int cutNu = 3;
  const char name[20] = "wgt_MaCCQE";

  // Get data files
  TFile *f     = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root", "READ");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/data/lownu/nue_output.root", "READ");

  // Get systematic uncertainties
  TFile *f_sys = new TFile(Form("../xS_covmtr/total_sigmtr%d_5sig.root",cutNu), "READ");
  TFile *f_fl  = new TFile(Form("../flux_covmtr/flux_covmtr%d_10000.root",cutNu),"READ");

  TH2D *hsys   = (TH2D*)f_sys->Get( "hcv" );
  TH2D *hfl    = (TH2D*)f_fl->Get("hcv");

  TH1D *CCm_nom = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma3",name,cutNu));
  TH1D *CCe_nom = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma3",name,cutNu));
  TH1D *nue_nom  = (TH1D*)f_nue->Get("hElep0");
  

  // only need the ND FHC part, which is the first 52 bins probably
  for( int x = 0; x < nbins; ++x ) {
    for( int y = 0; y < nbins; ++y ) {
      sysmx[x][y]  = hfl->GetBinContent(x+1, y+1) + hsys->GetBinContent(x+1, y+1);		//total sys = flux unc + cross section unc
      statmx[x][y] = 0.;

      if(x==y){
        if(x<nbins_CC)                  statmx[x][y] = CCm_nom->GetBinContent( x+1 );
        if(x>=nbins_CC && x<2*nbins_CC) statmx[x][y] = CCe_nom->GetBinContent( x-nbins_CC+1 );
        if(x>=2*nbins_CC)               statmx[x][y] = nue_nom->GetBinContent( x-2*nbins_CC+1 );
      }

    }
  }
  


  // MAKE FRACTIONAL COVARIANCE MATRIX
  for( int i = 0; i < nbins; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      double sum = statmx[i][j] + sysmx[i][j];
        if(i<nbins_CC){
          if(j<nbins_CC)                  covmx[i][j] = sum/(CCm_nom->GetBinContent(i+1) * CCm_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) covmx[i][j] = sum/(CCm_nom->GetBinContent(i+1) * CCe_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               covmx[i][j] = sum/(CCm_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}
      

        if(i>=nbins_CC && i<2*nbins_CC){ 
          if(j<nbins_CC)                  covmx[i][j] = sum/(CCe_nom->GetBinContent(i-nbins_CC+1) * CCm_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) covmx[i][j] = sum/(CCe_nom->GetBinContent(i-nbins_CC+1) * CCe_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               covmx[i][j] = sum/(CCe_nom->GetBinContent(i-nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}


        if(i>=2*nbins_CC){
          if(j<nbins_CC)                  covmx[i][j] = sum/(nue_nom->GetBinContent(i-2*nbins_CC+1) * CCm_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) covmx[i][j] = sum/(nue_nom->GetBinContent(i-2*nbins_CC+1) * CCe_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               covmx[i][j] = sum/(nue_nom->GetBinContent(i-2*nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}
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
  
  
  double mean;
  double val, old;
  for( int i = 0; i < N; ++i ) {
    mean = 0.;
    for( int j = 0; j < nbins; ++j ) {
      val = rand->Gaus( 0., 1. );
      scales[i][j] = val;
      mean += val;
    }

    mean /= N;
    // force each column mean to be exactly 0, this just eliminates tiny statistical fluctuations in the mean weight
    for( int j = 0; j < nbins; ++j ) {
      old = scales[i][j];
      scales[i][j] = old - mean;
    }
  }


  double covar;
  for( int i = 0; i < nbins; ++i ) { // columns
    for( int j = 0; j < nbins; ++j ) { // columns
      // compute column covariance
      covar = 0.;
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
  scales = scales*inverse;
  scales *= chol;


  TFile *out = new TFile(Form("FC%d_%d.root",cutNu,N),"RECREATE");
  scales.Write("scales");
  out->Close();




 

}

