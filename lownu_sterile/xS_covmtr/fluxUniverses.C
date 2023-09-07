//static const int N = 10000; // number of universes
static const int N = 60; // number of universes
static const int nbins = 52;

const int n_mu = 19; // number of muon bins in covariance
const int n_e = 7; // number of electron bins in covariance

// bin edges -- fix these to be whatever they actually are
const double mubins[20] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,12.,16.,20.,40.,100.};
const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};

void fluxUniverses()
{

  // covariance matrix
  TFile * covfile = new TFile( "/nashome/q/qvuong/cov_matrix/total_covariance_DUNE_opt.root", "OLD" );
  TH2D * hcovmx = (TH2D*) covfile->Get( "total_covariance" );

  // only need the ND FHC part, which is the first 52 bins probably
  TMatrixD covmx( nbins, nbins );
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
  TMatrixD scales( N, nbins );
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

  TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
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


  for( int u = 0; u < N; ++u ) {

  // Fill the varied histograms
  // declare 1000 histograms of lepton energy here
  TChain * tree = new TChain( "cafTree", "cafTree" );
  for(int i = 150; i<151; i++){
  tree->Add( Form("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%d.root",i) );
  std::cout << "File number:" << i << "\n";
  }

  const Int_t nbinsX = 480; const Int_t nbinsY = 100;
  Double_t xEdges[nbinsX+1], yEdges[nbinsY+1];
  xEdges[0]=yEdges[0]=0;
  for(int i=0; i<nbinsX+1; i++)
  {
    if(i<200)                xEdges[i+1] = xEdges[i] + 0.02;
    else if(i>=200 && i<240) xEdges[i+1] = xEdges[i] + 0.1;
    else if(i>=240 && i<400) xEdges[i+1] = xEdges[i] + 0.2;
    else                     xEdges[i+1] = xEdges[i] + 1.0;
  }
  for(int i=0; i<nbinsY+1; i++)
  {
    yEdges[i+1]=yEdges[i]+0.16;
  }

  TH1D *e_hElep_w0 = new TH1D("e_hElep_w0","",100,0,16);
  TH1D *e_hElep_w3 = new TH1D("e_hElep_w3","",100,0,16);
  TH1D *m_hElep_w0 = new TH1D("m_hElep_w0","",100,0,16);
  TH1D *m_hElep_w3 = new TH1D("m_hElep_w3","",100,0,16);
  TH2D *m_hElepVsEv0 = new TH2D("m_hElepVsEv0","",n_mu,mubins,100,0,16);
  TH2D *m_hElepVsEv3 = new TH2D("m_hElepVsEv3","",n_mu,mubins,nbinsY,yEdges);
  TH2D *e_hElepVsEv0 = new TH2D("e_hElepVsEv0","",n_e,ebins,nbinsY,yEdges);
  TH2D *e_hElepVsEv3 = new TH2D("e_hElepVsEv3","",n_e,ebins,nbinsY,yEdges);

  std::cout << "test0" << "\n";

  double vtx_x, vtx_y, vtx_z; // the position where the neutrino interaction occurred, in cm
  int nuPDG; // PDG code of the neutrino, numu = 14, nue = 12, antineutrinos are negative
  double Ev; // the energy of the neutrino, in GeV
  int LepPDG; // PDG code of the final-state lepton, mu = 13, e = 11
  double LepE; // Total energy of the final-state lepton; note that true nu is not saved but nu = Ev - LepE
  double LepNuAngle; // angle between lepton and neutrino
  double eP, eN, ePip, ePim, ePi0, eOther; // energy in the final state due to different particles
  double Ev_reco, Elep_reco; // the measured neutrino energy and lepton energy. Reco nu = Ev_reco - Elep_reco
  int reco_numu, reco_nue; // = 1 if the reconstruction thinks it's a muon or electron
  int muon_contained, muon_tracker, muon_ecal; // different ways that the muon can be measured
  double Ehad_veto; // hadronic energy near the edge of the detector, which is a hint that we might not have measured all the energy

  tree->SetBranchAddress( "vtx_x", &vtx_x );
  tree->SetBranchAddress( "vtx_y", &vtx_y );
  tree->SetBranchAddress( "vtx_z", &vtx_z );
  tree->SetBranchAddress( "Ev", &Ev );
  tree->SetBranchAddress( "nuPDG", &nuPDG );
  tree->SetBranchAddress( "LepPDG", &LepPDG );
  tree->SetBranchAddress( "LepE", &LepE );
  tree->SetBranchAddress( "LepNuAngle", &LepNuAngle );
  tree->SetBranchAddress( "Ev_reco", &Ev_reco );
  tree->SetBranchAddress( "Elep_reco", &Elep_reco );
  tree->SetBranchAddress( "reco_numu", &reco_numu );
  tree->SetBranchAddress( "reco_nue", &reco_nue );
  tree->SetBranchAddress( "muon_contained", &muon_contained );
  tree->SetBranchAddress( "muon_tracker", &muon_tracker );
  tree->SetBranchAddress( "muon_ecal", &muon_ecal );
  tree->SetBranchAddress( "Ehad_veto", &Ehad_veto );
  tree->SetBranchAddress( "eP", &eP );
  tree->SetBranchAddress( "eN", &eN );
  tree->SetBranchAddress( "ePip", &ePip );
  tree->SetBranchAddress( "ePim", &ePim );
  tree->SetBranchAddress( "ePi0", &ePi0 );
  tree->SetBranchAddress( "eOther", &eOther );
  tree->SetBranchAddress( "LepNuAngle", &LepNuAngle );

  std::cout << "test1" << "\n";
  TRandom *r1 = new TRandom(8888);

  int Nevents = tree->GetEntries(); // replace with number events in TTree
  for( int evt = 0; evt < Nevents; ++ evt ) {

    if( evt % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", evt*100.0/Nevents, Nevents );
    tree->GetEntry(evt);

    // apply the reco evente selection here

    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;
    double nu = eP + eN + ePip + ePim + ePi0 + eOther;

    if(LepPDG == 13){
      if( reco_numu && (muon_contained || muon_tracker || muon_ecal)){
        int fluxbin = 0; 
        for(int b = 0; b < n_mu+1; b++){
          if(Ev>mubins[b]){
            if(Ev<mubins[b+1]) fluxbin = b+1;  // determine which bin the event is in from 0-51 based on the neutrino flavor and true neutrino energy
            else continue;}
        }
        std::cout << Ev << "\t" << fluxbin << "\n";
          //double lepton_energy = 0.; // find the reco lepton energy
          // loop over universes
        //for( int u = 0; u < N; ++u ) {
        double evtwgt = scales[u][fluxbin]; // blow[i] is the offset for this flux, i.e. 104 for FD FHC numu
        if(nu<10.0){ 
          m_hElepVsEv0->Fill(Ev,Elep_reco,evtwgt);
          m_hElep_w0->Fill(Elep_reco,evtwgt);}
        if(nu<0.3){  
          m_hElepVsEv3->Fill(Ev,Elep_reco,evtwgt); // fill your histogram with this weight
          m_hElep_w3->Fill(Elep_reco,evtwgt);}
      }
    }
    
/*
    if(LepPDG == 11){
      double LepE_sm = r1->Gaus(LepE,0.05);
      if(LepE_sm*LepNuAngle*LepNuAngle*1E3>3.){
      if( reco_nue ){
        int fluxbin = 0; 
        for(int b = 0; b < n_e+1; b++){
          if(Ev>mubins[b]){
            if(Ev<mubins[b+1]) fluxbin = 38+b+1;  // determine which bin the event is in from 0-51 based on the neutrino flavor and true neutrino energy
            else continue;}
        }
        std::cout << Ev << "\t" << fluxbin << "\n";
          //double lepton_energy = 0.; // find the reco lepton energy
          // loop over universes
          double evtwgt = scales[u][fluxbin]; // blow[i] is the offset for this flux, i.e. 104 for FD FHC numu
          if(nu<10.0){ 
            e_hElepVsEv0->Fill(Ev,Elep_reco,evtwgt);
            e_hElep_w0->Fill(Elep_reco,evtwgt);}
          if(nu<0.3){  
            e_hElepVsEv3->Fill(Ev,Elep_reco,evtwgt); // fill your histogram with this weight
            e_hElep_w3->Fill(Elep_reco,evtwgt);}
        }
      }
    }
*/  
    }

    m_hElep_w0->SetTitle("CC_mu Reconstructed lepton energy (nu<10.0)");
    m_hElep_w0->GetXaxis()->SetTitle("Elep_reco (GeV)");
    m_hElep_w3->SetTitle("CC_mu Reconstructed lepton energy (nu<0.3)");
    m_hElep_w3->GetXaxis()->SetTitle("Elep_reco (GeV)");
    e_hElep_w0->SetTitle("CC_e Reconstructed lepton energy (nu<10.0)");
    e_hElep_w0->GetXaxis()->SetTitle("Elep_reco (GeV)");
    e_hElep_w3->SetTitle("CC_e Reconstructed lepton energy (nu<0.3)");
    e_hElep_w3->GetXaxis()->SetTitle("Elep_reco (GeV)");
    TCanvas *c = new TCanvas("c","",1600,1200);
    c->Divide(2,2);
    c->cd(1);
    m_hElep_w0->Draw();
    c->cd(2);
    e_hElep_w0->Draw();
    c->cd(3);
    m_hElep_w3->Draw();
    c->cd(4);
    e_hElep_w3->Draw();
    c->SaveAs(Form("hElep_covmtr_%d.pdf",u));

    TFile *out = new TFile(Form("covmtr%d.root",u),"RECREATE");
    m_hElep_w0->Write();
    m_hElep_w3->Write();
    e_hElep_w0->Write();
    e_hElep_w3->Write();
    out->Close();
  
  }

}
