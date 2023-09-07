void flux()
{
  TH1D *h  = new TH1D("h","",50,0,15);
  TH1D *h1 = new TH1D("h1","",50,0,15);
  TH1D *h2 = new TH1D("h2","",50,0,15);

  double total_pot1 = 0., total_pot2 = 0., total_pot;

  for( int i = 0; i < 141; ++i ) { 
  TFile * tf = new TFile( Form("/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/Flux/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017/neutrino/gsimpleND/gsimple_DUNEND_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_%05d_%05d.root", i+1, i) );

  if(i%20==0) printf( "run %d\n", i );

  TTree * tree = (TTree*)tf->Get("flux");
  TTree * meta = (TTree*)tf->Get("meta");

  //std::cout << tree << std::endl;
  double pot;
  meta->SetBranchAddress("meta.protons",&pot);
  const int N = meta->GetEntries();
  //printf( "Got %d files\n", N );
  for(int ii=0; ii<N; ii++)
  {
    meta->GetEntry(ii);
    total_pot1 += pot;
  }

  int pdg;
  double E;
  double vtxx, vtxy;

  tree->SetBranchAddress("entry.E",&E);
  tree->SetBranchAddress("entry.vtxx",&vtxx);
  tree->SetBranchAddress("entry.vtxy",&vtxy);
  tree->SetBranchAddress("entry.pdg",&pdg);

  const int N2 = tree->GetEntries();
  //printf( "Got %d entries\n", N2 );

  for(int ii=0; ii<tree->GetEntries(); ii++)
  {
    tree->GetEntry(ii);
    //std::cout << pdg << "\t" << E << "\n";
    //if( pdg==12 && abs(vtxx)<3 && vtxy>-1.516 && vtxy<0.516 ) h_e->Fill(E);
    if( pdg==14 && abs(vtxx)<3 && vtxy>-1.516 && vtxy<0.516 ) h1->Fill(E);
  }

  tf->Close();
  delete tf;
  }


  for( int i = 141; i < 249; ++i ) { 
  //if(i<141)
  TFile * tf = new TFile( Form("/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/Flux/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017/neutrino/gsimpleND/gsimple_DUNEND_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_%05d_%05d.root", i+2, i) );
  //else
  //TFile * tf = new TFile( Form("/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/Flux/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017/neutrino/gsimpleND/gsimple_DUNEND_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_%05d_%05d.root", i+2, i) );

  if(i%20==0) printf( "run %d\n", i );
  //printf( "/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/Flux/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017/neutrino/gsimpleND/gsimple_DUNEND_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_%05d_%05d.root", i+1, i);

  TTree * tree = (TTree*)tf->Get("flux");
  TTree * meta = (TTree*)tf->Get("meta");

  //std::cout << tree << std::endl;
  double pot;
  meta->SetBranchAddress("meta.protons",&pot);
  const int N = meta->GetEntries();
  //printf( "Got %d files\n", N );
  for(int ii=0; ii<N; ii++)
  {
    meta->GetEntry(ii);
    total_pot2 += pot;
    //std::cout << pot << "\t" << std::endl;
  }

  int pdg;
  double E;
  double vtxx, vtxy;

  tree->SetBranchAddress("entry.E",&E);
  tree->SetBranchAddress("entry.vtxx",&vtxx);
  tree->SetBranchAddress("entry.vtxy",&vtxy);
  tree->SetBranchAddress("entry.pdg",&pdg);

  const int N2 = tree->GetEntries();
  //printf( "Got %d entries\n", N2 );

  for(int ii=0; ii<tree->GetEntries(); ii++)
  {
    tree->GetEntry(ii);
    //std::cout << pdg << "\t" << E << "\n";
    //if( pdg==12 && abs(vtxx)<3 && vtxy>-1.516 && vtxy<0.516 ) h_e->Fill(E);
    if( pdg==14 && abs(vtxx)<3 && vtxy>-1.516 && vtxy<0.516 ) h2->Fill(E);
  }

  tf->Close();
  delete tf;
  }

  total_pot = total_pot1 + total_pot2;
  double scalePOT = 1./total_pot;
  std::cout << total_pot << "\t" << scalePOT << "\n";

  h->Add(h1);
  h->Add(h2);
  h->Scale(scalePOT);

  TCanvas *c = new TCanvas("c","",800,600);
  h->SetTitle("Corrected numu flux");
  h->GetXaxis()->SetTitle("Ev (GeV)");
  h->Draw();
  c->SaveAs("numu_flux.png");

  TFile *out = new TFile("flux_50bins.root","RECREATE");
  h->Write();
  out->Close();


}
