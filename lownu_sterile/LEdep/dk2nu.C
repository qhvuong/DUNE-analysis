void dk2nu()
{
  TRandom3 *rando = new TRandom3(12345); // 12345 is the most random of all seeds

  TH2D *h = new TH2D("h","",200,0,5,200,300,600);

  TFile * fout = new TFile( "/dune/app/users/qvuong/data/lownu/LEdep/LEdep.root", "RECREATE" );
  TTree * nudir = new TTree( "nudir", "neutrinodirection" );
  double nu_px, nu_py, nu_pz, nuE; // the neutrino momentum 4-vector
  double vtx_x, vtx_y, vtx_z, creation_z; // the neutrino *interaction* vertex and creation point in NuMI-z
  double dx, dy, dz, disp, vx, vy, vz;
  double wght, fipr, fipz, parent_p;
  int nu_pdg, parent_pdg;
  int is_kdar, POT;
  int count=0;

  nudir->Branch( "nuPx", &nu_px, "nuPx/D" );
  nudir->Branch( "nuPy", &nu_py, "nuPy/D" );
  nudir->Branch( "nuPz", &nu_pz, "nuPz/D" );
  nudir->Branch( "nuE", &nuE, "nuE/D" );
  nudir->Branch( "parent_p", &parent_p, "parent_p/D" );
  nudir->Branch( "vtx_x", &vtx_x, "vtx_x/D" );
  nudir->Branch( "vtx_y", &vtx_y, "vtx_y/D" );
  nudir->Branch( "vtx_z", &vtx_z, "vtx_z/D" );
  nudir->Branch( "nu_pdg", &nu_pdg, "nu_pdg/I" );
  nudir->Branch( "parent_pdg", &parent_pdg, "parent_pdg/I" );
  nudir->Branch( "fipz", &fipz, "fipz/D" );
  nudir->Branch( "fipr", &fipr, "fipr/D" );
  nudir->Branch( "creation_z", &creation_z, "creation_z/D" );
  nudir->Branch( "wght", &wght, "wght/D" );
  nudir->Branch( "is_kdar", &is_kdar, "is_kdar/I" );
  nudir->Branch( "vx", &vx, "vx" );
  nudir->Branch( "vy", &vy, "vy" );
  nudir->Branch( "vz", &vz, "vz/D" );

  TChain *meta = new TChain("dkmetaTree", "dkmetaTree");
  TChain *tree = new TChain("dk2nuTree", "dk2nuTree");

  for(int i=1; i<11; i++) {
    meta->Add(Form("/pnfs/dune/persistent/stash/Flux/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017/neutrino/flux/g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_%05d.dk2nu.root",i));
    tree->Add(Form("/pnfs/dune/persistent/stash/Flux/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017/neutrino/flux/g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_%05d.dk2nu.root",i));
  }

  double total_POT = 0.;
  int NFiles = meta->GetEntries();
  for( int i = 0; i < NFiles; ++i ) {
    meta->GetEntry(i);
    double pot = meta->GetLeaf("pots")->GetValue();
    total_POT += pot;
  }

  printf( "The total POT is %3.3g\n", total_POT );

  const int N = tree->GetEntries();

  for(int ii=0; ii<N; ii++)
  {	
    tree->GetEntry(ii);
    if( ii % 10000 == 0 ) printf( "Event %d of %d...\n", ii, N );

    int nuPDG = tree->GetLeaf("decay.ntype")->GetValue(); // neutrino PDG
    int parentPDG = tree->GetLeaf("decay.ptype")->GetValue(); // neutrino PDG
    double nu_vtxX = tree->GetLeaf("decay.vx")->GetValue(); // decay vertex / neutrino birthplace
    double nu_vtxY = tree->GetLeaf("decay.vy")->GetValue();
    double nu_vtxZ = tree->GetLeaf("decay.vz")->GetValue();
    double parent_px = tree->GetLeaf("decay.pdpx")->GetValue(); // momentum of parent, should be zero for DAR
    double parent_py = tree->GetLeaf("decay.pdpy")->GetValue();
    double parent_pz = tree->GetLeaf("decay.pdpz")->GetValue();
    double nu_E = tree->GetLeaf("decay.necm")->GetValue(0); // decay center of mass energy, gets updated later for DIF events
    double FIPx = tree->GetLeaf("ancestor.startx")->GetValue(1);
    double FIPy = tree->GetLeaf("ancestor.starty")->GetValue(1);
    double FIPz = tree->GetLeaf("ancestor.startz")->GetValue(1);
    double wt_multiplier = tree->GetLeaf("decay.nimpwt")->GetValue(0);

    double parent_p3 = sqrt( parent_px*parent_px + parent_py*parent_py + parent_pz*parent_pz );

    nu_pdg = nuPDG;
    fipz = FIPz;
    fipr = sqrt(FIPx*FIPx + FIPy*FIPy);

    creation_z = nu_vtxZ;
    parent_pdg = parentPDG;
    parent_p = parent_p3;

    if( parentPDG != 321 || parent_p3 > 0. ) is_kdar = 0;
    else is_kdar = 1;

    const TVector3 beam_dir(0.000000,-0.100828,0.994904);
    vtx_x = rando->Rndm()*600. - 300.;
    vtx_y = rando->Rndm()*200. - 100.;
    vtx_z = rando->Rndm()*500.;

    dx = vtx_x - nu_vtxX;
    dy = vtx_y - nu_vtxY;
    dz = (vtx_z + 57400) - nu_vtxZ;

    vx = nu_vtxX;
    vy = nu_vtxY;
    vz = nu_vtxZ;

    if(nu_vtxZ > 20000) {
      count ++;
      std::cout << vtx_z << "\t" << nu_vtxZ << "\t" << dz << "\n";
    }

    disp = sqrt(dx*dx + dy*dy + dz*dz);
    TVector3 p_vec(dx, dy ,dz);

    double parent_m;
    if (abs(parentPDG) == 13) parent_m = .1056583745;
    else if (abs(parentPDG) == 211) parent_m = .13957061;
    else if (abs(parentPDG) == 130) parent_m = .497611;
    else if (abs(parentPDG) == 321) parent_m = .493677;
    else if (abs(parentPDG) == 2112) parent_m = .9395654133;
    else printf("\n \n \n \n \n %d \n \n \n \n \n \n", parentPDG);

    double parent_E = sqrt(parent_p3*parent_p3 + parent_m*parent_m);
    double gamma = parent_E / parent_m;
    double beta_mag = sqrt((gamma*gamma - 1) / (gamma*gamma));

    double costheta = (dx*parent_px + dy*parent_py + dz*parent_pz)/ (parent_p3*disp);
    if (costheta > 1) costheta = 1.;
    else if (costheta < -1) costheta = -1.;
    else if (TMath::IsNaN(costheta) == true) costheta=0;

    double E_ratio = 1. / (gamma * (1. - beta_mag * costheta));
    wght = wt_multiplier * E_ratio*E_ratio / (4 * TMath::Pi() * disp*disp);

    nuE = nu_E * E_ratio;
    p_vec *= (nuE / disp);
    double p_nu[3]; //* This unrotated vector will come in handy for muon decay in flight weighting.
    p_nu[0] = p_vec[0]; //*
    p_nu[1] = p_vec[1]; //*
    p_nu[2] = p_vec[2]; //*
    p_vec.RotateUz(beam_dir);
    nu_px = p_vec[0];
    nu_py = p_vec[1];
    nu_pz = p_vec[2];

    int ndecay;
    if (abs(parentPDG) == 13) ndecay = tree->GetLeaf("decay.ndecay")->GetValue();
    else ndecay = 0;

    if (ndecay == 11 || ndecay == 12) {
      double ppenergy = tree->GetLeaf("decay.ppenergy")->GetValue();
      double ppdxdz = tree->GetLeaf("decay.ppdxdz")->GetValue();
      double ppdydz = tree->GetLeaf("decay.ppdydz")->GetValue();
      double pppz = tree->GetLeaf("decay.pppz")->GetValue();
      double muparpx = tree->GetLeaf("decay.muparpx")->GetValue();
      double muparpy = tree->GetLeaf("decay.muparpy")->GetValue();
      double muparpz = tree->GetLeaf("decay.muparpz")->GetValue();
      double mupare = tree->GetLeaf("decay.mupare")->GetValue();

      double beta[3];
      beta[0]=parent_px / parent_E;
      beta[1]=parent_py / parent_E;
      beta[2]=parent_pz / parent_E;
      double partial = gamma*(beta[0]*p_nu[0]+beta[1]*p_nu[1]+beta[2]*p_nu[2]);
      partial = (E_ratio*nu_E)-partial / (gamma+1.);
      double p_dcm_nu[4];
      for (int i=0;i<3;i++) p_dcm_nu[i]=p_nu[i]-beta[i]*gamma*partial;
      p_dcm_nu[3]=0.;
      for (int i=0;i<3;i++) p_dcm_nu[3]+=p_dcm_nu[i]*p_dcm_nu[i];
      p_dcm_nu[3]=sqrt(p_dcm_nu[3]);

      gamma=ppenergy / parent_m;
      beta[0] = ppdxdz * pppz / ppenergy;
      beta[1] = ppdydz * pppz / ppenergy;
      beta[2] = pppz / ppenergy;
      partial = gamma*(beta[0]*muparpx+beta[1]*muparpy+beta[2]*muparpz);
      partial = mupare - partial / (gamma+1.);
      double p_pcm_mp[4];
      p_pcm_mp[0]=muparpx-beta[0]*gamma*partial;
      p_pcm_mp[1]=muparpy-beta[1]*gamma*partial;
      p_pcm_mp[2]=muparpz-beta[2]*gamma*partial;
      p_pcm_mp[3]=0.;
      for (int i=0;i<3;i++) p_pcm_mp[3]+=p_pcm_mp[i]*p_pcm_mp[i];
      p_pcm_mp[3]=sqrt(p_pcm_mp[3]);
      double wt_ratio = 1.;

      if (p_pcm_mp[3] != 0.) {
        double costh = (p_dcm_nu[0]*p_pcm_mp[0]+p_dcm_nu[1]*p_pcm_mp[1]+p_dcm_nu[2]*p_pcm_mp[2])/(p_dcm_nu[3]*p_pcm_mp[3]);
        if (costh>1.) costh = 1.;
        else if (costh<-1.) costh = -1.;
        if (abs(nuPDG) == 12) wt_ratio = 1.-costh;
        else if (abs(nuPDG) == 14) {
          double xnu = 2.* nu_E / parent_m;
          wt_ratio = ( (3.-2.*xnu) - (1.-2.*xnu)*costh ) / (3.-2.*xnu);
        }
        else {
          cout << "NuWeight:: Bad neutrino type"<<endl;
        }
      }

      if (TMath::IsNaN(wt_ratio) == true) wt_ratio = 1;
      wght *= wt_ratio;
    }

    if(nuPDG==14 && nuE<5) h->Fill(nuE, (57400-nu_vtxZ)/100., wght);

    nudir->Fill();
  }

  POT = total_POT;
  nudir->Branch( "POT", &POT, "POT/I" )->Fill();

  gStyle->SetPalette(kColorPrintableOnGrey);
  TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  h->GetXaxis()->SetTitle("E_nu (GeV)");
  h->GetYaxis()->SetTitle("L (m)");
  TCanvas *c = new TCanvas("c","",800,600);
  h->Draw("colz");
  c->SaveAs("LEdep.png");


  std::cout << count << "\n";
  nudir->Write();
  //nu_vtxZ->Write();
  fout->Close();

}
