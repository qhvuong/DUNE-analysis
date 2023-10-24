#include "TemplateFitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TGraph.h"
#include "TLegend.h"
#include <TRandom3.h>

TemplateFitter::TemplateFitter(TH1D * CC_templates_m[nbins_Ev], TH1D * CC_templates_m_nc[nbins_Ev], TH1D * CC_templates_e[nbins_Ev], TH1D * nue_templates_m[nbins_Ev], TH1D * nue_templates_m_w[nbins_Ev], TH1D * nue_templates_e[nbins_Ev], TH1D * nue_templates_e_w[nbins_Ev], TH1D * LEdep_m[29], TH1D * LEdep_e[29] )
{
  for( int i = 0; i < nbins_Ev; ++i ) {
  CC_m_templates[i] = CC_templates_m[i];
  CC_nc_m_templates[i] = CC_templates_m_nc[i];
  CC_e_templates[i] = CC_templates_e[i];
  nue_m_templates[i] = nue_templates_m[i];
  nue_w_m_templates[i] = nue_templates_m_w[i];
  nue_e_templates[i] = nue_templates_e[i];
  nue_w_e_templates[i] = nue_templates_e_w[i];
  }

  for( int i = 0; i < 29; i++ ) {
  LE_m[i] = LEdep_m[i];
  LE_e[i] = LEdep_e[i];
  }
}

void TemplateFitter::setEnergyBins( double bins[nbins_Ev+1] )
{
  for( int i = 0; i < nbins_Ev+1; ++ i ) {m_energy_bins[i] = bins[i];}
}

void TemplateFitter::setPara( char var[20], int par, double oscpar[3], int nuCut, double seed[3], double fitPara_m[29][7], double fitPara_e[29][7] )
{
  name  = var;
  para  = par;
  cutNu = nuCut;
  for(int i = 0; i < 3; i++){
    ospar[i] = oscpar[i];
    std::cout << ospar[i];
  }
  b0 = oscpar[0];
  b1 = oscpar[1];
  b2 = oscpar[2];

  s0 = seed[0];
  s1 = seed[1];
  s2 = seed[2];

  for(int i=1;i<29;i++){
    fitP_m[i][0] = fitPara_m[i][0];
    fitP_m[i][1] = fitPara_m[i][1];
    fitP_m[i][2] = fitPara_m[i][2];
    fitP_m[i][3] = fitPara_m[i][3];
    fitP_m[i][4] = fitPara_m[i][4];
    fitP_m[i][5] = fitPara_m[i][5];
    fitP_m[i][6] = fitPara_m[i][6];

    fitP_e[i][0] = fitPara_e[i][0];
    fitP_e[i][1] = fitPara_e[i][1];
    fitP_e[i][2] = fitPara_e[i][2];
    fitP_e[i][3] = fitPara_e[i][3];
    fitP_e[i][4] = fitPara_e[i][4];
    fitP_e[i][5] = fitPara_e[i][5];
    fitP_e[i][6] = fitPara_e[i][6];
  }
}

double TemplateFitter::getPmue( double energy, double Uee2, double Umm2, double dm2, double L )
{
  //double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2mue2 = 4 * Uee2 * Umm2;
  double prob = s2mue2  * pow(sin(del),2);
  return prob;
}
double TemplateFitter::getPee( double energy, double Uee2, double Umm2, double dm2, double L )
{
  //double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2ee2 = 4 * Uee2 * (1 - Uee2);
  double prob = 1.0 - s2ee2  * pow(sin(del),2);
  return prob;
}
double TemplateFitter::getPmm( double energy, double Uee2, double Umm2, double dm2, double L )
{
  //double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2mm2 = 4 * Umm2 * (1 - Umm2);
  double prob = 1.0 - s2mm2  * pow(sin(del),2);
  return prob;
}

double TemplateFitter::getAvgPmue( double energy, double Uee2, double Umm2, double dm2, double ft[7] )
{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * Uee2 * Umm2;

  L=L1;
  double prob_u = - avg1 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);
  L=L0;
  double prob_l = - avg1 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);

  double prob1 = prob_u - prob_l;


  L=L2;
  double term1_u = 2*pow(k,4) * L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_u = 6*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_u = 3*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );
  L=L1;
  double term1_l = 2*pow(k,4) * L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_l = 6*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_l = 3*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );

  prob_u = A/(48*pow(k,4)) * (term1_u - term2_u - term3_u);
  prob_l = A/(48*pow(k,4)) * (term1_l - term2_l - term3_l);

  double prob2 = prob_u - prob_l;


  L=L3;
  prob_u = - avg2 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);
  L=L2;
  prob_l = - avg2 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);

  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);

  return prob;
}

double TemplateFitter::getAvgPee( double energy, double Uee2, double Umm2, double dm2, double ft[7] )
{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * Uee2 * (1.-Uee2);

  L=L1;
  double prob_u = avg1 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L0;
  double prob_l = avg1 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  double prob1 = prob_u - prob_l;

  L=L2;
  double term1_u = -2*(A-2)*pow(k,4)*L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_u = 6*A*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_u = 3*A*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );
  L=L1;
  double term1_l = -2*(A-2)*pow(k,4)*L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_l = 6*A*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_l = 3*A*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );
  prob_u = 1/(48*pow(k,4)) * (term1_u + term2_u + term3_u);
  prob_l = 1/(48*pow(k,4)) * (term1_l + term2_l + term3_l);
  double prob2 = prob_u - prob_l;

  L=L3;
  prob_u = avg2* ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L2;
  prob_l = avg2* ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);
  return prob;
}

double TemplateFitter::getAvgPmm( double energy, double Uee2, double Umm2, double dm2, double ft[7] )
{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * Umm2 * (1.-Umm2);

  L=L1;
  double prob_u = avg1 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L0;
  double prob_l = avg1 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  double prob1 = prob_u - prob_l;

  L=L2;
  double term1_u = -2*(A-2)*pow(k,4)*L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_u = 6*A*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_u = 3*A*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );
  L=L1;
  double term1_l = -2*(A-2)*pow(k,4)*L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_l = 6*A*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_l = 3*A*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );
  prob_u = 1/(48*pow(k,4)) * (term1_u + term2_u + term3_u);
  prob_l = 1/(48*pow(k,4)) * (term1_l + term2_l + term3_l);
  double prob2 = prob_u - prob_l;

  L=L3;
  prob_u = avg2* ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L2;
  prob_l = avg2* ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);
  return prob;
}


TMatrixD covmtr(nbins,nbins);

void TemplateFitter::setCovmtr( double bincontent[nbins+1][nbins+1] )
{
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      covmtr[i][j] = bincontent[i][j];
    }
  }
}

void TemplateFitter::getTarget( double *par )
{
  CCe_tgt->Reset();
  CCm_tgt->Reset();
  nue_tgt->Reset();

  gRandom = new TRandom(12345); 
  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m  = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_em = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_mm = (TH1D*) CC_tp_e->Clone();

  TH1D * nue_tp = (TH1D*) nue_m_templates[0]->Clone();
  nue_tp->Reset();
  TH1D * nue_tp_em = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_ee = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_me = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_mm = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_os = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_unos = (TH1D*) nue_tp->Clone(); 

  TCanvas *c0 = new TCanvas("c0","",800,600);

  int iL;
/*
  TH2D *PmueE = new TH2D("PmueE","",100,0,20,100,0.,0.002);
  TH2D *PmueE1 = new TH2D("PmueE1","",100,0,20,100,0.,0.002);
  TH2D *PemuE = new TH2D("PemuE","",100,0,20,100,0.,0.002);
  TH2D *PemuE1 = new TH2D("PemuE1","",100,0,20,100,0.,0.002);
  TH2D *PeeE = new TH2D("PeeE","",100,0,20,100,0.8,1.);
  TH2D *PeeE1 = new TH2D("PeeE1","",100,0,20,100,0.8,1.);
  TH2D *PmmE = new TH2D("PmmE","",100,0,20,100,0.95,1.);
  TH2D *PmmE1 = new TH2D("PmmE1","",100,0,20,100,0.95,1.);
*/
  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < nbins_Ev; ++i ) {
    double mue = 0;
    double emu = 0;
    double ee  = 0;
    double mm  = 0;

    if(i<240) 		     iL = (int) i/10;
    else if(i>=240 && i<400) iL = (int) 24+(i-240)/40;
    else		     iL = 28;

    for(int k=0; k<7; k++) {
      ft_m[k] = fitP_m[iL][k];
      ft_e[k] = fitP_e[iL][k];
    }

    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getAvgPmue(e, par[0], par[1], par[2], ft_m); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      emu = emu + getAvgPmue(e, par[0], par[1], par[2], ft_e); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getAvgPee(e, par[0], par[1], par[2], ft_e);
      mm  = mm  + getAvgPmm(e, par[0], par[1], par[2], ft_m);
    }
    double Pmue = mue/1001.0;
    double Pemu = emu/1001.0;
    double Pee  = ee/1001.0;
    double Pmm  = mm/1001.0;

    CC_tp_me->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pemu);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pemu);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }


  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);

  CCe_tgt->Add(CC_tp_e);
  CCm_tgt->Add(CC_tp_m);
  nue_tgt->Add(nue_tp);

/*
  TRandom3 *rando = new TRandom3(12345);

  for(int bx = 1; bx <= CC_tp_e->GetNbinsX(); bx++){
    double mean = CC_tp_e->GetBinContent(bx);
    double fl_bc = rando->Poisson(mean);
    CCe_tgt->AddBinContent(bx, fl_bc);
  }

  for(int bx = 1; bx <= CC_tp_m->GetNbinsX(); bx++){
    double mean = CC_tp_m->GetBinContent(bx);
    double fl_bc = rando->Poisson(mean);
    CCm_tgt->AddBinContent(bx, fl_bc);
  }

  for(int bx =  1; bx <= nue_tp->GetNbinsX(); bx++){
    double mean = nue_tp->GetBinContent(bx);
    double fl_bc = rando->Poisson(mean);
    nue_tgt->AddBinContent(bx, fl_bc);
  }
*/
}


// function whose return Minuit mimizes, must take const double* and return double
double TemplateFitter::getChi2( double * par )
{
  gRandom = new TRandom(12345);
  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m  = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_em = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_mm = (TH1D*) CC_tp_e->Clone();

  TH1D * nue_tp = (TH1D*) nue_m_templates[0]->Clone();
  nue_tp->Reset();
  TH1D * nue_tp_em = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_ee = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_me = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_mm = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_os = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_unos = (TH1D*) nue_tp->Clone(); 

  //TRandom *rand = new TRandom(4444);
  //double L = rand->Uniform(0.37,0.57);

  int iL;

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < nbins_Ev; ++i ) {
    double mue = 0.;
    double emu = 0.;
    double ee = 0.;
    double mm = 0.;

    if(i<240)                iL = (int)i/10;
    else if(i>=240 && i<400) iL = (int) 24+(i-240)/40;
    //else if(i>=400 && i<440) iL = 28;
    else                     iL = 28;

    for(int k=0; k<7; k++) {
      ft_m[k] = fitP_m[iL][k];
      ft_e[k] = fitP_e[iL][k];
    }

    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getAvgPmue(e, par[0], par[1], par[2], ft_m); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      emu = emu + getAvgPmue(e, par[0], par[1], par[2], ft_e); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getAvgPee(e, par[0], par[1], par[2], ft_e);
      mm  = mm  + getAvgPmm(e, par[0], par[1], par[2], ft_m);
    }

    double Pmue = mue/1001.0;
    double Pemu = emu/1001.0;
    double Pee  = ee/1001.0;
    double Pmm  = mm/1001.0;

    CC_tp_me->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pemu);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pemu);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }

  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);

  // calculate the chi2 with the "data" target
  double chi2 = 0.0;

  TMatrixD target(nbins, 1);
  TMatrixD temp(nbins, 1);
  TMatrixD unc(nbins, nbins);

  for( int bx = 0; bx < nbins; bx++ ) {
    for( int by = 0; by < nbins; by++ ) {
      unc[bx][by] = 0.;
  }
  }

  for( int bx = 0; bx < nbins; bx++ ) {
    if(bx<nbins_CC)                   temp[bx][0] = CC_tp_m->GetBinContent(bx+1);
    if(bx>=nbins_CC && bx<2*nbins_CC) temp[bx][0] = CC_tp_e->GetBinContent(bx-nbins_CC+1);
    if(bx>=2*nbins_CC)                temp[bx][0] = nue_tp->GetBinContent(bx-2*nbins_CC+1);

    if(bx<nbins_CC)                   target[bx][0] = CCm_tgt->GetBinContent(bx+1);
    if(bx>=nbins_CC && bx<2*nbins_CC) target[bx][0] = CCe_tgt->GetBinContent(bx-nbins_CC+1);
    if(bx>=2*nbins_CC)                target[bx][0] = nue_tgt->GetBinContent(bx-2*nbins_CC+1);

    if(bx<nbins_CC)                   unc[bx][bx] = CCm_tgt->GetBinContent(bx+1);
    if(bx>=nbins_CC && bx<2*nbins_CC) unc[bx][bx] = CCe_tgt->GetBinContent(bx-nbins_CC+1);
    if(bx>=2*nbins_CC)                unc[bx][bx] = nue_tgt->GetBinContent(bx-2*nbins_CC+1);
  }

  //covmtr = covmtr + unc; 
  covmtr = unc; 

  TMatrixD diff = temp - target;
  TMatrixD diff_T(TMatrixD::kTransposed, diff);
  TDecompSVD svd(covmtr);
  TMatrixD inv = svd.Invert();

  TMatrixD diff_cov(1, nbins);
  for( int bx = 0; bx < nbins; bx++ ) {
    for( int by = 0; by < nbins; by++ ) {
      diff_cov[0][bx] += diff_T[0][by] * inv[by][bx];
    }
    chi2 += diff_cov[0][bx] * diff[bx][0];
  }

  std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";
  return chi2;

}


void TemplateFitter::Draw()
{
  int N = 50;

  double p0[1], p1[1], p2[1];
  p0[0] = b0;
  p1[0] = b1;
  p2[0] = b2;


  TH2D *h0 = new TH2D("h0","",N,0,0.1, N,0,12.0);
  TH2D *h1 = new TH2D("h1","",N,0,12.0,N,0,0.1);
  TH2D *h2 = new TH2D("h2","",N,0,0.1, N,0,0.1);
/*
  TH2D *h0L = new TH2D("h0L","",N,0,0.1, N,0,12.0);
  TH2D *h1L = new TH2D("h1L","",N,0,12.0,N,0,0.1);
  TH2D *h2L = new TH2D("h2L","",N,0,0.1, N,0,0.1);
*/
  TH2D *h0z = new TH2D("h0z","",N,0.0      ,2.0*p1[0],N,0.5*p2[0],1.5*p2[0]);
  TH2D *h1z = new TH2D("h1z","",N,0.5*p2[0],1.5*p2[0],N,0.0      ,2.0*p0[0]);
  TH2D *h2z = new TH2D("h2z","",N,0.0      ,2.0*p0[0],N,0.0      ,2.0*p1[0]);

  TH2D *h0d = new TH2D("h0d","",N,0,0.1, N,0,12.0);
  TH2D *h1d = new TH2D("h1d","",N,0,12.0,N,0,0.1);
  TH2D *h2d = new TH2D("h2d","",N,0,0.1, N,0,0.1);
/*
  TH2D *h0dL = new TH2D("h0dL","",N,0,0.1, N,0,12.0);
  TH2D *h1dL = new TH2D("h1dL","",N,0,12.0,N,0,0.1);
  TH2D *h2dL = new TH2D("h2dL","",N,0,0.1, N,0,0.1);
*/
  TH2D *h0dz = new TH2D("h0dz","",N,0.0      ,2.0*p1[0],N,0.5*p2[0],1.5*p2[0]);
  TH2D *h1dz = new TH2D("h1dz","",N,0.5*p2[0],1.5*p2[0],N,0.0      ,2.0*p0[0]);
  TH2D *h2dz = new TH2D("h2dz","",N,0.0      ,2.0*p0[0],N,0.0      ,2.0*p1[0]);

  double par[3], parz[3], bin[3];
  par[0] = b0;
  par[1] = b1;
  par[2] = b2;
  double chi2 = 0.0;
  double chi2z = 0.0;

  double nu, Ev;
  if(cutNu == 0)      nu = 10.0;
  else if(cutNu == 3) nu = 0.3;

  if(cutEv == 0)      Ev = 3.0;
  else if(cutEv == 1) Ev = 0.8;
  else if(cutEv == 2) Ev = 0.5;

  double chi2t = getChi2(par);
  double diff, diffz;

  par[0]  = b0;
  parz[0] = b0;
  for(int j=1; j<=N; j++) {
    par[1]  = h0->GetXaxis()->GetBinCenter(j);
    parz[1] = h0z->GetXaxis()->GetBinCenter(j);
    for(int k=1; k<=N; k++) {
      par[2]  = h1->GetXaxis()->GetBinCenter(k);
      parz[2] = h1z->GetXaxis()->GetBinCenter(k);
      chi2  = getChi2(par);
      chi2z = getChi2(parz);
      diff  = sqrt(std::fabs(chi2 - chi2t));
      diffz = sqrt(std::fabs(chi2z - chi2t));
      if( k==N ) std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\t" << diff << "\t" << j*k*100./(3*N*N) << "%" << "\n";
      h0->Fill(par[1], par[2], chi2);
      //h0L->Fill(par[1], par[2], chi2);
      h0z->Fill(parz[1], parz[2], chi2z);
      h0d->Fill(par[1], par[2], diff);
      //h0dL->Fill(par[1], par[2], diff);
      h0dz->Fill(parz[1], parz[2], diffz);
    }
  }

  par[1]  = b1;
  parz[1] = b1;
  for(int k=1; k<=N; k++) {
    par[2]  = h1->GetXaxis()->GetBinCenter(k);
    parz[2] = h1z->GetXaxis()->GetBinCenter(k);
    for(int i=1; i<=N; i++) {
      par[0]  = h2->GetXaxis()->GetBinCenter(i);
      parz[0] = h2z->GetXaxis()->GetBinCenter(i);
      chi2  = getChi2(par);
      chi2z = getChi2(parz);
      diff  = sqrt(std::fabs(chi2 - chi2t));
      diffz = sqrt(std::fabs(chi2z - chi2t));
      if( i==N ) std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\t" << diff << "\t" << (N*N+k*i)*100./(3*N*N) << "%" << "\n";
      h1->Fill(par[2], par[0], chi2);
      //h1L->Fill(par[2], par[0], chi2);
      h1z->Fill(parz[2], parz[0], chi2z);
      h1d->Fill(par[2], par[0], diff);
      //h1dL->Fill(par[2], par[0], diff);
      h1dz->Fill(parz[2], parz[0], diffz);
    }
  }

  par[2]  = b2;
  parz[2] = b2;
  for(int i=1; i<=N; i++) {
    par[0]  = h2->GetXaxis()->GetBinCenter(i);
    parz[0] = h2z->GetXaxis()->GetBinCenter(i);
    for(int j=1; j<=N; j++) {
      par[1]  = h0->GetXaxis()->GetBinCenter(j);
      parz[1] = h0z->GetXaxis()->GetBinCenter(j);
      chi2  = getChi2(par);
      chi2z = getChi2(parz);
      diff = sqrt(std::fabs(chi2 - chi2t));
      diffz = sqrt(std::fabs(chi2z - chi2t));
      if( j==N ) std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\t" << diff << "\t" << (2*N*N+i*j)*100./(3*N*N) << "%" << "\n";
      h2->Fill(par[0], par[1], chi2);
      //h2L->Fill(par[0], par[1], chi2);
      h2z->Fill(parz[0], parz[1], chi2z);
      h2d->Fill(par[0], par[1], diff);
      //h2dL->Fill(par[0], par[1], diff);
      h2dz->Fill(parz[0], parz[1], diffz);
    }
  }

  int n = 1;
  TGraph *g0 = new TGraph(n,p1,p2);
  TGraph *g1 = new TGraph(n,p2,p0);
  TGraph *g2 = new TGraph(n,p0,p1);

  g0->SetMarkerColor(kRed);
  g0->SetMarkerSize(5);
  g1->SetMarkerColor(kRed);
  g1->SetMarkerSize(5);
  g2->SetMarkerColor(kRed);
  g2->SetMarkerSize(5);

  h0->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0->GetXaxis()->SetTitle("Umm2");
  h0->GetYaxis()->SetTitle("dm2");
  h0->SetStats(0);
  h1->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1->GetXaxis()->SetTitle("dm2");
  h1->GetYaxis()->SetTitle("Uee2");
  h1->SetStats(0);
  h2->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2->GetXaxis()->SetTitle("Uee2");
  h2->GetYaxis()->SetTitle("Umm2");
  h2->SetStats(0);
  h0z->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0z->GetXaxis()->SetTitle("Umm2");
  h0z->GetYaxis()->SetTitle("dm2");
  h0z->SetStats(0);
  h1z->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1z->GetXaxis()->SetTitle("dm2");
  h1z->GetYaxis()->SetTitle("Uee2");
  h1z->SetStats(0);
  h2z->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2z->GetXaxis()->SetTitle("Uee2");
  h2z->GetYaxis()->SetTitle("Umm2");
  h2z->SetStats(0);
/*
  h0L->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0L->GetXaxis()->SetTitle("Umm2");
  h0L->GetYaxis()->SetTitle("dm2");
  h0L->SetStats(0);
  h1L->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1L->GetXaxis()->SetTitle("dm2");
  h1L->GetYaxis()->SetTitle("Uee2");
  h1L->SetStats(0);
  h2L->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2L->GetXaxis()->SetTitle("Uee2");
  h2L->GetYaxis()->SetTitle("Umm2");
  h2L->SetStats(0);
*/
  h0d->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0d->GetXaxis()->SetTitle("Umm2");
  h0d->GetYaxis()->SetTitle("dm2");
  h0d->SetStats(0);
  h1d->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1d->GetXaxis()->SetTitle("dm2");
  h1d->GetYaxis()->SetTitle("Uee2");
  h1d->SetStats(0);
  h2d->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2d->GetXaxis()->SetTitle("Uee2");
  h2d->GetYaxis()->SetTitle("Umm2");
  h2d->SetStats(0);
  h0dz->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0dz->GetXaxis()->SetTitle("Umm2");
  h0dz->GetYaxis()->SetTitle("dm2");
  h0dz->SetStats(0);
  h1dz->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1dz->GetXaxis()->SetTitle("dm2");
  h1dz->GetYaxis()->SetTitle("Uee2");
  h1dz->SetStats(0);
  h2dz->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2dz->GetXaxis()->SetTitle("Uee2");
  h2dz->GetYaxis()->SetTitle("Umm2");
  h2dz->SetStats(0);
/*
  h0dL->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0dL->GetXaxis()->SetTitle("Umm2");
  h0dL->GetYaxis()->SetTitle("dm2");
  h0dL->SetStats(0);
  h1dL->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1dL->GetXaxis()->SetTitle("dm2");
  h1dL->GetYaxis()->SetTitle("Uee2");
  h1dL->SetStats(0);
  h2dL->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2dL->GetXaxis()->SetTitle("Uee2");
  h2dL->GetYaxis()->SetTitle("Umm2");
  h2dL->SetStats(0);
*/
  double chi2L_max = 1E5;
  double chi2_max = 2E4;
  double chi2z_max = 400;
  double diffL_max = 800;
  double diff_max = 300;
  double diffz_max = 30;

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  h0->SetMaximum(chi2L_max);
  h1->SetMaximum(chi2L_max);
  h2->SetMaximum(chi2L_max);
  h0d->SetMaximum(diffL_max);
  h1d->SetMaximum(diffL_max);
  h2d->SetMaximum(diffL_max);

  TCanvas *cchi2 = new TCanvas("cchi2","",1200,300);
  cchi2->Divide(3,1);
  cchi2->cd(1);
  gPad->SetLogz();
  h0->Draw("colz");
  g0->Draw("same C*");
  cchi2->cd(2);
  gPad->SetLogz();
  h1->Draw("colz");
  g1->Draw("same C*");
  cchi2->cd(3);
  gPad->SetLogz();
  h2->Draw("colz");
  g2->Draw("same C*");
  cchi2->SaveAs(Form("%s_chi2Surface_stat_%d%d.png",name,para,cutNu)); 
  
  TCanvas *cchi2_z = new TCanvas("cchi2_z","",1200,300);
  cchi2_z->Divide(3,1);
  cchi2_z->cd(1);
  h0z->Draw("colz");
  g0->Draw("same C*");
  cchi2_z->cd(2);
  h1z->Draw("colz");
  g1->Draw("same C*");
  cchi2_z->cd(3);
  h2z->Draw("colz");
  g2->Draw("same C*");
  cchi2_z->SaveAs(Form("%s_chi2Surface_stat_%d%d_zoom.png",name,para,cutNu)); 
 
  TCanvas *cdiff = new TCanvas("cdiff","",1200,300);
  cdiff->Divide(3,1);
  cdiff->cd(1);
  gPad->SetLogz();
  h0d->Draw("colz");
  g0->Draw("same C*");
  cdiff->cd(2);
  gPad->SetLogz();
  h1d->Draw("colz");
  g1->Draw("same C*");
  cdiff->cd(3);
  gPad->SetLogz();
  h2d->Draw("colz");
  g2->Draw("same C*");
  cdiff->SaveAs(Form("%s_chi2Diff_stat_%d%d.png",name,para,cutNu)); 
  
  TCanvas *cdiff_z = new TCanvas("cdiff_z","",1200,300);
  cdiff_z->Divide(3,1);
  cdiff_z->cd(1);
  h0dz->Draw("colz");
  g0->Draw("same C*");
  cdiff_z->cd(2);
  h1dz->Draw("colz");
  g1->Draw("same C*");
  cdiff_z->cd(3);
  h2dz->Draw("colz");
  g2->Draw("same C*");
  cdiff_z->SaveAs(Form("%s_chi2Diff_stat_%d%d_zoom.png",name,para,cutNu)); 

  TFile *out = new TFile(Form("%s_chi2_stat_%d%d.root",name,para,cutNu),"RECREATE");
  h0->Write();
  h1->Write();
  h2->Write();
  h0z->Write();
  h1z->Write();
  h2z->Write();
  h0d->Write();
  h1d->Write();
  h2d->Write();
  h0dz->Write();
  h1dz->Write();
  h2dz->Write();
  out->Close();

/*
  TCanvas *cchi2_0L = new TCanvas("cchi2_0L","",800,600);
  cchi2_0L->SetLogz(1);
  h0L->SetMaximum(chi2L_max);
  h0L->Draw("colz");
  g0->Draw("same C*");
  cchi2_0L->SaveAs(Form("%s_chi2Surface_stat_%d%d_0_Log.png",name,para,cutNu));
  TCanvas *cchi2_0 = new TCanvas("cchi2_0","",800,600);
  cchi2_0->SetLogz(0);
  h0->SetMaximum(chi2_max);
  h0->Draw("colz");
  g0->Draw("same C*");
  cchi2_0->SaveAs(Form("%s_chi2Surface_stat_%d%d_0.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_0z = new TCanvas("cchi2_0z","",800,600);
  cchi2_0z->SetLogz(0);
  h0z->SetMaximum(chi2z_max);
  h0z->Draw("colz");
  g0->Draw("same C*");
  cchi2_0z->SaveAs(Form("%s_chi2Surface_stat_%d%d_0_zoom.png",name,para,cutNu,cutEv));

  TCanvas *cchi2_1L = new TCanvas("cchi2_1L","",800,600);
  cchi2_1L->SetLogz(1);
  h1L->SetMaximum(chi2L_max);
  h1L->Draw("colz");
  g1->Draw("same C*");
  cchi2_1L->SaveAs(Form("%s_chi2Surface_stat_%d%d_1_Log.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_1 = new TCanvas("cchi2_1","",800,600);
  cchi2_1->SetLogz(0);
  h1->SetMaximum(chi2_max);
  h1->Draw("colz");
  g1->Draw("same C*");
  cchi2_1->SaveAs(Form("%s_chi2Surface_stat_%d%d_1.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_1z = new TCanvas("cchi2_1z","",800,600);
  cchi2_1z->SetLogz(0);
  h1z->SetMaximum(chi2z_max);
  h1z->Draw("colz");
  g1->Draw("same C*");
  cchi2_1z->SaveAs(Form("%s_chi2Surface_stat_%d%d_1_zoom.png",name,para,cutNu,cutEv));

  TCanvas *cchi2_2L = new TCanvas("cchi2_2L","",800,600);
  cchi2_2L->SetLogz(1);
  h2L->SetMaximum(chi2L_max);
  h2L->Draw("colz");
  g2->Draw("same C*");
  cchi2_2L->SaveAs(Form("%s_chi2Surface_stat_%d%d_2_Log.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_2 = new TCanvas("cchi2_2","",800,600);
  cchi2_2->SetLogz(0);
  h2->SetMaximum(chi2_max);
  h2->Draw("colz");
  g2->Draw("same C*");
  cchi2_2->SaveAs(Form("%s_chi2Surface_stat_%d%d_2.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_2z = new TCanvas("cchi2_2z","",800,600);
  cchi2_2z->SetLogz(0);
  h2z->SetMaximum(chi2z_max);
  h2z->Draw("colz");
  g2->Draw("same C*");
  cchi2_2z->SaveAs(Form("%s_chi2Surface_stat_%d%d_2_zoom.png",name,para,cutNu,cutEv));

  TCanvas *cdiff_0L = new TCanvas("cdiff_0L","",800,600);
  cdiff_0L->SetLogz(1);
  h0dL->SetMaximum(diffL_max);
  h0dL->Draw("colz");
  g0->Draw("same C*");
  cdiff_0L->SaveAs(Form("%s_chi2Diff_stat_%d%d_0_Log.png",name,para,cutNu,cutEv));
  TCanvas *cdiff_0 = new TCanvas("cdiff_0","",800,600);
  cdiff_0->SetLogz(0);
  h0d->SetMaximum(diff_max);
  h0d->Draw("colz");
  g0->Draw("same C*");
  cdiff_0->SaveAs(Form("%s_chi2Diff_stat_%d%d_0.png",name,para,cutNu,cutEv));
  TCanvas *cdiff_0z = new TCanvas("cdiff_0z","",800,600);
  cdiff_0z->SetLogz(0);
  h0dz->SetMaximum(diffz_max);
  h0dz->Draw("colz");
  g0->Draw("same C*");
  cdiff_0z->SaveAs(Form("%s_chi2Diff_stat_%d%d_0_zoom.png",name,para,cutNu,cutEv));

  TCanvas *cdiff_1L = new TCanvas("cdiff_1L","",800,600);
  cdiff_1L->SetLogz(1);
  h1dL->SetMaximum(diffL_max);
  h1dL->Draw("colz");
  g1->Draw("same C*");
  cdiff_1L->SaveAs(Form("%s_chi2Diff_stat_%d%d_1_Log.png",name,para,cutNu,cutEv));
  TCanvas *cdiff_1 = new TCanvas("cdiff_1","",800,600);
  cdiff_1->SetLogz(0);
  h1d->SetMaximum(diff_max);
  h1d->Draw("colz");
  g1->Draw("same C*");
  cdiff_1->SaveAs(Form("%s_chi2Diff_stat_%d%d_1.png",name,para,cutNu,cutEv));
  TCanvas *cdiff_1z = new TCanvas("cdiff_1z","",800,600);
  cdiff_1z->SetLogz(0);
  h1dz->SetMaximum(diffz_max);
  h1dz->Draw("colz");
  g1->Draw("same C*");
  cdiff_1z->SaveAs(Form("%s_chi2Diff_stat_%d%d_1_zoom.png",name,para,cutNu,cutEv));

  TCanvas *cdiff_2L = new TCanvas("cdiff_2L","",800,600);
  cdiff_2L->SetLogz(1);
  h2dL->SetMaximum(diffL_max);
  h2dL->Draw("colz");
  g2->Draw("same C*");
  cdiff_2L->SaveAs(Form("%s_chi2Diff_stat_%d%d_2_Log.png",name,para,cutNu));
  TCanvas *cdiff_2 = new TCanvas("cdiff_2","",800,600);
  cdiff_2->SetLogz(0);
  h2d->SetMaximum(diff_max);
  h2d->Draw("colz");
  g2->Draw("same C*");
  cdiff_2->SaveAs(Form("%s_chi2Diff_stat_%d%d_2.png",name,para,cutNu));
  TCanvas *cdiff_2z = new TCanvas("cdiff_2z","",800,600);
  cdiff_2z->SetLogz(0);
  h2dz->SetMaximum(diffz_max);
  h2dz->Draw("colz");
  g2->Draw("same C*");
  cdiff_2z->SaveAs(Form("%s_chi2Diff_stat_%d%d_2_zoom.png",name,para,cutNu));

  TFile *out = new TFile(Form("%s_chi2_stat_%d%d%d.root",name,para,cutNu,cutEv),"RECREATE");
  h0->Write();
  h1->Write();
  h2->Write();
  h0L->Write();
  h1L->Write();
  h2L->Write();
  h0z->Write();
  h1z->Write();
  h2z->Write();
  h0d->Write();
  h1d->Write();
  h2d->Write();
  h0dL->Write();
  h1dL->Write();
  h2dL->Write();
  h0dz->Write();
  h1dz->Write();
  h2dz->Write();
  out->Close();

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
*/

}

