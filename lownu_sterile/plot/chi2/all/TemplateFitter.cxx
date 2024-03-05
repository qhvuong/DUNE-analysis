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
#include "TMatrixDEigen.h"
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

  CCm_tgt = (TH1D*)CC_m_templates[0]->Clone();
  CCe_tgt = (TH1D*)CC_m_templates[0]->Clone();
  nue_tgt = (TH1D*)nue_m_templates[0]->Clone();
  CCm_tgt->Reset();
  CCe_tgt->Reset();
  nue_tgt->Reset();
}

void TemplateFitter::setEnergyBins( double bins[nbins_Ev+1] )
{
  for( int i = 0; i < nbins_Ev+1; ++ i ) {m_energy_bins[i] = bins[i];}
}

void TemplateFitter::setPara( char var[20], int nuCut, double par[3], double fitPara_m[29][7], double fitPara_e[29][7] )
{
  name  = var;
  cutNu = nuCut;

  b0 = par[0];
  b1 = par[1];
  b2 = par[2];

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


TMatrixD covmx(nbins,nbins);
TMatrixD invmx(nbins,nbins);
TMatrixD flmx(nbins,nbins);
TMatrixD sigmx(nbins,nbins);
TMatrixD statmx(nbins,nbins);

void TemplateFitter::setCovmtr( double flmx_bct[nbins+1][nbins+1], double sigmx_bct[nbins+1][nbins+1] )
{
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      flmx[i][j]  = flmx_bct[i][j];
      sigmx[i][j] = sigmx_bct[i][j];
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
  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < nbins_Ev; ++i ) {
    if(i<240)                iL = (int)i/10;
    else if(i>=240 && i<400) iL = (int) 24+(i-240)/40;
    else                     iL = 28;
    
    for(int k=0; k<7; k++) {
    ft_m[k] = fitP_m[iL][k];
    ft_e[k] = fitP_e[iL][k];
    }

    double Pmue, Pemu, Pmm, Pee;
    if(par[2]==0.){
    Pemu = Pmue = 0.;
    Pmm = Pee = 1.;
    }
    else{
      if(par[0]==0. && par[1]==0.){
	Pemu = Pmue = 0.;
	Pmm = Pee = 1.;
      }

      if(par[0]==0. || par[1]==0.){
        Pmue = Pemu = 0.;

        if(par[0]==0.){
        Pee = 1.;
	double mm = 0;
	for(int j = 0; j<1001; j++){
	  double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
	  mm  = mm  + getAvgPmm(e, par[0], par[1], par[2], ft_m);
	}
	Pmm  = mm/1001.0;
        }
        
        if(par[1]==0.){
        Pmm = 1.;
	double ee = 0;
	for(int j = 0; j<1001; j++){
	  double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
	  ee  = ee  + getAvgPee(e, par[0], par[1], par[2], ft_e);
	}
	Pee  = ee/1001.0;
        }
      }

      else{
	double mue = 0;
	double emu = 0;
	double ee = 0;
	double mm = 0;
	
	for(int j = 0; j<1001; j++){
	  double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
	  mue = mue + getAvgPmue(e, par[0], par[1], par[2], ft_m);
	  emu = emu + getAvgPmue(e, par[0], par[1], par[2], ft_e);
	  ee  = ee  + getAvgPee(e, par[0], par[1], par[2], ft_e);
	  mm  = mm  + getAvgPmm(e, par[0], par[1], par[2], ft_m);
	}
	
	Pmue = mue/1001.0;
	Pemu = emu/1001.0;
	Pee  = ee/1001.0;
	Pmm  = mm/1001.0;
      }
    } 
   
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

  for( int bx = 0; bx < nbins; bx++ ) {
    for( int by = 0; by < nbins; by++ ) {
      statmx[bx][by] = 0.;
      if(bx==by){
        if(bx<nbins_CC)                   statmx[bx][bx] = CCm_tgt->GetBinContent(bx+1);
        if(bx>=nbins_CC && bx<2*nbins_CC) statmx[bx][bx] = CCe_tgt->GetBinContent(bx-nbins_CC+1);
        if(bx>=2*nbins_CC)                statmx[bx][bx] = nue_tgt->GetBinContent(bx-2*nbins_CC+1);
      }
    }
  }

  covmx = statmx + flmx + sigmx;

  TDecompSVD svd(covmx);

  if(!svd.Decompose()) {
    std::cerr << "SVD Decomposition failed!" << std::endl;}

  TVectorD singularValues = svd.GetSig();

  TMatrixDEigen eigen(covmx);
  TVectorD eiV = eigen.GetEigenValuesRe();

  TMatrixD Sinv(nbins, nbins);
  Sinv.Zero();

  for(int i=0; i<nbins; i++){
    if(eiV[i] < 0) continue;
    Sinv(i, i) = 1.0 / singularValues[i]; // Invert non-zero singular values
  }

  TMatrixD U = svd.GetU();
  TMatrixD V = svd.GetV();
  TMatrixD Ut(TMatrixD::kTransposed, U);

  invmx = V*Sinv*Ut;

}


// function whose return Minuit mimizes, must take const double* and return double
double TemplateFitter::getChi2( double * par )
{
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


  int iL;

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < nbins_Ev; ++i ) {
/*
  if(par[0]==0. && par[1]==0. && par[2]==0.){
    CC_tp_me->Add(CC_nc_m_templates[i], 0.0);
    CC_tp_ee->Add(CC_e_templates[i], 1.0);
    
    CC_tp_em->Add(CC_e_templates[i], 0.0);
    CC_tp_mm->Add(CC_m_templates[i], 1.0);
    
    nue_tp_me->Add(nue_w_m_templates[i], 0.0);
    nue_tp_mm->Add(nue_m_templates[i], 1.0);
    nue_tp_em->Add(nue_w_e_templates[i], 0.0);
    nue_tp_ee->Add(nue_e_templates[i], 1.0);
  } 
  
  else{
    double mue = 0;
    double emu = 0;
    double ee = 0;
    double mm = 0;
    double L  = 0.5;
    
    if(i<240)                iL = (int)i/10;
    else if(i>=240 && i<400) iL = (int) 24+(i-240)/40;
    else                     iL = 28;
    
    for(int k=0; k<7; k++) {
      ft_m[k] = fitP_m[iL][k];
      ft_e[k] = fitP_e[iL][k];
    } 
    
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getAvgPmue(e, par[0], par[1], par[2], ft_m);
      emu = emu + getAvgPmue(e, par[0], par[1], par[2], ft_e);
      ee  = ee  + getAvgPee(e, par[0], par[1], par[2], ft_e);
      mm  = mm  + getAvgPmm(e, par[0], par[1], par[2], ft_m);
    } 
    

    double Pmue = mue/1001.0;
    double Pemu = emu/1001.0;
    double Pee  = ee/1001.0;
    double Pmm  = mm/1001.0;
*/  
    if(i<240)                iL = (int)i/10;
    else if(i>=240 && i<400) iL = (int) 24+(i-240)/40;
    else                     iL = 28;
    
    for(int k=0; k<7; k++) {
    ft_m[k] = fitP_m[iL][k];
    ft_e[k] = fitP_e[iL][k];
    }

    double Pmue, Pemu, Pmm, Pee;
    if(par[2]==0.){
    Pemu = Pmue = 0.;
    Pmm = Pee = 1.;
    }
    else{
      if(par[0]==0. && par[1]==0.){
	Pemu = Pmue = 0.;
	Pmm = Pee = 1.;
      }

      if(par[0]==0. || par[1]==0.){
        Pmue = Pemu = 0.;

        if(par[0]==0.){
        Pee = 1.;
	double mm = 0;
	for(int j = 0; j<1001; j++){
	  double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
	  mm  = mm  + getAvgPmm(e, par[0], par[1], par[2], ft_m);
	}
	Pmm  = mm/1001.0;
        }
        
        if(par[1]==0.){
        Pmm = 1.;
	double ee = 0;
	for(int j = 0; j<1001; j++){
	  double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
	  ee  = ee  + getAvgPee(e, par[0], par[1], par[2], ft_e);
	}
	Pee  = ee/1001.0;
        }
      }

      else{
	double mue = 0;
	double emu = 0;
	double ee = 0;
	double mm = 0;
	
	for(int j = 0; j<1001; j++){
	  double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
	  mue = mue + getAvgPmue(e, par[0], par[1], par[2], ft_m);
	  emu = emu + getAvgPmue(e, par[0], par[1], par[2], ft_e);
	  ee  = ee  + getAvgPee(e, par[0], par[1], par[2], ft_e);
	  mm  = mm  + getAvgPmm(e, par[0], par[1], par[2], ft_m);
	}
	
	Pmue = mue/1001.0;
	Pemu = emu/1001.0;
	Pee  = ee/1001.0;
	Pmm  = mm/1001.0;
      }
    }  

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
  TMatrixD target(nbins, 1);
  TMatrixD temp(nbins, 1);
  TMatrixD unc(nbins, nbins);
  TMatrixD cov(nbins, nbins);

  for( int bx = 0; bx < nbins; bx++ ) {
    if(bx<nbins_CC)                   temp[bx][0] = CC_tp_m->GetBinContent(bx+1);
    if(bx>=nbins_CC && bx<2*nbins_CC) temp[bx][0] = CC_tp_e->GetBinContent(bx-nbins_CC+1);
    if(bx>=2*nbins_CC)                temp[bx][0] = nue_tp->GetBinContent(bx-2*nbins_CC+1);

    if(bx<nbins_CC)                   target[bx][0] = CCm_tgt->GetBinContent(bx+1);
    if(bx>=nbins_CC && bx<2*nbins_CC) target[bx][0] = CCe_tgt->GetBinContent(bx-nbins_CC+1);
    if(bx>=2*nbins_CC)                target[bx][0] = nue_tgt->GetBinContent(bx-2*nbins_CC+1);

  }
  TMatrixD diff = temp - target;
  TMatrixD diffT(TMatrixD::kTransposed, diff);

  TMatrixD chi2 = diffT*invmx*diff;

  //std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2[0][0] << "\n";
  return chi2[0][0];

}


void TemplateFitter::Draw()
{
  int N = 100;

  double p0[1], p1[1], p2[1];
  p0[0] = b0;
  p1[0] = b1;
  p2[0] = b2;

  double Umin = 1e-4, Umax = 0.5;
  double dmMin = 1e-2, dmMax = 100.0;
  double log_min = TMath::Log10(Umin);
  double log_max = TMath::Log10(Umax);
  double binWidth = (log_max-log_min)/N;
  double mlog_min = TMath::Log10(dmMin);
  double mlog_max = TMath::Log10(dmMax);
  double mbinWidth = (mlog_max-mlog_min)/N;

  std::vector<double> binEdges(N + 1, 0);
  std::vector<double> mbinEdges(N + 1, 0);
  for (int i = 0; i <= N; ++i) {
    binEdges[i] = TMath::Power(10, log_min + i * binWidth);
    mbinEdges[i] = TMath::Power(10, mlog_min + i * mbinWidth);
  }

  TH2D *h0 = new TH2D("h0","",N,&binEdges[0], N,&mbinEdges[0]);
  TH2D *h1 = new TH2D("h1","",N,&binEdges[0], N,&mbinEdges[0]);
  TH2D *h2 = new TH2D("h2","",N,&binEdges[0], N,&binEdges[0]);
  
  TH2D *h0d = new TH2D("h0d","",N,&binEdges[0], N,&mbinEdges[0]);
  TH2D *h1d = new TH2D("h1d","",N,&binEdges[0], N,&mbinEdges[0]);
  TH2D *h2d = new TH2D("h2d","",N,&binEdges[0], N,&binEdges[0]);
  
  double par[3], parL[3], bin[3];
  par[0] = b0;
  par[1] = b1;
  par[2] = b2;
  double chi2 = 0.0;
  double chi2L = 0.0;

  double nu=0.3, Ev=3.0;

  double chi2t = getChi2(par);
  double diff, diffz;

  par[0]  = b0;
  for(int i=1; i<=N; i++) {
    par[1]  = h0->GetXaxis()->GetBinCenter(i);
    for(int j=1; j<=N; j++) {
      par[2]  = h0->GetYaxis()->GetBinCenter(j);
      chi2  = getChi2(par);
      diff  = sqrt(std::fabs(chi2 - chi2t));
      if( j==N ) std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\t" << diff << "\t" << i*j*100./(3*N*N) << "%" << "\n";
      h0->Fill(par[1], par[2], chi2);
      h0d->Fill(par[1], par[2], diff);
    }
  }

  par[1]  = b1;
  for(int i=1; i<=N; i++) {
    par[0]  = h1->GetXaxis()->GetBinCenter(i);
    for(int j=1; j<=N; j++) {
      par[2]  = h1->GetYaxis()->GetBinCenter(j);
      chi2  = getChi2(par);
      diff  = sqrt(std::fabs(chi2 - chi2t));
      if( j==N ) std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\t" << diff << "\t" << (N*N+i*j)*100./(3*N*N) << "%" << "\n";
      h1->Fill(par[0], par[2], chi2);
      h1d->Fill(par[0], par[2], diff);
    }
  }

  par[2]  = b2;
  for(int i=1; i<=N; i++) {
    par[1]  = h2->GetXaxis()->GetBinCenter(i);
    for(int j=1; j<=N; j++) {
      par[0]  = h2->GetYaxis()->GetBinCenter(j);
      chi2  = getChi2(par);
      diff = sqrt(std::fabs(chi2 - chi2t));
      if( j==N ) std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\t" << diff << "\t" << (2*N*N+i*j)*100./(3*N*N) << "%" << "\n";
      h2->Fill(par[1], par[0], chi2);
      h2d->Fill(par[1], par[0], diff);
    }
  }

  int n = 1;
  TGraph *g0 = new TGraph(n,p1,p2);
  TGraph *g1 = new TGraph(n,p0,p2);
  TGraph *g2 = new TGraph(n,p1,p0);

  g0->SetMarkerColor(kGreen);
  g0->SetMarkerSize(2);
  g1->SetMarkerColor(kGreen);
  g1->SetMarkerSize(2);
  g2->SetMarkerColor(kGreen);
  g2->SetMarkerSize(2);

  h0->SetTitle(Form("#chi^{2} Surface (U_{e4}^{2} = %.2f)",b0));
  h0->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  h0->GetYaxis()->SetTitle("#Deltam^{2}");
  h0->SetStats(0);
  h1->SetTitle(Form("#chi^{2} Surface (U_{#mu4}^{2} = %.2f)",b1));
  h1->GetXaxis()->SetTitle("U_{e4}^{2}");
  h1->GetYaxis()->SetTitle("#Deltam^{2}");
  h1->SetStats(0);
  h2->SetTitle(Form("#chi^{2} Surface (#Deltam^{2} = %.1f)",b2));
  h2->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  h2->GetYaxis()->SetTitle("U_{e4}^{2}");
  h2->SetStats(0);

  h0d->SetTitle(Form("#sqrt{#Delta#chi^{2}} Surface (U_{e4}^{2} = %.2f)",b0));
  h0d->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  h0d->GetYaxis()->SetTitle("#Deltam^{2}");
  h0d->SetStats(0);
  h1d->SetTitle(Form("#sqrt{#Delta#chi^{2}} Surface (U_{#mu4}^{2} = %.2f)",b1));
  h1d->GetXaxis()->SetTitle("U_{e4}^{2}");
  h1d->GetYaxis()->SetTitle("#Deltam^{2}");
  h1d->SetStats(0);
  h2d->SetTitle(Form("#sqrt{#Delta#chi^{2}} Surface (#Deltam^{2} = %.1f)",b2));
  h2d->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  h2d->GetYaxis()->SetTitle("U_{e4}^{2}");
  h2d->SetStats(0);
  
  //double c2Min = 0.1, c2Max = 1e5;
  //double d2Min = 1.0, d2Max = 1e3;

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);
/*
  h0->SetMinimum(c2Min);
  h1->SetMinimum(c2Min);
  h2->SetMinimum(c2Min);
  h0->SetMaximum(c2Max);
  h1->SetMaximum(c2Max);
  h2->SetMaximum(c2Max);

  h0d->SetMinimum(d2Min);
  h1d->SetMinimum(d2Min);
  h2d->SetMinimum(d2Min);
  h0d->SetMaximum(d2Max);
  h1d->SetMaximum(d2Max);
  h2d->SetMaximum(d2Max);
*/

  const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu/chi2_surface";

  TCanvas *cchi2 = new TCanvas("cchi2","",1200,300);
  cchi2->Divide(3,1);
  cchi2->cd(1);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();
  h0->Draw("colz");
  g0->Draw("same C*");
  cchi2->cd(2);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();
  h1->Draw("colz");
  g1->Draw("same C*");
  cchi2->cd(3);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();
  h2->Draw("colz");
  g2->Draw("same C*");
  cchi2->SaveAs(Form("%s/chi2Surface_all%d_%d.png",data_path,cutNu,N)); 
  
  TCanvas *cdiff = new TCanvas("cdiff","",1200,300);
  cdiff->Divide(3,1);
  cdiff->cd(1);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();
  h0d->Draw("colz");
  g0->Draw("same C*");
  cdiff->cd(2);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();
  h1d->Draw("colz");
  g1->Draw("same C*");
  cdiff->cd(3);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();
  h2d->Draw("colz");
  g2->Draw("same C*");
  cdiff->SaveAs(Form("%s/chi2Diff_all%d_%d.png",data_path,cutNu,N)); 
  

  TFile *out = new TFile(Form("%s/chi2_all%d_%d.root",data_path,cutNu,N),"RECREATE");
  h0->Write();
  h1->Write();
  h2->Write();
  h0d->Write();
  h1d->Write();
  h2d->Write();
  out->Close();
}

