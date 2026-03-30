#include "TemplateFitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
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
#include <TMatrixDEigen.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPaveText.h>

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

void TemplateFitter::setPara( char var[20], int nuCut, double fitPara_m[29][7], double fitPara_e[29][7] )
{
  name  = var;
  cutNu = nuCut;

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


double TemplateFitter::getAvgPme( double energy, double Ue42, double Um42, double Ut42, double dm2, double ft[7] )
{
  if(dm2==0.) return 0.;
  else{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * Ue42 * Um42;

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

  return prob;}
}

double TemplateFitter::getAvgPmt( double energy, double Ue42, double Um42, double Ut42, double dm2, double ft[7] )
{
  if(dm2==0.) return 0.;
  else{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * Ut42 * Um42;

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

  return prob;}
}

double TemplateFitter::getAvgPee( double energy, double Ue42, double Um42, double Ut42, double dm2, double ft[7] )
{
  if(dm2==0.) return 1.;
  else{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * Ue42 * (1.-Ue42);

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
  prob_u = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L2;
  prob_l = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);

  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);

  return prob;}
}

double TemplateFitter::getAvgPmm( double energy, double Ue42, double Um42, double Ut42, double dm2, double ft[7] )
{
  if(dm2==0.) return 1.;
  else{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * Um42 * (1.-Um42);

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
  prob_u = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L2;
  prob_l = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);

  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);

  return prob;}
}




TMatrixD covmx(nbins,nbins);
TMatrixD invmx(nbins,nbins);
TMatrixD flmx(nbins,nbins);
TMatrixD sigmx(nbins,nbins);
TMatrixD statmx(nbins,nbins);
TMatrixD detmx(nbins,nbins);

void TemplateFitter::setCovmtr( double flmx_bct[nbins+1][nbins+1], double sigmx_bct[nbins+1][nbins+1], double detmx_bct[nbins+1][nbins+1] )
{
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      flmx[i][j]  = flmx_bct[i][j];
      sigmx[i][j] = sigmx_bct[i][j];
      detmx[i][j] = detmx_bct[i][j];
    }
  }
}



void TemplateFitter::getTarget( double *par )
{
  CCe_tgt->Reset();
  CCm_tgt->Reset();
  nue_tgt->Reset();

  //std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << par[3] << "\n";

  // Create a histogram "temp" from the templates
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

    if(i<240)                iL = (int)i/10;
    else if(i>=240 && i<400) iL = (int) 24+(i-240)/40;
    else                     iL = 28;
    
    for(int k=0; k<7; k++) {
    ft_m[k] = fitP_m[iL][k];
    ft_e[k] = fitP_e[iL][k];
    }

    double Pme, Pem, Pmm, Pee, Pmt;
    double me = 0;
    double em = 0;
    double ee = 0;
    double mm = 0;
    double mt = 0;
    
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      me = me + getAvgPme(e, par[0], par[1], par[2], par[3], ft_m);
      mt = mt + getAvgPmt(e, par[0], par[1], par[2], par[3], ft_m);
      em = em + getAvgPme(e, par[0], par[1], par[2], par[3], ft_e);
      ee = ee + getAvgPee(e, par[0], par[1], par[2], par[3], ft_e);
      mm = mm + getAvgPmm(e, par[0], par[1], par[2], par[3], ft_m);
    }
    
    Pme = me/1001.0;
    Pem = em/1001.0;
    Pee  = ee/1001.0;
    Pmm  = mm/1001.0;
    Pmt  = mt/1001.0;

    CC_tp_me->Add(CC_nc_m_templates[i], Pme);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pem);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pme);
    nue_tp_mm->Add(nue_m_templates[i], (Pmm+Pmt));
    nue_tp_em->Add(nue_w_e_templates[i], Pem);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }
  
  for(int i=0; i<nbins_CC; i++){
    if(i>=44) {
      CC_tp_mm->SetBinContent(i+1, 0.);
      CC_tp_em->SetBinContent(i+1, 0.);
    }
  }
  

  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me);    CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);
    

  CCe_tgt->Add(CC_tp_e);
  CCm_tgt->Add(CC_tp_m);
  nue_tgt->Add(nue_tp);

  // combining all 3 samples into 1 histograms
  TH1D * tgt = new TH1D("tgt","",nbins,0,nbins);

  for( int bx = 0; bx < nbins; bx++ ) {
    if(bx<nbins_CC)                   tgt->SetBinContent(bx+1, CCm_tgt->GetBinContent(bx+1));
    if(bx>=nbins_CC && bx<2*nbins_CC) tgt->SetBinContent(bx+1, CCe_tgt->GetBinContent(bx-nbins_CC+1));
    if(bx>=2*nbins_CC)                tgt->SetBinContent(bx+1, nue_tgt->GetBinContent(bx-2*nbins_CC+1));

    for( int by = 0; by < nbins; by++ ) {
      statmx[bx][by] = 0.;
      if(bx==by){
        if(bx<nbins_CC)                   statmx[bx][bx] = CCm_tgt->GetBinContent(bx+1);
        if(bx>=nbins_CC && bx<2*nbins_CC) statmx[bx][bx] = CCe_tgt->GetBinContent(bx-nbins_CC+1);
        if(bx>=2*nbins_CC)                statmx[bx][bx] = nue_tgt->GetBinContent(bx-2*nbins_CC+1);
      }
    }
  }

  //std::cout << statmx[1][1] << "\t" << detmx[1][1] << "\n";

  for(int bx=0; bx<nbins; bx++){
    for(int by=0; by<nbins; by++){
      detmx[bx][by] = detmx[bx][by] * tgt->GetBinContent(bx+1) * tgt->GetBinContent(by+1);
    }
  }
  //std::cout << statmx[1][1] << "\t" << flmx[1][1] << "\t" << sigmx[1][1] << "\t" << detmx[1][1] << "\n";

  covmx = statmx + flmx + sigmx + detmx;

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



double TemplateFitter::bfChi2( double Ue42, double Um42, double Ut42, double dm2 )
{
  double par[4];
  par[0] = Ue42;
  par[1] = Um42;
  par[2] = Ut42;
  par[3] = dm2;

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

    if(i<240)                iL = (int)i/10;
    else if(i>=240 && i<400) iL = (int) 24+(i-240)/40;
    else                     iL = 28;
    
    for(int k=0; k<7; k++) {
    ft_m[k] = fitP_m[iL][k];
    ft_e[k] = fitP_e[iL][k];
    }

    double Pme, Pem, Pmm, Pee, Pmt;
    double me = 0;
    double em = 0;
    double ee = 0;
    double mm = 0;
    double mt = 0;
    
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      me = me + getAvgPme(e, par[0], par[1], par[2], par[3], ft_m);
      mt = mt + getAvgPmt(e, par[0], par[1], par[2], par[3], ft_m);
      em = em + getAvgPme(e, par[0], par[1], par[2], par[3], ft_e);
      ee = ee + getAvgPee(e, par[0], par[1], par[2], par[3], ft_e);
      mm = mm + getAvgPmm(e, par[0], par[1], par[2], par[3], ft_m);
    }
    
    Pme = me/1001.0;
    Pem = em/1001.0;
    Pee  = ee/1001.0;
    Pmm  = mm/1001.0;
    Pmt  = mt/1001.0;

    CC_tp_me->Add(CC_nc_m_templates[i], Pme);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pem);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pme);
    nue_tp_mm->Add(nue_m_templates[i], Pmm + Pmt);
    nue_tp_em->Add(nue_w_e_templates[i], Pem);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }
  
  for(int i=0; i<nbins_CC; i++){
    if(i>=44) {
      CC_tp_mm->SetBinContent(i+1, 0.);
      CC_tp_em->SetBinContent(i+1, 0.);
    }
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

  for( int bx = 0; bx < nbins; bx++ ) {
    if(bx<nbins_CC)                   temp[bx][0] = CC_tp_m->GetBinContent(bx+1);
    if(bx>=nbins_CC && bx<2*nbins_CC) temp[bx][0] = CC_tp_e->GetBinContent(bx-nbins_CC+1);
    if(bx>=2*nbins_CC)                temp[bx][0] = nue_tp->GetBinContent(bx-2*nbins_CC+1);

    if(bx<nbins_CC)                   target[bx][0] = CCm_tgt->GetBinContent(bx+1);
    if(bx>=nbins_CC && bx<2*nbins_CC) target[bx][0] = CCe_tgt->GetBinContent(bx-nbins_CC+1);
    if(bx>=2*nbins_CC)                target[bx][0] = nue_tgt->GetBinContent(bx-2*nbins_CC+1);

    //std::cout << bx << "\t" << target[bx][0] << "\t" << temp[bx][0] << "\n";

  }

  TMatrixD diff = temp - target;
  TMatrixD diffT(TMatrixD::kTransposed, diff);

  TMatrixD chi2 = diffT*invmx*diff;

  //std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << par[3] << "\t" << chi2[0][0] << "\n";
  return chi2[0][0];
}

