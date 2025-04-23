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
#include "TMatrixT.h"
#include "TDecompSVD.h"
#include "TGraph.h"
#include "TLegend.h"
#include <TRandom3.h>
#include <TMatrixDEigen.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPaveText.h>

using namespace std;

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
  CCm_tgt_thr = (TH1D*)CCm_tgt->Clone();
  CCe_tgt_thr = (TH1D*)CCe_tgt->Clone();
  nue_tgt_thr = (TH1D*)nue_tgt->Clone();
}

void TemplateFitter::setEnergyBins( double bins[nbins_Ev+1] )
{
  for( int i = 0; i < nbins_Ev+1; ++ i ) {m_energy_bins[i] = bins[i];}
}


void TemplateFitter::setPara( char var[20], int nuCut, double fitPara_m[29][7], double fitPara_e[29][7], int universe )
{
  name  = var;
  cutNu = nuCut;
  uni = universe;

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
  double del = 1.27*L*dm2/energy;
  double s2mue2 = 4 * Uee2 * Umm2;
  double prob = s2mue2  * pow(sin(del),2);
  return prob;
}
double TemplateFitter::getPee( double energy, double Uee2, double Umm2, double dm2, double L )
{
  double del = 1.27*L*dm2/energy;
  double s2ee2 = 4 * Uee2 * (1 - Uee2);
  double prob = 1.0 - s2ee2  * pow(sin(del),2);
  return prob;
}
double TemplateFitter::getPmm( double energy, double Uee2, double Umm2, double dm2, double L )
{
  double del = 1.27*L*dm2/energy;
  double s2mm2 = 4 * Umm2 * (1 - Umm2);
  double prob = 1.0 - s2mm2  * pow(sin(del),2);
  return prob;
}

double TemplateFitter::getAvgPmue( double energy, double Uee2, double Umm2, double dm2, double ft[7] )
{
  if(dm2==0.) return 0.;

  else{

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

  return prob;}
}

double TemplateFitter::getAvgPee( double energy, double Uee2, double Umm2, double dm2, double ft[7] )
{
  if(dm2==0.) return 1.;

  else{

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
  prob_u = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L2;
  prob_l = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);

  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);

  return prob;}
}

double TemplateFitter::getAvgPmm( double energy, double Uee2, double Umm2, double dm2, double ft[7] )
{
  if(dm2==0.) return 1.;

  else{

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
  prob_u = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L2;
  prob_l = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);

  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);

  return prob;}
}


TMatrixD scales(1,nbins);
TMatrixD covmx(nbins,nbins);
TMatrixD invmx(nbins,nbins);
TMatrixD flmx_fr(nbins,nbins);
TMatrixD flmx(nbins,nbins);
TMatrixD sigmx_fr(nbins,nbins);
TMatrixD sigmx(nbins,nbins);
TMatrixD statmx(nbins,nbins);

void TemplateFitter::setCovmtr( double flmx_bct[nbins+1][nbins+1], double sigmx_bct[nbins+1][nbins+1], double wgt_bct[nbins+1] )
{
  for(int i=0; i<nbins; i++) {
    scales[0][i]  = wgt_bct[i];
    //cout << i << "\t" << scales[0][i] << "\n";
    for(int j=0; j<nbins; j++) {
      flmx_fr[i][j]  = flmx_bct[i][j];
      sigmx_fr[i][j] = sigmx_bct[i][j];
    }
  }
}



double TemplateFitter::getTarget( double *par )
{
  CCe_tgt->Reset();
  CCm_tgt->Reset();
  nue_tgt->Reset();
  CCe_tgt_thr->Reset();
  CCm_tgt_thr->Reset();
  nue_tgt_thr->Reset();

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

  int iL;

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < nbins_Ev; ++i ) {
    CC_tp_me->Add(CC_nc_m_templates[i], 0.0);
    CC_tp_ee->Add(CC_e_templates[i], 1.0);

    CC_tp_em->Add(CC_e_templates[i], 0.0);
    CC_tp_mm->Add(CC_m_templates[i], 1.0);

    nue_tp_me->Add(nue_w_m_templates[i], 0.0);
    nue_tp_mm->Add(nue_m_templates[i], 1.0);
    nue_tp_em->Add(nue_w_e_templates[i], 0.0);
    nue_tp_ee->Add(nue_e_templates[i], 1.0);
  }

  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em); CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);

  CCe_tgt->Add(CC_tp_e);
  CCm_tgt->Add(CC_tp_m);
  nue_tgt->Add(nue_tp);

  //std::cout << scales.GetNrows() << "\n";
  //std::cout << scales.GetNcols() << "\n";

  for(int i=0; i<nbins; i++){
    if(i<nbins_CC){
      int b = i;
      double bct = CC_tp_m->GetBinContent(b+1);
      double wgt = scales[0][i];
      CCm_tgt_thr->SetBinContent(b+1, bct*(1.+wgt));
      //cout << i << "\t" << bct << "\t" << wgt << "\t" << CCm_tgt->GetBinContent(b+1) << "\t" << CCm_tgt_thr->GetBinContent(b+1) << "\n";
    }

    if(i>=nbins_CC && i<2*nbins_CC){
      int b = i-nbins_CC;
      double bct = CC_tp_e->GetBinContent(b+1);
      double wgt = scales[0][i];
      CCe_tgt_thr->SetBinContent(b+1, bct*(1.+wgt));
      //cout << i << "\t" << bct << "\t" << wgt << "\t" << CCe_tgt->GetBinContent(b+1) << "\t" << CCe_tgt_thr->GetBinContent(b+1) << "\n";
    }

    if(i>=2*nbins_CC){
      int b = i-2*nbins_CC;
      double bct = nue_tp->GetBinContent(b+1);
      double wgt = scales[0][i];
      nue_tgt_thr->SetBinContent(b+1, bct*(1.+wgt));
      //cout << i << "\t" << bct << "\t" << wgt << "\t" << nue_tgt->GetBinContent(b+1) << "\t" << nue_tgt_thr->GetBinContent(b+1) << "\n";
    }
  }

  double throw_area = CCm_tgt_thr->Integral() + CCe_tgt_thr->Integral() + nue_tgt_thr->Integral();
  double area       = CCm_tgt->Integral()     + CCe_tgt->Integral()     + nue_tgt->Integral();

  double ratio = throw_area/area;


  for( int bx = 0; bx < nbins; bx++ ) {
    for( int by = 0; by < nbins; by++ ) {
      statmx[bx][by] = 0.;
      if(bx==by){
        if(bx<nbins_CC)                   statmx[bx][bx] = CCm_tgt_thr->GetBinContent(bx+1);
        if(bx>=nbins_CC && bx<2*nbins_CC) statmx[bx][bx] = CCe_tgt_thr->GetBinContent(bx-nbins_CC+1);
        if(bx>=2*nbins_CC)                statmx[bx][bx] = nue_tgt_thr->GetBinContent(bx-2*nbins_CC+1);
      }
    }
  }

  for(int bx=0; bx<nbins; bx++)  {
  for(int by=0; by<nbins; by++)  {
    if(bx < nbins_CC) {
      if(by<nbins_CC)                 flmx[bx][by] = flmx_fr[bx][by] * CCm_tgt_thr->GetBinContent(bx+1) * CCm_tgt_thr->GetBinContent(by+1); 
      if(nbins_CC <= by < 2*nbins_CC) flmx[bx][by] = flmx_fr[bx][by] * CCm_tgt_thr->GetBinContent(bx+1) * CCe_tgt_thr->GetBinContent(by-nbins_CC+1); 
      if(by >= 2*nbins_CC)            flmx[bx][by] = flmx_fr[bx][by] * CCm_tgt_thr->GetBinContent(bx+1) * nue_tgt_thr->GetBinContent(by-2*nbins_CC+1); 
    }
    if(nbins_CC <= bx < 2*nbins_CC) {
      if(by<nbins_CC)                 flmx[bx][by] = flmx_fr[bx][by] * CCe_tgt_thr->GetBinContent(bx-nbins_CC+1) * CCm_tgt_thr->GetBinContent(by+1); 
      if(nbins_CC <= by < 2*nbins_CC) flmx[bx][by] = flmx_fr[bx][by] * CCe_tgt_thr->GetBinContent(bx-nbins_CC+1) * CCe_tgt_thr->GetBinContent(by-nbins_CC+1); 
      if(by >= 2*nbins_CC)            flmx[bx][by] = flmx_fr[bx][by] * CCe_tgt_thr->GetBinContent(bx-nbins_CC+1) * nue_tgt_thr->GetBinContent(by-2*nbins_CC+1); 
    }
    if(bx >= 2*nbins_CC) {
      if(by<nbins_CC)                 flmx[bx][by] = flmx_fr[bx][by] * nue_tgt_thr->GetBinContent(bx-2*nbins_CC+1) * CCm_tgt_thr->GetBinContent(by+1); 
      if(nbins_CC <= by < 2*nbins_CC) flmx[bx][by] = flmx_fr[bx][by] * nue_tgt_thr->GetBinContent(bx-2*nbins_CC+1) * CCe_tgt_thr->GetBinContent(by-nbins_CC+1); 
      if(by >= 2*nbins_CC)            flmx[bx][by] = flmx_fr[bx][by] * nue_tgt_thr->GetBinContent(bx-2*nbins_CC+1) * nue_tgt_thr->GetBinContent(by-2*nbins_CC+1); 
    }
  }}


  //covmx = statmx + flmx + sigmx;
  //covmx = statmx + flmx;
  covmx = statmx + flmx;



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

  //invmx.Print();

  return ratio;
  

}








// function whose return Minuit mimizes, must take const double* and return double
double TemplateFitter::getChi2( const double * par )
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
    if(i<240)                iL = (int)i/10;
    else if(i>=240 && i<400) iL = (int) 24+(i-240)/40;
    else                     iL = 28;

    for(int k=0; k<7; k++) {
    ft_m[k] = fitP_m[iL][k];
    ft_e[k] = fitP_e[iL][k];
    }

    double Pmue, Pemu, Pmm, Pee;
/*
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
*/
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

    std::cout << i << "\t" << Pmue << "\t" << Pee << "\t" << Pmm << "\n";

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

    if(bx<nbins_CC)                   target[bx][0] = CCm_tgt_thr->GetBinContent(bx+1);
    if(bx>=nbins_CC && bx<2*nbins_CC) target[bx][0] = CCe_tgt_thr->GetBinContent(bx-nbins_CC+1);
    if(bx>=2*nbins_CC)                target[bx][0] = nue_tgt_thr->GetBinContent(bx-2*nbins_CC+1);

  }
  TMatrixD diff = temp - target;
  TMatrixD diffT(TMatrixD::kTransposed, diff);

  TMatrixD chi2 = diffT*invmx*diff;

  //std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2[0][0] << "\n";
  return chi2[0][0];
}



double TemplateFitter::bfChi2( double Uee2, double Umm2, double dm2 )
{
  double par[3];
  par[0] = Uee2;
  par[1] = Umm2;
  par[2] = dm2;

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

    double Pmue, Pemu, Pmm, Pee;
/*
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
*/
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
    
    std::cout << i << "\t" << Pmue << "\t" << Pee << "\t" << Pmm << "\n";

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


  TH1D *htgt = new TH1D("htgt","",nbins,0,nbins);
  TH1D *htgt_thr = new TH1D("htgt_thr","",nbins,0,nbins);
  for(int i=0; i<nbins; i++){
    if(i<nbins_CC)                  htgt->SetBinContent(i+1, CCm_tgt->GetBinContent(i+1));
    if(i>=nbins_CC && i<2*nbins_CC) htgt->SetBinContent(i+1, CCe_tgt->GetBinContent(i-nbins_CC+1));
    if(i>=2*nbins_CC)               htgt->SetBinContent(i+1, nue_tgt->GetBinContent(i-2*nbins_CC+1));
  }
  for(int i=0; i<nbins; i++){
    if(i<nbins_CC)                  htgt_thr->SetBinContent(i+1, CCm_tgt_thr->GetBinContent(i+1));
    if(i>=nbins_CC && i<2*nbins_CC) htgt_thr->SetBinContent(i+1, CCe_tgt_thr->GetBinContent(i-nbins_CC+1));
    if(i>=2*nbins_CC)               htgt_thr->SetBinContent(i+1, nue_tgt_thr->GetBinContent(i-2*nbins_CC+1));
  }

  htgt->SetLineColor(kBlack);
  htgt->SetLineWidth(2.0);
  htgt_thr->SetLineColor(kGreen);
  htgt_thr->SetLineWidth(2.0);
  double x1 = htgt->GetXaxis()->GetBinUpEdge(nbins_CC);
  double x2 = htgt->GetXaxis()->GetBinUpEdge(2*nbins_CC);


  TH1D *hOsc = new TH1D("hOsc","",nbins,0,nbins);
  TH1D *hnOsc = new TH1D("hnOsc","",nbins,0,nbins);
  for(int i=0; i<nbins; i++){
    if(i<nbins_CC)                  hOsc->SetBinContent(i+1, CC_tp_em->GetBinContent(i+1));
    if(i>=nbins_CC && i<2*nbins_CC) hOsc->SetBinContent(i+1, CC_tp_me->GetBinContent(i-nbins_CC+1));
    if(i>=2*nbins_CC)               hOsc->SetBinContent(i+1, nue_tp_os->GetBinContent(i-2*nbins_CC+1));
  }
  for(int i=0; i<nbins; i++){
    if(i<nbins_CC)                  hnOsc->SetBinContent(i+1, CC_tp_mm->GetBinContent(i+1));
    if(i>=nbins_CC && i<2*nbins_CC) hnOsc->SetBinContent(i+1, CC_tp_ee->GetBinContent(i-nbins_CC+1));
    if(i>=2*nbins_CC)               hnOsc->SetBinContent(i+1, nue_tp_unos->GetBinContent(i-2*nbins_CC+1));
  }

  hOsc->SetFillColor(kRed);
  hnOsc->SetFillColor(kBlue);


  THStack *h = new THStack("h","");
  h->Add(hnOsc);
  h->Add(hOsc);

  TCanvas *c = new TCanvas("c","",900,600);
  c->SetLogy();
  h->Draw("hist");
  htgt->Draw("same");
  htgt_thr->Draw("same");
  c->Update();
  TLatex latex;
  latex.SetNDC(false);
  latex.SetTextSize(0.03);      // Set text size
  latex.SetTextAlign(22);       // Center alignment
  double yPosition = pow(10, (gPad->GetUymin() + gPad->GetUymax())/2.);         // Align to the middle of the pad
  latex.DrawLatex(nbins_CC/2, yPosition, "CCm");        // Position for CCm
  latex.DrawLatex(nbins_CC + nbins_CC/2, yPosition, "CCe");     // Position for CCe
  latex.DrawLatex(2*nbins_CC + nbins_nue/2, yPosition, "nue");  // Position for nue
  h->GetXaxis()->SetTitle("bin number");
  h->GetYaxis()->SetTitle("entries (/yr.POT)");
  TLine *line1 = new TLine(x1, 0, x1, pow(10, gPad->GetUymax()));
  TLine *line2 = new TLine(x2, 0, x2, pow(10, gPad->GetUymax()));
  line1->SetLineColor(kBlack);
  line2->SetLineColor(kBlack);
  line1->SetLineWidth(2.0);
  line2->SetLineWidth(2.0);
  line1->Draw("same");
  line2->Draw("same");
  TLegend *lg = new TLegend(0.70,0.75,0.9,0.9);
  lg->AddEntry(htgt,"target before throw");
  lg->AddEntry(htgt_thr,"target after throw");
  lg->AddEntry(hOsc,"oscillated");
  lg->AddEntry(hnOsc,"unoscillated");
  lg->Draw();
  c->SaveAs(Form("fit_flux_alg_%d.png",uni));


  // calculate the chi2 with the "data" target
  TMatrixD target(nbins, 1);
  TMatrixD temp(nbins, 1);
  TMatrixD unc(nbins, nbins);
  TMatrixD cov(nbins, nbins);

  for( int bx = 0; bx < nbins; bx++ ) {
    if(bx<nbins_CC)                   temp[bx][0] = CC_tp_m->GetBinContent(bx+1);
    if(bx>=nbins_CC && bx<2*nbins_CC) temp[bx][0] = CC_tp_e->GetBinContent(bx-nbins_CC+1);
    if(bx>=2*nbins_CC)                temp[bx][0] = nue_tp->GetBinContent(bx-2*nbins_CC+1);

    if(bx<nbins_CC)                   target[bx][0] = CCm_tgt_thr->GetBinContent(bx+1);
    if(bx>=nbins_CC && bx<2*nbins_CC) target[bx][0] = CCe_tgt_thr->GetBinContent(bx-nbins_CC+1);
    if(bx>=2*nbins_CC)                target[bx][0] = nue_tgt_thr->GetBinContent(bx-2*nbins_CC+1);

  }
  TMatrixD diff = temp - target;
  diff.Print();

  TMatrixD diffT(TMatrixD::kTransposed, diff);

  TMatrixD chi2 = diffT*invmx*diff;

  //std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";
  return chi2[0][0];

}

