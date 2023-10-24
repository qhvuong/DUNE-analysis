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

double b0,b1,b2;
double s0,s1,s2;
double L_tg, L_tp;
void TemplateFitter::setPara( char var[20], double oscpar[3], int par, int nuCut, double seed[3], double fitPara_m[29][7], double fitPara_e[29][7] )
{
  name  = var;
  para = par;
  cutNu = nuCut;
  for(int i = 0; i < 3; i++){
    ospar[i] = oscpar[i];
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
  prob_u = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L2;
  prob_l = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);

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
  prob_u = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L2;
  prob_l = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);

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


const int N = 40;
TH2D *h0  = new TH2D("h0","",N,0,0.1, N,0,12.0);
TH2D *h1  = new TH2D("h1","",N,0,12.0,N,0,0.1);
TH2D *h2  = new TH2D("h2","",N,0,0.1, N,0,0.1);
TH2D *h0z = new TH2D("h0z","",N,0.0*b1,2.0*b1,N,0.5*b2,1.5*b2);
TH2D *h1z = new TH2D("h1z","",N,0.5*b2,1.5*b2,N,0.0*b0,2.0*b0);
TH2D *h2z = new TH2D("h2z","",N,0.0*b0,2.0*b0,N,0.0*b1,2.0*b1);


void TemplateFitter::getTarget( double *oscpar )
{
  CCe_tgt->Reset();
  CCm_tgt->Reset();
  nue_tgt->Reset();

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
    double mue = 0;
    double emu = 0;
    double ee = 0;
    double mm = 0;
    double L = 0.5;

    if(i<240)                iL = (int)i/10;
    else if(i>=240 && i<400) iL = (int) 24+(i-240)/40;
    else		     iL = 28;

    for(int k=0; k<7; k++) {
      ft_m[k] = fitP_m[iL][k];
      ft_e[k] = fitP_e[iL][k];
    }

    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getAvgPmue(e, ospar[0], ospar[1], ospar[2], ft_m); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      emu = emu + getAvgPmue(e, ospar[0], ospar[1], ospar[2], ft_e); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getAvgPee(e, ospar[0], ospar[1], ospar[2], ft_e);
      mm  = mm  + getAvgPmm(e, ospar[0], ospar[1], ospar[2], ft_m);
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


  THStack *he = new THStack("he","");
  CCe_tgt->SetMarkerStyle(kStar);
  CCe_tgt->SetMarkerSize(1);
  CCe_tgt->SetMarkerColor(8);
  CC_tp_me->SetFillColor(kRed);
  CC_tp_ee->SetFillColor(kBlue);
  he->Add(CC_tp_me);
  he->Add(CC_tp_ee);
  TCanvas *ce = new TCanvas("ce","",800,600);
  //he->SetMaximum(60000);
  he->Draw("hist");
  CCe_tgt->Draw("same");
  TLegend *legend_e = new TLegend(0.55,0.70,0.9,0.9);
  legend_e->AddEntry(CCe_tgt,"CCe data");
  legend_e->AddEntry(CC_tp_me,"oscillated #nu_{#mu}#rightarrow#nu_{e}");
  legend_e->AddEntry(CC_tp_ee,"unoscillated #nu_{e}#rightarrow#nu_{e}");
  legend_e->Draw();
  ce->SaveAs(Form("CCe_tgt_%d%d.png",para,cutNu));

  THStack *hm = new THStack("hm","");
  CCm_tgt->SetMarkerStyle(kStar);
  CCm_tgt->SetMarkerSize(1);
  CCm_tgt->SetMarkerColor(8);
  CC_tp_em->SetFillColor(kRed);
  CC_tp_mm->SetFillColor(kBlue);
  hm->Add(CC_tp_em);
  hm->Add(CC_tp_mm);
  TCanvas *cm = new TCanvas("cm","",800,600);
  hm->Draw("hist");
  CCm_tgt->Draw("same");
  TLegend *legend_m = new TLegend(0.55,0.70,0.9,0.9);
  legend_m->AddEntry(CCm_tgt,"CCm data");
  legend_m->AddEntry(CC_tp_me,"oscillated #nu_{e}#rightarrow#nu_{#mu}");
  legend_m->AddEntry(CC_tp_mm,"unoscillated #nu_{#mu}->#nu_{#mu}");
  legend_m->Draw();
  cm->SaveAs(Form("CCm_tgt_%d%d.png",para,cutNu));

  THStack *hnue = new THStack("hnue","");
  nue_tgt->SetMarkerStyle(kStar);
  nue_tgt->SetMarkerSize(1);
  nue_tgt->SetMarkerColor(8);
  nue_tp_os->SetFillColor(kRed);
  nue_tp_unos->SetFillColor(kBlue);
  hnue->Add(nue_tp_os);
  hnue->Add(nue_tp_unos);
  TCanvas *cnue = new TCanvas("cnue","",800,600);
  gPad->SetLogy();
  hnue->Draw("hist");
  nue_tgt->Draw("same");
  TLegend *legend_nue = new TLegend(0.55,0.70,0.9,0.9);
  legend_nue->AddEntry(nue_tgt,"#nu+e data");
  legend_nue->AddEntry(nue_tp_os,"oscillated #nu_{#mu}#rightarrow#nu_{e} & #nu_{e}#rightarrow#nu_{#mu}");
  legend_nue->AddEntry(nue_tp_unos,"unoscillated #nu_{e}#rightarrow#nu_{e} & #nu_{#mu}->#nu_{#mu}");
  legend_nue->Draw();
  cnue->SaveAs(Form("nue_tgt_%d%d.png",para,cutNu));


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

  //TRandom *rand = new TRandom(4444);
  //double L = rand->Uniform(0.37,0.57);

  int iL;

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < nbins_Ev; ++i ) {
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
  TMatrixD cov(nbins, nbins);

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

  cov = covmtr + unc;

  TMatrixD diff = temp - target;
  TMatrixD diff_T(TMatrixD::kTransposed, diff);
  TDecompSVD svd(cov);
  TMatrixD inv = svd.Invert();

  TMatrixD diff_cov(1, nbins);
  for( int bx = 0; bx < nbins; bx++ ) {
    for( int by = 0; by < nbins; by++ ) {
      diff_cov[0][bx] += diff_T[0][by] * inv[by][bx];
    }
    chi2 += diff_T[0][bx] * diff[bx][0];
  }

  std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";
  return chi2;

}

double bf0[1], bf1[1], bf2[1];

bool TemplateFitter::doFit( double &Uee2, double &Umm2, double &dm2 )
{
  // Make a Minuit fitter object
  ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2"); 
  fitter->SetMaxFunctionCalls(1000000); // maximum number of times to try to find the minimum before failing
  fitter->SetMaxIterations(1000000);
  fitter->SetTolerance(0.1); // You might have to play with this -- how close to the correct value do you need to be?

  // The variables will be normalizations of the templates, we will start the with seed values of 1.0
  // fourth argument is step size, i.e. how much to change the normalization by at each step
  fitter->SetVariable( 0, "Uee2", s0, 0.00001 );
  fitter->SetVariable( 1, "Umm2", s1, 0.00001 );
  fitter->SetVariable( 2, "dm2",  s2, 0.01 );
  fitter->SetVariableLowerLimit(0, 0.0);
  fitter->SetVariableLowerLimit(1, 0.0);
  fitter->SetVariableLowerLimit(2, 0.0);

  // 3 free parameters = theta, dm2
  ROOT::Math::Functor lf( this, &TemplateFitter::getChi2, 3 );
  ROOT::Math::Functor functor( lf, 3 );
  fitter->SetFunction( functor );

  // Go!
  fitter->Minimize();

/*
  if( fitter->Status() != 0 ) {
    std::cout << "Something bad happened" << std::endl;
    return false;
  }
*/

  const double *bestfit = fitter->X();
  Uee2 = bestfit[0];
  Umm2 = bestfit[1];
  dm2 = bestfit[2];
  bf0[0] = bestfit[0];
  bf1[0] = bestfit[1];
  bf2[0] = bestfit[2];
  
  double chi2 = fitter->MinValue();
  //h->Fill(Uee2, Umm2, dm2);
  return true;

}

double TemplateFitter::bfChi2( double *par )
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
/*
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getPmue(e, par[0], par[1], par[2], L);
      emu = emu + getPmue(e, par[0], par[1], par[2], L);
      ee  = ee  + getPee(e, par[0], par[1], par[2], L);
      mm  = mm  + getPmm(e, par[0], par[1], par[2], L);
    }
*/

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

  THStack *he = new THStack("he","");
  CCe_tgt->SetMarkerStyle(kStar);
  CCe_tgt->SetMarkerSize(1);
  CCe_tgt->SetMarkerColor(8);
  CC_tp_me->SetFillColor(kRed);
  CC_tp_ee->SetFillColor(kBlue);
  he->Add(CC_tp_me);
  he->Add(CC_tp_ee);
  TCanvas *ce = new TCanvas("ce","",800,600);
  he->Draw("hist");
  CCe_tgt->Draw("same");
  TLegend *legend_e = new TLegend(0.55,0.70,0.9,0.9);
  legend_e->AddEntry(CCe_tgt,"CCe data");
  legend_e->AddEntry(CC_tp_me,"oscillated #nu_{#mu}#rightarrow#nu_{e}");
  legend_e->AddEntry(CC_tp_ee,"unoscillated #nu_{e}#rightarrow#nu_{e}");
  legend_e->Draw();
  ce->SaveAs(Form("fit_stat_CCe_%d%d.png",para,cutNu));

  THStack *hm = new THStack("hm","");
  CCm_tgt->SetMarkerStyle(kStar);
  CCm_tgt->SetMarkerSize(1);
  CCm_tgt->SetMarkerColor(8);
  CC_tp_em->SetFillColor(kRed);
  CC_tp_mm->SetFillColor(kBlue);
  hm->Add(CC_tp_em);
  hm->Add(CC_tp_mm);
  TCanvas *cm = new TCanvas("cm","",800,600);
  hm->Draw("hist");
  CCm_tgt->Draw("same");
  TLegend *legend_m = new TLegend(0.55,0.70,0.9,0.9);
  legend_m->AddEntry(CCm_tgt,"CCm data");
  legend_m->AddEntry(CC_tp_me,"oscillated #nu_{e}#rightarrow#nu_{#mu}");
  legend_m->AddEntry(CC_tp_mm,"unoscillated #nu_{#mu}->#nu_{#mu}");
  legend_m->Draw();
  cm->SaveAs(Form("fit_stat_CCm_%d%d.png",para,cutNu));

  THStack *hnue = new THStack("hnue","");
  nue_tgt->SetMarkerStyle(kStar);
  nue_tgt->SetMarkerSize(1);
  nue_tgt->SetMarkerColor(8);
  nue_tp_os->SetFillColor(kRed);
  nue_tp_unos->SetFillColor(kBlue);
  hnue->Add(nue_tp_os);
  hnue->Add(nue_tp_unos);
  TCanvas *cnue = new TCanvas("cnue","",800,600);
  gPad->SetLogy();
  hnue->Draw("hist");
  nue_tgt->Draw("same");
  TLegend *legend_nue = new TLegend(0.55,0.70,0.9,0.9);
  legend_nue->AddEntry(nue_tgt,"#nu+e data");
  legend_nue->AddEntry(nue_tp_os,"oscillated #nu_{#mu}#rightarrow#nu_{e} & #nu_{e}#rightarrow#nu_{#mu}");
  legend_nue->AddEntry(nue_tp_unos,"unoscillated #nu_{e}#rightarrow#nu_{e} & #nu_{#mu}->#nu_{#mu}");
  legend_nue->Draw();
  cnue->SaveAs(Form("fit_stat_nue_%d%d.png",para,cutNu));

  // calculate the chi2 with the "data" target
  double chi2 = 0.0;

  TMatrixD target(nbins, 1);
  TMatrixD temp(nbins, 1);
  TMatrixD unc(nbins, nbins);
  TMatrixD cov(nbins, nbins);

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

  cov = covmtr + unc;

  TMatrixD diff = temp - target;
  TMatrixD diff_T(TMatrixD::kTransposed, diff);
  TDecompSVD svd(cov);
  TMatrixD inv = svd.Invert();

  TMatrixD diff_cov(1, nbins);
  for( int bx = 0; bx < nbins; bx++ ) {
    for( int by = 0; by < nbins; by++ ) {
      diff_cov[0][bx] += diff_T[0][by] * inv[by][bx];
    }
    chi2 += diff_cov[0][bx] * diff[bx][0];
  }
  return chi2;
}

/*
double TemplateFitter::bfChi2 ( double *par )
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

  covmtr = covmtr + unc;

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

  return chi2;
}

*/
