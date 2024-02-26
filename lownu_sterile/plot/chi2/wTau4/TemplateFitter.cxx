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

void TemplateFitter::setPara( char var[20], int par, double oscpar[4], int nuCut, double seed[4], double fitPara_m[29][7], double fitPara_e[29][7] )
{
  name  = var;
  para  = par;
  cutNu = nuCut;
  for(int i = 0; i < 4; i++){
    ospar[i] = oscpar[i];
  }
  b0 = oscpar[0];
  b1 = oscpar[1];
  b2 = oscpar[2];
  b3 = oscpar[3];

  s0 = seed[0];
  s1 = seed[1];
  s2 = seed[2];
  s3 = seed[3];

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

double TemplateFitter::getPme( double energy, double Uee2, double Umm2, double dm2, double L )
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

double TemplateFitter::getAvgPme( double energy, double Uee2, double Umm2, double dm2, double ft[7] )
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

void TemplateFitter::setCovmtr( double fl_bct[nbins+1][nbins+1], double sig_bct[nbins+1][nbins+1] )
{
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      flmx[i][j]  = fl_bct[i][j];
      sigmx[i][j] = sig_bct[i][j];
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
  TH1D * nue_tp_mt = (TH1D*) nue_tp->Clone();
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
    double me = 0., em = 0., mt = 0.;
    double ee = 0., mm = 0.;

    if(i<240) 		     iL = (int) i/10;
    else if(i>=240 && i<400) iL = (int) 24+(i-240)/40;
    else		     iL = 28;

    for(int k=0; k<7; k++) {
      ft_m[k] = fitP_m[iL][k];
      ft_e[k] = fitP_e[iL][k];
    }

    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      me = me + getAvgPme(e, par[0], par[1], par[2], ft_m); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2, par[3] = Utt2
      em = em + getAvgPme(e, par[0], par[1], par[2], ft_e); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2, par[3] = Utt2
      ee = ee + getAvgPee(e, par[0], par[1], par[2], ft_e);
      mm = mm + getAvgPmm(e, par[0], par[1], par[2], ft_m);
      mt = mt + getAvgPmm(e, par[0], par[1], par[3], ft_m);
    }
    double Pme = me/1001.0;
    double Pem = em/1001.0;
    double Pee = ee/1001.0;
    double Pmm = mm/1001.0;
    double Pmt = mt/1001.0;

    CC_tp_me->Add(CC_nc_m_templates[i], Pme);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pem);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pme);
    nue_tp_mm->Add(nue_m_templates[i], Pmm+Pmt);
    nue_tp_em->Add(nue_w_e_templates[i], Pem);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
    //nue_tp_mt->Add(nue_w_m_templates[i], Pmt);
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

  invmx = covmx;
  invmx.Invert();

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
  legend_nue->AddEntry(nue_tp_os,"oscillated #nu_{#mu}#rightarrow#nu_{e} & #nu_{e}#rightarrow#nu_{#mu} & #nu_{#mu}#rightarrow#nu_{#tau}");
  legend_nue->AddEntry(nue_tp_unos,"unoscillated #nu_{e}#rightarrow#nu_{e} & #nu_{#mu}->#nu_{#mu}");
  legend_nue->Draw();
  cnue->SaveAs(Form("nue_tgt_%d%d.png",para,cutNu));


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
  TH1D * nue_tp_mt = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_os = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_unos = (TH1D*) nue_tp->Clone(); 

  //TRandom *rand = new TRandom(4444);
  //double L = rand->Uniform(0.37,0.57);

  int iL;

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < nbins_Ev; ++i ) {
    double me = 0., em = 0., mt = 0.;
    double ee = 0., mm = 0.;

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
      me = me + getAvgPme(e, par[0], par[1], par[2], ft_m); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2, par[3] = Utt2
      em = em + getAvgPme(e, par[0], par[1], par[2], ft_e); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2, par[3] = Utt2
      ee = ee + getAvgPee(e, par[0], par[1], par[2], ft_e);
      mm = mm + getAvgPmm(e, par[0], par[1], par[2], ft_m);
      mt = mt + getAvgPmm(e, par[0], par[1], par[3], ft_m);
    }
    double Pme = me/1001.0;
    double Pem = em/1001.0;
    double Pee = ee/1001.0;
    double Pmm = mm/1001.0;
    double Pmt = mt/1001.0;

    CC_tp_me->Add(CC_nc_m_templates[i], Pme);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pem);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pme);
    nue_tp_mm->Add(nue_m_templates[i], Pmm+Pmt);
    nue_tp_em->Add(nue_w_e_templates[i], Pem);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
    //nue_tp_mt->Add(nue_w_m_templates[i], Pmt);
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
  }

  TMatrixD diff = temp - target;
  TMatrixD diffT(TMatrixD::kTransposed, diff);

  TMatrixD chi2 = diffT*invmx*diff; 

  std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2[0][0] << "\n";
  return chi2[0][0];

}


void TemplateFitter::Draw()
{
  int N = 20;

  double p0[1], p1[1], p2[1], p3[1];
  p0[0] = b0;
  p1[0] = b1;
  p2[0] = b2;
  p3[0] = b3;


  TH2D *h0 = new TH2D("h0","",N,0,0.1, N,0,12.0);
  TH2D *h1 = new TH2D("h1","",N,0,12.0,N,0,0.1);
  TH2D *h2 = new TH2D("h2","",N,0,0.1, N,0,0.1);
  
  TH2D *h0z = new TH2D("h0z","",N,0.0      ,2.0*p1[0],N,0.5*p2[0],1.5*p2[0]);
  TH2D *h1z = new TH2D("h1z","",N,0.5*p2[0],1.5*p2[0],N,0.0      ,2.0*p0[0]);
  TH2D *h2z = new TH2D("h2z","",N,0.0      ,2.0*p0[0],N,0.0      ,2.0*p1[0]);

  TH2D *h0d = new TH2D("h0d","",N,0,0.1, N,0,12.0);
  TH2D *h1d = new TH2D("h1d","",N,0,12.0,N,0,0.1);
  TH2D *h2d = new TH2D("h2d","",N,0,0.1, N,0,0.1);
  
  TH2D *h0dz = new TH2D("h0dz","",N,0.0      ,2.0*p1[0],N,0.5*p2[0],1.5*p2[0]);
  TH2D *h1dz = new TH2D("h1dz","",N,0.5*p2[0],1.5*p2[0],N,0.0      ,2.0*p0[0]);
  TH2D *h2dz = new TH2D("h2dz","",N,0.0      ,2.0*p0[0],N,0.0      ,2.0*p1[0]);

  double par[4], parz[4], bin[4];
  par[0] = b0;
  par[1] = b1;
  par[2] = b2;
  par[3] = parz[3] = bin[3] = 0.02;
  double chi2 = 0.0;
  double chi2z = 0.0;

  double nu=0.3, Ev=3.0;

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
  cchi2->SetLogz();
  h0->Draw("colz");
  g0->Draw("same C*");
  cchi2->cd(2);
  cchi2->SetLogz();
  h1->Draw("colz");
  g1->Draw("same C*");
  cchi2->cd(3);
  cchi2->SetLogz();
  h2->Draw("colz");
  g2->Draw("same C*");
  cchi2->SaveAs(Form("%s_chi2Surface_sys_%d%d.png",name,para,cutNu)); 
  
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
  cchi2_z->SaveAs(Form("%s_chi2Surface_sys_%d%d_zoom.png",name,para,cutNu)); 
 
  TCanvas *cdiff = new TCanvas("cdiff","",1200,300);
  cdiff->Divide(3,1);
  cdiff->cd(1);
  cdiff->SetLogz();
  h0d->Draw("colz");
  g0->Draw("same C*");
  cdiff->cd(2);
  cdiff->SetLogz();
  h1d->Draw("colz");
  g1->Draw("same C*");
  cdiff->cd(3);
  cdiff->SetLogz();
  h2d->Draw("colz");
  g2->Draw("same C*");
  cdiff->SaveAs(Form("%s_chi2Diff_sys_%d%d.png",name,para,cutNu)); 
  
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
  cdiff_z->SaveAs(Form("%s_chi2Diff_sys_%d%d_zoom.png",name,para,cutNu)); 

}

