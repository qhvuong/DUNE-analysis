#ifndef TEMPLATEFITTER_H 
#define TEMPLATEFITTER_H 

#include "Math/Minimizer.h"
#include "TH1.h"
#include "TMatrixD.h"

static const int nbins_Ev = 430;
static const int nbins_CC = 56;
static const int nbins_nue = 8;
static const int nbins = 2*nbins_CC + nbins_nue;

class TemplateFitter {

  public:

    TemplateFitter(TH1D * CC_templates_m[nbins_Ev], TH1D * CC_templates_m_nc[nbins_Ev], TH1D * CC_templates_e[nbins_Ev], TH1D * nue_templates_m[nbins_Ev], TH1D * nue_templates_m_w[nbins_Ev], TH1D * nue_templates_e[nbins_Ev], TH1D * nue_templates_e_w[nbins_Ev], TH1D * LEdep_m[29], TH1D * LEdep_e[29] );
    ~TemplateFitter(){};
    void setEnergyBins(double bins[nbins_Ev+1]);
    void setCovmtr(double flmx_bct[nbins+1][nbins+1], double sigmx_bct[nbins+1][nbins+1]);
    void setPara(char var[20], int nuCut, double fitPara_m[29][7], double fitPara_e[29][7]);
    bool doFitFine1(double seed[3], double &Uee2, double &Umm2, double &dm2);
    bool doFitFine2(double seed[3], double &Uee2, double &Umm2, double &dm2);
    bool doFitCoarse(double seed[3], double &Uee2, double &Umm2, double &dm2);
    void getTarget(double *par_tgt);
    double bfChi2(double Uee2, double Umm2, double dm2);
    double noChi2(double *par_no);


  private:

    double getPmue(double energy, double Uee2, double Umm2, double dm2, double L);
    double getPee(double energy, double Uee2, double Umm2, double dm2, double L);
    double getPmm(double energy, double Uee2, double Umm2, double dm2, double L);
    double getAvgPmue(double energy, double Uee2, double Umm2, double dm2, double ft[7]);
    double getAvgPee(double energy, double Uee2, double Umm2, double dm2, double ft[7]);
    double getAvgPmm(double energy, double Uee2, double Umm2, double dm2, double ft[7]);
    double getChi2(const double * par);

    // The templates are reconstructed lepton energy, in a slice of true neutrino energy
    TH1D * CC_m_templates[nbins_Ev];
    // The templates are reconstructed lepton energy, in a slice of true neutrino energy no cuts
    TH1D * CC_nc_m_templates[nbins_Ev];
    // The template for the intrinsic nue
    TH1D * CC_e_templates[nbins_Ev]; 
    TH1D * nue_m_templates[nbins_Ev];
    TH1D * nue_w_m_templates[nbins_Ev];
    TH1D * nue_e_templates[nbins_Ev];
    TH1D * nue_w_e_templates[nbins_Ev];
    TH1D * LE_m[29];
    TH1D * LE_e[29];
    // define the energy bins used in the template
    double m_energy_bins[nbins_Ev+1];
    double s[3], fitP_m[29][7], fitP_e[29][7], ft_m[7], ft_e[7];
    char *name;
    int cutNu;

/*
    double CCEdges[nbins_CC+1];
    void binning(){
      CCEdges[0]=0., CCEdges[1]=0.3;
      for(int i=1; i<nbins_CC+1; i++)
      {
      if(i<38)               CCEdges[i+1] = CCEdges[i] + 0.1;
      else if(i>=38 && i<43) CCEdges[i+1] = CCEdges[i] + 0.2;
      else if(i>=43 && i<48) CCEdges[i+1] = CCEdges[i] + 0.4;
      else if(i>=48 && i<53) CCEdges[i+1] = CCEdges[i] + 0.8;
      else if(i>=53 && i<55) CCEdges[i+1] = CCEdges[i] + 1.5;
      else                   CCEdges[i+1] = CCEdges[i] + 2.0;
      }
    }
    const double nueEdges[9] = {0., 0.3, 0.6, 0.92, 1.3, 1.75, 2.45, 3.9, 16.0};

    // this is the thing you are trying to fit to, i.e. the data distribution
    TH1D * CCe_tgt = new TH1D("CCe_tgt","",nbins_CC,CCEdges);
    TH1D * CCm_tgt = new TH1D("CCm_tgt","",nbins_CC,CCEdges);
    TH1D * nue_tgt = new TH1D("nue_tgt","",nbins_nue,nueEdges);
*/
    TH1D * CCe_tgt = new TH1D();
    TH1D * CCm_tgt = new TH1D();
    TH1D * nue_tgt = new TH1D();
};

#endif
