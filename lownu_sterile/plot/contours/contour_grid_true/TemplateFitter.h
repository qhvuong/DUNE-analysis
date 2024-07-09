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
    bool doFitFine1(double seed[4], double &Ue42, double &Um42, double &Ut42, double &dm2);
    bool doFitFine2(double seed[4], double &Ue42, double &Um42, double &Ut42, double &dm2);
    bool doFitCoarse(double seed[4], double &Ue42, double &Um42, double &Ut42, double &dm2);
    void getTarget(double *par_tgt);
    double bfChi2(double Ue42, double Um42, double Ut42, double dm2);


  private:

    double getAvgPme(double energy, double Ue42, double Um42, double Ut42, double dm2, double ft[7]);
    double getAvgPmt(double energy, double Ue42, double Um42, double Ut42, double dm2, double ft[7]);
    double getAvgPee(double energy, double Ue42, double Um42, double Ut42, double dm2, double ft[7]);
    double getAvgPmm(double energy, double Ue42, double Um42, double Ut42, double dm2, double ft[7]);
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


    // this is the thing you are trying to fit to, i.e. the data distribution
    TH1D * CCe_tgt = new TH1D();
    TH1D * CCm_tgt = new TH1D();
    TH1D * nue_tgt = new TH1D();
};

#endif
