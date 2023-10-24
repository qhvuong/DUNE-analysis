#ifndef TEMPLATEFITTER_H 
#define TEMPLATEFITTER_H 

#include "Math/Minimizer.h"
#include "TH1.h"
#include "TMatrixD.h"

static const int nbins_Ev = 430;
static const int nbins = 250;
static const int nbins_CC = 100;
static const int nbins_nue = 50;


class TemplateFitter {

  public:

    TemplateFitter(TH1D * CC_templates_m[nbins_Ev], TH1D * CC_templates_m_nc[nbins_Ev], TH1D * CC_templates_e[nbins_Ev], TH1D * nue_templates_m[nbins_Ev], TH1D * nue_templates_m_w[nbins_Ev], TH1D * nue_templates_e[nbins_Ev], TH1D * nue_templates_e_w[nbins_Ev], TH1D * LEdep_m[29], TH1D * LEdep_e[29] );
    ~TemplateFitter(){};
    void setEnergyBins(double bins[nbins_Ev+1]);
    void setCovmtr(double bincontent[nbins+1][nbins+1]);
    void setPara(char var[20], int par, double oscpar[3], int nuCut, double seed[3], double fitPara_m[29][7], double fitPara_e[29][7]);
    void Draw();
    void getTarget(double oscpar[3]);
    double getChi2(double * par);


  private:

    double getPmue(double energy, double Uee2, double Umm2, double dm2, double L);
    double getPee(double energy, double Uee2, double Umm2, double dm2, double L);
    double getPmm(double energy, double Uee2, double Umm2, double dm2, double L);
    double getAvgPmue(double energy, double Uee2, double Umm2, double dm2, double ft[7]);
    double getAvgPee(double energy, double Uee2, double Umm2, double dm2, double ft[7]);
    double getAvgPmm(double energy, double Uee2, double Umm2, double dm2, double ft[7]);

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
    double b0,b1,b2;
    double s0,s1,s2;
    double m_energy_bins[nbins_Ev+1];
    double ospar[3], s[3], fitP_m[29][7], fitP_e[29][7], ft_m[7], ft_e[7];
    char *name;;
    int para, cutNu, cutEv;
    // this is the thing you are trying to fit to, i.e. the data distribution

    TH1D * CCe_tgt = new TH1D("CCe_tgt","",nbins_CC,0,16);
    TH1D * CCm_tgt = new TH1D("CCm_tgt","",nbins_CC,0,16);
    TH1D * nue_tgt = new TH1D("nue_tgt","",nbins_nue,0,16);

/*
    TH1D * CCm_tgt;
    TH1D * CCe_tgt;
    TH1D * nue_tgt;
*/   
    //TMatrixD covmtr(300, 300); 
};
/*
class TemplateFitter_nue {

  public:

    TemplateFitter_nue(TH1D * templates_nue[100], TH1D * intrinsic_nue, TH1D * nue_tgt);
    ~TemplateFitter_nue(){};
    void setEnergyBins(double bins[101]);
    bool doFit(double &theta, double &dm2);

  private:

    double getOscProb(double energy, double theta, double dm2);
    double getChi2(const double * par);

    // The templates are reconstructed lepton energy, in a slice of true neutrino energy
    TH1D * m_templates_nue[100];
    // The template for the intrinsic nue
    TH1D * m_intrinsic_nue; 
    // define the energy bins used in the template
    double m_energy_bins_nue[101];
    // this is the thing you are trying to fit to, i.e. the data distribution
    TH1D * m_nue_tgt;
};*/

#endif
