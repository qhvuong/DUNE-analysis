#ifndef TEMPLATEFITTER_H
#define TEMPLATEFITTER_H

#include "TH1.h"
#include "TMatrixD.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <vector>
#include <utility>
#include "TCanvas.h"


//Oscillation variable names
constexpr const char* OscParNames[4] = {"ue42", "um42", "ut42", "dm2"};

//Oscillation parameter ranges
static const double par_min[4] = {0., 0., 0., 0.};
// static const double par_min[4] = {1E-5, 1E-5, 1E-5, 1E-3};
// static const double par_max[4] = {.394, .489, .718, 1000.};
static const double par_max[4] = {1., 1., 1., 1000.};

//Neutrino bins
static const int nbins_Ev = 430;

//Lepton bins
static const int nbinsCCm = 41;
static const int nbinsCCe = 54;
static const int nbinsCC = 54;
static const int nbins_nue = 8;
static const int nbins = 2*nbinsCC + nbins_nue;


struct CovarianceMatrix {
  double FluxMx[nbins + 1][nbins + 1];
  double CrossSectionMx[nbins + 1][nbins + 1];
  double DetectorMx[nbins + 1][nbins + 1];
  double FCWeights[nbins + 1];
};

struct LE_FitParameters {
  TH1D* LEdep_m[29];
  TH1D* LEdep_e[29];
  double LEfit_m[29][7];
  double LEfit_e[29][7];
};

// struct TemplateGroup_CC {
//     TH1D* mCC;
//     TH1D* mCC_NoCut;
//     TH1D* eCC;
//     TH1D* eCC_inv;
//     TH1D* eCC_inv_w;
// };

// struct TemplateGroup_nue {
//     TH1D* m_nue;
//     TH1D* m_nue_w;
//     TH1D* e_nue;
//     TH1D* e_nue_w;
//     TH1D* m_nue_inv;
//     TH1D* m_nue_inv_w;
//     TH1D* e_nue_inv;
//     TH1D* e_nue_inv_w;
// };

struct TemplateGroup_CC {
    TH1D* mCC_cc;
    TH1D* mCC_nc;
    TH1D* mCC_w;
    TH1D* eCC_cc;
    TH1D* eCC_nc;
    TH1D* eCC_BKGm;
    TH1D* eCC_BKGm_w;
    TH1D* eCC_BKGe;
    TH1D* eCC_BKGe_w;
};

struct TemplateGroup_nue {
    TH1D* m_nue;
    TH1D* m_nue_w;
    TH1D* e_nue;
    TH1D* e_nue_w;
    TH1D* nue_BKGe_cc;
    TH1D* nue_BKGe_nc;
};

class TemplateFitter {
  public:
    TemplateFitter(TemplateGroup_CC (&ccTemplates)[nbins_Ev], TemplateGroup_nue (&nueTemplates)[nbins_Ev], const std::vector<int>& dropBins = {});
    ~TemplateFitter() {}

    void setEnergyBins(double bins[nbins_Ev + 1]);
    void setCovmtr(const CovarianceMatrix& cov);
    void setLEParameters(LE_FitParameters LEPars);
    bool doFit(int nFunctionCalls, int nIterations, double Tolerance, const char* varName[4], double seed[4], double stepSize[4], std::vector<double>& parBf, std::vector<double>& parError, double& chi2, bool fixedUt42 = false);
    double getTarget(double par[4], bool originalTgt = false);
    double CalculateChi2(double par[4], bool plotBestfit = false, const char* outName = "");

    // === Public getters for safe access ===
    const TemplateGroup_CC* getCCTemplates() const { return CC; }
    const TemplateGroup_nue* getNueTemplates() const { return nue; }
    const LE_FitParameters& getLEFitParameters() const { return LEfit; }
    const CovarianceMatrix& getCovMatrix() const { return cov; }
    const double* getBinEdges() const { return binEdges; }
    // void plotOscillationComponents(TCanvas* c, const double par[4],
    //                            TH1D* CC_tp_mm, TH1D* CC_tp_me, TH1D* CC_tp_ee, TH1D* CC_tp_em,
    //                            TH1D* nue_tp_mm, TH1D* nue_tp_me, TH1D* nue_tp_ee, TH1D* nue_tp_em,
    //                            bool plotBestfit = false) const;
    void plotOscillationComponents(
        TCanvas* c, const double par[4],
        // CCνμ (μ-like sample)
        TH1D* CCm_sig_em, TH1D* CCm_sig_mm,
        // CCνe background components (from NC)
        TH1D* CCm_bkg_nc,
        // CCνe signal (e-like sample)
        TH1D* CCe_sig_me, TH1D* CCe_sig_ee,
        // CCνe background components (from NC and from ν–e leakage)
        TH1D* CCe_bkg_nc, TH1D* CCe_bkg_me, TH1D* CCe_bkg_ee, TH1D* CCe_bkg_em, TH1D* CCe_bkg_mm,
        // ν–e signal (inverse Eθ² selection)
        TH1D* nue_sig_mm, TH1D* nue_sig_me, TH1D* nue_sig_ee, TH1D* nue_sig_em,
        // ν–e background (from CCνe leakage into ν–e)
        TH1D* nue_bkg_ee, TH1D* nue_bkg_em,
        bool plotBestfit) const;


  private:
    double getAvgPme(double energy, double ue42, double um42, double ut42, double dm2, double ft[7]) const;
    double getAvgPmt(double energy, double ue42, double um42, double ut42, double dm2, double ft[7]) const;
    double getAvgPee(double energy, double ue42, double um42, double ut42, double dm2, double ft[7]) const;
    double getAvgPet(double energy, double ue42, double um42, double ut42, double dm2, double ft[7]) const;
    double getAvgPmm(double energy, double ue42, double um42, double ut42, double dm2, double ft[7]) const;
    TMatrixD BuildPredictionVector(const double par[4], bool plotBestfit = false, const char* plotName = "") const;
    double CalcChi2Core(const TMatrixD& prediction) const;
    double getChi2(const double* par);

    TemplateGroup_CC CC[nbins_Ev];
    TemplateGroup_nue nue[nbins_Ev];
    LE_FitParameters LEfit;
    CovarianceMatrix cov;
    double binEdges[nbins_Ev + 1];

    std::vector<int> dropBins_;

    TMatrixD invmx;
    TMatrixD myTarget;

    // Instance-local matrices for thread safety
    TMatrixD covmx = TMatrixD(nbins, nbins);
    TMatrixD statmx = TMatrixD(nbins, nbins);
    TMatrixD flmx = TMatrixD(nbins, nbins);
    TMatrixD sigmx = TMatrixD(nbins, nbins);
    TMatrixD detmx = TMatrixD(nbins, nbins);
    TMatrixD sysmx = TMatrixD(nbins, nbins);
    TMatrixD target = TMatrixD(nbins, 1);
    TMatrixD FCtarget = TMatrixD(nbins, 1);
};


namespace FitUtils {

inline std::vector<std::pair<double, double>> GetDM2SeedsAndSteps(int nPoints_dm2) {
    if (nPoints_dm2 == 4) {
        return {
            {0.1, 1.0}, {1.0, 1.0}, {10.0, 2.0}, {100.0, 10.0}
        };
    } else if (nPoints_dm2 == 8) {
        return {
            {0.1, 1.0}, {0.4, 1.0}, {1.0, 1.0}, {3.0, 1.0},
            {8.0, 1.0}, {20.0, 2.0}, {50.0, 5.0}, {100.0, 10.0}
        };
    } else if (nPoints_dm2 == 16) {
        return {
            {0.1, 1.0}, {0.4, 1.0}, {0.7, 1.0}, {1.0, 1.0}, 
            {2.0, 1.0}, {3.0, 1.0}, {4.0, 1.0}, {5.0, 1.0}, 
            {6.0, 2.0}, {7.0, 2.0}, {8.0, 2.0}, {9.0, 2.0},
            {10.0, 5.0}, {50.0, 5.0}, {100.0, 10.0}, {200.0, 20.0}
        };
    } else {
        throw std::invalid_argument("nPoints_dm2 must be one of {4, 8, 16}");
    }
}



inline std::vector<double> GetStage1StepSizes() {
    return {0.03, 0.03, 0.03, 10.0};
}

inline std::vector<double> GetStage3StepSizes() {
    return {0.001, 0.001, 0.001, 0.1};
}

inline std::vector<std::vector<double>> GetStage1SeedSets() {
    return {
        {par_min[0], par_min[1], par_min[2], par_min[3]},
        {par_max[0], par_max[1], par_max[2], 100.0},
        {par_min[0], par_min[1], par_min[2], 100.0},
        {par_max[0], par_max[1], par_min[2], par_min[3]}
    };
}

constexpr double Stage1Tolerance = 1.0;
constexpr double Stage2Tolerance = 0.1;
constexpr double Stage3Tolerance = 0.01;

constexpr int Stage1MaxCalls = 5000;
constexpr int Stage1MaxIters = 500;
constexpr int Stage2MaxCalls = 2000;
constexpr int Stage2MaxIters = 200;
constexpr int Stage3MaxCalls = 1000;
constexpr int Stage3MaxIters = 200;

} // namespace FitUtils

#endif // TEMPLATEFITTER_H
