//
//  CMS_TEST_ANALYSIS.cc 
//  
//
//  Created by Yaxian Mao on 3/9/15.
//
//

#include <stdio.h>

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
namespace Rivet {
    class CMS_TEST_ANALYSIS : public Analysis {
    public:
        /// @name Constructors etc.
        //@{
        /// Constructor
        CMS_TEST_ANALYSIS()
        : Analysis("CMS_TEST_ANALYSIS")
        //,
        //_filter(fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3))),
        //_trimmer(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03))),
        //_pruner(fastjet::Pruner(fastjet::cambridge_algorithm, 0.1, 0.5))
        {
            setBeams(PROTON, PROTON);
            setNeedsCrossSection(true);
        }
        //@}
    public:
        /// @name Analysis methods
        //@{
        /// Book histograms and initialise projections before the run
           const double etaBins[5] = {0., 0.5, 1.0, 1.5, 2.0};
        //   const int N_ETA_BINS = 4 ; //= sizeof(etaBins)/sizeof(double)-1 ;
        void init() {
           FinalState fs(-2.4, 2.4, 0*GeV);
         //  addProjection(fs, "FS");
            // addProjection(fs, "FS");
            // Jet collections
            addProjection(FastJets(fs, FastJets::ANTIKT, 0.3), "JetsAK3");
            addProjection(FastJets(fs, FastJets::ANTIKT, 0.5), "JetsAK5");
            addProjection(FastJets(fs, FastJets::ANTIKT, 0.7), "JetsAK7");

           //Histograms
        for (size_t i = 0; i < 4; ++i ) {
          _hist_sigma_AK3[i] = bookHistogram1D(i+1+0*4, 1, 1);
          _hist_sigma_AK5[i] = bookHistogram1D(i+1+1*4, 1, 1);
          _hist_sigma_AK7[i] = bookHistogram1D(i+1+2*4, 1, 1);
       }
    }
        // Find the pT histogram bin index for value pt (in GeV), to hack a 2D histogram equivalent
        /// @todo Use a YODA axis/finder alg when available
        size_t findBin(double etaJ) {
         for (size_t ibin = 0; ibin < 4; ++ibin) {
         if (inRange(etaJ, etaBins[ibin], etaBins[ibin+1])) return ibin;
         }
        return 4;
        }
        /// Perform the per-event analysis
        void analyze(const Event& event) {
            const double weight = event.weight();

             const PseudoJets& fjAK3 = applyProjection<FastJets>(event, "JetsAK3").pseudoJetsByPt( 30.0*GeV );
             const PseudoJets& fjAK5 = applyProjection<FastJets>(event, "JetsAK5").pseudoJetsByPt( 30.0*GeV );
             const PseudoJets& fjAK7 = applyProjection<FastJets>(event, "JetsAK7").pseudoJetsByPt( 30.0*GeV );


            // ... and fill the histograms
             for(size_t a = 0;a<fjAK3.size();++a){
               const fastjet::PseudoJet& jetsAK3 = fjAK3[a];
              double ieta = fabs(jetsAK3.eta());
            const size_t njetBin = findBin(ieta);
           if (njetBin >=4) vetoEvent;
               
               _hist_sigma_AK3[njetBin]->fill((jetsAK3.pt())/GeV, weight);
             }
             for(size_t a = 0;a<fjAK5.size();++a){
               const fastjet::PseudoJet& jetsAK5 = fjAK5[a];
              double ieta = fabs(jetsAK5.eta());
            const size_t njetBin = findBin(ieta);
           if (njetBin >= 4) vetoEvent;

               _hist_sigma_AK5[njetBin]->fill((jetsAK5.pt())/GeV, weight);
             }
              for(size_t a = 0;a<fjAK7.size();++a){
               const fastjet::PseudoJet& jetsAK7 = fjAK7[a];
              double ieta = fabs(jetsAK7.eta());
            const size_t njetBin = findBin(ieta);
           if (njetBin >= 4) vetoEvent;

               _hist_sigma_AK7[njetBin]->fill((jetsAK7.pt())/GeV, weight);
             }

        }
        /// Normalise histograms etc., after the run
        void finalize() {
          //  const double normalizationVal = crossSection()/sumOfWeights(); // the 2 is for absolute eta from -2 to +2
          //  const double normalizationVal = 1000; // the 2 is for absolute eta from -2 to +2
            const double normalizationVal = crossSection()/(5000*201); // total cross section/Nevents
            for (size_t i = 0; i < 4; ++i ) {
             normalize(_hist_sigma_AK3[i], normalizationVal);
            normalize(_hist_sigma_AK5[i], normalizationVal);
            normalize(_hist_sigma_AK7[i], normalizationVal);
        //   _hist_sigma_AK3[i].scale(crossSection()/sumOfWeights()/2, this);
       //    _hist_sigma_AK5[i].scale(crossSection()/sumOfWeights()/2, this);
       //    _hist_sigma_AK7[i].scale(crossSection()/sumOfWeights()/2, this);
          } 
        }

  private:
   //  const double etaBins[] = {0., 0.5, 1.0, 1.5, 2.0};
   //  const int N_ETA_BINS = sizeof(etaBins)/sizeof(double)-1 ;
  //  BinnedHistogram<double> _hist_sigma_AK3[4];
 //   BinnedHistogram<double> _hist_sigma_AK5[4];
 //   BinnedHistogram<double> _hist_sigma_AK7[4];
     AIDA::IHistogram1D*  _hist_sigma_AK3[4];
    AIDA::IHistogram1D* _hist_sigma_AK5[4];
    AIDA::IHistogram1D* _hist_sigma_AK7[4];

    };
    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(CMS_TEST_ANALYSIS);
}
