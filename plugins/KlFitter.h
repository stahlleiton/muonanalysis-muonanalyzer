//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// dimuon vertex fitter class

#ifndef KLFITTER_H
#define KLFITTER_H

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "NtupleContent.h"
#include "helper.h"

class KlFitter{
  public:
    KlFitter(std::vector<reco::TransientTrack> &);
    ~KlFitter();

    void fillNtuple( NtupleContent &nt);
    bool status(){ return  status_;}
    float prob(){ return prob_; }

  private:
    TransientVertex dimuvtx;
    std::vector<reco::TransientTrack> refited;
    bool status_=true;
    float prob_=0;    
};


#endif
