//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// dimuon vertex fitter class

#ifndef KLFITTER_H
#define KLFITTER_H

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "NtupleContent.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "helper.h"

class KlFitter {
 public:
  KlFitter(const std::vector<reco::TransientTrack> &);
  ~KlFitter();

  void fillNtuple(NtupleContent &nt);
  bool status() { return status_; }
  float prob() { return prob_; }

 private:
  TransientVertex dimuvtx;
  std::vector<reco::TransientTrack> refited;
  bool status_ = true;
  float prob_ = 0;
};

#endif
