#include "KlFitter.h"

KlFitter::KlFitter(const std::vector<reco::TransientTrack> &vecttrk) {
  KalmanVertexFitter vtxFitter(true);
  dimuvtx = vtxFitter.vertex(vecttrk);
  if (!dimuvtx.isValid())
    status_ = false;
  else {
    refited = dimuvtx.refittedTracks();
    prob_ = ChiSquaredProbability(dimuvtx.totalChiSquared(),
                                  dimuvtx.degreesOfFreedom());
  }
}

KlFitter::~KlFitter() {};

void KlFitter::fillNtuple(NtupleContent &nt) {
  if (status_) {
    nt.pair_svprob = prob_;
    nt.pair_fit_mass =
        DimuonMass(refited[0].track(), refited[1].track());
  } else {
    nt.pair_svprob = -1;
    nt.pair_fit_mass = -1;
  }
}
