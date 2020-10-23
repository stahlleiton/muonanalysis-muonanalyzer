//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// flat tree branches/var declaration

#ifndef NTUPLECONTENT_H
#define NTUPLECONTENT_H
#include <string>
#include <algorithm>
#include <vector>
#include "TTree.h"

class NtupleContent {
 public:
  NtupleContent();
  virtual ~NtupleContent();
  void SetTree(TTree *t1);
  void CreateBranches(const std::vector<std::string> &);
  void CreateExtraTrgBranches(const std::vector<std::string> &);
  void CreateProbeFlagBranches(const std::vector<std::string> &);
  void ClearBranches();

  // Standard stuff
  int run;
  int event;
  int ls;
  bool fromFullAOD;

  // Beamspot and vertex
  float BSpot_x;
  float BSpot_y;
  float BSpot_z;
  float pv_x;
  float pv_y;
  float pv_z;
  int nvertices;

  // Pileup
  float trueNumInteractions;
  int puNumInteractions;
  double Rho;

  // Number of muons
  int nmuons;
  int ntag;
  int iprobe;

  // Triggers
  std::map<std::string, bool> trigger;
  std::vector<float> trg_pt;
  std::vector<float> trg_eta;
  std::vector<float> trg_phi;

  std::vector<float> prb_pt;
  std::vector<float> prb_eta;
  std::vector<float> prb_phi;
  // Triggers - probe
  std::map<std::string, bool> probe_trg;

  // Flags - probe
  std::map<std::string, bool> probe_flag;

  float genmu1_pt;
  float genmu1_eta;
  float genmu1_phi;
  float genmu2_pt;
  float genmu2_eta;
  float genmu2_phi;

  // Tag properties
  float tag_pt;
  float tag_eta;
  float tag_phi;
  bool tag_isLoose;
  bool tag_isMedium;
  bool tag_isTight;
  bool tag_isSoft;
  bool tag_isHighPt;
  float tag_relIso04;
  bool tag_isMatchedGen;

  // Probe properties
  float probe_pt;
  float probe_eta;
  float probe_phi;
  bool probe_isLoose;
  bool probe_isMedium;
  bool probe_isTight;
  bool probe_isSoft;
  bool probe_isHighPt;
  bool probe_isMuMatched;
  bool probe_isPF;
  bool probe_isSA;
  bool probe_isTracker;
  bool probe_isGlobal;
  bool probe_isdSA;
  bool probe_isdGlobal;
  bool probe_isCosmic;
  bool probe_isGood;
  bool probe_isHighPurity;
  float probe_relIso04;
  bool probe_isMatchedGen;

  float probe_validFraction;
  float probe_trkChi2;
  float probe_positionChi2;
  float probe_trkKink;
  float probe_segmentCompatibility;
  float probe_trackerLayers;
  float probe_pixelLayers;
  float probe_muonStations;
  float probe_muonHits;
  float probe_DTHits;
  float probe_CSCHits;
  float probe_pterr;
  float probe_dxy;
  float probe_dz;

  float probe_dsa_pt;
  float probe_dsa_eta;
  float probe_dsa_phi;
  float probe_dsa_trkChi2;
  float probe_dsa_muonStations;
  float probe_dsa_muonHits;
  float probe_dsa_DTHits;
  float probe_dsa_CSCHits;
  float probe_dsa_pterr;
  float probe_dsa_dxy;
  float probe_dsa_dz;

  // Pair properties
  float pair_pt;
  float pair_mass;
  float pair_eta;
  float pair_phi;
  float pair_fit_mass;
  float pair_svprob;
  float pair_dz;
  float pair_first_pair;

 private:
  TTree *t1;
};
#endif
