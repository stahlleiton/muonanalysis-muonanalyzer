//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// flat tree branches/var declaration

#ifndef NTUPLECONTENT_H
#define NTUPLECONTENT_H
#include "TTree.h"
#include <vector>
#include <string>
#include "TString.h"

class NtupleContent {
public:
  NtupleContent();
  virtual ~NtupleContent();
  void SetTree(TTree *t1);
  void CreateBranches(const std::vector<std::string> &);
  void CreateExtraTrgBranches(const std::vector<std::string> &);
  void ClearBranches();

  // standard stuff
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

  //number of muons
  int nmuons;
  int ntag;
  int iprobe;

  //triggers
  bool trigger[10];
  std::vector<float> trg_pt;
  std::vector<float> trg_eta;
  std::vector<float> trg_phi;

  std::vector<float> prb_pt;
  std::vector<float> prb_eta;
  std::vector<float> prb_phi;
  //triggers - probe
  bool probe_trg[10] = {false};

  float genmu1_pt;
  float genmu1_eta;
  float genmu1_phi;
  float genmu2_pt;
  float genmu2_eta;
  float genmu2_phi;

  // tag properties
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

  // probe properties
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
  bool probe_isGlobal;
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
  float probe_dxy;
  float probe_dz;

  // pair properties
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
