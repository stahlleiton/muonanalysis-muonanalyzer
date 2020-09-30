#include "NtupleContent.h"

NtupleContent::NtupleContent() {}

NtupleContent::~NtupleContent() {}

void NtupleContent::SetTree(TTree *mytree) { t1 = mytree; }

void NtupleContent::CreateBranches(const std::vector<std::string> &HLTs) {
  // General
  t1->Branch("run", &run);
  t1->Branch("event", &event);
  t1->Branch("ls", &ls);
  t1->Branch("fromFullAOD", &fromFullAOD);
  t1->Branch("BSpot_x", &BSpot_x);
  t1->Branch("BSpot_y", &BSpot_y);
  t1->Branch("BSpot_z", &BSpot_z);
  t1->Branch("pv_x", &pv_x);
  t1->Branch("pv_y", &pv_y);
  t1->Branch("pv_z", &pv_z);
  t1->Branch("nVertices", &nvertices);
  t1->Branch("nTrueInteractions", &trueNumInteractions);
  t1->Branch("nPUInteractions", &puNumInteractions);
  t1->Branch("rho", &Rho);
  t1->Branch("nmuons", &nmuons);
  t1->Branch("ntag", &ntag);
  t1->Branch("genmu1_pt", &genmu1_pt);
  t1->Branch("genmu1_eta", &genmu1_eta);
  t1->Branch("genmu1_phi", &genmu1_phi);
  t1->Branch("genmu2_pt", &genmu2_pt);
  t1->Branch("genmu2_eta", &genmu2_eta);
  t1->Branch("genmu2_phi", &genmu2_phi);
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++)
    t1->Branch(TString(HLTs[ihlt]), &trigger[ihlt]);
  // Tag specific
  t1->Branch("tag_pt", &tag_pt);
  t1->Branch("tag_eta", &tag_eta);
  t1->Branch("tag_phi", &tag_phi);
  t1->Branch("tag_isLoose", &tag_isLoose);
  t1->Branch("tag_isMedium", &tag_isMedium);
  t1->Branch("tag_isTight", &tag_isTight);
  t1->Branch("tag_isSoft", &tag_isSoft);
  t1->Branch("tag_isHighPt", &tag_isHighPt);
  t1->Branch("tag_relIso04", &tag_relIso04);
  t1->Branch("tag_isMatchedGen", &tag_isMatchedGen);
  // Probe specific
  t1->Branch("iprobe", &iprobe);
  t1->Branch("probe_pt", &probe_pt);
  t1->Branch("probe_eta", &probe_eta);
  t1->Branch("probe_phi", &probe_phi);
  t1->Branch("probe_isLoose", &probe_isLoose);
  t1->Branch("probe_isMedium", &probe_isMedium);
  t1->Branch("probe_isTight", &probe_isTight);
  t1->Branch("probe_isSoft", &probe_isSoft);
  t1->Branch("probe_isHighPt", &probe_isHighPt);
  t1->Branch("probe_isMuMatched", &probe_isMuMatched);
  t1->Branch("probe_isPF", &probe_isPF);
  t1->Branch("probe_isSA", &probe_isSA);
  t1->Branch("probe_isTracker", &probe_isTracker);
  t1->Branch("probe_isGlobal", &probe_isGlobal);
  t1->Branch("probe_isdSA", &probe_isdSA);
  t1->Branch("probe_isdGlobal", &probe_isdGlobal);
  t1->Branch("probe_isCosmic", &probe_isCosmic);
  //  t1->Branch("probe_isGood", &probe_isGood);
  t1->Branch("probe_isHighPurity", &probe_isHighPurity);
  t1->Branch("probe_validFraction", &probe_validFraction);
  t1->Branch("probe_trkChi2", &probe_trkChi2);
  t1->Branch("probe_positionChi2", &probe_positionChi2);
  t1->Branch("probe_trkKink", &probe_trkKink);
  // t1->Branch("probe_segmentCompatibility", &probe_segmentCompatibility);
  t1->Branch("probe_trackerLayers", &probe_trackerLayers);
  t1->Branch("probe_pixelLayers", &probe_pixelLayers);
  t1->Branch("probe_muonStations", &probe_muonStations);
  t1->Branch("probe_muonHits", &probe_muonHits);
  t1->Branch("probe_DTHits", &probe_DTHits);
  t1->Branch("probe_CSCHits", &probe_CSCHits);
  t1->Branch("probe_pterr", &probe_pterr);
  t1->Branch("probe_dxy", &probe_dxy);
  t1->Branch("probe_dz", &probe_dz);
  t1->Branch("probe_relIso04", &probe_relIso04);
  t1->Branch("probe_isMatchedGen", &probe_isMatchedGen);
  t1->Branch("probe_dsa_muonStations", &probe_dsa_muonStations);
  t1->Branch("probe_dsa_muonHits", &probe_dsa_muonHits);
  t1->Branch("probe_dsa_DTHits", &probe_dsa_DTHits);
  t1->Branch("probe_dsa_CSCHits", &probe_dsa_CSCHits);
  t1->Branch("probe_dsa_pterr", &probe_dsa_pterr);
  t1->Branch("probe_dsa_dxy", &probe_dsa_dxy);
  t1->Branch("probe_dsa_dz", &probe_dsa_dz);
  t1->Branch("probe_dsa_trkChi2", &probe_dsa_trkChi2);
  t1->Branch("probe_dsa_pt", &probe_dsa_pt);
  t1->Branch("probe_dsa_eta", &probe_dsa_eta);
  t1->Branch("probe_dsa_phi", &probe_dsa_phi);
  // Pair specific
  t1->Branch("pair_pt", &pair_pt);
  t1->Branch("pair_eta", &pair_eta);
  t1->Branch("pair_phi", &pair_phi);
  t1->Branch("pair_mass", &pair_mass);
  t1->Branch("pair_fit_mass", &pair_fit_mass);
  t1->Branch("pair_svprob", &pair_svprob);
  t1->Branch("pair_dz", &pair_dz);
}

void NtupleContent::CreateExtraTrgBranches(
    const std::vector<std::string> &HLTs) {
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++)
    t1->Branch(TString("probe_" + HLTs[ihlt]), &probe_trg[ihlt]);
}

void NtupleContent::ClearBranches() {
  run = -1;
  event = -1;
  ls = -1;
  BSpot_x = -99;
  BSpot_y = -99;
  BSpot_z = -99;
  pv_x = -99;
  pv_y = -99;
  pv_z = -99;
  nvertices = 0;
  trueNumInteractions = -1.0;
  puNumInteractions = -1;
  Rho = -1;
  nmuons = 0;
  ntag = 0;

  for (unsigned int itrg = 0; itrg < 10; itrg++) trigger[itrg] = false;

  for (unsigned int itrg = 0; itrg < 10; itrg++) probe_trg[itrg] = false;
  
  // Gens
  genmu1_pt = 0;
  genmu1_eta = -99;
  genmu1_phi = -99;
  genmu2_pt = 0;
  genmu2_eta = -99;
  genmu2_phi = -99;

  trg_pt.clear();
  trg_eta.clear();
  trg_phi.clear();
  prb_pt.clear();
  prb_eta.clear();
  prb_phi.clear();

  tag_pt = 0;
  tag_eta = -99;
  tag_phi = -99;
  tag_isLoose = false;
  tag_isMedium = false;
  tag_isTight = false;
  tag_isSoft = false;
  tag_isHighPt = false;
  tag_relIso04 = -99;
  tag_isMatchedGen = false;

  iprobe = 0;
  probe_pt = 0;
  probe_eta = -99;
  probe_phi = -99;
  probe_isLoose = false;
  probe_isMedium = false;
  probe_isTight = false;
  probe_isSoft = false;
  probe_isHighPt = false;
  probe_isMuMatched = false;
  probe_isPF = false;
  probe_isSA = false;
  probe_isTracker = false;
  probe_isGlobal = false;
  probe_isdSA = false;
  probe_isdGlobal = false;
  probe_isCosmic = false;
  probe_isGood = false;
  probe_isHighPurity = false;
  probe_validFraction = -99;
  probe_trkChi2 = -99;
  probe_positionChi2 = -99;
  probe_trkKink = -99;
  probe_segmentCompatibility = -99;
  probe_trackerLayers = -99;
  probe_pixelLayers = -99;
  probe_muonStations = -99;
  probe_muonHits = -99;
  probe_DTHits = -99;
  probe_CSCHits = -99;
  probe_pterr = 0;
  probe_dxy = -99;
  probe_dz = -99;
  probe_relIso04 = -99;
  probe_isMatchedGen = false;
  probe_dsa_muonStations = -99;
  probe_dsa_muonHits = -99;
  probe_dsa_DTHits = -99;
  probe_dsa_CSCHits = -99;
  probe_dsa_pterr = 0;
  probe_dsa_dxy = -99;
  probe_dsa_dz = -99;
  probe_dsa_trkChi2 = -99;
  probe_dsa_pt = -99;
  probe_dsa_eta = -99;
  probe_dsa_phi = -99;

  pair_pt = 0;
  pair_mass = 0;
  pair_eta = -99;
  pair_phi = -99;
  pair_fit_mass = 0;
  pair_svprob = 0;
  pair_dz = -99;
}
