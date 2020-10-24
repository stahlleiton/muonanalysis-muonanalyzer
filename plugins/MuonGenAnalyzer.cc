#include "MuonGenAnalyzer.h"

MuonGenAnalyzer::MuonGenAnalyzer(){};

MuonGenAnalyzer::~MuonGenAnalyzer(){};

void MuonGenAnalyzer::SetInputs(
    const edm::Event& iEvent,
    const edm::EDGetTokenT<reco::GenParticleCollection>& gens_,
    const int& momPdg_) {
  iEvent.getByToken(gens_, gens);
  for (const auto& gen : *gens) {
    if (fabs(gen.pdgId()) != 13) continue;
    if (fabs(gen.mother()->pdgId()) != momPdg_) continue;
    gmuons.push_back(reco::GenParticleRef(gens, &gen - &gens->at(0)));
  }
  if (gmuons.size() != 2)
    std::cout << "Warning the decay " << momPdg_
              << " 2 muons not found. Gen branches will remain empty"
              << std::endl;
}

void MuonGenAnalyzer::FillNtuple(NtupleContent& nt) {
  if (gmuons.size() != 2) return;
  const bool isOrder = gmuons[0]->pt() > gmuons[1]->pt();
  const auto& gmuon1 = isOrder ? gmuons[0] : gmuons[1];
  const auto& gmuon2 = isOrder ? gmuons[1] : gmuons[0];
  nt.genmu1_pt = gmuon1->pt();
  nt.genmu1_eta = gmuon1->eta();
  nt.genmu1_phi = gmuon1->phi();
  nt.genmu2_pt = gmuon2->pt();
  nt.genmu2_eta = gmuon2->eta();
  nt.genmu2_phi = gmuon2->phi();
}
