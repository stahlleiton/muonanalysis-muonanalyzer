#include "MuonGenAnalyzer.h"

MuonGenAnalyzer::MuonGenAnalyzer(){};

MuonGenAnalyzer::~MuonGenAnalyzer(){};

void MuonGenAnalyzer::SetInputs(const edm::Event& iEvent,
                                const edm::EDGetTokenT<edm::View<reco::GenParticle> >& gens_,
                                const int& momPdg_) {
  iEvent.getByToken(gens_, gens);
  std::vector<TLorentzVector> gmuons;
  std::vector<int> gcharge;
  for (const auto& gen : *gens) {
    if (fabs(gen.pdgId()) != 13)
      continue;
    if (fabs(gen.mother()->pdgId()) != momPdg_)
      continue;
    TLorentzVector temp;
    temp.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), MU_MASS);
    gmuons.push_back(temp);
    gcharge.push_back(gen.charge());
  }
  if (gmuons.size() == 2) {
    if (gmuons[0].Pt() > gmuons[1].Pt()) {
      gmuon1 = gmuons[0];
      gmuon2 = gmuons[1];
      gcharge1 = gcharge[0];
      gcharge2 = gcharge[1];
    } else {
      gmuon1 = gmuons[1];
      gmuon2 = gmuons[0];
      gcharge1 = gcharge[1];
      gcharge2 = gcharge[0];
    }
  } else {
    std::cout << "Warning the decay " << momPdg_ << " 2 muons not found. Gen branches will remain empty" << std::endl;
  }
}

void MuonGenAnalyzer::FillNtuple(NtupleContent& nt) {
  if (success) {
    nt.genmu1_pt = gmuon1.Pt();
    nt.genmu1_eta = gmuon1.Eta();
    nt.genmu1_phi = gmuon1.Phi();
    nt.genmu2_pt = gmuon2.Pt();
    nt.genmu2_eta = gmuon2.Eta();
    nt.genmu2_phi = gmuon2.Phi();
  } else {
    nt.genmu1_pt = 0;
    nt.genmu1_eta = -99;
    nt.genmu1_phi = -99;
    nt.genmu2_pt = 0;
    nt.genmu2_eta = -99;
    nt.genmu2_phi = -99;
  }
}
