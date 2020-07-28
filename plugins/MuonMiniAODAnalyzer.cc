// Package:    HLTAnalysis/TriggerAnalyzer
// Class:      TriggerAnalyzer
// 
/**\class TriggerAnalyzer TriggerAnalyzer.cc HLTAnalysis/TriggerAnalyzer/plugins/TriggerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2017 17:40:23 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include <vector>
#include "TTree.h"
#include <string>
#include <iostream>
#include "DataFormats/Common/interface/Ref.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "NtupleContent.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "helper.h"
#include "KlFitter.h"
#include "MuonBranches.h"
#include "MuonGenAnalyzer.h"

/*namespace edm {
  class ConfigurationDescriptions;
  }*/



using namespace std;
//using namespace edm;


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
class MuonMiniAODAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  
public:
  typedef std::vector<std::pair<pat::Muon,reco::TransientTrack>> RecoTrkAndTransientTrkCollection;
  explicit MuonMiniAODAnalyzer(const edm::ParameterSet&);
  ~MuonMiniAODAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  bool HLTaccept( const edm::Event& ,NtupleContent & ,std::vector<std::string>&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
 

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken PFCands_;
  edm::EDGetToken LostTracks_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_; 
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genToken_;
  std::vector<std::string> HLTPaths_;
  std::vector<std::string> ProbePaths_;
  const unsigned int tagQual_;
  const StringCutObjectSelector<pat::Muon>  tagSelection_; //kinematic cuts for tag
  const bool HighPurity_;
  const StringCutObjectSelector<pat::PackedCandidate>  probeSelection_; //kinematic cuts for probe
  const double pairMassMin_;
  const double pairMassMax_;
  const double pairDz_;
  const bool RequireVtxCreation_; //if true skip pairs that do not create gthat do not have a vertex
  const double minSVtxProb_; //min probability of a vertex to be kept. If <0 inactive
  const double maxdz_trk_mu_;
  const double maxpt_relative_dif_trk_mu_;
  const double maxdr_trk_mu_;
  const unsigned momPdgId_;
  const double genRecoDrMatch_;


  edm::Service<TFileService> fs;
  TTree * t1;
  NtupleContent nt;
 
    // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonMiniAODAnalyzer::MuonMiniAODAnalyzer(const edm::ParameterSet& iConfig): 
 beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter <edm::InputTag>("beamSpot"))),
 vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
 muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  PFCands_(consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("PFCands"))),
 LostTracks_(consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("lostTracks"))),
trgresultsToken_(consumes<edm::TriggerResults >(iConfig.getParameter<edm::InputTag>("triggerResults"))),
 genToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("gen"))),
 HLTPaths_(iConfig.getParameter<std::vector<std::string>>("HLTPaths")),
 ProbePaths_(iConfig.getParameter<std::vector<std::string>>("ProbePaths")),
 tagQual_(iConfig.getParameter<unsigned>("tagQuality")),
 tagSelection_(iConfig.getParameter<std::string>("tagSelection")),
 HighPurity_(iConfig.getParameter<bool>("ProbeHPyrity")),
 probeSelection_(iConfig.getParameter<std::string>("probeSelection")),
 pairMassMin_(iConfig.getParameter<double>("pairMassMin")),
 pairMassMax_(iConfig.getParameter<double>("pairMassMax")),
 pairDz_(iConfig.getParameter<double>("pairDz")),
 RequireVtxCreation_(iConfig.getParameter<bool>("RequireVtxCreation")),
 minSVtxProb_(iConfig.getParameter<double>("minSVtxProb")),
 maxdz_trk_mu_(iConfig.getParameter<double>("maxDzProbeTrkMuon")),
 maxpt_relative_dif_trk_mu_(iConfig.getParameter<double>("maxRelPtProbeTrkMuon")),
 maxdr_trk_mu_(iConfig.getParameter<double>("maxDRProbeTrkMuon")),
 momPdgId_(iConfig.getParameter<unsigned>("momPdgId")),
 genRecoDrMatch_(iConfig.getParameter<double>("genRecoDrMatch"))

 
{
//  edm::ParameterSet runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");
}

MuonMiniAODAnalyzer::~MuonMiniAODAnalyzer()
{
  // cout<<"total "<<trg_counter<<" fires "<<fire_counter<<" l3 "<<l3_counter<<endl;
   // do anything here that needs to be done at desctruction time

}


//
// member functions
bool
MuonMiniAODAnalyzer::HLTaccept( const edm::Event& iEvent, NtupleContent &nt,std::vector<std::string>& HLTPaths){
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  edm::TriggerNames trigName;
  trigName = iEvent.triggerNames(*trigResults);
  bool EvtFire=false;
  unsigned int ipath=0;
  for (auto path:HLTPaths){
   bool TrgFire=false;
   for( unsigned int itrg = 0; itrg < trigResults->size(); ++itrg ) {
     TString TrigPath =trigName.triggerName(itrg);
     if (!trigResults->accept(itrg)) continue;
     if (!TrigPath.Contains(path)) continue;
     EvtFire=true;
     TrgFire=true;
   }
   nt.trigger[ipath]=TrgFire;
   ipath++;
  }
 return EvtFire;
}
// ------------ method called for each event  ------------

void
MuonMiniAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  //skip evts if there are no vertices
  if (vertices->size()==0) return;
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_,muons);
  edm::Handle<std::vector<pat::PackedCandidate>> pfcands;
  iEvent.getByToken(PFCands_,pfcands);
  edm::Handle<std::vector<pat::PackedCandidate>>lostTracks;
  iEvent.getByToken(LostTracks_,lostTracks);
  edm::ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);

  //information about run
  nt.ClearBranches();
  nt.run=iEvent.id().run();         nt.ls=iEvent.luminosityBlock();
  nt.fromFullAOD=false;
  nt.BSpot_x= theBeamSpot->x0();    nt.BSpot_y= theBeamSpot->y0();
  nt.BSpot_z= theBeamSpot->z0();

  bool goodVtx=false;
  for ( const reco::Vertex & vtx: *vertices){
    if ( vtx.isFake() || !vtx.isValid() ) continue;
    nt.pv_x=vtx.x();  nt.pv_y=vtx.y();    nt.pv_z=vtx.z();
    goodVtx=true;
    break;
  } 
  if (!goodVtx) return; 
  if (!HLTaccept( iEvent, nt, HLTPaths_) ) return;
//  HLTaccept(iEvent, nt.doublemu_trg, DoubleMuPaths_);

   //gen information
  MuonGenAnalyzer genmu;
  std::vector<unsigned> matched_muon_idx;
  if(! iEvent.isRealData() ) {
    genmu.SetInputs(iEvent,genToken_,momPdgId_);
    genmu.FillNtuple(nt);
    auto reco_match_genmu1= MatchReco<pat::Muon> (*muons,nt.genmu1_eta, nt.genmu1_phi, genRecoDrMatch_);
    auto reco_match_genmu2= MatchReco<pat::Muon> (*muons,nt.genmu2_eta, nt.genmu2_phi, genRecoDrMatch_);
    if (reco_match_genmu1.first) matched_muon_idx.push_back(reco_match_genmu1.second);
    if (reco_match_genmu2.first) matched_muon_idx.push_back(reco_match_genmu2.second);
  }


  // find triggering muon
  RecoTrkAndTransientTrkCollection tag_muon_ttrack; 
  std::vector<bool> genmatched_tag;
  for( const pat::Muon &mu : *muons){
    if ( mu.passed(pow(2,tagQual_)) ) continue;
    bool fired=false;
    for ( const std::string path: HLTPaths_){
      char cstr[ (path+"*").size() + 1];
      strcpy( cstr, (path+"*").c_str() );
      if ( !mu.triggered( cstr ) ) continue;
      fired=true;
      break;
    }
    if (!fired) continue;
    if ( !tagSelection_(mu)) continue;
    tag_muon_ttrack.emplace_back(
                     std::make_pair(
                       mu, reco::TransientTrack(*mu.bestTrack(),&(*bField) )
                     ) 
                   );
    if ( std::find(matched_muon_idx.begin(),matched_muon_idx.end(),&mu-&muons->at(0)) != matched_muon_idx.end() )
      genmatched_tag.push_back(true); 
    else
      genmatched_tag.push_back(false);
  }
  if ( tag_muon_ttrack.size()==0) return;
  nt.nmuons=muons->size();
  nt.ntag=tag_muon_ttrack.size();

   
  // add Lost Tracks to Packed cands
  std::vector<pat::PackedCandidate> tracks;
  for (const auto container : {pfcands,lostTracks} ){
    for ( const pat::PackedCandidate & trk: *container){
      if( !probeSelection_(trk) ) continue;
      if( !trk.hasTrackDetails()) continue;
      if ( HighPurity_ && !trk.trackHighPurity() ) continue;
      tracks.emplace_back(trk);
    }
  }
  std::vector<unsigned> matched_track_idx;
  if (!iEvent.isRealData()){
    auto reco_match_genmu1= MatchReco<pat::PackedCandidate> (tracks,nt.genmu1_eta, nt.genmu1_phi, genRecoDrMatch_);
    auto reco_match_genmu2= MatchReco<pat::PackedCandidate> (tracks,nt.genmu2_eta, nt.genmu2_phi, genRecoDrMatch_);
    if (reco_match_genmu1.first) matched_track_idx.push_back(reco_match_genmu1.second);
    if (reco_match_genmu2.first) matched_track_idx.push_back(reco_match_genmu2.second);
  }
  std::pair<std::vector<unsigned>,std::vector<unsigned>> trk_muon_map;
  for (const auto &mu: *muons){
    float minDR=1000;
    unsigned int idx_trk;
    for (const auto &trk: tracks){
      if (mu.charge() != trk.charge()) continue;
      if (fabs(mu.vz()-trk.vz())>maxdz_trk_mu_) continue;
      if (fabs(mu.pt()-trk.pt())/mu.pt()>maxpt_relative_dif_trk_mu_) continue;
      float DR=deltaR(mu.eta(),mu.phi(),trk.eta(),trk.phi());
      if (minDR<DR) continue;
      minDR=DR;
      idx_trk=&trk -&tracks[0];
    }
    if (minDR>maxdr_trk_mu_) continue;
    trk_muon_map.first.push_back(idx_trk); trk_muon_map.second.push_back(&mu-&muons->at(0));
  }
    

  // final pair selection
  for ( const auto &tag: tag_muon_ttrack ){
    for (const auto &probe: tracks ){
      if ( tag.first.charge() == probe.charge()) continue;
      if (fabs(tag.first.vz() - probe.vz())> pairDz_ ) continue;
      float mass = DimuonMass(tag.first.pt(), tag.first.eta(), tag.first.phi(),
                              probe.pt(), probe.eta(), probe.phi()
                              );
      
      if ( mass<pairMassMin_ || mass>pairMassMax_) continue;
      std::vector<reco::TransientTrack> trk_pair={tag.second,
                                      reco::TransientTrack(probe.pseudoTrack(),&(*bField)) };

      KlFitter vtx(trk_pair);
      if (RequireVtxCreation_ &&! vtx.status()) continue;
      if (minSVtxProb_>0 && vtx.prob()<minSVtxProb_) continue;

      math::PtEtaPhiMLorentzVector mu1(tag.first.pt(), tag.first.eta(), tag.first.phi(), MU_MASS);
      math::PtEtaPhiMLorentzVector mu2(probe.pt(), probe.eta(), probe.phi(), MU_MASS);

      FillTagBranches<pat::Muon,pat::PackedCandidate> (tag.first,tracks,nt);
      nt.tag_isMatchedGen=genmatched_tag[&tag-&tag_muon_ttrack[0]];

      std::vector<unsigned>::iterator it = std::find( trk_muon_map.first.begin(), trk_muon_map.first.end(), &probe - &tracks[0]);
      if ( it != trk_muon_map.first.end() ){
         unsigned idx=std::distance(trk_muon_map.first.begin(),it);
         FillProbeBranches<pat::Muon,pat::PackedCandidate> (muons->at(trk_muon_map.second[idx]),tracks,nt,true);
       /*  for ( const std::string path: ProbePaths_){
           char cstr[ (path+"*").size() + 1];
           strcpy( cstr, (path+"*").c_str() );
           nt.probe_trg[&path-&ProbePaths_[0]]=muons->at(trk_muon_map.second[idx]).triggered(cstr);
         }*/
         
      } else{
         reco::Muon fakeMuon;
         fakeMuon.setP4(mu2);
         FillProbeBranches<reco::Muon,pat::PackedCandidate> (fakeMuon,tracks,nt,false);
       //  for ( const std::string path: ProbePaths_) nt.probe_trg[&path-&ProbePaths_[0]]=false;
      }
      nt.iprobe++;
      nt.probe_isHighPurity=probe.trackHighPurity();
      FillPairBranches<pat::Muon,pat::PackedCandidate>(tag.first,probe,nt);
      if (std::find(matched_track_idx.begin(),matched_track_idx.end(),&probe-&tracks[0]) != matched_track_idx.end() )
        nt.probe_isMatchedGen=true;
     else nt.probe_isMatchedGen=false;
      t1->Fill();
    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonMiniAODAnalyzer::beginJob()
{
 t1=fs->make<TTree>("tree","tree");
 nt.SetTree(t1);
 nt.CreateBranches(HLTPaths_);
// if (ProbePaths_.size()>0) nt.CreateExtraTrgBranches(ProbePaths_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonMiniAODAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonMiniAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


///////////////////////
  
//define this as a plug-in
DEFINE_FWK_MODULE(MuonMiniAODAnalyzer);

