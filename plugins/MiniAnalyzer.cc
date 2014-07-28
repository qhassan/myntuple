// -*- C++ -*-
//
// Package:    UserCode/MiniAnalyzer
// Class:      MiniAnalyzer
// 
/**\class MiniAnalyzer MiniAnalyzer.cc Test/MiniAnalyzer/plugins/MiniAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Qamar Ul Hassan
//         Created:  Sun, 13 Jul 2014 06:22:18 GMT
//
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
 
#include <TTree.h>
#include <TClonesArray.h>
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>


using namespace edm;
using namespace std;
using namespace reco;
using namespace pat; 

//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAnalyzer(const edm::ParameterSet&);
  ~MiniAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
      
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;

  //TH1F* fHistnew_Histo;

  std::unordered_map<std::string,TH1F*> histContainer_;
  std::unordered_map<std::string,TH2F*> histContainer2d_; 

  //TH1F* fHistnew_Histo;

  bool muon_selection, electron_veto, jet_selection;

//std::vector<float> muonpt_;

TTree* tree_;

Int_t njets_;
Float_t number_of_vertices_[1000];
Float_t rho_[1000];
Float_t muonpt_[1000];
Float_t muoncharge_[1000];
Float_t muoneta_[1000];
Float_t muonphi_[1000];
Float_t jetpt_[1000];
Float_t jeteta_[1000];
Float_t jetphi_[1000];
Float_t jetcsv_[1000];
Float_t jetvtx_mass_[1000];
Float_t pileup_[1000];
Float_t metpt_[1000];
Float_t meteta_[1000];
Float_t metphi_[1000];
Float_t mt_[1000]; 
Float_t mt3_[1000];
Float_t ncsv_[1000];
Float_t nvtx_[1000];
Float_t mt41b_[1000];
Float_t mt42b_[1000];
Float_t met42b_pt_[1000];
Float_t met42b_eta_[1000];
Float_t met42b_phi_[1000];
Float_t mt4_[1000];

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
MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig) :
//  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),  
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
{
  //now do what ever initialization is needed

}


MiniAnalyzer::~MiniAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //   using namespace edm;

  muon_selection = false;
  electron_veto = false;
  jet_selection = false; 

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  /*
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    std::cout << "\n === TRIGGER PATHS === " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    std::cout << "Trigger " << names.triggerName(i) <<
    //                ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
    ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
    << std::endl;
    }
  */
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &primVtx = vertices->front();

number_of_vertices_[0] = vertices->size();

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  float rho=*rhoH;
		    
//  histContainer_["cutflow"]->Fill(0);


  //
  // MUONS
  // cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId  
  //
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  std::vector<const pat::Muon *> selectedMuons;        
  for (const pat::Muon &mu : *muons) { 
    if( mu.isPFMuon() 
	&& mu.isGlobalMuon() 
	&& mu.normChi2() < 10 
	&& mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 
	&& mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 
	&& mu.dB() < 0.2 
	&& mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 
	&& mu.numberOfMatchedStations() > 1 
	&& (mu.chargedHadronIso()+max(0.,mu.neutralHadronIso()+mu.photonIso()-0.50*mu.puChargedHadronIso()))/mu.pt() < 0.12 )
      {
	
	if( mu.pt() > 26 && fabs(mu.eta()) < 2.1 )
	  {
	    selectedMuons.push_back( &mu );
	    muonpt_[0] = mu.pt();
//	    muonpt_ .at(0)=(mu.pt());	    
//	    muonpt_.push_back(mu.pt());
	    muoncharge_[0] = mu.charge();
	    muoneta_[0] = mu.eta();	
	    muonphi_[0] = mu.phi();
	  }
      }
  }
/*
for(unsigned int i=0;i<muonpt_.size();i++)
{muonpt_ = muonpt_[i].pt();}
*/
  //require at least one muon
  if(selectedMuons.size()==1); //histContainer_["cutflow"]->Fill(1);
  else return;

  //
  // ELECTRONS
  // cf. https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification  
  //
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  std::vector<const pat::Electron *> selectedElectrons;
  for (const pat::Electron &el : *electrons) {        	

    //use a cut based id
    bool passVetoId( EgammaCutBasedEleId::PassWP( EgammaCutBasedEleId::VETO, el.isEB(), el.pt(), el.eta(),
                                                  el.deltaEtaSuperClusterTrackAtVtx(), 
						  el.deltaPhiSuperClusterTrackAtVtx(),
						  el.sigmaIetaIeta(),
						  el.hadronicOverEm(),
						  (1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy()),
						  fabs(el.gsfTrack()->dxy(primVtx.position())),
						  fabs(el.gsfTrack()->dz(primVtx.position())),
                                                  0., 0., 0., 
						  !(el.passConversionVeto()), 
						  el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(),
						  rho) );
    
    if( passVetoId 
	&& fabs(el.superCluster()->eta()) > 1.4442 && fabs(el.superCluster()->eta()) < 1.5660 
	&& el.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 
	&& el.dB() < 0.02 
	&& el.passConversionVeto() == true 
	&& (el.chargedHadronIso()+max(0.,el.neutralHadronIso()+el.photonIso()-0.50*el.puChargedHadronIso()))/el.pt() < 0.1 )
      {
	if(el.pt() > 20 && fabs(el.eta()) < 2.5)
	  {
	    selectedElectrons.push_back(&el);
	  }
      }
  }
  
  //require at least one electron
  if(selectedElectrons.size()==0); //histContainer_["cutflow"]->Fill(2);
  else return;

  //
  // JETS
  //
  uint32_t nCSVMtags(0);
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  std::vector<const pat::Jet *> selectedJets;
  for (const pat::Jet &j : *jets) {
	  
    float dR2muon=deltaR(j,*(selectedMuons[0]));
    float rawEnergy(j.energy()*j.jecFactor("Uncorrected"));
    if ( j.numberOfDaughters() > 1 
	 && (j.neutralHadronEnergy() + j.HFHadronEnergy())/rawEnergy < 0.99 
	 && j.neutralEmEnergyFraction() < 0.99 
	 && (j.chargedEmEnergyFraction() < 0.99 || fabs(j.eta()) >= 2.4)
	 && (j.chargedHadronEnergyFraction() > 0. || fabs(j.eta()) >= 2.4) 
	 && (j.chargedMultiplicity() > 0 || fabs(j.eta()) >= 2.4)
	 && dR2muon>0.5)
      {
	if( fabs(j.eta()) < 2.5 && j.pt() > 30)
	  {
	    selectedJets.push_back( &j );
	    float csv=j.bDiscriminator("combinedSecondaryVertexBJetTags");
	    if(csv>0.679) nCSVMtags++;
	    jetpt_[njets_] = j.pt();	    
	    jeteta_[njets_] = j.eta();
	    jetphi_[njets_] = j.phi();
	    jetcsv_[njets_] = j.bDiscriminator("combinedSecondaryVertexBJetTags");
	    jetvtx_mass_[njets_] = j.userFloat("vtxMass");
	    pileup_[njets_] = j.userFloat("pileupJetId:fullDiscriminant");	
	  }
      }
  }

  //
  // MET
  //
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
	metpt_[0] = mets->at(0).pt();
	meteta_[0] = mets->at(0).eta();
	metphi_[0] = mets->at(0).phi();
  float dphi_met_mu = deltaPhi(muonphi_[0], metphi_[0]); // use the function to restrict to the 0,pi range

mt_[njets_]=sqrt(2*muonpt_[0]*metpt_[0]*(1-cos(dphi_met_mu)));
//mt_=sqrt(2*mu.pt()*metpt_*(1-cos(dphi_met_mu)));

  //
  // FINAL SELECTION PLOTS
  //
  if(selectedJets.size()>=3) 
    {
// mt_[njets_] = mt3_;

    }
  if(selectedJets.size()>=4) 
    {
//mt_[njets_] = mt4_;
ncsv_[njets_] = nCSVMtags;
nvtx_[njets_] = number_of_vertices_[0];
      if(nCSVMtags>=1) 
	{
//mt_[njets_] = mt41b_;
	}
      if(nCSVMtags>=2) 
	{
//mt_[njets_] = mt42b_;
met42b_pt_[njets_] = metpt_[0];
met42b_eta_[njets_] = meteta_[0];
met42b_phi_[njets_] = metphi_[0];
	}
	tree_->Fill();      
    }
  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;

	tree_ = fs->make<TTree>("AnaTree", "AnaTree");

tree_->Branch("njets_", &njets_, "njets_/I"); 
tree_->Branch("number_of_vertices", &number_of_vertices_, "nvtx/F"); 
/*
tree_->SetBranchAddress("number_of_vertices", &number_of_vertices_); 
tree_->SetBranchAddress("rho", &rho_); 
tree_->SetBranchAddress("muonpt", &muonpt_); 
tree_->SetBranchAddress("muoncharge", &muoncharge_); 
tree_->SetBranchAddress("muoneta", &muoneta_); 
tree_->SetBranchAddress("muonphi", &muonphi_); 
tree_->SetBranchAddress("jetpt", &jetpt_); 
tree_->SetBranchAddress("jeteta", &jeteta_); 
tree_->SetBranchAddress("jetphi", &jetphi_); 
tree_->SetBranchAddress("jetcsv", &jetcsv_); 
tree_->SetBranchAddress("jetvtxmass", &jetvtx_mass_); 
tree_->SetBranchAddress("pileup", &pileup_); 
tree_->SetBranchAddress("metpt", &metpt_); 
tree_->SetBranchAddress("meteta", &meteta_); 
tree_->SetBranchAddress("metphi", &metphi_); 
tree_->SetBranchAddress("mt3jets", &mt3_); 
tree_->SetBranchAddress("ncsvjets", &ncsv_); 
tree_->SetBranchAddress("nvtx", &number_of_vertices_); 
tree_->SetBranchAddress("mt41b", &mt41b_); 
tree_->SetBranchAddress("mt42b", &mt42b_); 
tree_->SetBranchAddress("mt4jets", &mt4_); 
tree_->SetBranchAddress("met42b_pt", &met42b_pt_); 
tree_->SetBranchAddress("met42b_eta", &met42b_eta_); 
tree_->SetBranchAddress("met42b_phi", &met42b_phi_); 
*/
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  MiniAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  MiniAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  MiniAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  MiniAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
