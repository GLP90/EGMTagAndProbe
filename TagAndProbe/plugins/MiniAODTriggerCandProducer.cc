//#include <cmath>
// system include files
//#include <memory>
//#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"

//#include "DataFormats/Common/interface/Ptr.h"
//#include "DataFormats/Common/interface/PtrVector.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include <DataFormats/Math/interface/deltaR.h>

//#include <set>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
#include <string>

template <class T>
class MiniAODTriggerCandProducer : public edm::EDProducer {

  typedef std::vector<T> TCollection;
  typedef edm::Ref<TCollection> TRef;
  typedef edm::RefVector<TCollection> TRefVector;

public:
  explicit MiniAODTriggerCandProducer(const edm::ParameterSet&);
  ~MiniAODTriggerCandProducer();

  void init(const edm::TriggerResults &result, const edm::TriggerNames & triggerNames);
  bool onlineOfflineMatching(const edm::TriggerNames & triggerNames, 
			     edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, 
			     math::XYZTLorentzVector p4, std::string filterLabel, float dRmin);

  //bool L1Matching(edm::Handle<edm::View<l1extra::L1EmParticle>> l1H, math::XYZTLorentzVector cand, float ptThreshold); 
  std::vector<math::XYZTLorentzVector> getHLTP4(const edm::TriggerNames & triggerNames, 
						edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, 
						std::vector<std::string> filterLabels);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  
  std::string filterName_;
  edm::EDGetTokenT<TCollection> inputs_;
  //edm::EDGetTokenT<edm::View<l1extra::L1EmParticle>> l1iso_;
  //edm::EDGetTokenT<edm::View<l1extra::L1EmParticle>> l1noniso_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  //edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  float dRMatch_;
  edm::ParameterSetID triggerNamesID_;
  std::map<std::string, unsigned int> trigger_indices;
  //std::map<std::string, unsigned int> prescales;
  //std::map<std::string, unsigned int> prescale_counter;
};
  
template <class T>
MiniAODTriggerCandProducer<T>::MiniAODTriggerCandProducer(const edm::ParameterSet& iConfig ):
  filterName_(iConfig.getParameter<std::string>("filterName")),
  inputs_(consumes<TCollection>(iConfig.getParameter<edm::InputTag>("inputs"))),
  //l1iso_(consumes<edm::View<l1extra::L1EmParticle>>(iConfig.getParameter<edm::InputTag>("l1Iso"))),
  //l1noniso_(consumes<edm::View<l1extra::L1EmParticle>>(iConfig.getParameter<edm::InputTag>("l1NonIso"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  //triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))) {
  dRMatch_(iConfig.getParameter<double>("dR")) {
  
  produces<TRefVector>();
}

template <class T>
MiniAODTriggerCandProducer<T>::~MiniAODTriggerCandProducer()
{}

template <class T>
void MiniAODTriggerCandProducer<T>::init(const edm::TriggerResults &result, const edm::TriggerNames & triggerNames) {
  
  trigger_indices.clear();
  for (unsigned int i=0; i<triggerNames.triggerNames().size(); i++) {    
    std::string trimmedName = HLTConfigProvider::removeVersion(triggerNames.triggerName(i));
    trigger_indices[trimmedName] = triggerNames.triggerIndex(triggerNames.triggerName(i));    
  }
}

template <class T>
void MiniAODTriggerCandProducer<T>::produce(edm::Event &iEvent, const edm::EventSetup &eventSetup) {

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  
  edm::Handle<TCollection> inputs;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(inputs_, inputs);

  // Create the output collection
  std::auto_ptr<TRefVector> outColRef(new TRefVector);
  
  if (!triggerBits.isValid()) {
    LogDebug("") << "TriggerResults product not found - returning result=false!";
    return;
  }

  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerBits);
  if (triggerNamesID_ != triggerNames.parameterSetID()) {
    triggerNamesID_ = triggerNames.parameterSetID();
    init(*triggerBits, triggerNames);
  } 
  
  for (size_t i=0; i<inputs->size(); i++) {
    TRef ref(inputs, i);
    if (onlineOfflineMatching(triggerNames, triggerObjects, ref->p4(), filterName_, dRMatch_)) {
      outColRef->push_back(ref);
    }
  }	  

  iEvent.put(outColRef);
}

template <class T>
std::vector<math::XYZTLorentzVector> MiniAODTriggerCandProducer<T>::getHLTP4(const edm::TriggerNames & triggerNames, 
									     edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, 
									     std::vector<std::string> filterLabels) {
  
  std::vector< math::XYZTLorentzVector> triggerCandidates;
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackPathNames(triggerNames); 
    for (std::string filter : filterLabels) {
      if (obj.hasFilterLabel(filter)) 
	triggerCandidates.push_back(obj.p4());
    }
  }
  
  return triggerCandidates;
}

template <class T>
bool MiniAODTriggerCandProducer<T>::onlineOfflineMatching(const edm::TriggerNames & triggerNames, 
							  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, 
							  math::XYZTLorentzVector p4, std::string filterLabel, float dRmin) {
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackPathNames(triggerNames); 
    if (obj.hasFilterLabel(filterLabel)) {
      //for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];  
      float dR = deltaR(p4, obj.p4());
      if (dR < dRmin)
	return true;
    }
  }

  return false;
}

typedef MiniAODTriggerCandProducer<pat::Electron> PatElectronTriggerCandProducer;
DEFINE_FWK_MODULE(PatElectronTriggerCandProducer);

typedef MiniAODTriggerCandProducer<pat::Photon> PatPhotonTriggerCandProducer;
DEFINE_FWK_MODULE(PatPhotonTriggerCandProducer);
