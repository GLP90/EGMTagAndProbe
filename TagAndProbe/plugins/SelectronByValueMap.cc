#ifndef SelectorByValueMap_h
#define SelectorByValueMap_h

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include <vector>

template <typename C>
struct CompatibleConfigurationType {
  typedef C type;
};
    
// "float" is not allowed, as it conflicts with "double"
template <>
struct CompatibleConfigurationType<float> {
  typedef double type;
};
  
template <typename T, typename C>
class SelectorByValueMap : public edm::EDProducer {
public:
  explicit SelectorByValueMap(edm::ParameterSet const & config);
  
private:
  typedef T candidate_type;
  typedef C selection_type;
  typedef typename CompatibleConfigurationType<selection_type>::type cut_type;
  typedef std::vector<candidate_type> candidateCollection;
  typedef edm::Ref<candidateCollection> candidateRef;
  typedef edm::RefVector<candidateCollection> candidateRefVector;

  void produce(edm::Event & event, edm::EventSetup const & setup);
    
  edm::EDGetTokenT<candidateCollection>    token_inputs;
  edm::EDGetTokenT<edm::ValueMap<selection_type> >  token_selection;
  cut_type m_cut;

  //edm::EDGetTokenT<std::vector<pat::Electron> >    token_electrons;
  //edm::EDGetTokenT<edm::ValueMap<bool> >  token_selection;    
  //bool m_cut;
};
  
template <typename T, typename C>
SelectorByValueMap<T,C>::SelectorByValueMap(edm::ParameterSet const & config) :
  token_inputs(consumes<candidateCollection>(config.getParameter<edm::InputTag>("input"))),
  token_selection(consumes<edm::ValueMap<selection_type> >(config.getParameter<edm::InputTag>("selection"))),
  m_cut(config.getParameter<cut_type>("cut")) {
  //token_electrons(consumes<std::vector<pat::Electron> >(config.getParameter<edm::InputTag>("input"))),
  //token_selection(consumes<edm::ValueMap<bool> >(config.getParameter<edm::InputTag>("selection"))),
  //m_cut(config.getParameter<bool>("cut"))
  //{
  //produces<edm::RefToBaseVector<candidate_type> >();
  //produces<edm::RefVector<candidate_type> >();
  //produces<edm::RefVector<pat::Electron> >();
  produces<candidateRefVector>();
}
  
template <typename T, typename C>
void SelectorByValueMap<T, C>::produce(edm::Event & event, const edm::EventSetup & setup) {
//void SelectorByValueMap::produce(edm::Event & event, const edm::EventSetup & setup)
//{
  std::auto_ptr<candidateRefVector> candidates(new candidateRefVector());
//std::auto_ptr<pat::ElectronRefVector> candidates(new pat::ElectronRefVector());
  
//  // read the collection of GsfElectrons from the Event
//  edm::Handle<std::vector<candidate_type> > h_electrons;
//  event.getByToken(token_electrons, h_electrons);
//  std::vector<candidate_type> const & electrons = * h_electrons;
//  
//  // read the selection map from the Event
//  edm::Handle<edm::ValueMap<selection_type> > h_selection;
//  event.getByToken(token_selection, h_selection);
//  edm::ValueMap<selection_type> const & selectionMap = * h_selection;
//  
//  for (unsigned int i = 0; i < electrons.size(); ++i) {
//    edm::Ref<candidate_type> ptr = electrons[i];
//    if (selectionMap[ptr] > m_cut)
//      candidates->push_back(ptr);
//  }
//  
//  // put the product in the event
//  event.put(candidates);


  edm::Handle<candidateCollection> h_inputs;
  event.getByToken(token_inputs, h_inputs);
  //std::vector<pat::Electron> const & electrons = * h_electrons;
  
  // read the selection map from the Event
  edm::Handle<edm::ValueMap<selection_type> > h_selection;
  event.getByToken(token_selection, h_selection);
  edm::ValueMap<selection_type> const & selectionMap = * h_selection;
  
  for (unsigned int i = 0; i < h_inputs->size(); ++i) {
    candidateRef ptr(h_inputs, i);
    if (selectionMap[ptr] >= m_cut)
      candidates->push_back(ptr);
  }
  
  event.put(candidates);
}

//} // namespace edm;

#endif

typedef SelectorByValueMap<pat::Electron, bool> PatElectronSelectorByValueMap;
DEFINE_FWK_MODULE(PatElectronSelectorByValueMap);
