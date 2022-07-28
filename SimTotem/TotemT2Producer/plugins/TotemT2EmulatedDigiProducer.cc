/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Author:
 *   Arkadiusz Cwikla
 *
 ****************************************************************************/

#include <memory>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "DataFormats/CTPPSDetId/interface/TotemT2DetId.h"
#include "DataFormats/TotemReco/interface/TotemT2Digi.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

class TotemT2EmulatedDigiProducer : public edm::stream::EDProducer<> {
public:
  explicit TotemT2EmulatedDigiProducer(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // CONSTANTS
  // only channel 4 contains meaningful data
  const short channelNumber_ = 4;
  // dummy chars just to satisfy TotemT2Digi constructor
  const unsigned char dummyGeo_ = 'g';
  const unsigned char dummyId_ = 'i';
  const unsigned char dummyMarker_ = 'm';
  // conversion constant between time slice and absolute time (in ns) - but inverted
  const double timeSliceNsInverted_ = 1024.0 / 25.0;
  // conversion constant between seconds and nanoseconds
  const double sToNs_ = std::pow(10, 9);

  // variables to read data from branches
  int channelNo_;
  double leadingEdge_;
  double trailingEdge_;

  std::string t2DataFile_;
  std::vector<std::vector<unsigned int>> testCasesVec_;
  std::unique_ptr<TChain> inputTree_;
};

TotemT2EmulatedDigiProducer::TotemT2EmulatedDigiProducer(const edm::ParameterSet& iConfig)
    : t2DataFile_(iConfig.getParameter<std::string>("t2DataFile")) {
  for (const auto& idsAndLimits : iConfig.getParameter<std::vector<edm::ParameterSet>>("testCasesSet")) {
    std::vector<unsigned int> ids = idsAndLimits.getParameter<std::vector<unsigned int>>("detId");
    unsigned int limit = idsAndLimits.getParameter<unsigned int>("eventLimit");

    testCasesVec_.push_back(ids);
    for (auto id : ids) {
      if (event_counters_.find(id) == event_counters_.end())
        event_counters_[id] = 0;
      event_counters_[id] += limit;
    }
  }

  inputTree_ = std::make_unique<TChain>("tree");

  inputTree_->Add(t2DataFile_.c_str());

  inputTree_->SetBranchAddress("channelNo", &channelNo_);
  inputTree_->SetBranchAddress("LeadingEdge", &leadingEdge_);
  inputTree_->SetBranchAddress("TrailingEdge", &trailingEdge_);

  produces<edmNew::DetSetVector<TotemT2Digi>>("TotemT2");
}

void TotemT2EmulatedDigiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto digi = std::make_unique<edmNew::DetSetVector<TotemT2Digi>>();

  inputTree_->GetEntry(iEvent.id().event());

  if (channelNo_ == channelNumber_) {
    unsigned short le = (unsigned short)floor(leadingEdge_ * sToNs_ * timeSliceNsInverted_);
    unsigned short te = (unsigned short)floor(trailingEdge_ * sToNs_ * timeSliceNsInverted_);
    TotemT2Digi digiTmp(dummyGeo_, dummyId_, dummyMarker_, le, te);

    for (auto vec : testCasesVec_) {
      for (const auto& id : vec) {
        // check if limit was reached
        if (event_counters_[id] > 0) {
          TotemT2DetId detId(id);
          edmNew::DetSetVector<TotemT2Digi>::FastFiller(*digi, detId).push_back(digiTmp);
        }
      }
    }
    iEvent.put(std::move(digi), "TotemT2");
  }

  void TotemT2EmulatedDigiProducer::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("t2DataFile")->setComment("path to fake T2 root data");

    edm::ParameterSetDescription idmap_valid;
    idmap_valid.add<std::vector<unsigned int>>("detId")->setComment("mapped TotemT2DetId's for some test case");

    desc.addVPSet("testCasesSet", idmap_valid);

    descriptions.add("totemT2EmulatedDigis", desc);
  }

  DEFINE_FWK_MODULE(TotemT2EmulatedDigiProducer);
