/****************************************************************************
 *
 * This is a part of PPS offline software.
 * Authors:
 *   Laurent Forthomme (laurent.forthomme@cern.ch)
 *   Nicola Minafra (nicola.minafra@cern.ch)
 *
 ****************************************************************************/

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"

#include "DataFormats/CTPPSDetId/interface/TotemT2DetId.h"
#include "DataFormats/CTPPSDigi/interface/TotemT2Digi.h"
#include "DataFormats/CTPPSReco/interface/TotemT2RecHit.h"

#include "RecoPPS/Local/interface/TotemT2RecHitProducerAlgorithm.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "CondFormats/DataRecord/interface/PPSTimingCalibrationRcd.h"

class TotemT2RecHitProducer : public edm::stream::EDProducer<> {
public:
  explicit TotemT2RecHitProducer(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<edm::DetSetVector<TotemT2Digi> > digiToken_;

  /// Label to timing calibration tag
  edm::ESInputTag timingCalibrationTag_;
  /// A watcher to detect timing calibration changes.
  edm::ESWatcher<PPSTimingCalibrationRcd> calibWatcher_;

  bool applyCalib_;
  TotemT2RecHitProducerAlgorithm algo_;
};

TotemT2RecHitProducer::TotemT2RecHitProducer(const edm::ParameterSet& iConfig)
    : digiToken_(consumes<edm::DetSetVector<TotemT2Digi> >(iConfig.getParameter<edm::InputTag>("digiTag"))),
      timingCalibrationTag_(iConfig.getParameter<std::string>("timingCalibrationTag")),
      applyCalib_(iConfig.getParameter<bool>("applyCalibration")),
      algo_(iConfig) {
  produces<edm::DetSetVector<TotemT2RecHit> >();
}

void TotemT2RecHitProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto pOut = std::make_unique<edm::DetSetVector<TotemT2RecHit> >();

  // get the digi collection
  edm::Handle<edm::DetSetVector<TotemT2Digi> > digis;
  iEvent.getByToken(digiToken_, digis);

  if (!digis->empty()) {
    if (applyCalib_ && calibWatcher_.check(iSetup)) {
      edm::ESHandle<PPSTimingCalibration> hTimingCalib;
      iSetup.get<PPSTimingCalibrationRcd>().get(timingCalibrationTag_, hTimingCalib);
      algo_.setCalibration(*hTimingCalib);
    }
    // get the geometry
    edm::ESHandle<CTPPSGeometry> geometry;
    iSetup.get<VeryForwardRealGeometryRecord>().get(geometry);

    // produce the rechits collection
    algo_.build(*geometry, *digis, *pOut);
  }

  iEvent.put(std::move(pOut));
}

void TotemT2RecHitProducer::fillDescriptions(edm::ConfigurationDescriptions& descr) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("digiTag", edm::InputTag("totemT2Digis", "TotemT2"))
      ->setComment("input digis collection to retrieve");
  desc.add<std::string>("timingCalibrationTag", "GlobalTag:TotemT2TimingCalibration")
      ->setComment("input tag for timing calibrations retrieval");
  desc.add<double>("timeSliceNs", 25.0 / 1024.0)
      ->setComment("conversion constant between timing bin size and nanoseconds");
  desc.add<bool>("applyCalibration", true)->setComment("switch on/off the timing calibration");

  descr.add("totemT2RecHits", desc);
}

DEFINE_FWK_MODULE(TotemT2RecHitProducer);
