/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Author:
 *   Laurent Forthomme
 *
 ****************************************************************************/

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/CTPPSDetId/interface/TotemT2DetId.h"
#include "DataFormats/TotemReco/interface/TotemT2Digi.h"
#include "DataFormats/TotemReco/interface/TotemT2RecHit.h"

#include "Geometry/Records/interface/TotemGeometryRcd.h"
#include "DQM/CTPPS/interface/TotemT2Segmentation.h"

#include <string>

class TotemT2DQMSource : public DQMEDAnalyzer {
public:
  TotemT2DQMSource(const edm::ParameterSet&);
  ~TotemT2DQMSource() override;

protected:
  void dqmBeginRun(const edm::Run&, const edm::EventSetup&) override;
  void bookHistograms(DQMStore::IBooker&, const edm::Run&, const edm::EventSetup&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  edm::ESGetToken<TotemGeometry, TotemGeometryRcd> geometryToken_;
  edm::EDGetTokenT<edm::DetSetVector<TotemT2Digi>> digiToken_;
  edm::EDGetTokenT<edm::DetSetVector<TotemT2RecHit>> rechitToken_;

  std::unique_ptr<TotemT2Segmentation> segm_;
  MonitorElement* m_digis_mult_[TotemT2DetId::maxPlane / 2];
  MonitorElement* m_rechits_mult_[TotemT2DetId::maxPlane / 2];
};

TotemT2DQMSource::TotemT2DQMSource(const edm::ParameterSet& iConfig)
    : geometryToken_(esConsumes<TotemGeometry, TotemGeometryRcd, edm::Transition::BeginRun>()),
      digiToken_(consumes<edm::DetSetVector<TotemT2Digi>>(iConfig.getParameter<edm::InputTag>("digisTag"))),
      rechitToken_(consumes<edm::DetSetVector<TotemT2RecHit>>(iConfig.getParameter<edm::InputTag>("rechitsTag"))) {}

TotemT2DQMSource::~TotemT2DQMSource() {}

void TotemT2DQMSource::dqmBeginRun(const edm::Run&, const edm::EventSetup&) {}

void TotemT2DQMSource::bookHistograms(DQMStore::IBooker& ibooker, const edm::Run&, const edm::EventSetup& iSetup) {
  ibooker.cd();
  ibooker.setCurrentFolder("CTPPS/TotemT2");

  const size_t summary_nbinsx = 50, summary_nbinsy = 50;

  for (unsigned int pl = 0; pl < TotemT2DetId::maxPlane / 2; ++pl) {
    m_digis_mult_[pl] = ibooker.book2D("digis multiplicity (plane " + std::to_string(pl) + ")",
                                       "x;y",
                                       summary_nbinsx,
                                       0.,
                                       summary_nbinsx,
                                       summary_nbinsy,
                                       0.,
                                       summary_nbinsy);
    m_rechits_mult_[pl] = ibooker.book2D("rechits multiplicity (plane " + std::to_string(pl) + ")",
                                         "x;y",
                                         summary_nbinsx,
                                         0.,
                                         summary_nbinsx,
                                         summary_nbinsy,
                                         0.,
                                         summary_nbinsy);
  }

  // build a segmentation helper for the size of histograms previously booked
  segm_ = std::make_unique<TotemT2Segmentation>(iSetup.getData(geometryToken_), summary_nbinsx, summary_nbinsy);
}

void TotemT2DQMSource::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // fill digis information
  for (const auto& ds_digis : iEvent.get(digiToken_)) {
    const TotemT2DetId detid(ds_digis.detId());
    for (const auto& digi : ds_digis) {
      segm_->fill(m_digis_mult_[detid.plane()]->getTH2D(), detid);
      (void)digi;  //FIXME make use of them
    }
  }
  // fill rechits information
  for (const auto& ds_rechits : iEvent.get(rechitToken_)) {
    const TotemT2DetId detid(ds_rechits.detId());
    for (const auto& rechit : ds_rechits) {
      segm_->fill(m_rechits_mult_[detid.plane()]->getTH2D(), detid);
      (void)rechit;  //FIXME make use of them
    }
  }
}

DEFINE_FWK_MODULE(TotemT2DQMSource);
