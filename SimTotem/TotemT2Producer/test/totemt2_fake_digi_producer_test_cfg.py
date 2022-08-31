import FWCore.ParameterSet.Config as cms

# test to check if producer works correctly
process = cms.Process('TESTDIGIPRODUCER')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimTotem.TotemT2Producer.totemT2EmulatedDigis_cfi')

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(6000))
process.verbosity = cms.untracked.PSet(input = cms.untracked.int32(-1))

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

from SimTotem.TotemT2Producer.totemT2DetIds_cff import *

test_cases = cms.VPSet(
    cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_0), eventLimit = cms.uint32(100)),
    cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_1), eventLimit = cms.uint32(200)),
    cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_2), eventLimit = cms.uint32(500)),
    cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_3), eventLimit = cms.uint32(1000)),

    cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_0), eventLimit = cms.uint32(10)),
    cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_1), eventLimit = cms.uint32(50)),
    cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_2), eventLimit = cms.uint32(100)),
    cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_3), eventLimit = cms.uint32(500)),
)

# fake digi producer
process.totemT2Digis = process.totemT2EmulatedDigis.clone(
    t2DataFile = cms.string("/eos/cms/store/group/dpg_ctpps/comm_ctpps/TotemT2/RecoTest/nT2_Nino.root"),
    testCasesSet = test_cases,
    nthZeroEdge = cms.uint32(20)
)

# fake digi analyzer
process.totemT2DigiAnalyzer = cms.EDAnalyzer('TotemT2DigiAnalyzer',
    digisTag = cms.InputTag('totemT2Digis', 'TotemT2'),
    numberOfAllEvents = cms.uint32(process.maxEvents.input.value()),
    testCasesSet = test_cases
)

process.path = cms.Path(
    process.totemT2Digis *
    process.totemT2DigiAnalyzer
)

process.schedule = cms.Schedule(
    process.path
)
