import FWCore.ParameterSet.Config as cms

process = cms.Process('RECODQM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Geometry.ForwardCommonData.totemT22021V2XML_cfi')
process.load('Geometry.ForwardGeometry.totemGeometryESModule_cfi')
process.load('SimTotem.TotemT2Producer.totemT2EmulatedDigis_cfi')
process.load('RecoPPS.Local.totemT2RecHits_cfi')
process.load('DQM.CTPPS.totemT2DQMSource_cfi')

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

# load DQM framework
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "CTPPS"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "CTPPS"

from SimTotem.TotemT2Producer.totemT2DetIds_cff import *

ids_arm_0_plane_0_tile_0 = ids_arm_0_plane_all_tile_0[0::2][:1]
ids_arm_0_plane_02_tile_1 = ids_arm_0_plane_all_tile_1[0::2][:2]
ids_arm_0_plane_024_tile_2 = ids_arm_0_plane_all_tile_2[0::2][:3]
ids_arm_0_plane_0246_tile_3 = ids_arm_0_plane_all_tile_3[0::2]
print(ids_arm_0_plane_0_tile_0, ids_arm_0_plane_02_tile_1, ids_arm_0_plane_024_tile_2)

# [2113929216] [2113945600, 2114076672] [2113961984, 2114093056, 2114224128]
test_cases = cms.VPSet(
    #cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_0), eventLimit = cms.uint32(100)),
    cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_0_tile_0), eventLimit = cms.uint32(100)),
    #cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_1), eventLimit = cms.uint32(200)),
    cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_02_tile_1), eventLimit = cms.uint32(200)),
    #cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_2), eventLimit = cms.uint32(500)),
    cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_024_tile_2), eventLimit = cms.uint32(500)),
    cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_0246_tile_3), eventLimit = cms.uint32(1000)),

    cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_0), eventLimit = cms.uint32(10)),
    cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_1), eventLimit = cms.uint32(50)),
    cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_2), eventLimit = cms.uint32(100)),
    cms.PSet(detId = cms.vuint32(*ids_arm_1_plane_all_tile_3), eventLimit = cms.uint32(500)),
)

test_cases_1 = cms.VPSet(
    cms.PSet(detId = cms.vuint32(*ids_arm_0_plane_all_tile_0)),
)

# fake digi producer
process.totemT2Digis = process.totemT2EmulatedDigis.clone(
    t2DataFile = cms.string('/eos/cms/store/group/dpg_ctpps/comm_ctpps/TotemT2/RecoTest/nT2_Nino.root'),
    # set of test cases
    testCasesSet = test_cases
)

process.path = cms.Path(
    process.totemT2Digis *
    process.totemT2RecHits *
    process.totemT2DQMSource
)

process.end_path = cms.EndPath(
    process.dqmEnv +
    process.dqmSaver
)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
