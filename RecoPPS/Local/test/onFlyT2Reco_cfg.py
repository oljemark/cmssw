import FWCore.ParameterSet.Config as cms

process = cms.Process('T2RECHIT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = "DEBUG"
# enable LogDebug messages only for specific modules
process.MessageLogger.debugModules = ["Totem"]

process.load('Configuration.EventContent.EventContent_cff')
from RecoPPS.Configuration.RecoCTPPS_EventContent_cff import RecoCTPPSAOD
from HLTrigger.Configuration.HLTrigger_EventContent_cff import HLTriggerAOD
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

#dummy = cms.untracked.FileInPath('RecoPPS/Local/data/run364983_ls0001_streamA_StorageManager.dat'),

# raw data source
process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring(
#        '/store/t0streamer/Data/PhysicsZeroBias2/000/369/585/run369585_ls0042_streamPhysicsZeroBias2_StorageManager.dat',
#        '/store/data/Run2023C/ZeroBiasNonColliding/AOD/PromptReco-v4/000/369/585/00000/82e14c5e-2b94-49c9-9fad-1259bca1d2ae.root',
#        '/store/data/Run2023C/ZeroBiasNonColliding/AOD/PromptReco-v4/000/369/585/00000/d7171c2e-f906-4a7f-89bd-12ef54efbf0d.root',
        '/store/data/Run2023C/ZeroBias11/AOD/PromptReco-v4/000/369/585/00000/dbd5bb7f-2ab6-4613-aedc-c671ba4f467d.root',
#        '/store/data/Run2023C/SpecialRandom7/AOD/PromptReco-v4/000/369/592/00000/c73f0118-153b-44a2-93a2-7f45e566d629.root',
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)


#Latest PPS GT candidate
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_Prompt_Candidate_2023_06_06_21_34_08', '')



# raw-to-digi conversion
#process.load('CalibPPS.ESProducers.totemT2DAQMapping_cff')
#process.load('EventFilter.CTPPSRawToDigi.totemT2Digis_cfi')
#process.totemT2Digis.rawDataTag = cms.InputTag("rawDataCollector")

# rechits production
process.load('Geometry.ForwardCommonData.totemT22021V2XML_cfi')
process.load('Geometry.ForwardGeometry.totemGeometryESModule_cfi')
process.load('RecoPPS.Local.totemT2RecHits_cfi')

process.output = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("file:reducedTOTEM-ZeroBias11-dbd5bb7f-2ab6-4613-aedc-c671ba4f467d-100000ev-plusHLTAOD.root"),
        outputCommands = RecoCTPPSAOD.outputCommands + HLTriggerAOD.outputCommands,
)

# execution configuration
process.p = cms.Path(
  process.totemT2RecHits
)

process.outpath = cms.EndPath(process.output)
