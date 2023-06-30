import FWCore.ParameterSet.Config as cms
import string

from Configuration.Eras.Era_Run3_cff import Run3

# needed to have ctppsDQMCalibrationSource properly working
from Configuration.Eras.Modifier_ctpps_cff import ctpps

process = cms.Process('RECODQM', Run3)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(6000) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM framework
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "CTPPS"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "CTPPS"
# raw data source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_13_0_6/ZeroBias/RECO/PromptNew_PPS_RelVal-v1/2590000/15b616d6-a572-4648-a794-449b86263d0e.root',
    )
)

#Latest PPS GT candidate
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_Prompt_Candidate_2023_06_06_21_34_08', '')

#Raw-to-digi
process.load('EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff')

# local RP reconstruction chain with standard settings
process.load("RecoPPS.Configuration.recoCTPPS_cff")

# CTPPS DQM modules
process.load("DQM.CTPPS.ctppsDQM_cff")
process.totemDAQMappingESSourceXML_TotemT2.verbosity = 0
process.totemT2Digis.RawUnpacking.verbosity = 0
process.totemT2Digis.RawToDigi.verbosity = 0
process.totemT2Digis.RawToDigi.testCRC = 2
process.totemT2Digis.RawToDigi.useOlderT2TestFile = False
process.totemT2Digis.RawToDigi.printUnknownFrameSummary = True
process.totemT2Digis.RawToDigi.printErrorSummary = True

process.path = cms.Path(
#  process.ctppsRawToDigi *
#  process.recoCTPPS *
  process.ctppsDQMCalibrationSource *
  process.ctppsDQMCalibrationHarvest
)

process.end_path = cms.EndPath(
    process.dqmEnv +
    process.dqmSaver
)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
