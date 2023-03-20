import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = "DEBUG"
# enable LogDebug messages only for specific modules
process.MessageLogger.debugModules = ["*"]

process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# raw data source
process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring(
#        '/store/t0streamer/Data/Physics/000/286/591/run286591_ls0521_streamPhysics_StorageManager.dat',
#        'file:../../../../run364538_ls0001_streamA_StorageManager.dat',
        'file:../../../../t2v2p1/run364983_ls0001_streamA_StorageManager.dat',
        'file:../../../../t2v2p1/run364983_ls0002_streamA_StorageManager.dat',
#        '/store/t0streamer/Minidaq/A/000/303/982/run303982_ls0001_streamA_StorageManager.dat',
    )
)
#process.source = cms.Source('PoolSource',
#    fileNames = cms.untracked.vstring(
#        'file:../../../../T2Commisioning-run358786/run358786_ls0001_streamA_StorageManager.dat',
##        '/store/data/Run2018D/ZeroBias/RAW/v1/000/324/747/00000/97A72F4B-786F-5A48-B97E-C596DD73BD77.root',
#    ),
#)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

# raw-to-digi conversion
process.load('CalibPPS.ESProducers.totemT2DAQMapping_cff')
process.load('EventFilter.CTPPSRawToDigi.totemT2Digis_cfi')
process.totemT2Digis.rawDataTag = cms.InputTag("rawDataCollector")
process.totemT2Digis.RawUnpacking.verbosity = cms.untracked.uint32(1)
process.totemT2Digis.RawToDigi.verbosity = cms.untracked.uint32(3)
process.totemT2Digis.RawToDigi.printUnknownFrameSummary = cms.untracked.uint32(3)
process.totemT2Digis.RawToDigi.printErrorSummary = cms.untracked.uint32(3)

# rechits production
process.load('Geometry.ForwardCommonData.totemT22021V2XML_cfi')
process.load('Geometry.ForwardGeometry.totemGeometryESModule_cfi')
process.load('RecoPPS.Local.totemT2RecHits_cfi')

process.output = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("file:output-miniDaq2303-T2testFile1-Sun19Mar-sampic-printUnknownVFATs-T2mapFrom0-5ev.root"),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_totemT2*_*_*',
    ),
)

# execution configuration
process.p = cms.Path(
    process.totemT2Digis
    * process.totemT2RecHits
)

process.outpath = cms.EndPath(process.output)
