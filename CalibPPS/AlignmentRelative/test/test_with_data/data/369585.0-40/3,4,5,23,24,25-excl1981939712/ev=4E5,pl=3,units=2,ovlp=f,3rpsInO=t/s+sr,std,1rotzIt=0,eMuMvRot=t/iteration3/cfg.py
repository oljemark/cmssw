import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("trackBasedAlignment", eras.Run3)

# =================== GlobalTag ===================
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '130X_dataRun3_Prompt_v3'
# =================== GlobalTag ===================

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
  statistics = cms.untracked.vstring(),
  destinations = cms.untracked.vstring('cout'),
  cout = cms.untracked.PSet(
    threshold = cms.untracked.string('WARNING')
  )
)

# input data
process.source = cms.Source("PoolSource",
    skipBadFiles = cms.untracked.bool(True),
    fileNames = cms.untracked.vstring(
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/585/00000/1b9d5007-f599-419d-bb85-15ee78ee54fc.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/585/00000/3d42ecab-f3d1-44e8-ba97-a80ad031fca1.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/585/00000/db1017a0-672a-451c-90a8-16d78ad0b205.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/585/00000/863d0329-2be8-4ebf-aa17-24956cd34d2c.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/585/00000/d32ea43b-2198-4be6-8f67-32516fd2278c.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/585/00000/83b45783-d37f-4b7e-9927-bc7926736048.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/585/00000/68919351-ea7f-4a1e-8e51-94ccb7a7afc4.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/585/00000/929235c7-f6b0-4e13-b35a-779c1703617a.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/585/00000/1a4263a4-28d9-4172-8ff0-14a13d279d59.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/585/00000/574a3c79-5364-42f1-a163-bff297ff2a42.root'
    ),
    lumisToProcess = cms.untracked.VLuminosityBlockRange("369585:0-369585:40"),
    inputCommands = cms.untracked.vstring(
        "keep *",
      # "drop *",
      # "keep TotemRPRecHitedmDetSetVector_*_*_*",
      # "keep CTPPSPixelRecHitedmDetSetVector_*_*_*",
    )
)

# geometry
process.load("Geometry.VeryForwardGeometry.geometryRPFromDD_2022_cfi")
process.ctppsGeometryESModule.buildMisalignedGeometry = cms.bool(True)
del(process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles[-1])
process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles.append("CalibPPS/AlignmentRelative/test/test_with_data/RP_Dist_Beam_Cent.xml")

# initial alignments
process.load("CalibPPS.ESProducers.ctppsRPAlignmentCorrectionsDataESSourceXML_cfi")
process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = cms.vstring('data/369585.0-40/3,4,5,23,24,25-excl1981939712/ev=4E5,pl=3,units=2,ovlp=f,3rpsInO=t/s+sr,std,1rotzIt=0,eMuMvRot=t/iteration2/results_cumulative_factored_Jan.xml')
process.ctppsRPAlignmentCorrectionsDataESSourceXML.verbosity = 1

process.alignPref=cms.ESPrefer("CTPPSRPAlignmentCorrectionsDataESSourceXML","ctppsRPAlignmentCorrectionsDataESSourceXML",RPRealAlignmentRecord=cms.vstring("CTPPSRPAlignmentCorrectionsData"))



# reco modules
process.load("RecoPPS.Local.totemRPLocalReconstruction_cff")

process.load("RecoPPS.Local.ctppsPixelLocalReconstruction_cff")

process.load("RecoPPS.Local.ctppsLocalTrackLiteProducer_cff")
process.ctppsLocalTrackLiteProducer.includeDiamonds = False

# aligner
process.load("CalibPPS.AlignmentRelative.ppsStraightTrackAligner_cfi")

process.ppsStraightTrackAligner.verbosity = 5

process.ppsStraightTrackAligner.tagUVPatternsStrip = cms.InputTag("totemRPUVPatternFinder")
process.ppsStraightTrackAligner.tagDiamondHits = cms.InputTag("")
process.ppsStraightTrackAligner.tagPixelHits = cms.InputTag("")
process.ppsStraightTrackAligner.tagPixelLocalTracks = cms.InputTag("ctppsPixelLocalTracks")

process.ppsStraightTrackAligner.maxEvents = int(4E5)

process.ppsStraightTrackAligner.rpIds = [3,4,5,23,24,25]
process.ppsStraightTrackAligner.excludePlanes = cms.vuint32(1981939712)
process.ppsStraightTrackAligner.z0 = -217000
process.ppsStraightTrackAligner.horizontalOffsets = cms.vstring('3:-2.927', '23:-3.171')

process.ppsStraightTrackAligner.maxResidualToSigma = 10
process.ppsStraightTrackAligner.minimumHitsPerProjectionPerRP = 3

process.ppsStraightTrackAligner.removeImpossible = True
process.ppsStraightTrackAligner.requireNumberOfUnits = 2
process.ppsStraightTrackAligner.requireOverlap = False
process.ppsStraightTrackAligner.requireAtLeast3PotsInOverlap = True
process.ppsStraightTrackAligner.additionalAcceptedRPSets = ""

process.ppsStraightTrackAligner.cutOnChiSqPerNdf = True
process.ppsStraightTrackAligner.chiSqPerNdfCut = 50

process.ppsStraightTrackAligner.maxTrackAx = 0.5E-3
process.ppsStraightTrackAligner.maxTrackAy = 2.5E-3

optimize="sr"
process.ppsStraightTrackAligner.resolveShR = True
process.ppsStraightTrackAligner.resolveShZ = False
process.ppsStraightTrackAligner.resolveRotZ = True

process.ppsStraightTrackAligner.constraintsType = "standard"
process.ppsStraightTrackAligner.standardConstraints.units = cms.vuint32(1,21)
process.ppsStraightTrackAligner.oneRotZPerPot = False
process.ppsStraightTrackAligner.useEqualMeanUMeanVRotZConstraints = True

process.ppsStraightTrackAligner.algorithms = cms.vstring("Jan")

process.ppsStraightTrackAligner.JanAlignmentAlgorithm.stopOnSingularModes = False

results_dir="/afs/cern.ch/user/f/foljemar/tmptst/j2/CMSSW_13_0_7_TOTEM/src/CalibPPS/AlignmentRelative/test/test_with_data/data/369585.0-40/3,4,5,23,24,25-excl1981939712/ev=4E5,pl=3,units=2,ovlp=f,3rpsInO=t/s+sr,std,1rotzIt=0,eMuMvRot=t/iteration3"

process.ppsStraightTrackAligner.taskDataFileName = "" # results_dir + "/task_data.root"

process.ppsStraightTrackAligner.fileNamePrefix = results_dir + "/results_iteration_"
process.ppsStraightTrackAligner.expandedFileNamePrefix = results_dir + "/results_cumulative_expanded_"
process.ppsStraightTrackAligner.factoredFileNamePrefix = results_dir + "/results_cumulative_factored_"

process.ppsStraightTrackAligner.diagnosticsFile = results_dir + '/diagnostics.root'
process.ppsStraightTrackAligner.buildDiagnosticPlots = True
process.ppsStraightTrackAligner.JanAlignmentAlgorithm.buildDiagnosticPlots = True

# processing sequence
process.p = cms.Path(
  # it is important to re-run part of the reconstruction as it may influence
  # the choice of rec-hits used in the alignment
  process.totemRPUVPatternFinder
  * process.totemRPLocalTrackFitter
  * process.ctppsPixelLocalTracks
  * process.ctppsLocalTrackLiteProducer
  * process.ppsStraightTrackAligner
)
