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
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/668d5798-c82f-4d90-bdc1-b3e6cecbbaab.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/29d45a0e-a675-447c-92f2-e6a52f66f721.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/39f599c0-c2dd-4868-988e-a2aadd0b01b7.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/3657e923-972e-432a-9c6b-f015b4f858eb.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/46f83c12-04a2-4e08-b0dc-ed2bcdc6036e.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/e4f2e705-bfd9-4a0f-9088-1787efa507fd.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/70e1e556-557f-40ae-b401-eff2e3c4a686.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/4c54a650-daf8-404b-beab-b8a56ce9f7c7.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/e6b914fb-99c3-4df8-a077-9997fb2f42f4.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/aaea93e0-8095-48b0-8e13-94b92cbfd0d4.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/6af5a06a-9960-4d6b-ba75-2649539b4af1.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/4f047e8e-732a-48bd-a57a-ec7cee0eb955.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/9c963826-890f-46f6-b653-74bb07a4b9e4.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/5187ee91-af23-47df-83dd-83d213294ea1.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/77519d33-5f8f-4f63-beb2-62a223559ef5.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/1baee240-8da6-4f30-9de4-ad19074e2180.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/9591618f-7208-4fd5-a0b3-8f5cfa2f6b52.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/9ec595e7-badd-497e-9832-96a5839b410d.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/17ff0718-87c8-4b8c-a8d8-23f2d79ddec7.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/18a719f5-eb42-4d19-87fc-365c9131effe.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/64c3807b-8aa8-46bb-8bc8-ed7fc265efe0.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/8c57abd7-a158-4e27-beb4-a8908605a94b.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/277d3fc0-9d56-4699-8601-4195c517483c.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/494d21aa-a833-4a99-a792-61d1056872cb.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/3d726cbe-931d-485d-acab-4f69c2022a7a.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/0555ff67-7ec5-4c52-b2cc-40c76d78ef1f.root',
        '/store/data/Run2023C/ZeroBias1/AOD/PromptReco-v4/000/369/596/00000/72dee819-b191-4eae-b974-e13926767410.root'
    ),
    lumisToProcess = cms.untracked.VLuminosityBlockRange("369596:0-369596:40"),
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
process.ctppsRPAlignmentCorrectionsDataESSourceXML.RealFiles = cms.vstring('data/369596.0-40/103,104,105,123,124,125-excl2007236608,2007269376/ev=4E5,pl=3,units=2,ovlp=f,3rpsInO=t/s+sr,std,1rotzIt=0,eMuMvRot=t/iteration2/results_cumulative_factored_Jan.xml')
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

process.ppsStraightTrackAligner.rpIds = [103,104,105,123,124,125]
process.ppsStraightTrackAligner.excludePlanes = cms.vuint32(2007236608,2007269376)
process.ppsStraightTrackAligner.z0 = +217000
process.ppsStraightTrackAligner.horizontalOffsets = cms.vstring('103:1.295', '123:2.179')

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
process.ppsStraightTrackAligner.maxTrackAy = 0.5E-3

optimize="sr"
process.ppsStraightTrackAligner.resolveShR = True
process.ppsStraightTrackAligner.resolveShZ = False
process.ppsStraightTrackAligner.resolveRotZ = True

process.ppsStraightTrackAligner.constraintsType = "standard"
process.ppsStraightTrackAligner.standardConstraints.units = cms.vuint32(101,121)
process.ppsStraightTrackAligner.oneRotZPerPot = False
process.ppsStraightTrackAligner.useEqualMeanUMeanVRotZConstraints = True

process.ppsStraightTrackAligner.algorithms = cms.vstring("Jan")

process.ppsStraightTrackAligner.JanAlignmentAlgorithm.stopOnSingularModes = False

results_dir="/afs/cern.ch/user/f/foljemar/tmptst/J3/CMSSW_13_0_7_TOTEM/src/CalibPPS/AlignmentRelative/test/test_with_data/data/369596.0-40/103,104,105,123,124,125-excl2007236608,2007269376/ev=4E5,pl=3,units=2,ovlp=f,3rpsInO=t/s+sr,std,1rotzIt=0,eMuMvRot=t/iteration3"

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
