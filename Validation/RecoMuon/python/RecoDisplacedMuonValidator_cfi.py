import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy
from Validation.RecoMuon.selectors_cff import muonTPSet, displacedMuonTPSet

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
recoDisplacedMuonValidator = DQMEDAnalyzer('RecoDisplacedMuonValidator',
    MuonServiceProxy,
    tpSelector = displacedMuonTPSet,

    usePFMuon = cms.untracked.bool(False),

    simLabel = cms.InputTag("mix","MergedTrackTruth"),
    tpRefVector = cms.bool(False),
    muonLabel = cms.InputTag("displacedMuons"),

    muAssocLabel = cms.InputTag("muonAssociatorByHitsHelper"),

    doAssoc = cms.untracked.bool(True),

    outputFileName = cms.untracked.string(''),
    subDir = cms.untracked.string('Muons/RecoDisplacedMuonV/'),
    trackType = cms.string("global"),
    #string cut selection
    selection = cms.string(""),

    wantTightMuon = cms.bool(False),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    primaryVertex = cms.InputTag('offlinePrimaryVertices'),

    #
    # Histogram dimensions     #
    #
    nBinP = cms.untracked.uint32(100),
    minP = cms.untracked.double(0.0),
    maxP = cms.untracked.double(1500.0),

    nBinPt = cms.untracked.uint32(100),
    minPt = cms.untracked.double(0.0),
    maxPt = cms.untracked.double(1500.0),

    doAbsEta = cms.untracked.bool(False),

    nBinEta = cms.untracked.uint32(50),
    minEta = cms.untracked.double(-2.5),
    maxEta = cms.untracked.double(2.5),

    nBinDxy = cms.untracked.uint32(100),
    minDxy = cms.untracked.double(-60), # 350 for DSA, 60 for DGL/DTK
    maxDxy = cms.untracked.double(60),

    nBinDz = cms.untracked.uint32(100),
    minDz = cms.untracked.double(-60),
    maxDz = cms.untracked.double(60),

    nBinPhi = cms.untracked.uint32(25),

    # Pull width     #
    nBinPull = cms.untracked.uint32(50),
    wPull = cms.untracked.double(5.0),

    nBinErr = cms.untracked.uint32(50),

    # |p| resolution     #
    minErrP = cms.untracked.double(-0.3),
    maxErrP = cms.untracked.double(0.3),

    # pT resolution     #
    minErrPt = cms.untracked.double(-0.3),
    maxErrPt = cms.untracked.double(0.3),

    # q/pT resolution     #
    minErrQPt = cms.untracked.double(-0.1),
    maxErrQPt = cms.untracked.double(0.1),

    # Eta resolution     #
    minErrEta = cms.untracked.double(-0.01),
    maxErrEta = cms.untracked.double(0.01),

    # Phi resolution     #
    minErrPhi = cms.untracked.double(-0.05),
    maxErrPhi = cms.untracked.double(0.05),

    # Dxy resolution     #
    minErrDxy = cms.untracked.double(-0.1),
    maxErrDxy = cms.untracked.double(0.1),

    # Dz resolution     #
    minErrDz = cms.untracked.double(-0.1),
    maxErrDz = cms.untracked.double(0.1),

    # Number of sim-reco associations     #
    nAssoc = cms.untracked.uint32(10),

    # Number of sim,reco Tracks     #
    nTrks = cms.untracked.uint32(50)
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(recoDisplacedMuonValidator, simLabel = "mixData:MergedTrackTruth")
