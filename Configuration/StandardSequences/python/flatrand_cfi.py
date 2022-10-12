import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("FlatRandomEThetaGunProducer",
        PGunParameters = cms.PSet(
            PartID = cms.vint32(2212),
            MinTheta = cms.double(0),
            MaxTheta = cms.double(0.1),
            MinPhi = cms.double(-3.1415),
            MaxPhi = cms.double(3.1415),
#           AddAntiParticle = cms.bool(False),
#            ECMS   = cms.double(ecms),
            MinEta = cms.double(4.4),
            MaxEta = cms.double(6.8),
#            Mint   = cms.double(t_min),
#            Maxt   = cms.double(t_max),
            MinE  = cms.double(2000.),
            MaxE  = cms.double(6500.)
            ),
        Verbosity = cms.untracked.int32(0),
        psethack = cms.string('single protons'),
        AddAntiParticle = cms.bool(False),
#        FireBackward = cms.bool(True),
#        FireForward  = cms.bool(True),
        firstRun = cms.untracked.uint32(1),
        )


