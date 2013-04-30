version="538"
tracks="generalTracks"
sample="PythiaZ2_2760GeV"

import FWCore.ParameterSet.Config as cms
process = cms.Process('jecdb')
process.load('CondCore.DBCommon.CondDBCommon_cfi')
process.CondDBCommon.connect = 'sqlite_file:JEC_PP2760GEV_CMSSW538_2013.db'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source('EmptySource')
process.PoolDBOutputService = cms.Service('PoolDBOutputService',
                                          process.CondDBCommon,
                                          toPut = cms.VPSet(


    cms.PSet(
    record = cms.string('AK3PF'),
    tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_'+tracks+'_'+sample+'_'+version+'_AK3PF'),
    label  = cms.string('AK3PF')
    ),
    cms.PSet(
    record = cms.string('AK4PF'),
    tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_'+tracks+'_'+sample+'_'+version+'_AK4PF'),
    label  = cms.string('AK4PF')
    ),
    cms.PSet(
    record = cms.string('AK5PF'),
    tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_'+tracks+'_'+sample+'_'+version+'_AK5PF'),
    label  = cms.string('AK5PF')
    ),
                                          cms.PSet(
    record = cms.string('AK3Calo'),
    tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_'+sample+'_'+version+'_AK3Calo'),
    label  = cms.string('AK3Calo')
    ),
                                          cms.PSet(
    record = cms.string('AK4Calo'),
    tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_'+sample+'_'+version+'_AK4Calo'),
    label  = cms.string('AK4Calo')
    ),
                                          cms.PSet(
    record = cms.string('AK5Calo'),
    tag    = cms.string('JetCorrectorParametersCollection_HI_PFTowers_'+sample+'_'+version+'_AK5Calo'),
    label  = cms.string('AK5Calo')
    ),


                                          )
)
process.dbWriterAK3PFTowers = cms.EDAnalyzer('JetCorrectorDBWriter',
                                       era    = cms.untracked.string('JEC_dijet'),
                                       algo   = cms.untracked.string('AK3PF')
                                       )
process.dbWriterAK4PFTowers = cms.EDAnalyzer('JetCorrectorDBWriter',
                                       era    = cms.untracked.string('JEC_dijet'),
                                       algo   = cms.untracked.string('AK4PF')
                                       )
process.dbWriterAK5PFTowers = cms.EDAnalyzer('JetCorrectorDBWriter',
                                       era    = cms.untracked.string('JEC_dijet'),
                                       algo   = cms.untracked.string('AK5PF')
                                       )

process.dbWriterAK3Calo = cms.EDAnalyzer('JetCorrectorDBWriter',
                                       era    = cms.untracked.string('JEC_dijet'),
                                       algo   = cms.untracked.string('AK3Calo')
                                       )
process.dbWriterAK4Calo = cms.EDAnalyzer('JetCorrectorDBWriter',
                                       era    = cms.untracked.string('JEC_dijet'),
                                       algo   = cms.untracked.string('AK4Calo')
                                       )
process.dbWriterAK5Calo = cms.EDAnalyzer('JetCorrectorDBWriter',
                                       era    = cms.untracked.string('JEC_dijet'),
                                       algo   = cms.untracked.string('AK5Calo')
                                       )

process.p = cms.Path(

process.dbWriterAK3PFTowers +
process.dbWriterAK4PFTowers +
process.dbWriterAK5PFTowers +
process.dbWriterAK3Calo +
process.dbWriterAK4Calo +
process.dbWriterAK5Calo 

)
'''
process.p = cms.Path(

process.dbWriterAK1PFTowers +
process.dbWriterAK2PFTowers +
process.dbWriterAK3PFTowers +
process.dbWriterAK4PFTowers +
process.dbWriterAK5PFTowers +
process.dbWriterAK6PFTowers +
process.dbWriterIC5Calo +
process.dbWriterAK1Calo +
process.dbWriterAK2Calo +
process.dbWriterAK3Calo +
process.dbWriterAK4Calo +
process.dbWriterAK5Calo +
process.dbWriterAK6Calo 

) 
'''
