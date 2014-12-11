# Auto generated configuration file
# using: 
# Revision: 1.381.2.28 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step4 --conditions INSERT_GT_HERE -s RAW2DIGI,L1Reco,RECO --scenario HeavyIons --datatier GEN-SIM-RECO --himix --eventcontent RECODEBUG -n 100 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO3')

# Jet Reconstruction

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.HiEventMixing_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring(
                            'file:/afs/cern.ch/user/d/dgulhan/workDir/hiReco_RAW2DIGI_L1Reco_RECO_1000_1_S4K.root',
                     # "root://se1.accre.vanderbilt.edu:1095//store/hidata/HIRun2011/HIHighPt/RECO/14Mar2014-v2/00000/003B85C7-D8B9-E311-81DF-FA163EDC6F8F.root"
                                                              )
                            )


# Output definition
process.RECODEBUGoutput = cms.OutputModule("PoolOutputModule",
                                           splitLevel = cms.untracked.int32(0),
                                           eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
                                           outputCommands = process.RECODEBUGEventContent.outputCommands,
                                           fileName = cms.untracked.string('DATA_MinBias_RECO.root'),
                                           dataset = cms.untracked.PSet(
    filterName = cms.untracked.string('MinBiasCollEvtSel'),
    dataTier = cms.untracked.string('GEN-SIM-RECO')
    ),
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_R_53_LV6::All', '')
# Path and EndPath definitions

process.load('RecoHI.HiJetAlgos.HiRecoJets_cff')
process.load('RecoHI.HiJetAlgos.HiRecoPFJets_cff')

## background for HF/Voronoi-style subtraction
process.voronoiBackgroundCaloEqualizeR0p5= process.voronoiBackgroundCalo.clone(equalizeR=cms.double(0.5))
process.voronoiBackgroundCaloEqualizeR0p6= process.voronoiBackgroundCalo.clone(equalizeR=cms.double(0.6))
                            
process.akVs5CaloJetsEqualizeR0p6 = process.akVs5CaloJets.clone(bkg = cms.InputTag("voronoiBackgroundCaloEqualizeR0p6"))
process.akVs4CaloJetsEqualizeR0p5 = process.akVs4CaloJets.clone(bkg = cms.InputTag("voronoiBackgroundCaloEqualizeR0p5"))

process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECODEBUGoutput_step = cms.EndPath(process.RECODEBUGoutput)

process.jetsequence = cms.Sequence(process.voronoiBackgroundCaloEqualizeR0p5+process.voronoiBackgroundCaloEqualizeR0p6+process.akVs5CaloJetsEqualizeR0p6+process.akVs4CaloJetsEqualizeR0p5)

process.jetrecostep = cms.Path(process.jetsequence)

# Schedule definition
process.schedule = cms.Schedule(process.jetrecostep,process.endjob_step,process.RECODEBUGoutput_step)

process.SimpleMemoryCheck=cms.Service("SimpleMemoryCheck",
                                      oncePerEventMode=cms.untracked.bool(False))

process.Timing=cms.Service("Timing",
                           useJobReport = cms.untracked.bool(True)
                           )