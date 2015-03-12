# Auto generated configuration file
# using: 
# Revision: 1.381.2.28 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/SevenTeV/Hadronizer_TuneZ2_7TeV_generic_LHE_pythia_cff.py -s GEN --filein pwgevents.lhe --datatier=GEN-SIM-RAW --conditions auto:mc --eventcontent RAWSIM --no_exec -n 5000 --python_filename=rivet_cfg.py --customise=Configuration/GenProduction/rivet_customize.py
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
import os

ivars = VarParsing.VarParsing('python')

ivars.maxEvents = 5000

#ivars.register('inputFiles',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_0.lhe',
#                VarParsing.VarParsing.multiplicity.singleton,
#                VarParsing.VarParsing.varType.string,
#               "Input LHE file")

ivars.register('output',
               'TuneZ2_POWHEG_CT10NLO_7TeV_supfactr250_ktmin5.aida',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
               "analysis output file")

ivars.register('sqrtS',
7000.0,
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.float,
"Center-of-mass energy")

ivars.register('sigma',
3.526e+10,
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.float,
"Cross section value from powheg in pico Barns")

ivars.parseArguments()

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(ivars.maxEvents)
)

#process.options = cms.untracked.PSet(
#wantSummary = cms.untracked.bool(True)
#)

# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring(
     ivars.inputFiles
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_0.lhe',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_100.lhe',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_101.lhe',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_102.lhe',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_103.lhe',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_104.lhe',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_105.lhe',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_106.lhe',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_107.lhe',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_108.lhe',
#               'file:/cms/store/user/ymao/pA5TEV/Generator/POWHEG/powheg_Job_109.lhe'
                )
    )

#mylist = FileUtils.loadListFromFile ('filelistAll.dat')
#for fname in mylist:
# process.source.fileNames.append('file:%s' % (ivars.inputFiles))
# process.source.fileNames.append('%s' % (fname))

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('runs Z2 Pythia6'),
    name = cms.untracked.string('$Source: /cvs/CMSSW/CMSSW/Configuration/GenProduction/python/SevenTeV/Z2onLHE_7TeV_pythia6_cff.py,v $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('Hadronizer_TuneZ2_7TeV_generic_LHE_pythia_cff_py_GEN.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.generator = cms.EDFilter("Pythia6HadronizerFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    comEnergy = cms.double(7000.0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'PARP(82)=1.832 ! pt cutoff for multiparton interactions', 
            'PARP(89)=1800. ! sqrts for which PARP82 is set', 
            'PARP(90)=0.275 ! Multiple interactions: rescaling power', 
            'MSTP(95)=6     ! CR (color reconnection parameters)', 
            'PARP(77)=1.016 ! CR', 
            'PARP(78)=0.538 ! CR', 
            'PARP(80)=0.1   ! Prob. colored parton from BBR', 
            'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
            'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
            'PARP(62)=1.025 ! ISR cutoff', 
            'MSTP(91)=1     ! Gaussian primordial kT', 
            'PARP(93)=10.0  ! primordial kT-max', 
            'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model'),
        processParameters = cms.vstring('MSEL=0          ! User defined processes', 
            'PMAS(5,1)=4.4   ! b quark mass', 
            'PMAS(6,1)=172.4 ! t quark mass'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from Configuration.GenProduction.rivet_customize
#from Configuration.GenProduction.rivet_customize import customise 

def customise(process):
        process.load('GeneratorInterface.RivetInterface.rivetAnalyzer_cfi')
        process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2015_HIN14001R5')
        process.rivetAnalyzer.CrossSection = cms.double(ivars.sigma)
#        process.rivetAnalyzer.CrossSection = cms.double(3.526e+10)
        process.rivetAnalyzer.OutputFile = cms.string(ivars.output)
        process.generation_step+=process.rivetAnalyzer
        process.schedule.remove(process.RAWSIMoutput_step)
        return(process)

#call to customisation function customise imported from Configuration.GenProduction.rivet_customize
process = customise(process)

# End of customisation functions
