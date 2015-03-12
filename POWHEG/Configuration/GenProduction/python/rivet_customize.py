import FWCore.ParameterSet.Config as cms

def customise(process):
        process.load('GeneratorInterface.RivetInterface.rivetAnalyzer_cfi')
        process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_TEST_ANALYSIS')
        process.rivetAnalyzer.CrossSection = cms.double(1.347e+10) 
        process.rivetAnalyzer.OutputFile = cms.string('TuneZ2_POWHEG_CT10NLO_7TeV_supfactr250_ktmin5.aida')
        process.generation_step+=process.rivetAnalyzer
        process.schedule.remove(process.RAWSIMoutput_step)
        return(process)

