import FWCore.ParameterSet.Config as cms

process = cms.Process("ggll")

#runOnMC = False
runOnMC = True
useAOD = True # AOD or MiniAOD?

leptonsType = 'Muon'
#leptonsType = 'Electron'
fetchProtons = False
if not runOnMC: fetchProtons = True

triggerList = [
    'HLT_DoubleMu43NoFiltersNoVtx_*'
#    'HLT_DoubleEle33_CaloIdL_MW_v*',
#    'HLT_Ele27_HighEta_Ele20_Mass55_v*',
#    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*'
]

#########################
#    General options    #
#########################

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound'),
    allowUnscheduled = cms.untracked.bool(True),
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Data
        #'/store/data/Run2017C/DoubleMuon/AOD/17Nov2017-v1/30001/EAA10504-C1D8-E711-9062-001E67396897.root'
        # MC
	'/store/mc/RunIIFall17DRPremix/GGToMuMu_Pt-50_Elastic_13TeV-lpair/AODSIM/94X_mc2017_realistic_v11-v1/40000/527A7A11-2A1B-E811-B9DD-A45D36FC89C4.root'
        #'/store/mc/RunIIFall17DRPremix/GGToEE_Pt-50_Elastic_13TeV-lpair/AODSIM/94X_mc2017_realistic_v11-v1/00000/CE748075-4E1B-E811-BBE4-3C4A92F8FC10.root'
        #'file:/tmp/antoniov/GGToMuMu_E0CAC014-2A1B-E811-BD0C-1458D04923EC.root'
    ),
    #firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltFilter.HLTPaths = cms.vstring(
#    'HLT_DoubleMu43NoFiltersNoVtx_*',
##    'HLT_DoubleEle33_CaloIdL_MW_v*',
##    'HLT_Ele27_HighEta_Ele20_Mass55_v*',
##    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*',
#)
process.hltFilter.HLTPaths = triggerList

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if not runOnMC: process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
else: process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


#########################
#     PAT-ification     #
#########################
## Look at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Core_Tools for more information

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('PATuple.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_offline*PrimaryVertices*_*_*',
        'keep *_selectedPatMuons*_*_*',
        'keep *_*lectron*_*_*',
        'keep *_selectedPatElectrons*_*_*',
        'keep *_selectedPat*Photons*_*_*',
        'keep *_selectedPatJets*_*_*',
        'keep *_*MET*_*_*',
        'keep *_*particleFlow*_*_*',
    ),
)
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
patAlgosToolsTask.add(process.patCandidatesTask)

process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
patAlgosToolsTask.add(process.selectedPatCandidatesTask)

from PhysicsTools.PatAlgos.tools.coreTools import runOnData
if not runOnMC:
    runOnData( process )

#########################
#      Electron ID      #
#########################

#from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer, setupVIDElectronSelection, setupAllVIDIdsInModule, DataFormat
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

switchOnVIDElectronIdProducer(process, DataFormat.AOD)
#setupAllVIDIdsInModule(process, 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff', setupVIDElectronSelection)
#setupAllVIDIdsInModule(process, 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff', setupVIDElectronSelection)
setupAllVIDIdsInModule(process, 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff', setupVIDElectronSelection)

#########################
#       Photon ID       #
#########################

switchOnVIDPhotonIdProducer(process, DataFormat.AOD)
#setupAllVIDIdsInModule(process, 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff', setupVIDPhotonSelection)
setupAllVIDIdsInModule(process, 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff', setupVIDPhotonSelection)

#########################
#       Analysis        #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.GammaGammaLL_cfi")

process.ggll_aod.triggersList = process.hltFilter.HLTPaths
#process.ggll_aod.leptonsType = cms.string('Muon')
#process.ggll_aod.leptonsType = cms.string('ElectronMuon')
#process.ggll_aod.leptonsType = cms.string('Electron')
process.ggll_aod.leptonsType = cms.string( leptonsType )
process.ggll_aod.muonTag = cms.InputTag("selectedPatMuons")
process.ggll_aod.electronTag = cms.InputTag("selectedPatElectrons")
process.ggll_aod.jetTag = cms.InputTag('selectedPatJets')
process.ggll_aod.runOnMC = cms.bool(runOnMC)
process.ggll_aod.fetchProtons = cms.bool( fetchProtons )
process.ggll_aod.mcpufile = 'pileup_MC.root'
process.ggll_aod.mcpupath = 'pileup'
process.ggll_aod.datapufile = 'pileup_Data.root'
process.ggll_aod.datapupath = 'pileup'

# E/gamma identification
process.ggll_aod.eleIdLabels = cms.PSet(
    #mediumLabel = cms.InputTag('mvaEleID-Spring16-GeneralPurpose-V1-wp90'),
    #tightLabel = cms.InputTag('mvaEleID-Spring16-GeneralPurpose-V1-wp80'),
    mediumLabel = cms.InputTag('mvaEleID-Fall17-iso-V1-wp90'),
    tightLabel = cms.InputTag('mvaEleID-Fall17-iso-V1-wp80'),
)
process.ggll_aod.phoIdLabels = cms.PSet(
    #mediumLabel = cms.InputTag('mvaPhoID-Spring16-nonTrig-V1-wp90'),
    #tightLabel = cms.InputTag('mvaPhoID-Spring16-nonTrig-V1-wp80'),
    mediumLabel = cms.InputTag('mvaPhoID-RunIIFall17-v1-wp80'),
    tightLabel = cms.InputTag('mvaPhoID-RunIIFall17-v1-wp90'),
)
#process.ggll_aod.eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90")
#process.ggll_aod.eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80")
#process.ggll_aod.phoMediumIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90")
#process.ggll_aod.phoTightIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp80")

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('output.root'),
    closeFileFast = cms.untracked.bool(True)
)

if not runOnMC:
    process.p = cms.Path(
	process.hltFilter*
	process.egmPhotonIDSequence*
	process.egmGsfElectronIDSequence*
	process.ggll_aod
    )
else:
    process.p = cms.Path(
	process.egmPhotonIDSequence*
	process.egmGsfElectronIDSequence*
	process.ggll_aod
    )

#process.outpath = cms.EndPath(process.out, patAlgosToolsTask)
process.outpath = cms.EndPath(patAlgosToolsTask)

