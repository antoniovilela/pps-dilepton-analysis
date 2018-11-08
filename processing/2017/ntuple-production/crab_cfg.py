from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

#config.General.requestName = 'GGToMuMu_Pt-50_Elastic_13TeV-GammaGammaLL-test-v4'
config.General.requestName = 'DoubleMuon_Run2017C-17Nov2017-test-v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RunGammaGammaLeptonLepton_cfg.py'
config.JobType.outputFiles = ['output.root']
#config.JobType.inputFiles = ['pileup_MC.root', 'pileup_Data.root']

#config.Data.outputPrimaryDataset = 'FPMC_gamgamZZ_anom_A0Z_1E-5_Lambda_2TeV_13TeV'
#inputFiles = open('files_FPMC_gamgamZZ_anom_A0Z_1E-5_Lambda_2TeV_13TeV-v1.txt').readlines()
#print "Adding files", inputFiles
#config.Data.userInputFiles = inputFiles

#config.Data.inputDataset = '/GGToMuMu_Pt-50_Elastic_13TeV-lpair/RunIIFall17DRPremix-94X_mc2017_realistic_v11-v1/AODSIM'
config.Data.inputDataset = '/DoubleMuon/Run2017C-17Nov2017-v1/AOD'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1
#config.Data.totalUnits = 5
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.totalUnits = 2000
config.Data.lumiMask = 'combined_RPIN_CMS.json'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/group/phys_pps/dilepton/MC/test/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName

config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_BR_UERJ'
#config.Site.whitelist = ['T2_CH_CERN']
#config.Site.whitelist = ['T2_BR_UERJ']

