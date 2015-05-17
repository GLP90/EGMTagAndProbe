import FWCore.ParameterSet.Config as cms

process = cms.Process("tnp")

MC_flag = True

HLTProcessName = "HLT"
INPUT_FILE_NAME  = "/store/relval/CMSSW_7_4_1/RelValZEE_13/MINIAODSIM/MCRUN2_74_V9_gensim_740pre7-v1/00000/1E35CCF8-32EC-E411-8F29-0025905A48D0.root"
OUTPUT_FILE_NAME = "TnPTree.root"

ELECTRON_ET_CUT_MIN = 10.0
ELECTRON_COLL = "slimmedElectrons"
ELECTRON_CUTS = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"

SUPERCLUSTER_COLL = "reducedEgamma:reducedSuperClusters"
SUPERCLUSTER_CUTS = "abs(eta)<2.5 && !(1.4442< abs(eta) <1.566) && et>" + str(ELECTRON_ET_CUT_MIN)   

#process.pileupReweightingProducer = cms.EDProducer("PileupWeightProducer",
#                                                   FirstTime = cms.untracked.bool(True)
#                                                   )

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.throw = cms.bool(True)
process.hltHighLevel.HLTPaths = ["HLT_Ele20WP60_Ele8_Mass55_v*",
                                 "HLT_Ele25WP60_SC4_Mass55_v*"]

tpHLTTagFilter   = "hltEle20WP60Ele8TrackIsoFilter"
tpHLTProbeFilter = "hltEle20WP60Ele8PixelMatchUnseededFilter" #hltEle20WP60Ele8Mass55Filter"

tpHLTTagFilter2   = "hltEle25WP60SC4TrackIsoFilter"
tpHLTProbeFilter2 = "hltEle25WP60SC4Mass55Filter"

####process.electronTriggerMatchHLTEle20WP60Ele8Mass55 = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
####                                                                    src     = cms.InputTag(ELECTRON_COLL),
####                                                                    matched = cms.InputTag( "selectedPatTrigger" ),
####                                                                    matchedCuts = cms.string('path( "HLT_Ele20WP60_Ele8_Mass55_v*")'),
####                                                                    maxDPtRel = cms.double( 0.5 ),
####                                                                    maxDeltaR = cms.double( 0.5 ),
####                                                                    resolveAmbiguities    = cms.bool( True ),
####                                                                    resolveByMatchQuality = cms.bool( True )
####                                                                    )
####
####from PhysicsTools.PatAlgos.tools.trigTools import SwitchOnTriggerMatchingStandAlone
####switchOnTriggerMatchingStandAlone.triggerMatchers = [process.electronTriggerMatchHLTEle20WP60Ele8Mass55]
#####, triggerProducer = 'patTrigger', path = '', hltProcess = 'HLT', outputModule = 'out', postfix = '' 

#######################################################################
#TRIGGER TO BE TESTED                                                 #
#######################################################################
####HLTPathToMeasure = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15"
####
####hltLeadingLegFilter  = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter"
####hltTrailingLegFilter = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter"
####
####hltTagsPassingLeadingLegHLT= cms.VInputTag(
####    cms.InputTag(HLTPathToMeasure, hltLeadingLegFilter, HLTProcessName),
####    )
####hltTagsPassingTrailingLegHLT =cms.VInputTag(
####    cms.InputTag(HLTPathToMeasure, hltTrailingLegFilter, HLTProcessName),
####    )


########################
##    ___            _           _      
##   |_ _|_ __   ___| |_   _  __| | ___ 
##    | || '_ \ / __| | | | |/ _` |/ _ \
##    | || | | | (__| | |_| | (_| |  __/
##   |___|_| |_|\___|_|\__,_|\__,_|\___|
##
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MCRUN2_74_V9'

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##   ____             _ ____ 
##  |  _ \ ___   ___ | / ___|  ___  _   _ _ __ ___ ___ 
##  | |_) / _ \ / _ \| \___ \ / _ \| | | | '__/ __/ _ \
##  |  __/ (_) | (_) | |___) | (_) | |_| | | | (_|  __/
##  |_|   \___/ \___/|_|____/ \___/ \__,_|_|  \___\___|
##  

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(INPUT_FILE_NAME),
                            #eventsToProcess=cms.untracked.VEventRange('1:28065077','1:28065427','1:28033364','1:54509262')
                            #eventsToProcess=cms.untracked.VEventRange('1:28065077','1:28065427','1:28033364','1:54509262','1:28033351','1:28033341','1:28033347')
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.goodVertexFilter = cms.EDFilter("VertexSelector",
                                        src = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
                                        filter = cms.bool(True),
                                        )
process.noScraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

process.fastFilter = cms.Sequence(process.goodVertexFilter)# + process.noScraping)

##    ____      __ _____ _           _                   
##   / ___|___ / _| ____| | ___  ___| |_ _ __ ___  _ __  
##  | |  _/ __| |_|  _| | |/ _ \/ __| __| '__/ _ \| '_ \ 
##  | |_| \__ \  _| |___| |  __/ (__| |_| | | (_) | | | |
##   \____|___/_| |_____|_|\___|\___|\__|_|  \___/|_| |_|
##  

process.goodElectrons = cms.EDFilter("PATElectronRefSelector",
                                     src = cms.InputTag( ELECTRON_COLL ),
                                     cut = cms.string( ELECTRON_CUTS )    
                                     )

process.goodElectronsTAG = cms.EDFilter("PATElectronRefSelector",
                                        src = cms.InputTag( ELECTRON_COLL ),
                                        cut = cms.string( "(abs(superCluster.eta)<=2.5) && !(1.4442<=abs(superCluster.eta)<=1.566) && pt >= 25.0" )
                                        )

##   ____                         ____ _           _            
##  / ___| _   _ _ __   ___ _ __ / ___| |_   _ ___| |_ ___ _ __ 
##  \___ \| | | | '_ \ / _ \ '__| |   | | | | / __| __/ _ \ '__|
##   ___) | |_| | |_) |  __/ |  | |___| | |_| \__ \ ||  __/ |   
##  |____/ \__,_| .__/ \___|_|   \____|_|\__,_|___/\__\___|_|   
##  

process.superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
                                           src = cms.InputTag(SUPERCLUSTER_COLL),
                                           particleType = cms.int32(11),
                                           )

process.goodSuperClusters = cms.EDFilter("CandViewSelector",
                                         src = cms.InputTag("superClusterCands"),
                                         cut = cms.string( SUPERCLUSTER_CUTS ),
                                         filter = cms.bool(True)
                                         )                                         


##    _____ _           _                     ___    _ 
##   | ____| | ___  ___| |_ _ __ ___  _ __   |_ _|__| |
##   |  _| | |/ _ \/ __| __| '__/ _ \| '_ \   | |/ _` |
##   | |___| |  __/ (__| |_| | | (_) | | | |  | | (_| |
##   |_____|_|\___|\___|\__|_|  \___/|_| |_| |___\__,_|
##   

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

#if useAOD == True :
#    dataFormat = DataFormat.AOD
#else :
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',]
#'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

process.passingCutBasedVeto = cms.EDProducer("PatElectronSelectorByValueMap",
                                             input     = cms.InputTag(ELECTRON_COLL),
                                             selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
                                             cut       = cms.bool(True)
                                             )

process.passingCutBasedLoose = process.passingCutBasedVeto.clone()
process.passingCutBasedLoose.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose")
process.passingCutBasedMedium = process.passingCutBasedVeto.clone()
process.passingCutBasedMedium.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium")
process.passingCutBasedTight = process.passingCutBasedVeto.clone()
process.passingCutBasedTight.selection = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight")

##    _____     _                         __  __       _       _     _             
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                |___/ |___/                                                |___/ 
##   

process.goodElectronsTagHLT = cms.EDProducer("PatElectronTriggerCandProducer",
                                             filterName = cms.vstring(tpHLTTagFilter),
                                             inputs     = cms.InputTag( ELECTRON_COLL ),
                                             bits       = cms.InputTag('TriggerResults::HLT'),
                                             objects    = cms.InputTag('selectedPatTrigger'),
                                             dR         = cms.double(0.3)
                                             )

process.goodElectronsProbeHLT = cms.EDProducer("PatElectronTriggerCandProducer",
                                               filterName = cms.vstring(tpHLTTagFilter),
                                               inputs     = cms.InputTag( ELECTRON_COLL ),
                                               bits       = cms.InputTag('TriggerResults::HLT'),
                                               objects    = cms.InputTag('selectedPatTrigger'),
                                               dR         = cms.double(0.3)
                                               )
  
#process.goodSuperClustersHLT = cms.EDProducer("MiniAODTriggerCandProducer",
#                                              filterName = cms.vstring(tpHLTProbeFilter),
#                                              electrons = cms.InputTag(goodSuperClusters),
#                                              bits            = cms.InputTag('TriggerResults::HLT'),
#                                              objects = cms.InputTag('selectedPatTrigger'),
#                                              dR = cms.double(0.3)
#                                          )
  
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag(ELECTRON_COLL)
process.ele_sequence = cms.Sequence(
    process.goodElectrons +
    process.goodElectronsTAG +
    process.goodElectronsProbeHLT +
    process.goodElectronsTagHLT +
    process.egmGsfElectronIDSequence +
    process.passingCutBasedVeto +
    process.passingCutBasedLoose +
    process.passingCutBasedMedium +
    process.passingCutBasedTight 
    )

process.sc_sequence = cms.Sequence(
    process.superClusterCands +
    process.goodSuperClusters 
    #process.goodSuperClustersHLT
    )

##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   

process.tagTightSC = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("goodElectronsTAG goodSuperClusters"), 
                                    checkCharge = cms.bool(False),
                                    cut = cms.string("40<mass<1000"),
                                    )

process.tagTightRECO = cms.EDProducer("CandViewShallowCloneCombiner",
                                      decay = cms.string("goodElectronsTAG@+ goodElectrons@-"), 
                                      checkCharge = cms.bool(True),
                                      cut = cms.string("40<mass<1000"),
                                    )


process.allTagsAndProbes = cms.Sequence(
    process.tagTightSC +
    process.tagTightRECO
    )

##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                   

process.McMatchSC = cms.EDProducer("MCTruthDeltaRMatcherNew",
                                   matchPDGId = cms.vint32(11),
                                   src = cms.InputTag("goodSuperClusters"),
                                   distMin = cms.double(0.3),
                                   matched = cms.InputTag("prunedGenParticles")
                                   )
                     
process.McMatchTag = cms.EDProducer("MCTruthDeltaRMatcherNew",
                                    matchPDGId = cms.vint32(11),
                                    src = cms.InputTag("goodElectronsTAG"),
                                    distMin = cms.double(0.2),
                                    matched = cms.InputTag("prunedGenParticles"),
                                    checkCharge = cms.bool(True)
                                    )

process.McMatchRECO = cms.EDProducer("MCTruthDeltaRMatcherNew",
                                     matchPDGId = cms.vint32(11),
                                     src = cms.InputTag("goodElectrons"), #("PassingTagLegHLTTightTAG"),
                                     distMin = cms.double(0.2),
                                     matched = cms.InputTag("prunedGenParticles"),
                                     checkCharge = cms.bool(True)
                                    )

process.mc_sequence = cms.Sequence(
    process.McMatchTag  +
    process.McMatchSC +
    process.McMatchRECO
    )

############################################################################
##    _____           _       _ ____            _            _   _  ____  ##
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |/ ___| ##
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | |  _  ##
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | |_| | ##
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_|\____| ##
##              |___/                                                     ##
##                                                                        ##
############################################################################
##    ____                      _     _           
##   |  _ \ ___ _   _ ___  __ _| |__ | | ___  ___ 
##   | |_) / _ \ | | / __|/ _` | '_ \| |/ _ \/ __|
##   |  _ <  __/ |_| \__ \ (_| | |_) | |  __/\__ \
##   |_| \_\___|\__,_|___/\__,_|_.__/|_|\___||___/
##
## I define some common variables for re-use later.
## This will save us repeating the same code for each efficiency category

ZVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    pt  = cms.string("pt"),
    mass  = cms.string("mass"),
    )   

SCProbeVariablesToStore = cms.PSet(
    probe_eta    = cms.string("eta"),
    probe_abseta = cms.string("abs(eta)"),
    probe_pt     = cms.string("pt"),
    probe_et     = cms.string("et"),
    probe_e      = cms.string("energy"),
)

ProbeVariablesToStore = cms.PSet(
    probe_Ele_eta    = cms.string("eta"),
    probe_Ele_abseta = cms.string("abs(eta)"),
    probe_Ele_pt     = cms.string("pt"),
    probe_Ele_et     = cms.string("et"),
    probe_Ele_e      = cms.string("energy"),
    probe_Ele_q      = cms.string("charge"),
    #probe_Ele_trackiso = cms.string("dr03TkSumPt"),
    #probe_Ele_reltrackiso = cms.string("dr03TkSumPt/pt"),
## super cluster quantities
    probe_sc_energy = cms.string("superCluster.energy"),
    probe_sc_et     = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    probe_sc_eta    = cms.string("superCluster.eta"),
    probe_sc_abseta = cms.string("abs(superCluster.eta)"),

#id based
    probe_Ele_dEtaIn        = cms.string("deltaEtaSuperClusterTrackAtVtx"),
    probe_Ele_dPhiIn        = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    probe_Ele_sigmaIEtaIEta = cms.string("sigmaIetaIeta"),
    probe_Ele_hoe           = cms.string("hadronicOverEm"),
    probe_Ele_ooemoop       = cms.string("(1.0/ecalEnergy - eSuperClusterOverP/ecalEnergy)"),
    #probe_Ele_mHits         = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfHits")
)

TagVariablesToStore = cms.PSet(
    Ele_eta    = cms.string("eta"),
    Ele_abseta = cms.string("abs(eta)"),
    Ele_pt     = cms.string("pt"),
    Ele_et     = cms.string("et"),
    Ele_e      = cms.string("energy"),
    Ele_q      = cms.string("charge"),
    
    ## super cluster quantities
    sc_energy = cms.string("superCluster.energy"),
    sc_et     = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    sc_eta    = cms.string("superCluster.eta"),
    sc_abseta = cms.string("abs(superCluster.eta)"),
)

CommonStuffForGsfElectronProbe = cms.PSet(
    variables = cms.PSet(ProbeVariablesToStore),
    ignoreExceptions =  cms.bool (True),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool(True),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
        mass60to120 = cms.string("60 < mass < 120")
        ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags       =  cms.PSet(
        passingHLT    = cms.InputTag("goodElectronsTagHLT"),
        passingVeto   = cms.InputTag("passingCutBasedVeto"),
        passingLoose  = cms.InputTag("passingCutBasedLoose"),
        passingMedium = cms.InputTag("passingCutBasedMedium"),
        passingTight  = cms.InputTag("passingCutBasedTight"),
        ),    
    )

CommonStuffForSuperClusterProbe = CommonStuffForGsfElectronProbe.clone()
CommonStuffForSuperClusterProbe.variables = cms.PSet(SCProbeVariablesToStore)

if MC_flag:
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(True),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(False),
        checkMotherInUnbiasEff = cms.bool(False),
        mcVariables = cms.PSet(
            probe_eta = cms.string("eta"),
            probe_abseta = cms.string("abs(eta)"),
            probe_pt  = cms.string("pt"),
            probe_et  = cms.string("et"),
            probe_e  = cms.string("energy"),
            probe_mass  = cms.string("mass"),
            ),
        mcFlags     =  cms.PSet(
            probe_flag = cms.string("pt>0")
            ),      
        )
else:
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(False)
        )
    
##   ____      __       __    ___                 ___    _ 
##  / ___|___ / _|      \ \  |_ _|___  ___       |_ _|__| |
## | |  _/ __| |_   _____\ \  | |/ __|/ _ \       | |/ _` |
## | |_| \__ \  _| |_____/ /  | |\__ \ (_) |  _   | | (_| |
##  \____|___/_|        /_/  |___|___/\___/  ( ) |___\__,_|
##                                           |/            

process.GsfElectronToSC = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                         CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
                                         tagProbePairs = cms.InputTag("tagTightSC"),
                                         arbitration   = cms.string("None"),
                                         flags = cms.PSet(),
                                         probeMatches  = cms.InputTag("McMatchSC"),
                                         allProbes     = cms.InputTag("goodSuperClusters"),
                                         )

process.GsfElectronToRECO = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                           mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
                                           tagProbePairs = cms.InputTag("tagTightRECO"),
                                           arbitration   = cms.string("None"),
                                           flags         = cms.PSet(passingHLT    = cms.InputTag("goodElectronsProbeHLT"),
                                                                    passingVeto   = cms.InputTag("passingCutBasedVeto"),
                                                                    passingLoose  = cms.InputTag("passingCutBasedLoose"),
                                                                    passingMedium = cms.InputTag("passingCutBasedMedium"),
                                                                    passingTight  = cms.InputTag("passingCutBasedTight"),
                                                                    ),                                               
                                           probeMatches  = cms.InputTag("McMatchRECO"),
                                           allProbes     = cms.InputTag("goodElectrons"),
                                           )

process.tree_sequence = cms.Sequence(
    #process.GsfElectronToSC +
    process.GsfElectronToRECO 
    )    


##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##
process.out = cms.OutputModule("PoolOutputModule", 
       fileName = cms.untracked.string("PFEE.root"),
       SelectEvents = cms.untracked.PSet( 
       SelectEvents = cms.vstring("p")
       )
    )
process.outpath = cms.EndPath(process.out)
# REMOVE THIS FOR DEBUGGING PROCESS
process.outpath.remove(process.out)

process.p = cms.Path(
    process.hltHighLevel +
    process.fastFilter +
    process.ele_sequence + 
    process.sc_sequence +
    ####process.GsfDRToNearestTau+
    ####process.ElecPt+
    process.allTagsAndProbes+ 
    ####process.pileupReweightingProducer +
    process.mc_sequence +
    process.tree_sequence
    ####* process.printTree * process.printDecay
    )

process.TFileService = cms.Service(
    "TFileService", fileName = cms.string( OUTPUT_FILE_NAME )
    )
