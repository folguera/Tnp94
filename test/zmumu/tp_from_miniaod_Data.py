import FWCore.ParameterSet.Config as cms

import subprocess

process = cms.Process("TagProbe")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(),
)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

import os
if "CMSSW_7_4_" in os.environ['CMSSW_VERSION']:

    #run 251168
    process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v1')
    sourcefilesfolder = "/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/168/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames = [ sourcefilesfolder+"/"+f for f in files.split() ]

    #run 251244
    sourcefilesfolder = "/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/244/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames.extend( [ sourcefilesfolder+"/"+f for f in files.split() ] )

    #run 251251
    sourcefilesfolder = "/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/251/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames.extend( [ sourcefilesfolder+"/"+f for f in files.split() ] )

    #run 251252
    sourcefilesfolder = "/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/252/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames.extend( [ sourcefilesfolder+"/"+f for f in files.split() ] )

    # to add following runs: 251491, 251493, 251496, ..., 251500 
    print process.source.fileNames
    #print process.source.fileNames, dataSummary
elif "CMSSW_7_6_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('76X_dataRun2_v15')
    process.source.fileNames = [
            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/00A3E567-75A8-E511-AD0D-0CC47A4D769E.root',
            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/06CC1B3A-FDA7-E511-B02B-00259073E388.root',
            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/0A9FEDA2-6DA8-E511-A451-002590596490.root',
            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/0AEF074D-EBA7-E511-B229-0002C94CDAF4.root',
            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/12998942-7BA8-E511-B1AA-003048FFCB84.root',
            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/145E4DB2-EFA7-E511-8E21-00266CF3DFE0.root',
            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/148E0F6C-EEA7-E511-A70E-0090FAA588B4.root',
            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/149A16F7-6DA8-E511-8A40-003048FFCC0A.root',
            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/18D542EB-FAA7-E511-A011-00259073E4E8.root',
            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/24537A2D-0BA8-E511-8D7C-20CF300E9ECF.root',
    ]
elif "CMSSW_8_0_"in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('80X_dataRun2_Prompt_v9')

    process.source.fileNames = [
        '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/276/283/00000/0001E5C0-AE44-E611-9F88-02163E014235.root',
        '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/276/283/00000/007E4250-AE44-E611-867E-02163E011AB6.root',
        '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/276/283/00000/00997A4B-B044-E611-9FBB-02163E011EDE.root',
        '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/276/283/00000/02BB51AA-B044-E611-8DB0-02163E014168.root',
        '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/276/283/00000/0466BA91-AE44-E611-825B-02163E0136EF.root',
        '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/276/283/00000/0485506E-AE44-E611-A24B-02163E0140ED.root',
        '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/276/283/00000/0494A580-B044-E611-993A-02163E012944.root',
        '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/276/283/00000/04C4B374-B044-E611-97D0-02163E011ECD.root',
        '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/276/283/00000/06056373-B044-E611-B41D-02163E0137AA.root',
        '/store/data/Run2016C/SingleMuon/AOD/PromptReco-v2/000/276/283/00000/064D926A-B044-E611-9CAA-02163E011FCC.root',
        ]
elif "CMSSW_9_2_"in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('92X_dataRun2_Express_v2')

    process.source.fileNames = [
        '/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/101/00000/0C01D9CD-D253-E711-9D2F-02163E013511.root'
    ]  
elif "CMSSW_9_4_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('91X_mcRun2_asymptotic_v3')

    process.source.fileNames = [
            'file:/afs/cern.ch/user/f/folguera/workdir//trees/SingleMuon_Run2017B_MiniAOD_17Nov2017-v1.root'
    ] 
else: raise RuntimeError, "Unknown CMSSW version %s" % os.environ['CMSSW_VERSION']

## SELECT WHAT DATASET YOU'RE RUNNING ON
TRIGGER="SingleMu"
#TRIGGER="DoubleMu"

## ==== Fast Filters ====
process.goodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlineSlimmedPrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
    filter = cms.bool(True),
)
process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")


if TRIGGER == "SingleMu":
    process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_Mu50_v*','HLT_IsoMu27_v*', 'HLT_IsoMu24_v*','HLT_IsoMu20_v*')
elif TRIGGER == "DoubleMu":
    process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_Mu8_v*', 'HLT_Mu17_v*',
                                                                  'HLT_Mu8_TrkIsoVVL_v*', 'HLT_Mu17_TrkIsoVVL_v*',
                                                                  'HLT_Mu17_TkMu8_v*', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*' )
else:
    raise RuntimeError, "TRIGGER must be 'SingleMu' or 'DoubleMu'"

process.triggerResultsFilter.l1tResults = "gtStage2Digis"
process.triggerResultsFilter.throw = False
process.triggerResultsFilter.hltResults = cms.InputTag("TriggerResults","","HLT")

#decomment when you have it
#process.triggerResultsFilterFake = process.triggerResultsFilter.clone(
#    triggerConditions = cms.vstring( 'HLT_Mu40_v*', 'HLT_Mu5_v*', 'HLT_Mu12_v*', 'HLT_Mu24_v*')
#)

process.fastFilter     = cms.Sequence(process.goodVertexFilter + process.triggerResultsFilter)

##    __  __                       
##   |  \/  |_   _  ___  _ __  ___ 
##   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|\___/|_| |_|___/
##                                 
## ==== Merge CaloMuons and Tracks into the collection of reco::Muons  ====
from MuonAnalysis.TagAndProbe.common_variables_miniAOD_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_miniAOD_cff")

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
## with some customization
process.muonMatchHLTL2.maxDeltaR = 0.3 # Zoltan tuning - it was 0.5
process.muonMatchHLTL3.maxDeltaR = 0.1
process.patMuonsWithTriggerSequence.remove(process.patTriggerFull)
process.patMuonsWithTriggerSequence.remove(process.patTrigger)
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "slimmedMuons")
useL1Stage2Candidates(process)
appendL1MatchingAlgo(process)
useExistingPATMuons(process,"slimmedMuons")
#addHLTL1Passthrough(process)
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
massSearchReplaceAnyInputTag(process.patMuonsWithTriggerSequence, cms.InputTag('patTrigger'), 'slimmedPatTrigger')
## end of customization for patMuonsWithTrigger

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > 15 && "+MuonIDFlags.Tight2012.value()+
                     " && !triggerObjectMatchesByCollection('hltIterL3MuonCandidates').empty()"+
                     " && pfIsolationR04().sumChargedHadronPt/pt < 0.2"),
)

##SFprocess.tagMuons = cms.EDFilter("PATMuonSelector",
##SF    src = cms.InputTag("slimmedMuons"),
##SF    cut = cms.string("pt > 15 && "+MuonIDFlags.Tight2012.value()+
##SF                     " && pfIsolationR04().sumChargedHadronPt/pt < 0.2"),
##SF)

process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(""),  # no real cut now
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('60 < mass && abs(daughter(0).vz - daughter(1).vz) < 4'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("None"),
    # probe variables: all useful ones
    variables = cms.PSet(
        AllVariables,
    ),
    flags = cms.PSet(
       TrackQualityFlags,
       MuonIDFlags,
       HighPtTriggerFlags,
       HighPtTriggerFlagsDebug,
    ),
    tagVariables = cms.PSet(
        AllVariables,
        nVertices   = cms.InputTag("nverticesModule"),
        isoTrk03Abs = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsIsoFromDepsTk"),
        isoTrk03Rel = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsRelIsoFromDepsTk"),
        instLumi = cms.InputTag("addEventInfo", "instLumi"),
    ),
    tagFlags = cms.PSet(HighPtTriggerFlags,HighPtTriggerFlagsDebug),
    pairVariables = cms.PSet(
        nJets30 = cms.InputTag("njets30Module"),
        dz      = cms.string("daughter(0).vz - daughter(1).vz"),
        pt      = cms.string("pt"), 
        rapidity = cms.string("rapidity"),
        deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
        probeMultiplicity = cms.InputTag("probeMultiplicity"),
        probeMultiplicity_TMGM = cms.InputTag("probeMultiplicityTMGM"),
        probeMultiplicity_Pt10_M60140 = cms.InputTag("probeMultiplicityPt10M60140"),
        ## New TuneP variables
        newTuneP_probe_pt            = cms.InputTag("newTunePVals", "pt"),
        newTuneP_probe_sigmaPtOverPt = cms.InputTag("newTunePVals", "ptRelError"),
        newTuneP_probe_trackType     = cms.InputTag("newTunePVals", "trackType"),
        newTuneP_mass                = cms.InputTag("newTunePVals", "mass"),
    ),
    pairFlags = cms.PSet(
        BestZ = cms.InputTag("bestPairByZMass"),
    ),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
    ignoreExceptions = cms.bool(True),                            
)
if TRIGGER == "DoubleMu":
    for K,F in MuonIDFlags.parameters_().iteritems():
        setattr(process.tpTree.tagFlags, K, F)


process.load("MuonAnalysis.TagAndProbe.muon.tag_probe_muon_extraIso_cfi")
process.load("PhysicsTools.PatAlgos.recoLayer0.pfParticleSelectionForIso_cff")

process.miniIsoSeq = cms.Sequence(
    process.pfParticleSelectionForIsoSequence +
    process.muonMiniIsoCharged + 
    process.muonMiniIsoPUCharged + 
    process.muonMiniIsoNeutrals + 
    process.muonMiniIsoPhotons 
)

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuons +
    process.tpPairs    +
    process.onePair    +
    process.nverticesModule +
    process.njets30Module +
    process.probeMultiplicities + 
    process.addEventInfo +
    process.bestPairByZMass + 
    process.newTunePVals +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.fastFilter +
#    process.mergedMuons                 *
    process.patMuonsWithTriggerSequence +
    process.tnpSimpleSequence
)

##    _____               _    _             
##   |_   _| __ __ _  ___| | _(_)_ __   __ _ 
##     | || '__/ _` |/ __| |/ / | '_ \ / _` |
##     | || | | (_| | (__|   <| | | | | (_| |
##     |_||_|  \__,_|\___|_|\_\_|_| |_|\__, |
##                                     |___/ 

###SF ## Then make another collection for standalone muons, using standalone track to define the 4-momentum
###SF process.muonsSta = cms.EDProducer("RedefineMuonP4FromTrack",
###SF     src   = cms.InputTag("muons"),
###SF     track = cms.string("outer"),
###SF )
###SF ## Match to trigger, to measure the efficiency of HLT tracking
###SF from PhysicsTools.PatAlgos.tools.helpers import *
###SF process.patMuonsWithTriggerSequenceSta = cloneProcessingSnippet(process, process.patMuonsWithTriggerSequence, "Sta")
###SF process.patMuonsWithTriggerSequenceSta.replace(process.patTriggerFullSta, process.patTriggerFull)
###SF process.patTriggerSta.src = 'patTriggerFull'
###SF process.muonMatchHLTL2Sta.maxDeltaR = 0.5
###SF process.muonMatchHLTL3Sta.maxDeltaR = 0.5
###SF massSearchReplaceAnyInputTag(process.patMuonsWithTriggerSequenceSta, "mergedMuons", "muonsSta")
###SF 
###SF ## Define probes and T&P pairs
###SF process.probeMuonsSta = cms.EDFilter("PATMuonSelector",
###SF     src = cms.InputTag("patMuonsWithTriggerSta"),
###SF     cut = cms.string("outerTrack.isNonnull"), # no real cut now
###SF )
###SF process.pseudoProbeSta = cms.EDFilter("MuonSelector",
###SF     src = cms.InputTag("muonsSta"),
###SF     cut = cms.string("outerTrack.isNonnull"),
###SF )
###SF 
###SF 
###SF process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta@-", cut = '40 < mass < 150')
###SF 
###SF process.onePairSta = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairsSta"), minNumber = cms.uint32(1))
###SF 
###SF process.pseudoPairsSta = process.tpPairsSta.clone(decay = "pseudoTag@+ pseudoProbeSta@-")
###SF process.onePseudoPairSta = process.onePairSta.clone(src = 'pseudoPairsSta')
###SF process.fastPseudoTnPSta = cms.Sequence(process.pseudoTag + process.muonsSta + process.pseudoProbeSta + process.pseudoPairsSta + process.onePseudoPairSta)
###SF 
###SF process.staToTkMatch.maxDeltaR     = 0.3
###SF process.staToTkMatch.maxDeltaPtRel = 2.
###SF process.staToTkMatchNoZ.maxDeltaR     = 0.3
###SF process.staToTkMatchNoZ.maxDeltaPtRel = 2.
###SF 
###SF process.load("MuonAnalysis.TagAndProbe.tracking_reco_info_cff")
###SF 
###SF process.tpTreeSta = process.tpTree.clone(
###SF     tagProbePairs = "tpPairsSta",
###SF     arbitration   = "OneProbe",
###SF     variables = cms.PSet(
###SF         KinematicVariables, 
###SF         StaOnlyVariables,
###SF         ## track matching variables
###SF         tk_deltaR     = cms.InputTag("staToTkMatch","deltaR"),
###SF         tk_deltaEta   = cms.InputTag("staToTkMatch","deltaEta"),
###SF         tk_deltaR_NoZ   = cms.InputTag("staToTkMatchNoZ","deltaR"),
###SF         tk_deltaEta_NoZ = cms.InputTag("staToTkMatchNoZ","deltaEta"),
###SF     ),
###SF     flags = cms.PSet(
###SF         outerValidHits = cms.string("outerTrack.numberOfValidHits > 0"),
###SF         TM  = cms.string("isTrackerMuon"),
###SF         Glb = cms.string("isGlobalMuon"),
###SF         Tk  = cms.string("track.isNonnull"),
###SF         StaTkSameCharge = cms.string("outerTrack.isNonnull && innerTrack.isNonnull && (outerTrack.charge == innerTrack.charge)"),
###SF     ),
###SF     tagVariables = cms.PSet(
###SF         pt = cms.string("pt"),
###SF         eta = cms.string("eta"),
###SF         phi = cms.string("phi"),
###SF         nVertices = cms.InputTag("nverticesModule"),
###SF         combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
###SF         chargedHadIso04 = cms.string("pfIsolationR04().sumChargedHadronPt"),
###SF         neutralHadIso04 = cms.string("pfIsolationR04().sumNeutralHadronEt"),
###SF         photonIso04 = cms.string("pfIsolationR04().sumPhotonEt"),
###SF         combRelIsoPF04dBeta = IsolationVariables.combRelIsoPF04dBeta,
###SF         l1rate = cms.InputTag("l1rate"),
###SF         bx     = cms.InputTag("l1rate","bx"),
###SF         instLumi = cms.InputTag("addEventInfo", "instLumi"),
###SF     ),
###SF     pairVariables = cms.PSet(
###SF         nJets30 = cms.InputTag("njets30ModuleSta"),
###SF         dz      = cms.string("daughter(0).vz - daughter(1).vz"),
###SF         pt      = cms.string("pt"), 
###SF         rapidity = cms.string("rapidity"),
###SF         deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
###SF     ),
###SF     pairFlags = cms.PSet(),
###SF )
###SF process.njets30ModuleSta = process.njets30Module.clone(pairs = "tpPairsSta")
###SF 
###SF process.tnpSimpleSequenceSta = cms.Sequence(
###SF     process.tagMuons +
###SF     process.oneTag     +
###SF     process.probeMuonsSta   +
###SF     process.tpPairsSta      +
###SF     process.onePairSta      +
###SF     process.nverticesModule +
###SF     process.staToTkMatchSequenceZ +
###SF     process.njets30ModuleSta +
###SF     process.addEventInfo +
###SF     process.l1rate #+
###SF   #  process.tpTreeSta
###SF )
###SF 
###SF ## Add extra RECO-level info
###SF if False:
###SF     process.tnpSimpleSequenceSta.replace(process.tpTreeSta, process.tkClusterInfo+process.tpTreeSta)
###SF     process.tpTreeSta.tagVariables.nClustersStrip = cms.InputTag("tkClusterInfo","siStripClusterCount")
###SF     process.tpTreeSta.tagVariables.nClustersPixel = cms.InputTag("tkClusterInfo","siPixelClusterCount")
###SF     process.tnpSimpleSequenceSta.replace(process.tpTreeSta, process.tkLogErrors+process.tpTreeSta)
###SF     process.tpTreeSta.tagVariables.nLogErrFirst = cms.InputTag("tkLogErrors","firstStep")
###SF     process.tpTreeSta.tagVariables.nLogErrPix   = cms.InputTag("tkLogErrors","pixelSteps")
###SF     process.tpTreeSta.tagVariables.nLogErrAny   = cms.InputTag("tkLogErrors","anyStep")
###SF 
###SF if True: 
###SF     process.tracksNoMuonSeeded = cms.EDFilter("TrackSelector",
###SF               src = cms.InputTag("packedPFCandidates"),
###SF               cut = cms.string("")  #|| ".join("isAlgoInMask('%s')" % a for a in [
###SF #              'initialStep', 'lowPtTripletStep', 'pixelPairStep', 'detachedTripletStep',
###SF #              'mixedTripletStep', 'pixelLessStep', 'tobTecStep', 'jetCoreRegionalStep',
###SF #              'lowPtQuadStep', 'highPtTripletStep', 'detachedQuadStep' ] ) )
###SF     )
###SF     process.pCutTracks0 = process.pCutTracks.clone(src = 'packedPFCandidates')
###SF     process.tkTracks0 = process.tkTracks.clone(src = 'pCutTracks0')
###SF     process.tkTracksNoZ0 = process.tkTracksNoZ.clone(src = 'tkTracks0')
###SF     process.preTkMatchSequenceZ.replace(
###SF             process.tkTracksNoZ, process.tkTracksNoZ +
###SF             process.tracksNoMuonSeeded + process.pCutTracks0 + process.tkTracks0 + process.tkTracksNoZ0)
###SF     process.staToTkMatch0 = process.staToTkMatch.clone(matched = 'tkTracks0')
###SF     process.staToTkMatchNoZ0 = process.staToTkMatchNoZ.clone(matched = 'tkTracksNoZ0')
###SF     process.staToTkMatchSequenceZ.replace( process.staToTkMatch, process.staToTkMatch + process.staToTkMatch0 )
###SF     process.staToTkMatchSequenceZ.replace( process.staToTkMatchNoZ, process.staToTkMatchNoZ + process.staToTkMatchNoZ0 )
###SF     process.tpTreeSta.variables.tk0_deltaR     = cms.InputTag("staToTkMatch0","deltaR")
###SF     process.tpTreeSta.variables.tk0_deltaEta   = cms.InputTag("staToTkMatch0","deltaEta")
###SF     process.tpTreeSta.variables.tk0_deltaR_NoZ   = cms.InputTag("staToTkMatchNoZ0","deltaR")
###SF     process.tpTreeSta.variables.tk0_deltaEta_NoZ = cms.InputTag("staToTkMatchNoZ0","deltaEta")
###SF 
###SF process.tagAndProbeSta = cms.Path( 
###SF     process.fastFilter +
###SF     process.fastPseudoTnPSta +
###SF     process.mergedMuons * process.patMuonsWithTriggerSequence +
###SF     process.muonsSta                       +
###SF     process.patMuonsWithTriggerSequenceSta +
###SF     process.tnpSimpleSequenceSta
###SF )
###SF 
###SF 
###SF if True: # turn on for tracking efficiency using L1 seeds
###SF     process.probeL1 = cms.EDFilter("CandViewSelector",
###SF         src = cms.InputTag("l1extraParticles"),
###SF         cut = cms.string("pt >= 5 && abs(eta) < 2.4"),
###SF     )
###SF     process.tpPairsTkL1 = process.tpPairs.clone(decay = "tagMuons@+ probeL1@-", cut = 'mass > 30')
###SF     process.onePairTkL1 = process.onePair.clone(src = 'tpPairsTkL1')
###SF     process.l1ToTkMatch    = process.staToTkMatch.clone(src = "probeL1", srcTrack="none")
###SF     process.l1ToTkMatchNoZ = process.staToTkMatchNoZ.clone(src = "probeL1", srcTrack="none")
###SF     process.l1ToTkMatch0    = process.staToTkMatch0.clone(src = "probeL1", srcTrack="none")
###SF     process.l1ToTkMatchNoZ0 = process.staToTkMatchNoZ0.clone(src = "probeL1", srcTrack="none")
###SF     process.tpTreeL1 = process.tpTreeSta.clone(
###SF         tagProbePairs = "tpPairsTkL1",
###SF         arbitration   = "OneProbe",
###SF         variables = cms.PSet(
###SF             KinematicVariables, 
###SF             bx = cms.string("bx"),
###SF             quality = cms.string("gmtMuonCand.quality"),
###SF             ## track matching variables
###SF             tk_deltaR     = cms.InputTag("l1ToTkMatch","deltaR"),
###SF             tk_deltaEta   = cms.InputTag("l1ToTkMatch","deltaEta"),
###SF             tk_deltaR_NoZ   = cms.InputTag("l1ToTkMatchNoZ","deltaR"),
###SF             tk_deltaEta_NoZ = cms.InputTag("l1ToTkMatchNoZ","deltaEta"),
###SF             ## track matching variables (early general tracks)
###SF             tk0_deltaR     = cms.InputTag("l1ToTkMatch0","deltaR"),
###SF             tk0_deltaEta   = cms.InputTag("l1ToTkMatch0","deltaEta"),
###SF             tk0_deltaR_NoZ   = cms.InputTag("l1ToTkMatchNoZ0","deltaR"),
###SF             tk0_deltaEta_NoZ = cms.InputTag("l1ToTkMatchNoZ0","deltaEta"),
###SF         ),
###SF         flags = cms.PSet(
###SF         ),
###SF         tagVariables = cms.PSet(
###SF             pt = cms.string("pt"),
###SF             eta = cms.string("eta"),
###SF             phi = cms.string("phi"),
###SF             nVertices   = cms.InputTag("nverticesModule"),
###SF             combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
###SF             chargedHadIso04 = cms.string("pfIsolationR04().sumChargedHadronPt"),
###SF             neutralHadIso04 = cms.string("pfIsolationR04().sumNeutralHadronEt"),
###SF             photonIso04 = cms.string("pfIsolationR04().sumPhotonEt"),
###SF             combRelIsoPF04dBeta = IsolationVariables.combRelIsoPF04dBeta,
###SF         ),
###SF         pairVariables = cms.PSet(
###SF             #nJets30 = cms.InputTag("njets30ModuleSta"),
###SF             pt      = cms.string("pt"),
###SF             rapidity = cms.string("rapidity"),
###SF             deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
###SF         ),
###SF         pairFlags = cms.PSet(),
###SF         allProbes     = cms.InputTag("probeL1"),
###SF     )
###SF     process.pseudoPairsTkL1 = process.tpPairsSta.clone(decay = "pseudoTag@+ probeL1@-", cut = 'mass > 30')
###SF     process.onePseudoPairTkL1 = process.onePairSta.clone(src = 'pseudoPairsTkL1')
###SF     process.fastPseudoTnPTkL1 = cms.Sequence(process.pseudoTag + process.probeL1 + process.pseudoPairsTkL1 + process.onePseudoPairTkL1)
###SF     process.tagAndProbeTkL1 = cms.Path(
###SF         process.fastFilter +
###SF         process.fastPseudoTnPTkL1 +
###SF         process.mergedMuons * process.patMuonsWithTriggerSequence +
###SF         process.tagMuons + 
###SF         process.oneTag     +
###SF         process.probeL1 +
###SF         process.tpPairsTkL1 +
###SF         process.onePairTkL1    +
###SF         process.preTkMatchSequenceZ +
###SF         process.l1ToTkMatch + process.l1ToTkMatchNoZ +
###SF         process.l1ToTkMatch0 + process.l1ToTkMatchNoZ0 +
###SF         process.nverticesModule + process.l1rate +
###SF         process.tpTreeL1
###SF     )
###SF 
###SF 
###SF ##    _____     _          ____       _            
###SF ##   |  ___|_ _| | _____  |  _ \ __ _| |_ ___  ___ 
###SF ##   | |_ / _` | |/ / _ \ | |_) / _` | __/ _ \/ __|
###SF ##   |  _| (_| |   <  __/ |  _ < (_| | ||  __/\__ \
###SF ##   |_|  \__,_|_|\_\___| |_| \_\__,_|\__\___||___/
###SF ##                                                 
###SF ##   
###SF process.load("MuonAnalysis.TagAndProbe.fakerate_all_cff")
###SF 
###SF process.fakeRateJetPlusProbeTree = process.tpTree.clone(
###SF     tagProbePairs = 'jetPlusProbe',
###SF     arbitration   = 'None', 
###SF     tagVariables = process.JetPlusProbeTagVariables,
###SF     tagFlags = cms.PSet(),
###SF     pairVariables = cms.PSet(deltaPhi = cms.string("deltaPhi(daughter(0).phi, daughter(1).phi)")), 
###SF     pairFlags     = cms.PSet(), 
###SF )
###SF process.fakeRateWPlusProbeTree = process.tpTree.clone(
###SF     tagProbePairs = 'wPlusProbe',
###SF     arbitration   = 'None', 
###SF     tagVariables = process.WPlusProbeTagVariables,
###SF     tagFlags = cms.PSet(),
###SF     pairVariables = cms.PSet(), 
###SF     pairFlags     = cms.PSet(SameSign = cms.string('daughter(0).daughter(0).charge == daughter(1).charge')), 
###SF )
###SF process.fakeRateZPlusProbeTree = process.tpTree.clone(
###SF     tagProbePairs = 'zPlusProbe',
###SF     arbitration   = 'None', 
###SF     tagVariables  = process.ZPlusProbeTagVariables,
###SF     tagFlags      = cms.PSet(),
###SF     pairVariables = cms.PSet(), 
###SF     pairFlags     = cms.PSet(), 
###SF )
###SF 
###SF process.fakeRateJetPlusProbe = cms.Path(
###SF     #process.fastFilterFake +
###SF     process.jetPlusProbeSequenceFast +
###SF     process.mergedMuons * process.patMuonsWithTriggerSequence +
###SF     process.tagMuons + process.probeMuons + 
###SF     process.jetPlusProbeSequence +
###SF     process.extraProbeVariablesSeq + 
###SF     process.fakeRateJetPlusProbeTree
###SF )
###SF process.fakeRateWPlusProbe = cms.Path(
###SF     process.fastFilter +
###SF     process.mergedMuons * process.patMuonsWithTriggerSequence +
###SF     process.tagMuons + process.probeMuons + 
###SF     process.wPlusProbeSequence +
###SF     process.extraProbeVariablesSeq + 
###SF     process.fakeRateWPlusProbeTree
###SF )
###SF process.fakeRateZPlusProbe = cms.Path(
###SF     process.fastFilter +
###SF     process.mergedMuons * process.patMuonsWithTriggerSequence +
###SF     process.tagMuons + process.probeMuons + 
###SF     process.zPlusProbeSequence +
###SF     process.extraProbeVariablesSeq + 
###SF     process.fakeRateZPlusProbeTree
###SF )

process.schedule = cms.Schedule(
   process.tagAndProbe, 
#   process.tagAndProbeSta, 
#   process.tagAndProbeTkL1
)

process.RandomNumberGeneratorService.tkTracksNoZ = cms.PSet( initialSeed = cms.untracked.uint32(81) )
process.RandomNumberGeneratorService.tkTracksNoZ0 = cms.PSet( initialSeed = cms.untracked.uint32(81) )


####SF if TRIGGER == "SingleMu": 
####SF     process.schedule.extend([
####SF        process.fakeRateJetPlusProbe,
####SF        process.fakeRateWPlusProbe,
####SF        process.fakeRateZPlusProbe,
####SF     ])

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ_Data.root"))
