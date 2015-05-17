import FWCore.ParameterSet.Config as cms

myVariables = cms.PSet(
    mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
    probe_pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
    probe_abseta = cms.vstring("Probe #eta", "0", "2.5", ""),
    probe_nJets = cms.vstring("Jet Mult.", "0", "7", ""),
    run = cms.vstring("Run", "160000", "200000", ""),
    event_nPV = cms.vstring("NVtx", "0", "200", ""),
    probe_sc_et = cms.vstring("Probe E_{T}", "0", "200", "GeV/c"),
    probe_sc_abseta = cms.vstring("Probe |#eta|", "0.", "2.5", ""),    
    probe_gsfEle_pt =  cms.vstring("Probe E_{T}", "0", "200", "GeV/c"),
    probe_gsfEle_et =  cms.vstring("Probe E_{T}", "0", "200", "GeV/c"),
    probe_et =  cms.vstring("Probe E_{T}", "0", "200", "GeV/c"),
    probe_gsfEle_abseta = cms.vstring("Probe |#eta|", "0.", "2.5", "")
    )

ptBins = cms.PSet(
    EfficiencyCategoryAndState = cms.vstring("probe_isWPMedium", "pass"),
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(
        probe_gsfEle_pt = cms.vstring("Probe p_{T}", "0", "200", "GeV/c"),
        probe_sc_abseta = cms.vstring("Probe |#eta|", "0.", "2.5", ""),
        #BinToPDFmap = cms.vstring("fittingFunction")
        )
    )

ptEtaBins = ptBins.clone()
ptEtaBins.BinnedVariables = cms.PSet(
    #probe_nJets = cms.vdouble(0,1),
    #probe_et = cms.vdouble(10,15,20,30,40,50,200),
    #probe_abseta = cms.vdouble(0.0,0.8,1.4442,1.556,2.0,2.5),
    probe_gsfEle_pt = cms.vdouble(10,15,20,30,40,50,200),
    probe_sc_abseta = cms.vdouble(0.0,0.8,1.4442,1.566,2.0,2.5),
    mcTrue = cms.vstring("true"),
    #probe_sc_et = cms.vdouble(10,15,20,30,40,50,200),
    #probe_sc_abseta = cms.vdouble(0.0,0.8,1.4442,1.556,2.0,2.5)
    )

effiVar = cms.PSet(binning = ptBins)

categories =  cms.PSet(
    probe_isWPMedium = cms.vstring("probe_passingSomething", "dummy[pass=1,fail=0]"),
    mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
)

CommonDummyTree = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                 InputFileNames = cms.vstring(""),
                                 InputDirectoryName = cms.string("ClassOneToClassTwo"),
                                 InputTreeName = cms.string("mcUnbias_tree"),#fitter_tree"),
                                 OutputFileName = cms.string(""),
                                 NumCPU = cms.uint32(1),
                                 SaveWorkspace = cms.bool(True),
                                 floatShapeParameters = cms.bool(True),
                                 Variables = myVariables,
                                 Categories = categories,
                                 PDFs = cms.PSet(gaussPlusLinear = 
                                                 cms.vstring("Voigtian::signal1(mass, mean1[90,80,100], width[2.495,-5,5],sigma1[2,-3,3])",
                                                             "Voigtian::signal2(mass, mean2[90,80,100], width[2.495,-5,5],sigma2[4,-10,10])",
                                                             "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
                                                             "Exponential::backgroundPass(mass, lp[-0.1,-1,1])",
                                                             "Exponential::backgroundFail(mass, lf[-0.1,-1,1])",
                                                             "efficiency[0.9,0,1]",
                                                             "signalFractionInPassing[0.9]"
                                                             ),
                                                 ),
                                 
                                 Efficiencies = effiVar
                                 )


def ConfigureTPAnalyzer(Analyzer, fileName, pred, binScheme, tree="fitter_tree"):
    """
    This one configures the TnP Analyzer
    """
    Analyzer.InputFileNames = cms.vstring("/hdfs/store/user/kaur/"+fileName+".root")
    print "/hdfs/store/user/kaur/"+fileName+".root"
    Analyzer.Categories.probe_isWPMedium = cms.vstring((pred, "dummy[pass=1,fail=0]"))
    Analyzer.Efficiencies.binning = binScheme
    Analyzer.Efficiencies.binning.EfficiencyCategoryAndState = cms.vstring(pred, "pass")
    Analyzer.InputTreeName = cms.string(tree)

    
def ActivateAndConfigureFit(Analyzer, pdfPSet, doBinnedFit, NBinsForFit=0):
    """
    This one configures the fit efficiency
    calculation.
    """
    Analyzer.PDFs = pdfPSet
    Analyzer.Efficiencies.binning.BinToPDFmap = cms.vstring("fittingFunction")
    Analyzer.binnedFit = cms.bool(doBinnedFit)
    Analyzer.binsForFit = cms.uint32(NBinsForFit)

            


