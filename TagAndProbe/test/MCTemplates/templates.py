import ROOT

def findBins(array, var):
    size = len(array)
    bin = "dump";
    for i in xrange(size-1):
        low = array[i]
        hi  = array[i+1]
        if (low <= var and hi > var):
            bin = str(low)+"To"+str(hi)
  
    return bin;

if __name__ == "__main__":

    #pts  = (10,15,20,30,40,50,200)
    #etas = (0.0,0.8,1.4442,1.566,2.0,2.5)

    pts  = (10,200)
    etas = (0.0,1.4442,1.566,2.5)

    inFile = ROOT.TFile("../TnPTree_MC.root")
    inFile.cd("GsfElectronToRECO")
    fDir = inFile.Get("GsfElectronToRECO")
    fChain = fDir.Get("fitter_tree")
    fChain.SetBranchStatus("*", 0)
    fChain.SetBranchStatus("pair_mass60to120", 1)
    fChain.SetBranchStatus("probe_Ele_pt", 1)
    fChain.SetBranchStatus("probe_sc_abseta", 1)
    fChain.SetBranchStatus("tag_passingTight", 1)
    fChain.SetBranchStatus("passingMedium", 1)
    fChain.SetBranchStatus("passingHLT", 1)
    fChain.SetBranchStatus("tag_passingHLT", 1)
    fChain.SetBranchStatus("mcTrue", 1)
    fChain.SetBranchStatus("mass", 1)
    
    outFile = ROOT.TFile("mcTruemadgraph_Medium.root", "RECREATE")
    #hCounter = ROOT.TH1D("hcounter", "hcounter", 4, 0, 4) # TH1D needed ???
    histos = dict()

    for binPt in xrange(len(pts)-1):
        for binEta in xrange(len(etas)-1):
            histNameSt = "hMass_"+str(pts[binPt])+"To"+str(pts[binPt+1])+"_"+str(etas[binEta])+"To"+str(etas[binEta+1])
            hp = histNameSt+"_Pass"
            hf = histNameSt+"_Fail"
            histos[hp] = ROOT.TH1D(hp, hp, 30, 60, 120)
            histos[hf] = ROOT.TH1D(hf, hf, 30, 60, 120)

    histos["dump"] = ROOT.TH1D("dump", "dump", 1, 0, 1000)
  
    nentries = fChain.GetEntriesFast()
    nbytes = 0
    nb = 0

    for jentry in xrange(nentries):
        fChain.GetEntry(jentry)

        if (jentry+1) % 10000 == 0:
            print "\rprocessing event %d/%d" % (jentry+1, nentries),"(%5.1f%%)  " % ((jentry+1) / float(nentries) * 100),
            import sys
            sys.stdout.flush()

        probePt = fChain.probe_Ele_pt # perche` ???
        probeEta = fChain.probe_sc_abseta

        ptadd = findBins(pts, probePt)
        etaadd = findBins(etas, probeEta)
        key = "hMass_"+ptadd+"_"+etaadd
        if (fChain.passingMedium):
            key = key + "_Pass"
        else:
            key = key + "_Fail"
        if ("dump" in key):
            key = "dump"

        if(fChain.mcTrue and fChain.pair_mass60to120 and fChain.tag_passingTight and fChain.tag_passingHLT and fChain.passingHLT): # FIXME add Tau removal 
            histos[key].Fill(fChain.mass) # FIXME add PUweight
  
    outFile.cd()
    for k in histos:
        histos[k].Write()
    outFile.Close()
