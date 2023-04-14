#Input directory where the files produced at the pre-selection level are
inputDir  = "stage2"

#Input directory where the files produced at the pre-selection level are
outputDir  = "final"

#Integrated luminosity for scaling number of events (required only if setting doScale to true)
intLumi = 5.0e+06 #in pb-1

#Scale event yields by intLumi and cross section (optional)
doScale = True

processList = {
    'p8_ee_ZZ_ecm240':{},#Run over the full statistics from stage2 input file <inputDir>/p8_ee_ZZ_ecm240.root. Keep the same output name as input
    'p8_ee_WW_ecm240':{}, #Run over the statistics from stage2 input files <inputDir>/p8_ee_WW_ecm240_out/*.root. Keep the same output name as input
    'p8_ee_ZH_ecm240':{}, #Run over the full statistics from stage2 input file <inputDir>/p8_ee_ZH_ecm240_out.root. Change the output name to MySample_p8_ee_ZH_ecm240
    'Signal_ecm240':{},
    'p8_ee_Zll_ecm240':{},
    'p8_ee_Zqq_ecm240':{},
    'wzp6_ee_eeH_ecm240':{},
    'wzp6_ee_mumuH_ecm240':{},
    'wzp6_ee_mumu_ecm240':{},
    'wzp6_ee_nuenueZ_ecm240':{},
    'wzp6_ee_nunuH_ecm240':{},
    'wzp6_ee_qqH_ecm240':{},
    'wzp6_ee_tautau_ecm240':{},
     'kkmcp8_ee_mumu_noFSR_ecm240':{}
}

#Link to the dictonary that contains all the cross section informations etc...
procDict = "FCCee_procDict_spring2021_IDEA.json"

#Add MySample_p8_ee_ZH_ecm240 as it is not an offical process
procDictAdd={"Signal_ecm240":{"numberOfEvents": 100000, "sumOfWeights": 100000, "crossSection": 0.001, "kfactor": 1.0, "matchingEfficiency": 1.0}}

#Number of CPUs to use
nCPUS = 4

#produces ROOT TTrees, default is False
#doTree = False

###Dictionnay of the list of cuts. The key is the name of the selection that will be added to the output file
cutList = {#"sel0":"Zcand_q == 0",
            #"sel1":"Zcand_q == -1 || Zcand_q == 1",
            "sel0":"Zcand_m > 0 && Zcand_m < 240",
            #"sel3":"MyFilter==true && (Zcand_m < 80 || Zcand_m > 100)"
            }


#Dictionary for the ouput variable/hitograms. The key is the name of the variable in the output files. "name" is the name of the variable in the input file, "title" is the x-axis label of the histogram, "bin" the number of bins of the histogram, "xmin" the minimum x-axis value and "xmax" the maximum x-axis value.
histoList = {
    "mz":{"name":"Zcand_m","title":"m_{#mu#mu} [GeV]","bin":125,"xmin":0,"xmax":240},
    "mz_zoom":{"name":"Zcand_m","title":"m_{#mu#mu} [GeV]","bin":100,"xmin":0,"xmax":20},
    "leptonic_recoil_m":{"name":"Zcand_recoil_m","title":"Z leptonic recoil [GeV]","bin":100,"xmin":0,"xmax":200},
    "leptonic_recoil_m_zoom":{"name":"Zcand_recoil_m","title":"Z leptonic recoil [GeV]","bin":200,"xmin":80,"xmax":160},
    "leptonic_recoil_m_zoom1":{"name":"Zcand_recoil_m","title":"Z leptonic recoil [GeV]","bin":100,"xmin":120,"xmax":140},
    "leptonic_recoil_m_zoom2":{"name":"Zcand_recoil_m","title":"Z leptonic recoil [GeV]","bin":200,"xmin":120,"xmax":140},
    "leptonic_recoil_m_zoom3":{"name":"Zcand_recoil_m","title":"Z leptonic recoil [GeV]","bin":400,"xmin":120,"xmax":140},
    "leptonic_recoil_m_zoom4":{"name":"Zcand_recoil_m","title":"Z leptonic recoil [GeV]","bin":800,"xmin":120,"xmax":140},
    "leptonic_recoil_m_zoom5":{"name":"Zcand_recoil_m","title":"Z leptonic recoil [GeV]","bin":2000,"xmin":120,"xmax":140},
    "leptonic_recoil_m_zoom6":{"name":"Zcand_recoil_m","title":"Z leptonic recoil [GeV]","bin":100,"xmin":130.3,"xmax":132.5},
    "mz_1D":{"cols":["Zcand_m"],"title":"m_{Z} [GeV]", "bins": [(40,80,100)]}, # 1D histogram (alternative syntax)
    "mz_recoil_2D":{"cols":["Zcand_m", "Zcand_recoil_m"],"title":"m_{Z} - leptonic recoil [GeV]", "bins": [(40,80,100), (100,120,140)]}, # 2D histogram
    "mz_recoil_3D":{"cols":["Zcand_m", "Zcand_recoil_m", "Zcand_recoil_m"],"title":"m_{Z} - leptonic recoil - leptonic recoil [GeV]", "bins": [(40,80,100), (100,120,140), (100,120,140)]}, # 3D histogram

    "missingET_x":{"name":"missingET_px","title":"Px","bin":100,"xmin":0,"xmax":150},
    "missingET_y":{"name":"missingET_py","title":"Py","bin":100,"xmin":0,"xmax":150},
    "missingET_z":{"name":"missingET_pz","title":"Pz","bin":100,"xmin":0,"xmax":150},
    "missingET_e":{"name":"missingET_e","title":"Missing Energy [GeV]","bin":100,"xmin":0,"xmax":150},
    "MCAngleBW":{"name":"MCdeltaR","title":"Monte Carlo delta R","bin":150,"xmin":0,"xmax":5},
    "RCAngleBW":{"name":"RCdeltaR","title":"Reconstructed deltaR","bin":150,"xmin":0,"xmax":5},
#    "FS_eta":{"name":"FSGenElectron_eta","title":"Eta","bin":100,"xmin":-4,"xmax":4},
#    "FS_theta":{"name":"FSGenElectron_theta","title":"Theta","bin":100,"xmin":-4,"xmax":4},
#    "FS_phi":{"name":"FSGenElectron_phi","title":"Phi","bin":100,"xmin":-4,"xmax":4},
    "muon_eta":{"name":"muon_eta","title":"Eta","bin":100,"xmin":0,"xmax":5},
    "resoAngleBW":{"name":"reso_deltaR","title":"deltaR with Resonance Matching Muons","bin":150,"xmin":0,"xmax":5},

    "acoplanarity":{"name":"acoplanarity","title":"Acoplanarity between selected muons","bin":150,"xmin":0,"xmax":3.5},
    "missingCosTheta":{"name":"cosTheta_miss","title":"cosTheta of missing energy vector","bin":150,"xmin":0,"xmax":1.5},
    "RecoMissingEnergy":{"name":"RecoMissingEnergy_e","title":"Missing Energy [GeV]","bin":100,"xmin":0,"xmax":150},
    # 2D histograms 
    "eta_MCdeltaR":{"cols":["muon_eta", "MCdeltaR"],"title":"eta - MCdeltaR", "bins": [(100,0,4), (100,0,4)]}, # 2D histogram
    "eta_RCdeltaR":{"cols":["muon_eta", "RCdeltaR"],"title":"eta - RCdeltaR", "bins": [(100,0,4), (100,0,4)]},
    "massVdR":{"cols":["Zcand_m", "reso_deltaR"],"title":"mass - RCdeltaR", "bins": [(100,0,240), (100,0,5)]}, # 2D histogram  


}
