import ROOT

# global parameters
intLumi        = 5.0e+06 #in pb-1
ana_tex        = 'e^{+}e^{-} #rightarrow ZH #rightarrow #mu^{+}#mu^{-} + X'
delphesVersion = '3.4.2'
energy         = 240.0
collider       = 'FCC-ee'
inputDir       = "/afs/cern.ch/user/h/hshaddix/LocallyRun/FCCAnalyses/final/"
formats        = ['png','pdf']
yaxis          = ['lin','log']
stacksig       = ['stack','nostack']
outdir      = "/afs/cern.ch/user/h/hshaddix/LocallyRun/FCCAnalyses/plots"

variables = ['mz','mz_zoom','leptonic_recoil_m','leptonic_recoil_m_zoom','missingET_x','missingET_y','missingET_z','missingET_e','MCAngleBW','RCAngleBW','muon_eta','eta_MCdeltaR','eta_RCdeltaR']#,"FS_eta","FS_theta","FS_phi"]#'muons','selected_muons','selected_muons_pt'] #add variables to include muons?

###Dictionary with the analysis name as a key, and the list of selections to be plotted for this analysis. The name of the selections should be the same than in the final selection
selections = {}
selections['HA3']   = ["sel0","sel1"]
selections['HA3_2'] = ["sel0","sel1"]

extralabel = {}
extralabel['sel0'] = "Selection: N_{Z} = 1"
extralabel['sel1'] = "Selection: N_{Z} = 1; 0 GeV < m_{Z} < 150 GeV" #adding similar selection lines could allow for varied plotting 

colors = {}
colors['HA3'] = ROOT.kRed
colors['WW'] = ROOT.kBlue+1
colors['ZZ'] = ROOT.kBlue+3
colors['ZH'] = ROOT.kBlue-3
colors['Zll'] = ROOT.kAzure+7
colors['Zqq'] = ROOT.kAzure+10
colors['eeH'] = ROOT.kGreen
colors['mumuH'] = ROOT.kGreen+3
colors['mumu_wzp6'] = ROOT.kMagenta
colors['nunuZ'] = ROOT.kYellow-5
colors['nunuH'] = ROOT.kYellow-8
colors['qqH'] = ROOT.kYellow-1
colors['tautau'] = ROOT.kCyan+3
colors['mumu_noFSR'] = ROOT.kMagenta+3
colors['VV'] = ROOT.kYellow 

colors['Zl/q'] = ROOT.kBlue+1
colors['nunuX'] = ROOT.kYellow-5
colors['mumu'] = ROOT.kMagenta
colors['xxH'] = ROOT.kGreen
colors['ll'] = ROOT.kCyan+3

plots = {}
plots['HA3'] = {'signal':{'HA3':['Signal_ecm240']},
               'backgrounds':{'WW':['p8_ee_WW_ecm240'],
                              'ZZ':['p8_ee_ZZ_ecm240'],
                              'ZH':['p8_ee_ZH_ecm240'],
                              'Zll':['p8_ee_Zll_ecm240'],
                              'Zqq':['p8_ee_Zqq_ecm240'],
                              'eeH':['wzp6_ee_eeH_ecm240'],
                              'mumuH':['wzp6_ee_mumuH_ecm240'],
                              'mumu_wzp6':['wzp6_ee_mumu_ecm240'],
                              'nunuZ':['wzp6_ee_nuenueZ_ecm240'],
                              'nunuH':['wzp6_ee_nunuH_ecm240'],
                              'qqH':['wzp6_ee_qqH_ecm240'],
                              'tautau':['wzp6_ee_tautau_ecm240'],
                              'mumu_noFSR':['kkmcp8_ee_mumu_noFSR_ecm240']
}
           }


plots['HA3_2'] = {'signal':{'HA3':['Signal_ecm240']},
                  'backgrounds':{'VV':['p8_ee_WW_ecm240','p8_ee_ZZ_ecm240','p8_ee_ZH_ecm240'],
                                 'Zl/q':['p8_ee_Zll_ecm240','p8_ee_Zqq_ecm240'],
                                 'nunuX':['wzp6_ee_nuenueZ_ecm240','wzp6_ee_nunuH_ecm240'],
                                 'mumu':['wzp6_ee_mumu_ecm240','kkmcp8_ee_mumu_noFSR_ecm240','wzp6_ee_mumuH_ecm240'],
                                 'xxH':['wzp6_ee_eeH_ecm240','wzp6_ee_qqH_ecm240'],
                                 'll':['wzp6_ee_tautau_ecm240']
                                 }
             }

legend = {}
legend['HA3'] = 'Higgs-A3'
legend['WW'] = 'WW'
legend['ZZ'] = 'ZZ'
legend['VV'] = 'VV boson'
legend['ZH'] = 'ZH'
legend['Zll'] = 'Zll'
legend['Zqq'] = 'Zqq'
legend['eeH'] = 'eeH'
legend['mumuH'] = 'mumuH'
legend['mumu_wzp6'] = 'mumu'
legend['nunuZ'] = 'nunuZ'
legend['nunuH'] = 'nunuH'
legend['qqH'] = 'qqH' 
legend['tautau'] = 'tautau' 
legend['mumu_noFSR'] = 'mumu_noFSR' 
legend['Zl/q'] = 'Zll & Zqq'
legend['nunuX'] = 'nunuZ & nunuH'
legend['mumu'] = 'mumu' 
legend['xxH'] = 'xxH'
legend['ll'] = 'tautau'
