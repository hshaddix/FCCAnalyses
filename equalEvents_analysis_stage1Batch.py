#Mandatory: List of processes
processList = {
    'p8_ee_ZZ_ecm240':{'fraction':.1,'chunks':10},#Run the full statistics in 10 jobs in output dir <outputDir>/p8_ee_ZZ_ecm240/chunk<N>.root
    'p8_ee_WW_ecm240':{'fraction':.1,'chunks':10},#Run the full statistics in 10 jobs in output dir <outputDir>/p8_ee_WW_ecm240/chunk<N>.root
    'p8_ee_ZH_ecm240':{'fraction':.1,'chunks':10}, #Run the full statistics in 10 jobs in output dir <outputDir>/p8_ee_ZH_ecm240/chunk<N>.root
    'p8_ee_Zqq_ecm240':{'fraction':.1,'chunks':10},
    'p8_ee_Zll_ecm240':{'fraction':.1,'chunks':10},
    'kkmcp8_ee_mumu_noFSR_ecm240':{'fraction':.1,'chunks':5},
    'wzp6_ee_eeH_ecm240':{'chunks':5},
    'wzp6_ee_mumuH_ecm240':{'chunks':10},
    'wzp6_ee_nuenueZ_ecm240':{'chunks':10},
    'wzp6_ee_nunuH_ecm240':{'chunks':10},
    'wzp6_ee_qqH_ecm240':{'chunks':10},
    'wzp6_ee_tautau_ecm240':{'chunks':10},
    'wzp6_ee_mumu_ecm240':{'fraction':.1,'chunks':10},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/spring2021/IDEA/"

#Optional: output directory, default is local dir
outputDir   = "stage1"

#Optional: ncpus, default is 4
nCPUS       = 8

#Optional running on HTCondor, default is False
runBatch    = True

#Optional batch queue name when running on HTCondor, default is workday
batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
compGroup = "group_u_FCC.local_gen"

#Optional output directory on eos, if specified files will be copied there once the batch job is done, default is empty
outputDirEos = "/eos/user/h/hshaddix/MainTest"

#Optional type for eos, needed when <outputDirEos> is specified. The default is FCC eos which is eospublic
eosType = "eosuser"


PDGID = 36
#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df

            #POTENTIAL ISSUE
            .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")                  
            .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")                  
            .Alias("Particle0", "Particle#0.index")                                       
            .Alias("Particle1", "Particle#1.index")

            # define an alias for muon index collection
            .Alias("AllMuon0", "AllMuon#0.index")
            # define the muon collection
            .Define("muons",  "ReconstructedParticle::get(AllMuon0, ReconstructedParticles)")
            #select muons on pT
            .Define("selected_muons", "ReconstructedParticle::sel_pt(0)(muons)")
            # create branch with muon transverse momentum
            .Define("selected_muons_pt", "ReconstructedParticle::get_pt(selected_muons)")
            # create branch with muon rapidity
            .Define("selected_muons_y",  "ReconstructedParticle::get_y(selected_muons)")
            # create branch with muon total momentum
            .Define("selected_muons_p",     "ReconstructedParticle::get_p(selected_muons)")
            # create branch with muon energy
            .Define("selected_muons_e",     "ReconstructedParticle::get_e(selected_muons)")
            # find zed candidates from  di-muon resonances
            .Define("zed_leptonic",         "ReconstructedParticle::resonanceBuilder(5)(selected_muons)")
            # create branch with zed mass
            .Define("zed_leptonic_m",       "ReconstructedParticle::get_mass(zed_leptonic)")
            # create branch with zed transverse momenta
            .Define("zed_leptonic_pt",      "ReconstructedParticle::get_pt(zed_leptonic)")
            # calculate recoil of zed_leptonic
            .Define("zed_leptonic_recoil",  "ReconstructedParticle::recoilBuilder(240)(zed_leptonic)")
            # create branch with recoil mass
            .Define("zed_leptonic_recoil_m","ReconstructedParticle::get_mass(zed_leptonic_recoil)")
            # create branch with leptonic charge
            .Define("zed_leptonic_charge","ReconstructedParticle::get_charge(zed_leptonic)")
            # Filter at least one candidate
            .Filter("zed_leptonic_recoil_m.size()>0")

            .Define("selected_muons_px", "ReconstructedParticle::get_px(selected_muons)")
            .Define("selected_muons_py", "ReconstructedParticle::get_py(selected_muons)")
            .Define("selected_muons_pz", "ReconstructedParticle::get_pz(selected_muons)")
            .Define('missingET_px', 'MissingET.momentum.x')
            .Define('missingET_py', 'MissingET.momentum.y')
            .Define('missingET_pz', 'MissingET.momentum.z')
            .Define('missingET_e', 'MissingET.energy')
            
            .Define("A2mumu_indices","FCCAnalyses::MCParticle::get_indices_ExclusiveDecay( %s, { -13, 13 }, true, false)( Particle, Particle1)"%(PDGID))
            .Define("ARecoParticles",  "if (A2mumu_indices.size()>0) return ReconstructedParticle2MC::selRP_matched_to_list( A2mumu_indices, MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle); else return ReconstructedParticle2MC::selRP_matched_to_list( A2mumu_indices, MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle);")
#            .Define("As",  "selMC_leg(0) ( A2mumu_indices, Particle )")                           
#            .Define("Mus",  "selMC_leg(1) ( A2mumu_indices, Particle )")                          

            .Define("deltaAlpha_ave","ReconstructedParticle::angular_separationBuilder(2)( ARecoParticles )")
            
            #Gen Level
            .Define("stable",  "MCParticle::sel_genStatus(1) ( Particle )")
            .Define("MC_Muminus",  "MCParticle::sel_pdgID( 13, false) ( stable )")
            .Define("MC_Muplus",  "MCParticle::sel_pdgID( -13, false) ( stable )")
            .Define("MC_Muminus_tlv", "MCParticle::get_tlv( MC_Muminus ) ")
            .Define("MC_Muplus_tlv", "MCParticle::get_tlv( MC_Muplus ) ")

            #.Define("deltaR", "return Muminus_tlv[0].DeltaR( Muplus_tlv[0] ) ; ")
            .Define("MCdeltaR", " if ( MC_Muminus_tlv.size() > 0 && MC_Muplus_tlv.size() > 0) return MC_Muminus_tlv[0].DeltaR( MC_Muplus_tlv[0] ) ; else return -9999. ;  ")

            #Reco Level
            .Define("RC_Muminus_q", "ReconstructedParticle::sel_charge(-1.0,false)(selected_muons)")
            .Define("RC_Muplus_q", "ReconstructedParticle::sel_charge(1.0,false)(selected_muons) ")
            .Define("RC_Muminus_tlv", "ReconstructedParticle::get_tlv( RC_Muminus_q  ) ")
            .Define("RC_Muplus_tlv", "ReconstructedParticle::get_tlv( RC_Muplus_q  ) ")
            .Define("RCdeltaR", " if ( RC_Muminus_tlv.size() > 0 && RC_Muplus_tlv.size() > 0) return RC_Muminus_tlv[0].DeltaR( RC_Muplus_tlv[0] ) ; else return -9999. ;  ")
            
            
            .Define("muon_eta", "ReconstructedParticle::get_eta(selected_muons)")

            .Define("selected_muons_no", "FCCAnalyses::ReconstructedParticle::get_n(selected_muons)")
            .Filter("selected_muons_no >= 2")
            .Define("zbuilder_result", "ReconstructedParticle::resonanceBuilder_mass_recoil(5,125,0,240, false)(selected_muons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
            .Define("zll_muons", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result[1],zbuilder_result[2]}")
            .Define("reso_deltaR", "ReconstructedParticle::deltaR(zll_muons)")

            # Checking acoplanarity for signal and WW                                                                                                                            
            .Define("acoplanarity", "ReconstructedParticle::acoplanarity(selected_muons)")

            # Reducing Z/gamma to mumu backgrounds                                                                                                                               
            #.Define("missingEnergy", "FCCAnalyses::missingEnergy(240., ReconstructedParticles)")                                                                                
            #.Define("cosTheta_miss", "FCCAnalyses::get_cosTheta_miss(missingEnergy)")                                                                                           
            .Define("cosTheta_miss", "ReconstructedParticle::get_cosTheta_miss(MissingET)")


            #FSGenParticle for e+e-
            .Define("GenElectron_PID", "MCParticle::sel_pdgID(11, true)(Particle)")

            .Define("FSGenElectron", "MCParticle::sel_genStatus(1)(GenElectron_PID)") #gen status==1 means final state particle (FS)
            .Define("n_FSGenElectron", "MCParticle::get_n(FSGenElectron)")

            .Define("FSGenElectron_eta", "if (n_FSGenElectron>0) return MCParticle::get_eta(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_theta", "if (n_FSGenElectron>0) return MCParticle::get_theta(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")
            .Define("FSGenElectron_phi", "if (n_FSGenElectron>0) return MCParticle::get_phi(FSGenElectron); else return MCParticle::get_genStatus(GenElectron_PID);")


        )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
            "selected_muons_pt",
            "selected_muons_y",
            "selected_muons_p",
            "selected_muons_e",
            "zed_leptonic_pt",
            "zed_leptonic_m",
            "zed_leptonic_charge",
            "zed_leptonic_recoil_m",
            "missingET_px",
            "missingET_py",
            "missingET_pz",
            "missingET_e",
            
            "deltaAlpha_ave",
            "MCdeltaR",
            "RCdeltaR",
            "muon_eta",
            "reso_deltaR",

            "acoplanarity",
            "cosTheta_miss"
        ]
        return branchList
