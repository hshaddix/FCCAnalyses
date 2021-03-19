#ifndef MYUTILS_ANALYZERS_H
#define MYUTILS_ANALYZERS_H
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/VertexData.h"

#include "TLorentzVector.h"

namespace myUtils{

  struct FCCAnalysesComposite{
    TLorentzVector particle;
    ROOT::VecOps::RVec<int> index;//index in the RP
    edm4hep::VertexData vertex;
    int charge;
    int mc_index;
  };

  struct filter_PV{
    filter_PV(bool arg_pv);
    bool m_pv=true;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, 
								      ROOT::VecOps::RVec<int> index);
  };


  ROOT::VecOps::RVec<FCCAnalysesComposite> add_truthmatched(ROOT::VecOps::RVec<FCCAnalysesComposite> comp,
							    ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
							    ROOT::VecOps::RVec<int> rp2mc,
							    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							    ROOT::VecOps::RVec<int> ind);
  //ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> rp2mc);

  bool isPV(edm4hep::ReconstructedParticleData recop, 
	    ROOT::VecOps::RVec<int> pvindex);

  ROOT::VecOps::RVec<int> getMC_parent(int parentindex, 
				       ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  
				       ROOT::VecOps::RVec<int> ind);

  int getMC_parent(int parentindex, 
		   edm4hep::MCParticleData in,  
		   ROOT::VecOps::RVec<int> ind);
  
  ROOT::VecOps::RVec<int> get_compmc(ROOT::VecOps::RVec<FCCAnalysesComposite> in);
  ROOT::VecOps::RVec<TLorentzVector> getFCCAnalysesComposite_particle(ROOT::VecOps::RVec<FCCAnalysesComposite> in);
  ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> getFCCAnalysesComposite_index(ROOT::VecOps::RVec<FCCAnalysesComposite> in);
  ROOT::VecOps::RVec<edm4hep::VertexData> getFCCAnalysesComposite_vertex(ROOT::VecOps::RVec<FCCAnalysesComposite> in);
  int getFCCAnalysesComposite_N(ROOT::VecOps::RVec<FCCAnalysesComposite> in);
  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_mass(ROOT::VecOps::RVec<FCCAnalysesComposite> in);
  
  ROOT::VecOps::RVec<FCCAnalysesComposite> build_Bu2D0Pi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							 ROOT::VecOps::RVec<FCCAnalysesComposite> d0, 
							 ROOT::VecOps::RVec<int> pions);

  struct build_D0 {
    build_D0(float arg_mass, float arg_p, bool arg_filterPV);
    float m_mass=0.05;    
    float m_p=1.;
    bool m_filterPV=true;
    ROOT::VecOps::RVec<FCCAnalysesComposite> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							 ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							 ROOT::VecOps::RVec<int> pions,
							 ROOT::VecOps::RVec<int> kaons, 
							 ROOT::VecOps::RVec<int> pvindex);
  };


  struct build_composite_vertex {
    build_composite_vertex(int arg_n, int arg_charge, float arg_masslow, float arg_masshigh, float arg_p, bool arg_cc, bool arg_filterPV);
    int m_n=3;
    int m_charge=0;
    float m_masslow=0.05;
    float m_masshigh=0.05;
    float m_p=1.;
    bool m_cc=true;
    bool m_filterPV=true;

    ROOT::VecOps::RVec<FCCAnalysesComposite> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							 ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							 ROOT::VecOps::RVec<int> in,
							 ROOT::VecOps::RVec<int> pvindex);
  };

  struct sel_PID {
    sel_PID(int arg_PDG);
    int m_PDG=211;    
    ROOT::VecOps::RVec<int> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);
  };

  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> PID(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							     ROOT::VecOps::RVec<int> recind, 
							     ROOT::VecOps::RVec<int> mcind, 
							     ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);
  
  ROOT::VecOps::RVec<float> awkwardtest(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,  
					ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
					ROOT::VecOps::RVec<int> recind, 
					ROOT::VecOps::RVec<int> mcind, 
					ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);
  
  float build_invmass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop, ROOT::VecOps::RVec<int> index);
  TLorentzVector build_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop, ROOT::VecOps::RVec<int> index);
  
}
#endif
