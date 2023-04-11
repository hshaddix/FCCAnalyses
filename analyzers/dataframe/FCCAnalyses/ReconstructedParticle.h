#ifndef  RECONSTRUCTEDPARTICLE_ANALYZERS_H
#define  RECONSTRUCTEDPARTICLE_ANALYZERS_H

#include "TROOT.h"
#include <cmath>
#include <cmath>
#include <vector>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "Math/Vector4D.h"

#include <iostream>
#include "defines.h"

namespace FCCAnalyses{

namespace ReconstructedParticle{

  // build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
  // technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair 
  
  struct resonanceBuilder_mass_recoil {
    float m_resonance_mass;
    float m_recoil_mass;
    float chi2_recoil_frac;
    float ecm;
    bool m_use_MC_Kinematics;
    resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs, ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> parents, ROOT::VecOps::RVec<int> daugthers) ;
  };
  
    
  inline float deltaR(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() != 2) return -1;
    
    ROOT::Math::PxPyPzEVector tlv1;
    tlv1.SetPxPyPzE(in.at(0).momentum.x, in.at(0).momentum.y, in.at(0).momentum.z, in.at(0).energy);

    ROOT::Math::PxPyPzEVector tlv2;
    tlv2.SetPxPyPzE(in.at(1).momentum.x, in.at(1).momentum.y, in.at(1).momentum.z, in.at(1).energy);
    
    return std::sqrt(std::pow(tlv1.Eta()-tlv2.Eta(), 2) + std::pow(tlv1.Phi()-tlv2.Phi(), 2));
  }
  
  // Missing energy selection as defined to create a filter 
  inline float get_cosTheta_miss(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> met){
    
    float costheta = 0.;
    if(met.size() > 0) {
        
      TLorentzVector lv_met;
      lv_met.SetPxPyPzE(met[0].momentum.x, met[0].momentum.y, met[0].momentum.z, met[0].energy);
      costheta = fabs(std::cos(lv_met.Theta()));

    }
    return costheta;
  }
  
  // calculate the cosine(theta) of the missing energy vector
  inline bool has_forward_photon(float cut, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    
    for (size_t i = 0; i < in.size(); ++i) {
      auto & p = in[i];
      TLorentzVector lv;
      lv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
      if(lv.Theta() > (M_PI-cut) || lv.Theta() < cut) return true;
    }
    return false;
  }
  
  // acoplanarity between two reco particles
  inline float acoplanarity(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() != 2) return -1;

    TLorentzVector p1;
    p1.SetXYZM(in[0].momentum.x, in[0].momentum.y, in[0].momentum.z, in[0].mass);

    TLorentzVector p2;
    p2.SetXYZM(in[1].momentum.x, in[1].momentum.y, in[1].momentum.z, in[1].mass);

    float acop = abs(p1.Phi() - p2.Phi());
    if(acop > M_PI) acop = 2 * M_PI - acop;
    acop = M_PI - acop;

    return acop;
  }

  /*
  resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_recoil_mass = arg_recoil_mass, chi2_recoil_frac = arg_chi2_recoil_frac, ecm = arg_ecm, m_use_MC_Kinematics = arg_use_MC_Kinematics;}

  Vec_rp resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil::operator()(Vec_rp legs, Vec_i recind, Vec_i mcind, Vec_rp reco, Vec_mc mc, Vec_i parents, Vec_i daugthers) {

    Vec_rp result;
    result.reserve(3);
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = legs.size();
  
    if(n > 1) {
      ROOT::VecOps::RVec<bool> v(n);
      std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
      do {
	std::vector<int> pair;
	rp reso;
	reso.charge = 0;
	TLorentzVector reso_lv; 
	for(int i = 0; i < n; ++i) {
	  if(v[i]) {
	    pair.push_back(i);
	    reso.charge += legs[i].charge;
	    TLorentzVector leg_lv;

	    if(m_use_MC_Kinematics) { // MC kinematics
	      int track_index = legs[i].tracks_begin;   // index in the Track array
	      int mc_index = ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
	      if (mc_index >= 0 && mc_index < mc.size()) {
		leg_lv.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y, mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
	      }
	    }
	    else { // reco kinematics
	      leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
	    }

	    reso_lv += leg_lv;
	  }
	}

	if(reso.charge != 0) continue; // neglect non-zero charge pairs
	reso.momentum.x = reso_lv.Px();
	reso.momentum.y = reso_lv.Py();
	reso.momentum.z = reso_lv.Pz();
	reso.mass = reso_lv.M();
	result.emplace_back(reso);
	pairs.push_back(pair);

      } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
      std::cout << "ERROR: resonanceBuilder_mass_recoil, at least two leptons required." << std::endl;
      exit(1);
    }
  
    if(result.size() > 1) {
  
      Vec_rp bestReso;
        
      int idx_min = -1;
      float d_min = 9e9;
      for (int i = 0; i < result.size(); ++i) {
            
	// calculate recoil
	auto recoil_p4 = TLorentzVector(0, 0, 0, ecm);
	TLorentzVector tv1;
	tv1.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
	recoil_p4 -= tv1;
      
	auto recoil_fcc = edm4hep::ReconstructedParticleData();
	recoil_fcc.momentum.x = recoil_p4.Px();
	recoil_fcc.momentum.y = recoil_p4.Py();
	recoil_fcc.momentum.z = recoil_p4.Pz();
	recoil_fcc.mass = recoil_p4.M();
            
	TLorentzVector tg;
	tg.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
        
	float boost = tg.P();
	float mass = std::pow(result.at(i).mass - m_resonance_mass, 2); // mass
	float rec = std::pow(recoil_fcc.mass - m_recoil_mass, 2); // recoil
	float d = (1.0-chi2_recoil_frac)*mass + chi2_recoil_frac*rec;
            
	if(d < d_min) {
	  d_min = d;
	  idx_min = i;
	}

     
      }
      if(idx_min > -1) { 
	bestReso.push_back(result.at(idx_min));
	auto & l1 = legs[pairs[idx_min][0]];
	auto & l2 = legs[pairs[idx_min][1]];
	bestReso.emplace_back(l1);
	bestReso.emplace_back(l2);
      }
      else {
	std::cout << "ERROR: resonanceBuilder_mass_recoil, no mininum found." << std::endl;
	exit(1);
      }
      return bestReso;
    }
    else {
      auto & l1 = legs[0];
      auto & l2 = legs[1];
      result.emplace_back(l1);
      result.emplace_back(l2);
      return result;
    }
  }    
  */


/*  
  struct get_legs_lv_first(ROOT::VecOps::RVec<FCCAnalyses::ReconstructedParticle::ResoWithLegs>& resoWithLegs) {
    return ROOT::VecOps::RVec<TLorentzVector>(resoWithLegs.size(), [&](size_t i) {
	return resoWithLegs[i].legs_lv.first;
      });
  }

  struct get_legs_lv_second(ROOT::VecOps::RVec<FCCAnalyses::ReconstructedParticle::ResoWithLegs>& resoWithLegs) {
    return ROOT::VecOps::RVec<TLorentzVector>(resoWithLegs.size(), [&](size_t i) {
	return resoWithLegs[i].legs_lv.second;
      });
  }

  struct get_deltaR_legs(ROOT::VecOps::RVec<FCCAnalyses::ReconstructedParticle::ResoWithLegs>& resoWithLegs) {
    return ROOT::VecOps::RVec<double>(resoWithLegs.size(), [&](size_t i) {
	const auto& legs = resoWithLegs[i].legs_lv;
	return ROOT::Math::VectorUtil::DeltaR(legs.first, legs.second);
      });
  }
  */

 /*
  auto lambda28 = [&](ROOT::VecOps::RVec<FCCAnalyses::ReconstructedParticle::ResoWithLegs>& resoWithLegs){return ReconstructedParticle::get_legs_lv_first(resoWithLegs);};
  auto lambda29 = [&](ROOT::VecOps::RVec<FCCAnalyses::ReconstructedParticle::ResoWithLegs>& resoWithLegs){return ReconstructedParticle::get_legs_lv_second(resoWithLegs);};

  auto legs_lv1 = fccutil::apply(eval(lambda28), resoWithLegs);
  auto legs_lv2 = fccutil::apply(eval(lambda29), resoWithLegs);
  */

  
  struct resonanceBuilder {
    float m_resonance_mass;
    resonanceBuilder(float arg_resonance_mass);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs);
  };
  
  /*
  struct resonanceBuilder {
    float m_resonance_mass;
    resonanceBuilder(float arg_resonance_mass);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs);
  };
  /// build the resonance from 2 particles from an arbitrary list of input ReconstructedPartilces. Keep the closest to the mass given as input
  //  struct resonanceBuilder {
  // float m_resonance_mass;
  // resonanceBuilder(float arg_resonance_mass);
  // ROOT::VecOps::RVec<ResoWithLegs> resonanceBuilder::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs);
  // };

  /// Legs used for resonance reconstruction 
   struct ResoWithLegs {
   edm4hep::ReconstructedParticleData reso;
   std::pair<TLorentzVector, TLorentzVector> legs_lv;
   };
  */
   /// build the resonance from 2 particles from an arbitrary list of input ReconstructedPartilces. Keep the closest to the mass given as input                                    
   //struct resonanceBuilder {
   // float m_resonance_mass;
   // resonanceBuilder(float arg_resonance_mass);
   // ROOT::VecOps::RVec<ResoWithLegs> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs);
   // };

  /// build the recoil from an arbitrary list of input ReconstructedPartilces and the center of mass energy
  struct recoilBuilder {
    recoilBuilder(float arg_sqrts);
    float m_sqrts = 240.0;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) ;
  };

  /// return the angular separations (min / max / average) between a collection of particles
  struct angular_separationBuilder {
    angular_separationBuilder( int arg_delta); //  0, 1, 2 = max, min, average
    int m_delta = 0;
    float operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) ;
  };

  /// select ReconstructedParticles with transverse momentum greater than a minimum value [GeV]
  struct sel_pt {
    sel_pt(float arg_min_pt);
    float m_min_pt = 1.; //> transverse momentum threshold [GeV]
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  };

  /// select ReconstructedParticles with momentum greater than a minimum value [GeV]
  struct sel_p {
    sel_p(float arg_min_p, float arg_max_p = 1e10);
    float m_min_p = 1.; //> momentum threshold [GeV]
    float m_max_p = 1e10; //< momentum threshold [GeV]
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  };

  /// select ReconstructedParticles with charge equal or in asolute value
  struct sel_charge {
    sel_charge(int arg_charge, bool arg_abs);
    float m_charge; //> charge condition
    bool  m_abs;//> absolute value of the charge
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  };

  /// select a list of reconstructed particles depending on the angle cosTheta axis
  struct sel_axis{
    bool m_pos = 0; //> Which hemisphere to select, false/0=cosTheta<0 true/1=cosTheta>0
    sel_axis(bool arg_pos);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<float> angle, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  };

  /// select a list of reconstructed particles depending on the status of a certain boolean flag
  struct sel_tag {
    bool m_pass; // if pass is true, select tagged jets. Otherwise select anti-tagged ones
    sel_tag(bool arg_pass);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<bool> tags, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  };




  /// return reconstructed particles
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> get(ROOT::VecOps::RVec<int> index, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the transverse momenta of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_pt(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the momenta of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_p(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the momenta of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_px(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the momenta of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_py(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the momenta of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_pz(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the pseudo-rapidity of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_eta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the rapidity of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_y(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the theta of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_theta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the phi of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_phi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the energy of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_e(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the masses of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_mass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the charges of the input ReconstructedParticles
  ROOT::VecOps::RVec<float> get_charge(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the type of the input ReconstructedParticles
  ROOT::VecOps::RVec<int> get_type(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the TlorentzVector of the input ReconstructedParticles
  ROOT::VecOps::RVec<TLorentzVector> get_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// return the TlorentzVector of the indexed input ReconstructedParticles
  TLorentzVector get_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, int index);

  /// return the TlorentzVector of the one input ReconstructedParticle
  TLorentzVector get_tlv(edm4hep::ReconstructedParticleData in);

  /// concatenate both input vectors and return the resulting vector
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> merge(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> y);

  /// remove elements of vector y from vector x
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> remove( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> y);

  /// return the size of the input collection
  int get_n(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// returns the bjet flavour
  ROOT::VecOps::RVec<bool> getJet_btag(ROOT::VecOps::RVec<int> index, ROOT::VecOps::RVec<edm4hep::ParticleIDData> pid, ROOT::VecOps::RVec<float> values);

  /// get number of b-jets
  int getJet_ntags(ROOT::VecOps::RVec<bool> in);

  
}//end NS ReconstructedParticle

}//end NS FCCAnalyses
#endif
