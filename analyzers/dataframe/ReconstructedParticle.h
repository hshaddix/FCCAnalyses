
#ifndef  RECONSTRUCTEDPARTICLE_ANALYZERS_H
#define  RECONSTRUCTEDPARTICLE_ANALYZERS_H

#include <cmath>
#include <vector>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"


/// select ReconstructedParticles with transverse momentum greater than a minimum value [GeV]
struct selRP_pT {
  selRP_pT(float arg_min_pt);
  float m_min_pt = 20; //> transverse momentum threshold [GeV]
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
};

/// return the transverse momenta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_pt(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the momenta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_p(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the momenta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_px(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
std::vector<float> getRP_px_std(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the momenta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_py(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
std::vector<float> getRP_py_std(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the momenta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_pz(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
std::vector<float> getRP_pz_std(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the pseudo-rapidity of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_eta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the rapidity of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_y(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the theta of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_theta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the phi of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_phi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the energy of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_e(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// return the masses of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_mass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in); 
std::vector<float> getRP_mass_std(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in); 

/// return the charges of the input ReconstructedParticles
ROOT::VecOps::RVec<float> getRP_charge(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in); 
std::vector<float> getRP_charge_std(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in); 

/// return the TlorentzVector of the input ReconstructedParticles
ROOT::VecOps::RVec<TLorentzVector> getRP_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
std::vector<TLorentzVector> getRP_tlv_std(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

/// concatenate both input vectors and return the resulting vector
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> mergeParticles(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> y);

/// return the size of the input collection
int getRP_n(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

#endif
