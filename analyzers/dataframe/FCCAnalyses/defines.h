#ifndef FCCANALYZER_DEFINES_H
#define FCCANALYZER_DEFINES_H

#include <cmath>
#include <vector>
#include <math.h>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "ReconstructedParticle2MC.h"




namespace FCCAnalyses {
    
using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_d = ROOT::VecOps::RVec<double>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;

using rp = edm4hep::ReconstructedParticleData;
using Vec_rp = ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>;
using Vec_mc = ROOT::VecOps::RVec<edm4hep::MCParticleData>;
using Vec_tlv = ROOT::VecOps::RVec<TLorentzVector>;
    
 
inline float sumScalar(Vec_f in){
    
    float tot = std::accumulate(in.begin(), in.end(), 0);
    return tot;
}

inline Vec_i indices_(const int& size, const int& start = 0) {
    Vec_i res(size, 0);
    std::iota(std::begin(res), std::end(res), start);
    return res;
}

    
}

#endif
