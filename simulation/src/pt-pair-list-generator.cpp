/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 25, 2019, 1:25 PM
 */

#include "simploce/simulation/pt-pair-list-generator.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/continuous.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/types/cvector_t.hpp"
#include "simploce/conf/s-conf.hpp"

namespace simploce {

    ProtonTransferPairListGenerator::ProtonTransferPairListGenerator(const bc_ptr_t& bc) :
            rcutoffPT_{conf::CUTOFF_DISTANCE_PT}, bc_{bc}
    {        
    }
    
    ProtonTransferPairListGenerator::prot_pair_list_t 
    ProtonTransferPairListGenerator::generate(const prot_cg_mod_ptr_t& cg) const
    {
        using prot_pair_list_t = ProtonTransferPairListGenerator::prot_pair_list_t;

        real_t rmax2 = rcutoffPT_() * rcutoffPT_();
        bc_ptr_t bc = bc_;

        return cg->doWithProtonatableBeads<prot_pair_list_t>([bc, rmax2] (const std::vector<prot_bead_ptr_t>& beads) {
            prot_pair_list_t pairlist{};
            for (auto i = beads.begin(); i != (beads.end() - 1); ++i) {
                auto pi = *i;
                const auto& ri = pi->position();
                for (auto j = i + 1; j != beads.end(); ++j) {
                    auto pj = *j;
                    const auto& rj = pj->position();
                    auto R = bc->apply(ri, rj);
                    real_t R2 = norm_square<real_t>(R);
                    if ( R2 < rmax2 ) {
                        prot_pair_t pair = std::make_pair(pi, pj);
                        pairlist.push_back(pair);
                    }
                }
            }
            return pairlist;
        });
    }
}
