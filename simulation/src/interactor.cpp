/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 6, 2019, 12:49 PM
 */

#include "simploce/simulation/interactor.hpp"
#include "simploce/simulation/pair-list-generator.hpp"
// #include "simploce/simulation/at-forcefield.hpp"
// #include "simploce/simulation/cg-forcefield.hpp"
#include "simploce/particle/atomistic.hpp"
// #include "simploce/particle/coarse-grained.hpp"
#include <memory>
#include <utility>

namespace simploce {
    
    /**
     * Result of energy/force calculation.
     */
    //using tdm_result_t = std::pair<energy_t, energy_t>;
        
    // Atom pair type
    //using atom_pair_t = std::pair<atom_ptr_t, atom_ptr_t>;
        
    //Bead pair lists type (single list).
    //using atom_pair_list_t = std::vector<atom_pair_t>;
        
    
    // Bead pair type
    //using bead_pair_t = std::pair<bead_ptr_t, bead_ptr_t>;
        
    //Bead pair lists type (single list).
    //using bead_pair_list_t = std::vector<bead_pair_t>;
    
    // Bead pair lists.
    //static std::vector<atom_pair_list_t> atomPairLists_{};

    // Bead pair lists.
    //static std::vector<bead_pair_list_t> beadPairLists_{};

    template <typename P>
    static PairLists<P> updatePairList_(const std::shared_ptr<pair_lists_generator<P>> &pairListGenerator,
                                        const std::shared_ptr<ParticleSystem<P, ParticleGroup<P>>> &particleSystem) {
        using p_ptr_t = typename Interactor<P>::p_ptr_t;
        using pg_ptr_t = typename Interactor<P>::pg_ptr_t;

        return std::move(particleSystem->template doWithAllFreeGroups<PairLists<P>>([pairListGenerator] (
                const std::vector<p_ptr_t> &all,
                const std::vector<p_ptr_t> &free,
                const std::vector<pg_ptr_t> &groups) {
            return pairListGenerator->generate(all, free, groups);
        }));
    }


    PairLists<Atom>
    updatePairLists(const std::shared_ptr<pair_lists_generator<Atom>> &pairListGenerator,
                    const std::shared_ptr<ParticleSystem<Atom, ParticleGroup<Atom>>> &atomistic) {
        return std::move(updatePairList_<Atom>(pairListGenerator, atomistic));
    }

    PairLists<Bead>
    updatePairLists(const std::shared_ptr<pair_lists_generator<Bead>> &pairListGenerator,
                    const std::shared_ptr<ParticleSystem<Bead, ParticleGroup<Bead>>> &coarseGrained) {
        return std::move(updatePairList_<Bead>(pairListGenerator, coarseGrained));
    }


    
    /*
    Interactor<Atom>::Interactor(const at_ff_ptr_t& forceField,
                                 const atom_pair_lists_gen_ptr_t& pairListGenerator) :
            forceField_{forceField}, pairListGenerator_{pairListGenerator}
    {        
    }
        
    std::pair<energy_t, energy_t> 
    Interactor<Atom>::interact(const sim_param_t& param, 
                               const at_mod_ptr_t& at)
    {
        static bool setup = false;
        static std::size_t npairlist = 0;
        static std::size_t counter = 0;
        
        if ( !setup ) {
            npairlist = param.get<std::size_t>("npairlists");
            setup = true;
        }
        if ( counter % npairlist == 0 || counter == 0) {
            this->updatePairLists_(at);
        } else {
            pairLists_.modified_(false);
        }
        
        tdm_result_t result = 
            at->doWithAllFreeGroups<tdm_result_t>([this] (const std::vector<atom_ptr_t>& all,
                                                      const std::vector<atom_ptr_t>& free,
                                                      const std::vector<atom_group_ptr_t>& groups) {
            for (auto p: all) {
                p->resetForce();
            }
            return this->forceField_->interact(all, free, groups, atomPairLists_);
        });
        
        counter += 1;
        return result;
    }
    
    std::pair<energy_t, energy_t> 
    Interactor<Atom>::interact(const atom_ptr_t& atom,
                               const sim_param_t& param, 
                               const at_mod_ptr_t& at)
    {
        throw std::domain_error(
            "Atomistic force field: No interaction energy calculator available."
        );
    }
    
    void 
    Interactor<Atom>::updatePairLists_(const at_mod_ptr_t& at)
    {
        pairLists_ = 
            at->doWithAllFreeGroups<PairLists<Atom>>([this] (const std::vector<atom_ptr_t>& all,
                                                             const std::vector<atom_ptr_t>& free,
                                                             const std::vector<atom_group_ptr_t>& groups) {
                return this->pairListGenerator_->generate(all, free, groups);
            });
        pairLists_.modified_(true);
    }
    
    std::string 
    Interactor<Atom>::id() const
    {
        return forceField_->id();
    }

    
    Interactor<Bead>::Interactor(const cg_ff_ptr_t& forcefield,
                                 const bead_pair_lists_gen_ptr_t& pairListGenerator) :
            forceField_{forcefield}, pairListGenerator_{pairListGenerator}
    {      
    }
        
    std::pair<energy_t, energy_t>
    Interactor<Bead>::interact(const sim_param_t& param, 
                               const cg_mod_ptr_t& cg)
    {
        static bool setup = false;
        static std::size_t npairlist = 0;
        static std::size_t counter = 0;
        
        if ( !setup ) {
            npairlist = param.get<std::size_t>("npairlists");
            setup = true;
        }
        if ( counter % npairlist == 0 || counter == 0) {
            this->updatePairLists_(cg);
        } else {
            pairLists_.modified_(false);
        }
        
        tdm_result_t result = 
            cg->doWithAllFreeGroups<tdm_result_t>([this] (const std::vector<bead_ptr_t>& all,
                                                      const std::vector<bead_ptr_t>& free,
                                                      const std::vector<bead_group_ptr_t>& groups) {
                for (auto p: all) {
                    p->resetForce();
                }
                return this->forceField_->interact(all, free, groups, pairLists_);
            });
        
        counter += 1;
        return result;
    }
    
    std::pair<energy_t, energy_t>
    Interactor<Bead>::interact(const bead_ptr_t& bead,
                               const sim_param_t& param, 
                               const cg_mod_ptr_t& cg)
    {
        return cg->doWithAllFreeGroups<tdm_result_t>([this, bead] (const std::vector<bead_ptr_t>& all,
                                                               const std::vector<bead_ptr_t>& free,
                                                               const std::vector<bead_group_ptr_t>& groups) {
            return this->forceField_->interact(bead, all, free, groups);
        });
    }
    
    void
    Interactor<Bead>::updatePairLists_(const cg_mod_ptr_t& cg)
    {
        pairLists_ = 
            cg->doWithAllFreeGroups<PairLists<Bead>>([this] (const std::vector<bead_ptr_t>& all,
                                                             const std::vector<bead_ptr_t>& free,
                                                             const std::vector<bead_group_ptr_t>& groups) {
                return this->pairListGenerator_->generate(all, free, groups);
            });
        pairLists_.modified_(true);
    }
    
    std::string
    Interactor<Bead>::id() const
    {
        return forceField_->id();
    }
*/
}