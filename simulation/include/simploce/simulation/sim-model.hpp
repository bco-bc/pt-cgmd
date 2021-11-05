/*
 * File:   sim-model.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 3, 2019, 11:59 AM
 */

#ifndef SIM_MODEL_HPP
#define SIM_MODEL_HPP

#include "simploce/particle/coarse-grained.hpp"
#include "sim-data.hpp"
#include "s-types.hpp"
#include <iostream>

namespace simploce {
    
    /**
     * Specifies a model for use in a simulation of particles.
     * @param P Particle type
     */
    template <typename P>
    class SimulationModel;
    
    /**
     * Specialization for beads.
     */
    template <>
    class SimulationModel<Bead> {
    public:
        
        // Noncopyable
        SimulationModel(const SimulationModel&) = delete;
        SimulationModel& operator = (const SimulationModel&) = delete;                
               
        /**
         * Constructor
         * @param cg Coarse grained particle model.
         * @param displacer Displacer.
         * @param forcefield Coarse grained force field.
         * @param box Simulation box.
         * @param bc Boundary condition.
         */
        SimulationModel(const cg_mod_ptr_t& cg,
                        const cg_displacer_ptr_t& displacer,                        
                        const cg_interactor_ptr_t& interactor,
                        const box_ptr_t& box,
                        const bc_ptr_t& bc);
        
        /**
         * Returns number of particles.
         * @return Number of particles.
         */
        std::size_t 
        size() const;
        
        /**
         * Lets all particles interact.
         * @param param Simulation parameters. 
         * @return Data (e.g. potential energy).
         */
        SimulationData 
        interact(const sim_param_t& param);
        
        /**
         * Lets given bead interact with all other particles.
         * @param bead Bead.
         * @param param Simulation parameters. 
         * @return Data (e.g. potential energy).
         */
        SimulationData 
        interact(const bead_ptr_t& bead, 
                 const sim_param_t& param);
        
        /**
         * Displaces the particles.
         * @param param Simulation parameters. Must provide "npairlists" parameter, the
         * number of steps between updating the particle pair list.
         * @return Simulation data (e.g. kinetic energy, temperature, etc).
         */
        SimulationData 
        displace(const sim_param_t& param);
        
        /**
         * Saves state.
         * @param beads Output beads into which state is saved.
         */
        void
        saveState(std::ostream& beads);
        
        /**
         * Reads state.
         * @param stream Input stream holding states.
         */
        void
        readState(std::istream& stream);
        
        /**
         * Reads simulation model from input stream.
         * @param stream Input stream.
         * @param catalog Particle specification catalog.
         * @return Input stream.
         */
        void 
        readFrom(std::istream& stream, const spec_catalog_ptr_t& catalog);
        
        /**
         * Returns simulation box.
         * @return Box.
         */
        box_ptr_t 
        box() const { return box_; }
        
        /**
         * Returns boundary condition.
         * @return Boundary condition.
         */
        bc_ptr_t 
        boundaryCondition() const { return bc_; }
        
        /**
         * Performs a 'task' with all particles. The given task must expose the 
         * operator 
         * <code>
         * R operator () (const std::vector<bead_ptr_t>& all,
         *                const std::vector<bead_ptr_t>& free,
         *                const std::vector<bead_group_ptr_t>& groups);
         * </code>
         * where 'all' represents all the particles in this model.
         * @param task Task of type TASK. May be a lambda expression.
         * @return Result of type R.
         */
        template <typename R, typename TASK>
        R 
        doWithAllFreeGroups(const TASK task) {
            return cg_->doWithAllFreeGroups<R>(task); 
        }
        
        /**
         * Assigns new displacer.
         * @param displacer Displacer.
         */
        void 
        displacer(const cg_displacer_ptr_t& displacer);
        
        /**
         * Returns displacer.
         * @return Displacer.
         */
        cg_displacer_ptr_t 
        displacer() const;
        
        /**
         * Returns interactor.
         * @return Interactor.
         */
        cg_interactor_ptr_t 
        interactor() const;
        
    private:
        
        friend std::ostream& operator << (std::ostream& stream, const SimulationModel<Bead>& sm);
        
        friend class SimulationModelFactory;
        
        SimulationModel();
                        
        cg_mod_ptr_t cg_;
        cg_displacer_ptr_t displacer_;
        cg_interactor_ptr_t interactor_;
        box_ptr_t box_;
        bc_ptr_t bc_;
        
    };
    
    std::ostream& operator << (std::ostream& stream, const SimulationModel<Bead>& sm);
        
}

#endif /* SIM_MODEL_HPP */

