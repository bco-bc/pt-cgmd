/*
 * Author: André H. Juffer.
 * Created on 21/10/2021, 21:09.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef PARTICLES_PROTONATABLE_BEAD_HPP
#define PARTICLES_PROTONATABLE_BEAD_HPP

#include "simploce/particle/p-types.hpp"
#include "bead.hpp"
#include <type_traits>
#include <iostream>

namespace simploce {

    /**
     * A bead capable of binding and releasing protons.
     * @tparam S Protonation state type. Must be <code>CopyAssignable</code> and
     * implement the <code>protonatable</code> interface.
     */
     template <typename S>
     class ProtonatableBead: public Bead {
     public:

         /**
          * Protonatable bead pointer type.
          */
         using prot_bead_ptr_t = std::shared_ptr<ProtonatableBead<S>>;

         /**
           * Creates new protonatable bead. All arguments are required.
           * @param id Unique bead identifier.
           * @param index Bead sequential index.
           * @param name Bead name (does not need to be unique).
           * @param spec Protonatable Bead specification.
           * @return Protonatable bead.
           */
         static prot_bead_ptr_t create(const id_t& id,
                                       std::size_t index,
                                       const std::string& name,
                                       const spec_ptr_t& spec);

         /**
          * Constructor. All arguments are required.
          * @param id Unique bead identifier.
          * @param index Bead sequential index.
          * @param name Bead name (does not need to be unique).
          * @param spec  Protonatable bead specification.
          */
         ProtonatableBead(const id_t& id,
                          std::size_t index,
                          const std::string& name,
                          const spec_ptr_t& spec);

         // Noncopyable.
         ProtonatableBead(const ProtonatableBead<S>& bead) = delete;
         ProtonatableBead<S>& operator = (const ProtonatableBead<S>& bead) = delete;

         ~ProtonatableBead() noexcept override = default;

         /**
          * Returns protonation state.
          * @return Protonation state.
          */
         S protonationState() const;

         /**
          * Sets protonation protonationState.
          * @param protonationState Protonation state.
          */
         void protonationState(const S& protonationState);

         charge_t charge() const override;

         /**
          * Is this bead protonated.
          * @return Result.
          */
         bool isProtonated() const;

         /**
          * Writes this bead to an output stream.
          * @param stream Output stream.
          */
         void write(std::ostream& stream) const override;

         /**
          * Writes bead state (protonation state) to an output stream.
          * @param stream Output stream.
          */
         void writeState(std::ostream& stream) const override;

         /**
          * Reads bead state (protonation state) from an input stream.
          * @param stream Input stream.
          */
         void readState(std::istream& stream) override;

     private:

         S protonationState_;

     };

     template <typename S>
     ProtonatableBead<S>::ProtonatableBead(const id_t& id,
                                           std::size_t index,
                                           const std::string& name,
                                           const spec_ptr_t& spec) :
        Bead(id, index, name, spec), protonationState_{} {
            if ( !std::is_copy_assignable<S>() ) {
                throw std::domain_error("Protonation state type 'S' is not CopyAssignable.");
            }
     }

     template <typename S>
     S
     ProtonatableBead<S>::protonationState() const {
        return protonationState_;
     }

     template <typename S>
     void
     ProtonatableBead<S>::protonationState(const S& protonationState) {
        protonationState_ = protonationState;
     }

     template <typename S>
     charge_t
     ProtonatableBead<S>::charge() const {
        return Bead::charge() + protonationState_.charge();
     }

     template <typename S>
     bool
     ProtonatableBead<S>::isProtonated() const {
         return protonationState_.isProtonated();
     }

     template <typename S>
     void
     ProtonatableBead<S>::write(std::ostream& stream) const {
         Particle::write(stream);
         stream << conf::SPACE << conf::PROTONATABLE << conf::SPACE << protonationState_;
     }

     template <typename S>
     void
     ProtonatableBead<S>::writeState(std::ostream& stream) const {
         Particle::writeState(stream);
         stream << protonationState_;
     }

     template <typename S>
     void
     ProtonatableBead<S>::readState(std::istream& stream) {
         Particle::readState(stream);
         stream >> protonationState_;
     }

     template <typename S>
     typename ProtonatableBead<S>::prot_bead_ptr_t
     ProtonatableBead<S>::create(const id_t& id,
                                 std::size_t index,
                                 const std::string& name,
                                 const spec_ptr_t& spec) {
        return prot_bead_ptr_t(new ProtonatableBead(id, index, name, spec));
    }

    template <typename S>
    std::ostream& operator << (std::ostream& stream, const ProtonatableBead<S>& protonatableBead) {
        protonatableBead.write(stream);
        return stream;
    }
}


#endif //PARTICLES_PROTONATABLE_BEAD_HPP
