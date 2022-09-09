/*
 * Author: André H. Juffer.
 * Created on 21/10/2021, 22:49.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/protonatable.hpp"
#include "simploce/units/units-mu.hpp"
#include <cstdlib>
#include <string>

using namespace simploce;

class Discrete : public protonatable {
    friend std::istream& operator >> (std::istream& stream, Discrete& discrete);
public:
    Discrete() : v_{0} {}

    void protonate() override { v_ = 1; }
    void deprotonate() override { v_ = 0; }

    int operator()() const {
        return v_;
    }

    charge_t charge() const override{
        return v_ * units::mu<real_t>::PROTON_CHARGE;
    }

    mass_t mass() const override {
        return v_ * units::mu<real_t>::PROTON_MASS;
    }

    bool isProtonated() const override {
        return v_ == 1;
    }

private:
    int v_;
};

std::ostream& operator << (std::ostream& stream, const Discrete& discrete) {
    stream << discrete();
    return stream;
}

std::istream& operator >> (std::istream& stream, Discrete& discrete) {
    stream >> discrete.v_;
    return stream;
}

class Continuous : public protonatable {
    friend std::istream& operator >> (std::istream& stream, Continuous& continuous);
public:
    Continuous() : x_{0.0}, I_{0.0} {}

    charge_t charge() const override {
        return x_ * units::mu<real_t>::PROTON_CHARGE;
    }

    mass_t mass() const override {
        return x_ * units::mu<real_t>::PROTON_MASS;
    }

    real_t x() const { return x_; }
    real_t I() const { return I_; }

    void protonate(const std::map<std::string, real_t>& values) override {
        static std::string x("x");
        static std::string I("I");
        auto iter = values.find(x);
        x_ = iter->second;
        iter = values.find(I);
        I_ = iter->second;
    }

    bool isProtonated() const override {
        return x_ > 0.0;
    }

private:
    real_t x_;
    real_t I_;
};

std::ostream& operator << (std::ostream& stream, const Continuous& continuous) {
    stream << continuous.x() << conf::SPACE << continuous.I();
    return stream;
}

std::istream& operator >> (std::istream& stream, Continuous& continuous) {
    stream >> continuous.x_ >> continuous.I_;
    return stream;
}

void test() {
    std::cout << "Discrete protonation state" << std::endl;

    using cg_t = ProtonatableCoarseGrained<Discrete>;

    auto spec = ParticleSpec::create("AP",
                                     1.0,
                                     1.0,
                                     1.0,
                                     7.4,
                                     false,
                                     "dp");
    Discrete discrete;
    cg_t cg;
    cg.box(factory::box(0.0));
    auto protonatableBead = cg.addProtonatableBead("OAP", spec);
    protonatableBead->protonationState(discrete);
    std::cout << "Protonatable Bead: " << std::endl << *protonatableBead << std::endl;
    std::cout << "Protonatable Bead charge: " << protonatableBead->charge() << std::endl;
    discrete.protonate();
    protonatableBead->protonationState(discrete);
    std::cout << "Protonatable Bead: " << std::endl << *protonatableBead << std::endl;
    std::cout << "Protonatable Bead charge: " << protonatableBead->charge() << std::endl;
    discrete.deprotonate();
    protonatableBead->protonationState(discrete);
    std::cout << "Protonatable Bead: " << std::endl << *protonatableBead << std::endl;
    std::cout << "Protonatable Bead charge: " << protonatableBead->charge() << std::endl;

    auto bead = cg.addBead("AP", spec);
    std::cout << "None-protonatable bead:" << std::endl << *bead << std::endl;

    std::cout << "Number of beads: " << cg.numberOfBeads() << std::endl;
    std::cout << "CoarseGrained:" << std::endl << cg << std::endl;

    std::cout << std::endl;
}

void test2() {
    std::cout << "Continuous protonation state:" << std::endl;

    using prot_cg_t = ProtonatableCoarseGrained<Continuous>;

    prot_cg_t cg;
    cg.box(factory::box(0.0));
    auto spec = ParticleSpec::create("AP", 1.0, 1.0, 1.0, 7.0, false, "dp");

    std::map<std::string, real_t> state{};

    Continuous continuous;
    auto protonatableBead = cg.addProtonatableBead("OAP", spec);
    protonatableBead->protonationState(continuous);
    std::cout << "Protonatable Bead: " << std::endl << *protonatableBead << std::endl;
    std::cout << "Protonatable Bead charge: " << protonatableBead->charge() << std::endl;
    std::cout << "Protonatable Bead isProtonated: " << protonatableBead->isProtonated() << std::endl;
    std::cout << std::endl;

    auto x = std::make_pair("x", 0.1);
    state.insert(x);
    auto I = std::make_pair("I", -1.01234);
    state.insert(I);
    continuous.protonate(state);
    protonatableBead->protonationState(continuous);
    std::cout << "Protonatable Bead: " << std::endl << *protonatableBead << std::endl;
    std::cout << "Protonatable Bead charge: " << protonatableBead->charge() << std::endl;

    std::cout << std::endl;
}


int main() {
    test();
    test2();
    return (EXIT_SUCCESS);
}

