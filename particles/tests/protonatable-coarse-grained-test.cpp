//
// Created by André H Juffer, Biocenter Oulu, University of Oulu, Finland
//
// Created on 10/25/21.
//

#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/particle/protonatable.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

class Continuous : public protonatable {
    friend std::istream& operator >> (std::istream& stream, Continuous& continuous);
public:
    Continuous() : protonatable{}, x_{}, I_{} {}

    void protonate(const std::map<std::string, real_t>& values) override {
        auto p = values.find("x");
        x_ = p->second;
        p = values.find("I");
        I_ = p->second;
    }

    bool isProtonated() const override {
        return x_ > 0.5;
    }

    charge_t charge() const override {
        return x_ * units::mu<real_t>::PROTON_CHARGE;
    }

    mass_t mass() const override {
        return x_ * units::mu<real_t>::PROTON_MASS;
    }

    real_t x() const { return x_; }
    void x(real_t x) { x_ = x; }
    real_t I() const { return I_; }
    void I(real_t I) { I_ = I; }

private:
    real_t x_;
    real_t I_;
};

std::ostream& operator << (std::ostream& stream, const Continuous& continuous) {
    stream << std::setw(conf::REAL_WIDTH) << continuous.x();
    stream << std::setw(conf::REAL_WIDTH) << continuous.I();
    return stream;
}

std::istream& operator >> (std::istream& stream, Continuous& continuous) {
    stream >> continuous.x_ >> continuous.I_;
    return stream;
}

void test1(const spec_catalog_ptr_t& catalog) {
    ProtonatableCoarseGrained<Continuous> cg;

    auto spec = catalog->lookup("HCOOH");
    auto bead = cg.addProtonatableBead("HCOOH", spec);

    std::cout << "Protonatable CG model: " << std::endl;
    std::cout << cg << std::endl;

    Continuous continuous;
    continuous.x(0.4);
    continuous.I(0.1);
    bead->protonationState(continuous);
    std::cout << "Protonatable CG model: " << std::endl;
    std::cout << cg << std::endl;
}

void test2(const spec_catalog_ptr_t& catalog) {
    std::cout << "Reading protonatable CG particle system:" << std::endl;
    std::ifstream input;
    std::string fileName = "/localdisk/resources/protonatable-coarse-grained-system.dat";
    util::open_input_file(input, fileName);
    auto cg = ProtonatableCoarseGrained<Continuous>::obtainFrom(input, catalog);
    std::cout << "Read protonatable CG model: " << std::endl;
    std::cout << *cg << std::endl;

    using prot_bead_ptr_t = ProtonatableCoarseGrained<Continuous>::prot_bead_ptr_t;
    cg->doWithProtonatableBeads<void>([] (const std::vector<prot_bead_ptr_t>& beads) {
        for (const auto& bead : beads) {
            std::cout << "Bead id: " << bead->id();
            std::cout << ", protonation state: " << bead->protonationState() << std::endl;
        }
    });
    std::cout << "Total charge: " << cg->charge() << std::endl;
    std::cout << "Protonation state: " << cg->protonationState() << std::endl;
    auto Na = catalog->lookup("Na+");
    cg->replaceGroupsByParticles(Na, 1);
    std::cout << "After Na+ insertion: " << std::endl << *cg << std::endl;
}

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fileName = "/localdisk/resources/particles-specs.dat";
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fileName);
    std::cout << *catalog << std::endl << std::endl;

    test1(catalog);
    test2(catalog);
    return (EXIT_SUCCESS);
}

