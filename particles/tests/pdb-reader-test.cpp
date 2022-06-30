//
// Created by ajuffer on 6/1/22.
//

#include "simploce/chem/pdb-reader.hpp"
#include "simploce/chem/input-source.hpp"
#include "simploce/chem/base-content-handler.hpp"
#include "simploce/util/logger.hpp"
#include <iostream>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(simploce::util::Logger::LOGINFO);

    PDBReader reader;
    auto source = std::make_shared<InputSource>("/wrk3/tests/1CUS.pdb");
    cont_handler_ptr_t handler{new BaseContentHandler};

    reader.parse(handler, source);


}