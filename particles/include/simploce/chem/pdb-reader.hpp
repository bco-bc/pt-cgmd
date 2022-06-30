/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 31, 2022, 14:33 PM
 */

#ifndef PARTICLES_PDB_READER_HPP
#define PARTICLES_PDB_READER_HPP

#include "chem-reader.hpp"

namespace simploce {


    /**
     * Parses PDB (protein) structures. Residues are reported to the content_handler as atom groups.
     * The PDB structure should at least contain a TITLE and/or a HEADER record and also one or more ATOM records.
     */
    class PDBReader : public chem_reader {

    public:

      PDBReader();

      ~PDBReader() override;

      void parse(const cont_handler_ptr_t& handler,
                 const input_source_ptr_t& source) override;

    };

}

#endif //PARTICLES_PDB_READER_HPP
