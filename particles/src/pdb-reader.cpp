/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 31, 2022, 14:33 PM
 */

#include "simploce/chem/pdb-reader.hpp"
#include "simploce/chem/input-source.hpp"
#include "simploce/chem/content-handler.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include "simploce/types/u-types.hpp"
#include "simploce/conf/p-conf.hpp"
#include "simploce/units/units-mu.hpp"
#include <algorithm>

namespace simploce {

    class PDBReaderHelper {

        using iter_t = std::vector<std::string>::const_iterator;

        friend class PDBReader;

        explicit PDBReaderHelper(util::Logger& logger);

        std::string findTitle(const std::vector<std::string>& content);
        bool newMolecule(const iter_t& current,
                         const std::vector<std::string>& content);

        void parseMolecule(const iter_t& current,
                           const std::vector<std::string>& content,
                           const cont_handler_ptr_t& handler);

        std::string compose(iter_t& iter, const std::vector<std::string> &content);

        void parseResidue(const iter_t& current,
                          const std::vector<std::string>& content,
                          const cont_handler_ptr_t& handler);

        void parseAtom(const iter_t& current,
                       const std::vector<std::string>& content,
                       const cont_handler_ptr_t& handler);

        util::Logger logger_;
        std::string chainId_;
        std::string previousChainId_;
        std::string prevResidueName_;
        std::size_t prevResidueIndex_;

        std::string TITLE{"TITLE"};
        std::string HEADER{"HEADER"};
        std::string COMPND{"COMPND"};
        std::string ATOM{"ATOM"};
        std::string HETATM{"HETATM"};
        std::string MODEL{"MODEL"};
        std::vector<std::string> RECORD_NAMES = {
                "HEADER", "TITLE", "OBSLTE", "CAVEAT", "COMPND", "SOURCE", "KEYWDS",
                "EXPDTA","AUTHOR", "REVDAT", "SPRSDE", "JRNL", "REMARK",
                "DBREF", "SEQADV", "SEQRES", "MODRES", "HET", "HETNAM", "HETSYN",
                "FORMUL","HELIX", "SHEET", "TURN","SSBOND", "LINK",
                "HYDBND", "SLTBRG", "CISPEP","SITE", "CRYST1", "CRYST2", "CRYST3",
                "ORIGXn", "ORIGX1", "ORIGX2", "ORIGX3", "SCALEn", "SCALE1", "SCALE2", "SCALE3",
                "MTRIXn", "TVECT",
                "MODEL", "ATOM", "SIGATM", "ANISOU", "SIGUIJ", "TER",
                "HETATM", "CONECT", "MASTER", "END"
        };
        std::vector<char> CHAIN_IDS = {'A', 'B', 'C', 'D', 'E',
                                       'F', 'G', 'H', 'I', 'J',
                                       'K', 'L', 'M', 'N', 'O',
                                       'P', 'Q', 'R', 'S', 'T',
                                       'U', 'V', 'X', 'Y', 'Z'};

    };


    PDBReader::PDBReader() = default;

    PDBReader::~PDBReader() = default;

    void
    PDBReader::parse(const cont_handler_ptr_t &handler,
                     const input_source_ptr_t &source) {
        util::Logger logger("simploce::PDBReader::parseIt()");
        PDBReaderHelper helper(logger);

        logger.trace(source->sourceId() + ": Parsing PDB content.");

        auto content = source->content();
        if (content.empty()) {
            util::logAndThrow(logger, source->sourceId() + ": No content in input source.");
        }
        if (!handler) {
            util::logAndThrow(logger, source->sourceId() + ": No content handler provided.");
        }

        auto title = helper.findTitle(content);
        handler->start(title);

        std::size_t lineCounter = 0;
        bool atoms{false};
        // Pointer to word in source input currently being processed.
        auto current = content.begin();
        std::vector<std::string> targets = { helper.ATOM, helper.MODEL, helper.HETATM };
        do {
            auto word = *current;
            logger.trace(word + ": Current word in input source.");

            // New record?
            auto first = helper.RECORD_NAMES.begin();
            auto last = helper.RECORD_NAMES.end();
            if ( std::find(first, last, word) != helper.RECORD_NAMES.end() )
                lineCounter += 1;

            if (std::find(targets.begin(), targets.end(), word) != targets.end()) {
                if (helper.newMolecule(current, content)) {
                    helper.parseMolecule(current, content, handler);
                }
                if (word == helper.ATOM || word == helper.HETATM) {
                    atoms = true;
                    helper.parseResidue(current, content, handler);
                    helper.parseAtom(current, content, handler);
                    // To next record.
                    auto iter = std::find(targets.begin(), targets.end(), *current);
                    while ( iter == targets.end() && current != content.end() ) {
                        ++current;
                        iter = std::find(targets.begin(), targets.end(), *current);
                    }
                }
            }
            ++current;
        } while (current != content.end());
        if ( !atoms ) {
            std::string msg = source->sourceId() + ": No '" + helper.ATOM + "' records were found";
            throw std::invalid_argument(msg);
        }

        handler->endAtomGroup();
        handler->endMolecule();
        handler->end();

        logger.debug("Estimated number of lines parsed: " + std::to_string(lineCounter));
        logger.trace("Parsing PDB completed.");
    }


    PDBReaderHelper::PDBReaderHelper(util::Logger& logger) :
        logger_{logger}, chainId_{}, previousChainId_{},  prevResidueName_{}, prevResidueIndex_{} {
    }

    std::string
    PDBReaderHelper::findTitle(const std::vector<std::string>& content) {
        logger_.trace("Looking for a title in input source");
        static std::string NA{"No title available."};
        // Find the HEADER record.
        auto iter = std::find(content.begin(), content.end(), HEADER);
        if (iter == content.end()) {
            // No HEADER, find a TITLE record.
            iter = std::find(content.begin(), content.end(), TITLE);
            if (iter == content.end()) {
                // No TITLE, find a COMPND record.
                iter = std::find(content.begin(), content.end(), COMPND);
                if (iter == content.end()) {
                    // No COMPND.
                    return NA;
                }
            }
        }

        // Found some record.
        ++iter;
        logger_.trace("Current word: " + *iter);
        if (*iter == ATOM || *iter == HETATM) {
            return "No title available";
        }

        // Construct title. Use the full record.
        auto title = compose(iter, content);
        logger_.trace(title + ": Found a title.");
        return title;
    }

    bool
    PDBReaderHelper::newMolecule(const iter_t &current, const std::vector<std::string> &content) {
        logger_.trace(*current + ": New molecule encountered?.");

        static std::vector<char> DIGITS = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
        static bool firstTime{true};

        if ( current + 1 >= content.end() || current + 5 >= content.end() ) {
            logger_.trace("No new molecule. Encountered end of input source.");
            //firstTime = false;
            return false;
        }

        // Check if there is new MODEL.
        if (*current == MODEL) {
            logger_.trace(*current + ": Encountered new " + MODEL + " record.");
            try {
                // Next word must be a positive integer.
                auto sn = boost::lexical_cast<int, std::string>(*(current + 1));
                if ( sn > 0) {
                    logger_.trace("Encountered new molecule.");
                    return true;
                } else {
                    logger_.trace("No new molecule.");
                    return false;
                }
            } catch (boost::bad_lexical_cast& exception) {
                // Not a serial number.
                logger_.trace("No new molecule.");
                return false;
            }
        }

        std::string chainId = *(current + 5);
        logger_.trace(chainId + ": Validating candidate Chain ID.");
        logger_.debug(chainId + ": Candidate Chain ID.");

        // Chain ID should be a single character.
        if (chainId.size() > 1) {
            logger_.debug(chainId + ": Not a Chain ID.");
            chainId = *(current + 4);
            if ( chainId.size() > 1) {
                logger_.debug(chainId + ": Not a Chain ID.");
                return false;
            }
        }

        auto ch = chainId[0];
        if (std::find(CHAIN_IDS.begin(), CHAIN_IDS.end(), ch) == CHAIN_IDS.end()) {
            // Not a chain id.
            return false;
        }

        // Found a chain ID.
        logger_.debug(chainId + ": Identified Chain ID.");
        if (!previousChainId_.empty()) {
            if (previousChainId_ != chainId) {
                // Update chainIds
                previousChainId_ = chainId;
                chainId_ = chainId;
                firstTime = false;
                logger_.trace("New molecule: Change of Chain ID.");
                return true;
            }
            previousChainId_ = chainId;
            chainId_ = chainId;
            return false;
        } else {
            previousChainId_ = chainId;
            chainId_ = chainId;
            if (firstTime) {
                logger_.trace("New molecule: First molecule in input source encountered.");
                firstTime = false;
                return true;
            } else {
                return false;
            }
        }
    }

    void
    PDBReaderHelper::parseMolecule(const iter_t &current,
                                   const std::vector<std::string> &content,
                                   const cont_handler_ptr_t &handler) {
        static bool firstTime = true;
        static std::size_t index = 0;

        logger_.trace("Parsing molecule.");

        if ( current + 1 >= content.end() ) {
            logger_.trace("Reached end of input source.");
            return;
        }

        if ( *current == MODEL ) {
            try {
                if ( !firstTime )
                    handler->endMolecule();
                auto serial = boost::lexical_cast<int, std::string>(*(current + 1));  // May refer to MODEL
                                                                                          // index number.
                if (serial < 0) {
                    logger_.debug(*(current + 1) + ": Not a suitable MODEL index number.");
                    auto moleculeName = MODEL + util::to_string(index);
                    handler->startMolecule(moleculeName);
                    index += 1;
                } else {
                    logger_.debug(*(current + 1) + ": Suitable MODEL index number.");
                    std::string moleculeName = *current + *(current + 1);
                    handler->startMolecule(moleculeName);
                }
            } catch (boost::bad_lexical_cast& exception) {
                logger_.debug(*(current + 1) + ": Not a suitable MODEL index number.");
                auto moleculeName = MODEL + util::to_string(index);
                handler->startMolecule(moleculeName);
                index += 1;
            }
        } else {
            auto moleculeName = chainId_;
            if ( !firstTime )
                handler->endMolecule();
            handler->startMolecule(moleculeName);
        }
        firstTime = false;
        logger_.trace("Parsing molecule complete.");
    }

    void
    PDBReaderHelper::parseResidue(const iter_t &current,
                                  const std::vector<std::string> &content,
                                  const cont_handler_ptr_t &handler) {
        logger_.trace(*current + ": Parsing residue.");

        if ( current + 4 >= content.end() ) {
            logger_.trace("Reached end of input source.");
            return;
        }

        auto candidate = current + 4;
        auto residueName = *candidate;
        logger_.trace(residueName + ": Candidate residue name.");
        if (residueName.size() < 3 || residueName.size() > 6) {  // Avoid short and long residue names.
            logger_.trace(residueName + ": Not a residue name.");
            candidate = current + 3;
            residueName = *candidate;
            logger_.trace(residueName + ": Candidate residue name.");
            if (residueName.size() < 3 || residueName.size() > 6) {
                logger_.trace(residueName + ": Not a residue name.");
                logger_.trace("Done parsing residue.");
                return;
            }
        }

        // May have found a residue name. Then there should a residue sequence number
        int residueIndex;
        try {
            // Assume there is no chain ID.
            residueIndex = boost::lexical_cast<int, std::string>(*(candidate + 1));
        } catch (boost::bad_lexical_cast& exception) {
            // There was possibly a chain ID.
            try {
                residueIndex = boost::lexical_cast<int, std::string>(*(candidate + 2));
            } catch (boost::bad_lexical_cast& exception) {
                logger_.trace(residueName + ": Not a residue name.");
                logger_.trace("Done parsing residue.");
                return;
            }
        }

        logger_.debug(residueName+ ": Accepted Residue Name.");
        logger_.debug(util::to_string(residueIndex) + ": Accepted Residue Index.");

        if (!prevResidueName_.empty() > 0) {
            if (residueName != prevResidueName_ || residueIndex != prevResidueIndex_) {
                prevResidueName_ = residueName;
                prevResidueIndex_ = residueIndex;
                handler->endAtomGroup();
                handler->startAtomGroup(residueName);
                handler->index(residueIndex);
            }
        } else {
            prevResidueName_ = residueName;
            prevResidueIndex_ = residueIndex;
            handler->startAtomGroup(residueName);
            handler->index(residueIndex);
        }

        logger_.trace("Done parsing residue.");
    }

    void
    PDBReaderHelper::parseAtom(const iter_t &current,
                               const std::vector<std::string> &content,
                               const cont_handler_ptr_t &handler) {
        logger_.trace(*current + ": Parsing atom.");

        // End of input source?
        if ( current + 10 >= content.end() ) {
            logger_.trace("Reached end of input source.");
            logger_.trace("Done parsing atom");
            return;
        }

        std::vector<std::string> targets = { ATOM, HETATM };
        if ( std::find(targets.begin(), targets.end(), *current) == targets.end() ) {
            logger_.trace("Done parsing atom");
            return;
        }

        try {
            logger_.trace(*(current + 1) + ": Candidate atom index.");
            logger_.trace(*(current + 2) + ": Candidate atom name.");
            auto atomIndex = boost::lexical_cast<int, std::string>(*(current + 1));
            auto atomName = *(current + 2);
            handler->startAtom(atomName);
            handler->index(atomIndex);
        } catch (boost::bad_lexical_cast& exception) {
            // No ATOM or HETATM
            logger_.trace("Done parsing atom");
            return;
        }

        // Determine x-coordinate.
        const auto small = std::numeric_limits<real_t>::epsilon();
        int index = 4; // Assume chain identifier, code for insertion of residues and alternate location indicator are not there.
        real_t x, y, z;
        bool found = false;
        do {
            try {
                auto coordinate = *(current + index);
                auto length = coordinate.length();
                logger_.trace(coordinate + ": Validating x-coordinate.");
                logger_.debug(coordinate + ": Candidate x-coordinate.");
                x = boost::lexical_cast<real_t, std::string>(coordinate);
                real_t fraction = std::modf(x, &fraction);
                if (std::fabs(fraction) > small)
                    found = true;
                else {
                    if (length >= 5) {
                        std::string t = coordinate.substr(length - 3, 4);
                        if (t == "000")
                            found = true;
                    }
                }
            }
            catch (boost::bad_lexical_cast&) {
                // Deliberately ignored!
            }
            ++index;
        } while (!found && index < 8);
        if (!found)
            util::logAndThrow(logger_, "'" + *(current + index) + "': Cannot find x-coordinate in ATOM record");
        logger_.debug(util::to_string(x) + ": Identified x-coordinate");

        // Determine y- and z-coordinate.
        index -= 1;
        try {
            y = boost::lexical_cast<real_t, std::string>(*(current + index + 1));
            z = boost::lexical_cast<real_t, std::string>(*(current + index + 2));
        } catch (boost::bad_lexical_cast &exception) {
            util::logAndThrow(logger_, "'" + *(current + index) + "': Cannot find y- and/or z-coordinate in ATOM record.");
        }
        x *= units::mu<real_t>::Angstrom_to_nm;
        y *= units::mu<real_t>::Angstrom_to_nm;
        z *= units::mu<real_t>::Angstrom_to_nm;
        handler->atomCoordinates(position_t{x, y, z});
        handler->endAtom();

        logger_.trace("Done parsing atom");
    }

    std::string
    PDBReaderHelper::compose(iter_t &iter, const std::vector<std::string> &content) {
        bool firstTime = true;
        std::string composition;
        bool nextRecord;
        do {
            if ( firstTime )
                composition += *iter;
            else {
                composition += conf::SPACE;
                composition += *iter;
            }
            firstTime = false;
            ++iter;
            auto first = RECORD_NAMES.begin();
            auto last = RECORD_NAMES.end();
            nextRecord = !(std::find(first, last, *iter) == RECORD_NAMES.end() && iter < content.end());
        } while (!nextRecord);
        return composition;
    }

}