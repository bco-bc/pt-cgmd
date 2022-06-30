/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 31, 2022, 14:44 PM
 */

#include "simploce/chem/input-source.hpp"
#include "simploce/util/file.hpp"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iterator>

namespace simploce {

    static void populate_(std::vector<std::string> &source, std::istream &stream) {
        std::istream_iterator<std::string> begin(stream), end;
        std::copy(begin, end, std::back_inserter(source));
    }

    InputSource::InputSource(const std::string &fileName) :
        sourceId_{fileName}, source_{} {
        std::ifstream stream;
        util::open_input_file(stream, fileName);
        populate_(source_, stream);
        stream.close();
    }

    std::vector<std::string>
    InputSource::content() const {
        return source_;
    }

    std::string
    InputSource::sourceId() const {
        return sourceId_;
    }
}