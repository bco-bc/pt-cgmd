/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 31, 2022, 14:44 PM
 */

#ifndef PARTICLES_INPUT_SOURCE_HPP
#define PARTICLES_INPUT_SOURCE_HPP

#include <vector>
#include <string>

namespace simploce {

    /**
     * Interface for a single input source for molecular content.
     */
    class InputSource {
    public:

        explicit InputSource(const std::string& fileName);

        /**
         * Get the content as a collection of string items. Each word in the original content is placed in a
         * separate item.
         * @return String items
        */
        std::vector<std::string> content() const;

        /**
         * Get a source identifier, if available.
         */
        std::string sourceId() const;

    private:

        std::string sourceId_;
        std::vector<std::string> source_;

    };

}

#endif //PARTICLES_INPUT_SOURCE_HPP
