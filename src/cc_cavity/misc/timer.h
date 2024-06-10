/*
 *  @BEGIN LICENSE
 *
 *  Hilbert: a space for quantum chemistry plugins to Psi4
 *
 *  Copyright (c) 2020 by its authors (LICENSE).
 *
 *  The copyrights for code used from other parties are included in
 *  the corresponding files.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 *  @END LICENSE
 */

#ifndef CC_CAVITY_TIMER_H
#define CC_CAVITY_TIMER_H
#include <string>

namespace hilbert {

    class Timer {

    private:
        /***** TIMERS *****/
        long double start_time_; // current start time
        long double end_time_; // current end time
        long double runtime_ = 0.0; // timer runtime_
        bool running_ = false; // is the timer running_?
        size_t n_calls_ = 0; // number of times the timer has been called

    public:


        Timer() = default;
        ~Timer() = default;

        void start(); // start timer
        void stop(bool increment = true); // stop timer
        void reset(); // reset timer

        static inline int precision_ = 3; // precision_ of the timer output

        /**
         * @brief return the time as a human readable string
         * @param time time in seconds
         * @return human readable string
         */
        std::string elapsed() const;

        /**
         * @brief return the time in seconds
         * @return time in seconds
         */
        std::string current_time() const;

        /**
         * format time as string
         */
        static std::string format_time(long double time);


        /**
         * @brief return the average time as a human readable string
         * @param time time in seconds
         * @return human readable string
         */
        std::string average_time() const;

        /**
         * Get runtime_ as double
         */
        long double get_runtime() const { return runtime_; }

        /**
         * Get number of calls
         */
        size_t num_calls() const { return n_calls_; }
    };

}


#endif //CC_CAVITY_TIMER_H
