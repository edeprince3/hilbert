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

#include "timer.h"
#include <omp.h>
#include <iomanip>
#include <cmath>

using namespace std;

namespace hilbert {

    void Timer::start(){
        start_time_ = omp_get_wtime();
        running_ = true;
    }

    void Timer::stop(bool increment) {
        end_time_ = omp_get_wtime();
        runtime_ += end_time_ - start_time_;
        running_ = false;
        if (increment) n_calls_++;
    }

    void Timer::reset() {
        start_time_ = 0.0;
        end_time_ = 0.0;
        runtime_ = 0.0;
        n_calls_ = 0;
        running_ = false;
    }

    std::string Timer::format_time(long double time) {
        time = std::fabs(time); // ignore negative time (shouldn't happen unless you're a time traveler)
        std::string unit = "s"; // default unit is seconds

        if (time < 1.0 && time != 0.0) {
            unit = "ms"; time *= 1.0e3; // milliseconds
            if (time < 1.0) {
                unit = "us"; time *= 1.0e3; // microseconds
                if (time < 1.0) {
                    unit = "ns"; time *= 1.0e3; // nanoseconds
                    if (time < 1.0) {
                        // faster than a clock cycle on modern CPUs... for now :D
                        unit = "ps"; time *= 1.0e3; // picoseconds
                        if (time < 1.0) {
                            unit = "fs"; time *= 1.0e3; // femtoseconds
                            if (time < 1.0) {
                                unit = "as"; time *= 1.0e3; // attoseconds
                                if (time < 1.0) {
                                    // I don't think you can measure time this small
                                    unit = "zs"; time *= 1.0e3; // zeptoseconds
                                    if (time < 1.0) {
                                        unit = "ys"; time *= 1.0e3; // yoctoseconds
                                    }}}}}}}} // I know, I know...
        else if (time >= 60.0) {
            unit = "m"; time /= 60.0; // minutes
            if (time >= 60.0) {
                unit = "h"; time /= 60.0; // hours
                if (time >= 24.0) {
                    unit = "d"; time /= 24.0; // days
                    if (time >= 7.0) {
                        unit = "w"; time /= 7.0; // weeks
                        if (time >= 4.0) {
                            unit = "mo"; time /= 4.0; // months
                            if (time >= 12.0) {
                                // hey, you never know how long a computer would run this code.
                                unit = "y"; time /= 12.0; // years
                                if (time >= 10.0) {
                                    // maybe you should consider buying a new computer
                                    unit = "dec"; time /= 10.0; // decades
                                    if (time >= 10.0) {
                                        // you should REALLY consider buying a new computer
                                        unit = "cen"; time /= 10.0; // centuries
                                        if (time >= 10.0) {
                                            // ok, fine. Don't take my advice.
                                            unit = "mill"; time /= 10.0; // millennia
                                            if (time >= 1000.0) {
                                                unit = "bill"; time /= 1.0e3; // billions of years
                                                if (time >= 1000.0) {
                                                    // have you ever heard of the big bang?
                                                    unit = "trill"; time /= 1.0e3; // trillions of years
                                                    if (time >= 1000) {
                                                        // have you ever heard of the heat death of the universe?
                                                        unit = "quad"; time /= 1.0e3; // quadrillions of years
                                                        if (time >= 1000) {
                                                            // you're most likely dead by now... maybe... hopefully?
                                                            unit = "quin"; time /= 1.0e3; // quintillions of years
                                                        }}}}}}}}}}}}} // I'm not going to write this out any further...
                                                                      // You get the idea... Does it ever end?

        // return the formatted time
        string time_str;
        stringstream ss;
        ss << std::right << std::setfill(' ') << std::setw(8) << std::fixed << std::setprecision(precision_)
           << time;
        ss >> time_str;

        return time_str + " " + unit;
    }

    string Timer::elapsed() const {
        return format_time(runtime_);
    }

    string Timer::average_time() const{
        return format_time(runtime_ / (double) n_calls_);
    }

    string Timer::current_time() const {
        return format_time(end_time_ - start_time_);
    }


}

