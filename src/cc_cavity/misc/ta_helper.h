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

#ifndef TA_HELPER_H
#define TA_HELPER_H
#include <tiledarray.h>
#include <cstdio>

namespace TA_Helper {

    using std::vector, std::initializer_list;
    using namespace TA;

    static inline size_t tile_size_ = static_cast<size_t>(-1); // default value to use all elements in a tile

    /**
     * create a new tiled range with the dimensions specified in the vector, N
     * @param N the target dimensions of the tiled range
     * @return the tiled range
     */
    TiledRange makeRange(const initializer_list<size_t> &N);
    TiledRange makeRange(const vector<size_t> &N);

    /**
     * create a new tiled array with the dimensions specified in the vector, N
     * @param world the world object
     * @param N the target dimensions of the tiled array
     * @param fillZero if true, the array will be filled with zeros
     * @return the tiled array
     */
    TArrayD makeTensor(World &world, const initializer_list<size_t> &N, bool fillZero = true);

    /**
     * create a new tiled array with the dimensions specified in the vector, N, and fill it with values from data
     * @param world the world object
     * @param N the target dimensions of the tiled array
     * @param data the data to fill the array with
     * @param Off the offsets of the data in each dimension (Default: {}, i.e. no offset)
     * @return the tiled array
     */
    TArrayD makeTensor(World &world, const initializer_list<size_t> &N, const double *data,
               initializer_list<size_t> Off = {});

    /**
     * create a new tiled array with the dimensions specified in the vector, N, and fill it with values from data
     * where data is a 2D array and NL and NR are the dimensions within the 2D array
     * @param world  the world object
     * @param NL  the left dimensions of the 2D array
     * @param NR the right dimensions of the 2D array
     * @param data the data to fill the array with
     * @param Off the offsets of the data in each dimension (Default: {}, i.e. no offset)
     * @return the tiled array
     */
    TArrayD makeTensor(World &world,
               const initializer_list<size_t> &NL,
               const initializer_list<size_t> &NR,
               const double *const *data,
               initializer_list<size_t> Off = {});


}

#endif