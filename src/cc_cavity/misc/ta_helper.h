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

using namespace std;
using namespace TA;
namespace Helper {
    template<typename T = double>
    struct TA_Helper {

        static inline size_t tile_size_ = -1; // default value to use all elements in a tile

        /**
         * create a new tiled range with the dimensions specified in the vector, N
         * @param N the target dimensions of the tiled range
         * @return the tiled range
         */
        static inline TiledRange makeRange(const initializer_list<size_t> &N);
        static inline TiledRange makeRange(const vector<size_t> &N);

        /**
         * create a new tiled array with the dimensions specified in the vector, N
         * @param world the world object
         * @param N the target dimensions of the tiled array
         * @param fillZero if true, the array will be filled with zeros
         * @return the tiled array
         */
        static inline TArray<T>
        makeTensor(World &world, const initializer_list<size_t> &N, bool fillZero);

        /**
         * create a new tiled array with the dimensions specified in the vector, N, and fill it with values from data
         * @param world the world object
         * @param N the target dimensions of the tiled array
         * @param data the data to fill the array with
         * @param Off the offsets of the data in each dimension (Default: {}, i.e. no offset)
         * @return the tiled array
         */
        static inline TArray<T>
        makeTensor(World &world, const initializer_list<size_t> &N, const T *data, initializer_list<size_t> Off = {});

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
        static inline TArray<T>
        makeTensor(World &world,
                   const initializer_list<size_t> &NL,
                   const initializer_list<size_t> &NR,
                   const T * const *data,
                   initializer_list<size_t> Off = {});

        /**
         * loop over all elements in the tiled array and apply a function, op, using each element
         * @param tensor the tiled array to iterate over
         * @param op the function to perform an operation using each element
         *           the function must take a Tensor<T> (the tile) and a vector<size_t> (the index) as arguments
         */
        static inline void
        forall(TArray<T> &tensor, function<void(Tensor<T> &, vector<size_t> &)> op);
        
        /**
         * extract elements from the tiled array, tensor, and store them in the array, data
         * @param N the dimensions of the tiled array
         * @param Off the offsets of the data in each dimension (Default: {}, i.e. no offset)
         * @result data the array to store the data in
         */
        static inline T*
        arrayFromTensor(TArray<T> A, const initializer_list<size_t> &N, initializer_list<size_t> Off = {});

    };
// make typedefs for TA_Helper
template <typename T>
using HelperT = TA_Helper<T>;
typedef HelperT<double> HelperD;
typedef HelperT<complex<double>> HelperC;
}
#endif
#include "ta_helper.cc"
