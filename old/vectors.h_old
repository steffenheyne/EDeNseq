// -*- C++ -*-
// Little library of vectors and sparse vectors
// Copyright (C) 2007- Leon Bottou

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA



// $Id: vectors.h,v 1.13 2007/10/02 20:40:06 cvs Exp $

#ifndef VECTORS_H
#define VECTORS_H 1

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <cstring>
#include <iostream>

#include <sstream>
#include <fstream>

#include <cassert>

#include <vector>
#include <cmath>

//#include "wrapper.h"

using namespace Eigen;

typedef float VFloat;
typedef SparseVector<VFloat> SVector;
typedef Matrix<VFloat, Dynamic, 1> FVector;

VFloat setOne(const VFloat x);

SVector binarize(const SVector &v);

SVector elementwise_product(const SVector &v1, const SVector &v2);

SVector elementwise_product_clamped(const SVector &v1, const SVector &v2, const VFloat hardlimit=1e9);

bool loadSparse(std::istream &f, SVector &v);

bool loadSparseBinary(std::istream &f, SVector &v);

std::string printSparse(const SVector &sv);

std::string printSparseBinary(const SVector &sv);

std::string serialize(const SVector &sv);

struct Pair
{
	int i;
	VFloat v;
};

#endif

/* -------------------------------------------------------------
   Local Variables:
   c++-font-lock-extra-types: ("\\sw+_t" "[A-Z]\\sw*[a-z]\\sw*" "std::\\sw+")
   End:
   ------------------------------------------------------------- */
