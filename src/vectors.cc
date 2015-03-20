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

// $Id: vectors.cpp,v 1.22 2009/01/31 12:37:33 cvs Exp $

#include "vectors.h"
#include <ios>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <cctype>
#include <stdexcept>

using namespace Eigen;

namespace {
	template<typename T> inline T min(T a, T b) {
		return (a < b) ? a : b;
	}
	template<typename T> inline T max(T a, T b) {
		return (a < b) ? b : a;
	}
}

VFloat setOne(const VFloat x) {
	return 1;
}

SVector binarize(const SVector &v) {
	return v.unaryExpr(std::ptr_fun(setOne));
}

SVector elementwise_product(const SVector &v1, const SVector &v2) {
	return v1.cwiseProduct(v2);
}

// clamp operation
template<typename Scalar>
struct CwiseClampOp {
		CwiseClampOp(const Scalar& hardlimit) :
				m_hardlimit(hardlimit) {
		}
		const Scalar operator()(const Scalar& x) const {
			return x < -m_hardlimit ? -m_hardlimit : (x > m_hardlimit ? m_hardlimit : x);
		}
		Scalar m_hardlimit;
};

SVector elementwise_product_clamped(const SVector &v1, const SVector &v2, const VFloat hardlimit) {
	return v1.cwiseProduct(v2).unaryExpr(CwiseClampOp<VFloat>(hardlimit));
}

bool loadSparse(std::istream &f, SVector &v) {
	v.setZero();
	for (;;) {
		int c = f.get();
		if (!f.good() || (c == '\n' || c == '\r'))
			break;
		if (::isspace(c))
			continue;
		int i;
		f.unget();
		f >> std::skipws >> i >> std::ws;
		if (f.get() != ':') {
			f.unget();
			f.setstate(std::ios::badbit);
			break;
		}
		VFloat x;
		f >> std::skipws >> x;
		if (!f.good())
			break;
		v.insert(i) = x;
	}
	return f.good();
}

bool loadSparseBinary(std::istream &f, SVector &v) {
	v.setZero();
	int npairs = 0;
	f.read((char*) &npairs, sizeof(int));
	if (npairs < 0)
		f.setstate(std::ios::badbit);
	if (!f.good())
		return false;
	for (int i = 0; i < npairs; i++) {
		Pair pair;
		f.read((char*) &pair, sizeof(Pair));
		if (f.good())
			v.insert(pair.i) = pair.v;
	}
	return f.good();
}

std::string printSparse(const SVector &sv) {
	// create sparse vector representation
	// coefficient format: "index:value"
	// separator: " "
	std::stringstream s;
	s << std::scientific << std::setprecision(sizeof(VFloat) == 4 ? 7 : 16);
	for (SVector::InnerIterator it(sv); it; ++it) {
		s << " " << it.index() << ":"<< std::setprecision(9) << it.value();
	}
	s << std::endl;
	return s.str();
}

std::string serialize(const SVector &sv) {
	// create sparse vector representation
	// coefficient format: "index:value"
	// separator: " "
	std::stringstream s;
	s << std::scientific << std::setprecision(sizeof(VFloat) == 4 ? 7 : 16);
	for (SVector::InnerIterator it(sv); it; ++it) {
		s << " " << it.index() << ":"<< std::setprecision(9) << it.value();
	}
	return s.str();
}

std::string printSparseBinary(const SVector &sv) {
	// create binary sparse vector representation
	// format: number of pairs, all pairs of index and value
	std::stringstream s;
	int npairs = sv.nonZeros();
	s.write((const char*) &npairs, sizeof(int));
	for (SVector::InnerIterator it(sv); it; ++it) {
		int index = it.index();
		s.write((const char*) &index, sizeof(int));
		VFloat value = it.value();
		s.write((const char*) &value, sizeof(VFloat));
	}
	return s.str();
}

/* -------------------------------------------------------------
 Local Variables:
 c++-font-lock-extra-types: ( "\\sw+_t" "[A-Z]\\sw*[a-z]\\sw*" "std::\\sw+")
 End:
 ------------------------------------------------------------- */

