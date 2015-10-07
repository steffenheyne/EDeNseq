# SVM with stochastic gradient

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA

#define HINGELOSS 1
#define SMOOTHHINGELOSS 2
#define SQUAREDHINGELOSS 3
#define LOGLOSS 10
#define LOGLOSSMARGIN 11

#CHECKLIMITS activates the limit check for graph operations

CXX=g++
OPTS=-g -O3 -Wno-deprecated -Wno-unused-local-typedefs -DLOSS=1 -std=c++11 -pthread -fopenmp -DUSEMULTITHREAD -DEIGEN_DONT_PARALLELIZE # -DDEBUGON # -DCHECKLIMITS
#OPTS=-g -Wno-deprecated -DLOSS=1 -static -DEIGEN_DONT_PARALLELIZE # -DDEBUGON # -DCHECKLIMITS
#OPTS=-g -Wno-deprecated -DLOSS=1 #-DDEBUGON
#OPTS=-g -Wno-deprecated -static -pg

CXXFLAGS=${OPTS} -Wall -I eigen-eigen-3.20 -Isparsehash-2.0.2 # -DNDEBUG # -fprofile-use -fprofile-correction # -fprofile-generate #-mtune=native -march=native
#CXXFLAGS=${OPTS} -Wall -I eigen-eigen-3.20 -Isparsehash-2.0.2 -fprofile-use -fprofile-correction # -fprofile-generate #-mtune=native -march=native 
#CXXFLAGS=${OPTS} -Wall -I eigen-eigen-3.20 -Isparsehash-2.0.2 -fprofile-generate #-mtune=native -march=native
#CXXFLAGS=${OPTS} -Wall -I eigen-eigen-3.20 -Isparsehash-2.0.2 -mtune=native -march=native # -ftree-vectorizer-verbose=1 -fopt-info-vec-missed
#CXXFLAGS= ${OPTS} 

#LIBS=-lm -lz -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
#LIBS=-lm -lz -lpthread -Wl,-u,pthread
LIBS=-lm -lz -lpthread

PROGRAMS=EDeNseq	

all: ${PROGRAMS}

clean:
	-rm ${PROGRAMS}
	-rm *.o
	
objects := $(patsubst %.cc,%.o,$(wildcard *.cc))
objects_mains := $(patsubst %,%.o,$(PROGRAMS))


EDeNseq : $(objects)
	${CXX} ${CXXFLAGS} ${OBABEL} $(filter-out $(objects_mains),$(objects)) $@.o ${LIBS} -o $@
	 
EDeNseq.o: EDeNseq.cc gzstream.h MinHashEncoder.h
	 ${CXX} ${CXXFLAGS} -c EDeNseq.cc -o EDeNseq.o

SeqClassifyManager.o:SeqClassifyManager.cc SeqClassifyManager.h MinHashEncoder.h

SeqClusterManager.o:SeqClusterManager.h MinHashEncoder.h MinHashEncoder.cc

TestManager.o:TestManager.cc TestManager.h MinHashEncoder.h

MinHashEncoder.o:MinHashEncoder.h Data.h

BaseManager.o:BaseManager.h

Data.o:Data.h

Parameters.o:Parameters.h Utility.h

Utility.o:Utility.h

gzstream.o: gzstream.cc gzstream.h
