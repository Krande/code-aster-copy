/**
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "ParallelUtilities/TemplateVectorTools.h"

int getSize( const MedVector< double >::ElementValue &in ) { return in.getComponentNumber(); };

int getTotalSize( const MedVector< double > &toCopy ) { return toCopy.totalSize(); };

void allocate( MedVector< double > &in, const int &size1, const int &size2 ) {
    in.setSize( size1 );
    in.setTotalSize( size2 );
};

void update( MedVector< double >::ElementValue in ) {};

void allocateOccurence( MedVector< double > &in, const int &pos, const int &size ) {
    in.setElement( pos, size );
};

int getSize( const MedVector< long int >::ElementValue &in ) { return in.getComponentNumber(); };

int getTotalSize( const MedVector< long int > &toCopy ) { return toCopy.totalSize(); };

void allocate( MedVector< long int > &in, const int &size1, const int &size2 ) {
    in.setSize( size1 );
    in.setTotalSize( size2 );
};

void update( MedVector< long int >::ElementValue in ) {};

void allocateOccurence( MedVector< long int > &in, const int &pos, const int &size ) {
    in.setElement( pos, size );
};

int getSize( const ArrayWrapper< JeveuxVectorReal >::ElementValue &in ) {
    return in.getComponentNumber();
};

void update( typename ArrayWrapper< JeveuxVectorReal >::ElementValue in ) {};

int getSize( const ArrayWrapper< JeveuxVectorLogical >::ElementValue &in ) {
    return in.getComponentNumber();
};

void update( typename ArrayWrapper< JeveuxVectorLogical >::ElementValue in ) {};

int getSize( const ArrayWrapper< JeveuxVectorLong >::ElementValue &in ) {
    return in.getComponentNumber();
};

void update( typename ArrayWrapper< JeveuxVectorLong >::ElementValue in ) {};

int getSize( const ArrayWrapper< VectorReal >::ElementValue &in ) {
    return in.getComponentNumber();
};

void update( typename ArrayWrapper< VectorReal >::ElementValue in ) {};

int getSize( const ArrayWrapper< VectorLong >::ElementValue &in ) {
    return in.getComponentNumber();
};

void update( typename ArrayWrapper< VectorLong >::ElementValue in ) {};
