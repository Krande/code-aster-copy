/**
 * @file MedFilter.cxx
 * @brief Implementation de MedFilter
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "IOManager/MedFilter.h"

#ifdef ASTER_HAVE_MED
MedFilter::MedFilter( const MedFilePointer &filePtr, med_int nentity, med_int nvaluesperentity,
                      med_int nconstituentpervalue, med_int constituentselect,
                      med_switch_mode switchmode, med_storage_mode storagemode, med_size start,
                      med_size stride, med_size count, med_size blocksize, med_size lastblocksize,
                      MedProfilePtr prof )
    : _filePtr( filePtr ),
      _nentity( nentity ),
      _nvaluesperentity( nvaluesperentity ),
      _nconstituentpervalue( nconstituentpervalue ),
      _constituentselect( constituentselect ),
      _switchmode( switchmode ),
      _storagemode( storagemode ),
      _start( start ),
      _stride( stride ),
      _count( count ),
      _blocksize( blocksize ),
      _lastblocksize( lastblocksize ),
      _profile( prof ) {
    _filter = MEDfilterAllocate( 1 );
    std::string profilename( "" );
    if ( _profile != nullptr )
        profilename = _profile->getName();
    MEDfilterBlockOfEntityCr( _filePtr.getFileId(), _nentity, _nvaluesperentity,
                              _nconstituentpervalue, _constituentselect, _switchmode, _storagemode,
                              profilename.c_str(), _start, _stride, _count, _blocksize,
                              _lastblocksize, _filter );
    _exists = true;
};

MedFilter::~MedFilter() {
    if ( _exists )
        MEDfilterDeAllocate( 1, _filter );
};

#endif
