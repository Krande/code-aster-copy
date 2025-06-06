#ifndef MEDFILTER_H_
#define MEDFILTER_H_

/**
 * @file MedFilter.h
 * @brief Fichier entete de la classe MedFilter
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "IOManager/MedFilePointer.h"
#include "IOManager/MedProfile.h"

#include <iostream>
#include <memory>
#include <string>

#ifdef ASTER_HAVE_MED
#include "med.h"
/**
 * @class MedFilter
 * @brief Med filter interface
 *        usefull in parallel to read only parts of mesh, field, ...
 * @author Nicolas Sellenet
 */
class MedFilter {
  private:
    const MedFilePointer &_filePtr;
    med_int _nentity;
    med_int _nvaluesperentity;
    med_int _nconstituentpervalue;
    med_int _constituentselect;
    med_switch_mode _switchmode;
    med_storage_mode _storagemode;
    MedProfilePtr _profile;
    med_size _start;
    med_size _stride;
    med_size _count;
    med_size _blocksize;
    med_size _lastblocksize;
    bool _exists = false;
    med_filter *_filter;

  public:
    /**
     * @typedef MedFilterPtr
     * @brief Pointeur intelligent vers un MedFilter
     */
    typedef std::shared_ptr< MedFilter > MedFilterPtr;

    /**
     * @brief med filter constructor
     * @param filePtr med file pointer
     * @param nentity total number of entities (over all processes)
     * @param nvaluesperentity number of values (eg: point Gauss number)
     * @param nconstituentpervalue  number of consituent (eg: component number)
     * @param constituentselect ??
     * @param switchmode ??
     * @param storagemode ??
     * @param start id of first read entity (for current process)
     * @param stride offset between entity id (often 1)
     * @param count number of read block
     * @param blocksize block size (local number of entities)
     * @param lastblocksize ??
     * @param
     */
    MedFilter( const MedFilePointer &filePtr, med_int nentity, med_int nvaluesperentity,
               med_int nconstituentpervalue, med_int constituentselect, med_switch_mode switchmode,
               med_storage_mode storagemode, med_size start, med_size stride, med_size count,
               med_size blocksize, med_size lastblocksize,
               MedProfilePtr prof = MedProfilePtr( nullptr ) );

    ~MedFilter();

    /** @brief get c pointer on filter */
    const med_filter *getPointer() const { return _filter; }
};

/**
 * @typedef MedFilterPtr
 * @brief Pointeur intelligent vers un MedFilter
 */
typedef std::shared_ptr< MedFilter > MedFilterPtr;

#endif
#endif /* MEDFILTER_H_ */
