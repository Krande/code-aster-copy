#ifndef MEDCALCULATIONSEQUENCE_H_
#define MEDCALCULATIONSEQUENCE_H_

/**
 * @file MedCalculationSequence.h
 * @brief Fichier entete de la classe MedCalculationSequence
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

#include "astercxx.h"

#ifdef ASTER_HAVE_MED
#include "med.h"

#include <iostream>
#include <memory>
#include <utility>

/**
 * @class MedCalculationSequence
 * @brief Med profile interface
 * @author Nicolas Sellenet
 */
class MedCalculationSequence {
  private:
    /** @brief step id */
    med_int _numdt = 0;
    /** @brief iteration id */
    med_int _numit = 0;
    /** @brief time step value */
    med_float _dt = 0.;

  public:
    /**
     * @typedef MedCalculationSequencePtr
     * @brief Pointeur intelligent vers un MedCalculationSequence
     */
    typedef std::shared_ptr< MedCalculationSequence > MedCalculationSequencePtr;

    /**
     * @brief Constructor
     * @param numdt step id
     * @param numit iteration id
     * @param dt time step value
     */
    MedCalculationSequence( int numdt, int numit, float dt )
        : _numdt( numdt ), _numit( numit ), _dt( dt ) {};

    /** @brief Get time step value */
    med_float getDt() const { return _dt; };

    /** @brief Get pair {step id, iteration id} */
    std::pair< med_int, med_int > getNumDtNumIt() const {
        return std::make_pair( _numdt, _numit );
    };
};

/**
 * @typedef MedCalculationSequencePtr
 * @brief Pointeur intelligent vers un MedCalculationSequence
 */
typedef std::shared_ptr< MedCalculationSequence > MedCalculationSequencePtr;

#endif
#endif /* MEDCALCULATIONSEQUENCE_H_ */
