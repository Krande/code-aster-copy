#ifndef MEDTYPES_H_
#define MEDTYPES_H_

/**
 * @file MedTypes.h
 * @brief Fichier entete de la classe MedProfile
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

#include "astercxx.h"

#ifdef ASTER_HAVE_MED
#include "med.h"

/** @brief med types used in code_aster */
constexpr std::array< med_geometry_type, 20 > medTypes = { 1,   102, 103, 104, 203, 206, 207,
                                                           204, 208, 209, 304, 310, 306, 315,
                                                           318, 305, 313, 308, 320, 327 };

#endif
#endif /* MEDTYPES_H_ */
