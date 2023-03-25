#ifndef VARIANTMODALBASIS_H_
#define VARIANTMODALBASIS_H_

/**
 * @file VariantModalBasis.h
 * @brief Fichier entete de la classe VariantModalBasis
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

#include "astercxx.h"

#include "aster_pybind.h"

#include "Results/GeneralizedModeResult.h"
#include "Results/ModeResult.h"

typedef std::variant< ModeResultPtr, GeneralizedModeResultPtr > ModalBasisVariant;

template < typename ObjectPointer >
ModalBasisVariant getModalBasis( ObjectPointer self ) {
    auto mat1 = self->getModalBasisFromGeneralizedModeResult();
    if ( mat1 != nullptr )
        return ModalBasisVariant( mat1 );
    auto mat2 = self->getModalBasisFromModeResult();
    return ModalBasisVariant( mat2 );
};

#endif /* VARIANTMODALBASIS_H_ */
