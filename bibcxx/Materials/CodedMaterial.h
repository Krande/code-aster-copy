#ifndef CODEDMATERIAL_H_
#define CODEDMATERIAL_H_

/**
 * @file CodedMaterial.h
 * @brief Fichier entete de la classe CodedMaterial
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

#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Materials/MaterialField.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

/**
 * @class CodedMaterial
 * @brief Coded material
 * @author Nicolas Sellenet
 */
class CodedMaterial : public DataStructure {
  private:
    MaterialFieldPtr _mater;
    ModelPtr _model;
    ConstantFieldOnCellsLongPtr _field;
    JeveuxVectorChar8 _grp;
    JeveuxVectorLong _nGrp;
    std::vector< JeveuxVectorLong > _vecOfCodiVectors;
    std::vector< JeveuxVectorReal > _vecOfR8;
    std::vector< JeveuxVectorLong > _vecOfIa;

  public:
    /**
     * @typedef CodedMaterialPtr
     * @brief Pointeur intelligent vers un CodedMaterial
     */
    typedef std::shared_ptr< CodedMaterial > CodedMaterialPtr;

    /**
     * @brief Constructeur
     */
    CodedMaterial( void ) = delete;

    CodedMaterial( const std::string &name, const MaterialFieldPtr &mater, const ModelPtr &model );

    CodedMaterial( const MaterialFieldPtr &mater, const ModelPtr &model )
        : CodedMaterial( ResultNaming::getNewResultName(), mater, model ) {};

    /**
     * @brief Function to allocate the coded material
     * @return return false if coded material already exists
     */
    bool allocate( bool force = false );

    /**
     * @brief Function to know if elastic properties depend of a function
     */
    bool constant() const;

    /**
     * @brief Get the .MATE_CODE
     */
    ConstantFieldOnCellsLongPtr getCodedMaterialField() const { return _field; };

    /**
     * @brief Function membre
     * @return material field
     */
    MaterialFieldPtr getMaterialField() const { return _mater; };

    bool exists() const;
};

/**
 * @typedef CodedMaterialPtr
 * @brief Pointeur intelligent vers un CodedMaterial
 */
typedef std::shared_ptr< CodedMaterial > CodedMaterialPtr;

#endif /* CODEDMATERIAL_H_ */
