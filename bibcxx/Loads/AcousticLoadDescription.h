#ifndef ACOUSTICLOADDESCRIPTION_H_
#define ACOUSTICLOADDESCRIPTION_H_

/**
 * @file AcousticLoad.h
 * @brief Fichier entete de la classe AcousticLoad
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

#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

template < typename ConstantFieldOnCellsType >
class AcousticLoadDescription : public DataStructure {

  public:
    typedef std::shared_ptr< ConstantFieldOnCellsType > ConstantFieldOnCellsTypePtr;

  private:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Vecteur Jeveux '.MODEL.NOMO' */
    JeveuxVectorChar8 _modelName;
    /** @brief FiniteElementDescriptor of load '.LIGRE' */
    FiniteElementDescriptorPtr _FEDesc;
    /** @brief Carte '.CIMPO' */
    ConstantFieldOnCellsTypePtr _imposedValues;
    /** @brief Carte '.CMULT' */
    ConstantFieldOnCellsComplexPtr _multiplier;
    /** @brief Carte '.VITFA' */
    ConstantFieldOnCellsTypePtr _speedValues;

  public:
    /**
     * @typedef AcousticLoadPtr
     * @brief Pointeur intelligent vers un AcousticLoad
     */
    typedef std::shared_ptr< AcousticLoadDescription > AcousticLoadDescriptionPtr;

    /**
     * @brief Constructeur
     */
    AcousticLoadDescription( void ) = delete;

    /**
     * @brief Constructeur
     */
    AcousticLoadDescription( const ModelPtr &model )
        : AcousticLoadDescription( ResultNaming::getNewResultName(), model ) {};

    /**
     * @brief Constructeur
     */
    AcousticLoadDescription( const std::string name, const ModelPtr &model )
        : DataStructure( name, 13, "CHAR_CHAC" ),
          _model( model ),
          _FEDesc( std::make_shared< FiniteElementDescriptor >( getName() + ".LIGRE",
                                                                model->getMesh() ) ),
          _modelName( JeveuxVectorChar8( getName() + ".MODEL.NOMO" ) ),
          _imposedValues(
              std::make_shared< ConstantFieldOnCellsType >( getName() + ".CIMPO", _FEDesc ) ),
          _multiplier(
              std::make_shared< ConstantFieldOnCellsComplex >( getName() + ".CMULT", _FEDesc ) ),
          _speedValues(
              std::make_shared< ConstantFieldOnCellsType >( getName() + ".VFACE", _FEDesc ) ) {};

    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    ConstantFieldOnCellsComplexPtr getMultiplicativeField() const { return _multiplier; };

    ConstantFieldOnCellsTypePtr getImposedField() const { return _imposedValues; }

    /**
     * @brief Get the model
     */
    ModelPtr getModel() { return _model; };

    bool hasLoadField( const std::string load_name ) const {
        if ( load_name == "VFACE" ) {
            return ( _speedValues && _speedValues->exists() );
        } else {
            AS_ASSERT( false );
        }

        return false;
    };

    ConstantFieldOnCellsTypePtr getConstantLoadField( const std::string name ) const {
        if ( name == "VFACE" ) {
            return _speedValues;
        } else {
            AS_ASSERT( false );
        }

        return nullptr;
    };

    /**
     * @brief Get the model
     */
    BaseMeshPtr getMesh() { return _model->getMesh(); };

    bool build() {
        _FEDesc->build();

        return true;
    };
};

/**********************************************************
 *  Explicit instantiation of template classes
 **********************************************************/

/** @typedef AcousticLoadDescriptionFunc Class d'une charge mécanique de fonctions */
typedef AcousticLoadDescription< ConstantFieldOnCellsComplex > AcousticLoadDescriptionComplex;

template < typename ConstantFieldOnCellsType >
using AcousticLoadDescriptionPtr =
    std::shared_ptr< AcousticLoadDescription< ConstantFieldOnCellsType > >;

#endif /* ACOUSTICLOADDESCRIPTION_H_ */
