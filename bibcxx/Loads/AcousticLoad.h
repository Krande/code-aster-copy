#ifndef ACOUSTICLOAD_H_
#define ACOUSTICLOAD_H_

/**
 * @file AcousticLoad.h
 * @brief Fichier entete de la classe AcousticLoad
 * @author Nicolas Sellenet
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
#include "Loads/AcousticLoadDescription.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

template < class ConstantFieldOnCellsType >
class AcousticLoad : public DataStructure {
  private:
    /** @brief Vecteur Jeveux '.TYPE' */
    JeveuxVectorChar8 _type;
    /** @brief sd_char_chth '.CHAC' */
    AcousticLoadDescriptionPtr< ConstantFieldOnCellsType > _acouLoadDesc;

  public:
    /**
     * @typedef AcousticLoadPtr
     * @brief Pointeur intelligent vers un AcousticLoad
     */
    typedef std::shared_ptr< AcousticLoad > AcousticLoadPtr;

    /**
     * @brief Constructeur
     */
    AcousticLoad( void ) = delete;

    /**
     * @brief Constructeur
     */
    AcousticLoad( const ModelPtr &model )
        : AcousticLoad( ResultNaming::getNewResultName(), model ) {};

    /**
     * @brief Constructeur
     */
    AcousticLoad( const std::string name, const ModelPtr &model )
        : DataStructure( name, 8, "CHAR_ACOU" ),
          _acouLoadDesc( std::make_shared< AcousticLoadDescription< ConstantFieldOnCellsType > >(
              getName() + ".CHAC", model ) ),
          _type( getName() + ".TYPE" ) {};

    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const {
        return _acouLoadDesc->getFiniteElementDescriptor();
    };

    AcousticLoadDescriptionPtr< ConstantFieldOnCellsType > getAcousticLoadDescription() const {
        return _acouLoadDesc;
    };

    ConstantFieldOnCellsComplexPtr getMultiplicativeField() const {
        return _acouLoadDesc->getMultiplicativeField();
    };

    auto getImposedField() const { return _acouLoadDesc->getImposedField(); }

    bool hasLoadField( const std::string &load_name ) const {
        return _acouLoadDesc->hasLoadField( load_name );
    };

    auto getConstantLoadField( const std::string name ) const {
        return _acouLoadDesc->getConstantLoadField( name );
    };

    /**
     * @brief Get the model
     */
    ModelPtr getModel() { return _acouLoadDesc->getModel(); };

    /**
     * @brief Get the mesh
     */
    BaseMeshPtr getMesh() { return _acouLoadDesc->getMesh(); };

    bool build() { return _acouLoadDesc->build(); };
};

/**********************************************************
 *  Explicit instantiation of template classes
 **********************************************************/

/** @typedef AcousticLoadComplex Class d'une charge mécanique réelle */
typedef AcousticLoad< ConstantFieldOnCellsComplex > AcousticLoadComplex;

/** @typedef AcousticLoad  */
template < class ConstantFieldOnCellsType >
using AcousticLoadPtr = std::shared_ptr< AcousticLoad< ConstantFieldOnCellsType > >;

typedef std::shared_ptr< AcousticLoadComplex > AcousticLoadComplexPtr;

/** @typedef std::list de AcousticLoad */
typedef std::list< AcousticLoadComplexPtr > ListAcouLoadComplex;
/** @typedef Iterateur sur une std::list de AcousticLoad */
typedef ListAcouLoadComplex::iterator ListAcouLoadComplexIter;
/** @typedef Iterateur constant sur une std::list de AcousticLoad */
typedef ListAcouLoadComplex::const_iterator ListAcouLoadComplexCIter;

#endif /* ACOUSTICLOAD_H_ */
