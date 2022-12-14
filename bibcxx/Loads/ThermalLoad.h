#ifndef THERMALLOAD_H_
#define THERMALLOAD_H_

/**
 * @file ThermalLoad.h
 * @brief Fichier entete de la classe ThermalLoad
 * @author Jean-Pierre Lefebvre
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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
#include "DataFields/FieldOnCells.h"
#include "Loads/ThermalLoadDescription.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ThermalLoad
 * @brief Classe definissant une charge thermique (issue d'AFFE_CHAR_THER)
 * @author Jean-Pierre Lefebvre
 */
template < class ConstantFieldOnCellsType >
class ThermalLoad : public DataStructure {
  private:
    /** @brief Vecteur Jeveux '.TYPE' */
    JeveuxVectorChar8 _type;
    /** @brief sd_char_chth '.CHTH' */
    ThermalLoadDescriptionPtr< ConstantFieldOnCellsType > _therLoadDesc;

  public:
    using ConstantFieldOnCellsTypePtr = std::shared_ptr< ConstantFieldOnCellsType >;

    /**
     * @brief Constructeur
     */
    ThermalLoad( void ) = delete;

    /**
     * @brief Constructeur
     */
    ThermalLoad( const ModelPtr &currentModel )
        : ThermalLoad( ResultNaming::getNewResultName(), currentModel ){};

    /**
     * @brief Constructeur
     */
    ThermalLoad( const std::string name, const ModelPtr &currentModel )
        : DataStructure( name, 8, "CHAR_THER" ),
          _therLoadDesc( std::make_shared< ThermalLoadDescription< ConstantFieldOnCellsType > >(
              getName() + ".CHTH", currentModel ) ),
          _type( getName() + ".TYPE" ){};

    ThermalLoadDescriptionPtr< ConstantFieldOnCellsType > getThermalLoadDescription() const {
        return _therLoadDesc;
    };

    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const {
        return _therLoadDesc->getFiniteElementDescriptor();
    };

    /**
     * @brief Get the model
     */
    ModelPtr getModel() const { return _therLoadDesc->getModel(); };

    /**
     * @brief Get the mesh
     */
    BaseMeshPtr getMesh() const { return _therLoadDesc->getMesh(); };

    ConstantFieldOnCellsTypePtr getConstantLoadField( const std::string name ) const {
        return _therLoadDesc->getConstantLoadField( name );
    }

    FieldOnCellsRealPtr getLoadField( const std::string name ) const {
        return _therLoadDesc->getLoadField( name );
    }
    bool hasLoadField( const std::string name ) const {
        return _therLoadDesc->hasLoadField( name );
    }

    bool hasLoadResult() const { return _therLoadDesc->hasLoadResult(); }

    std::string getLoadResultName() const { return _therLoadDesc->getLoadResultName(); }

    ConstantFieldOnCellsTypePtr getImposedField() const { return _therLoadDesc->getImposedField(); }

    ConstantFieldOnCellsRealPtr getMultiplicativeField() const {
        return _therLoadDesc->getMultiplicativeField();
    }

    JeveuxVectorChar8 getType() const { return _type; }

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() {
        _therLoadDesc->updateValuePointers();
        _type->updateValuePointer();
    };

    bool build() { return _therLoadDesc->build(); };
};

/**********************************************************/
/*  Explicit instantiation of template classes
/**********************************************************/

/** @typedef ThermalLoadReal Class d'une charge m??canique r??elle */
typedef ThermalLoad< ConstantFieldOnCellsReal > ThermalLoadReal;
/** @typedef ThermalLoadFunc Class d'une charge m??canique de fonctions */
typedef ThermalLoad< ConstantFieldOnCellsChar24 > ThermalLoadFunction;

typedef std::shared_ptr< ThermalLoadReal > ThermalLoadRealPtr;
typedef std::shared_ptr< ThermalLoadFunction > ThermalLoadFunctionPtr;

/** @typedef ThermalLoad  */
template < class ConstantFieldOnCellsType >
using ThermalLoadPtr = std::shared_ptr< ThermalLoad< ConstantFieldOnCellsType > >;

/** @typedef std::list de ThermalLoad */
typedef std::list< ThermalLoadRealPtr > ListTherLoadReal;
/** @typedef Iterateur sur une std::list de ThermalLoad */
typedef ListTherLoadReal::iterator ListTherLoadRealIter;
/** @typedef Iterateur constant sur une std::list de ThermalLoad */
typedef ListTherLoadReal::const_iterator ListTherLoadRealCIter;

/** @typedef std::list de ThermalLoad */
typedef std::list< ThermalLoadFunctionPtr > ListTherLoadFunction;
/** @typedef Iterateur sur une std::list de ThermalLoad */
typedef ListTherLoadFunction::iterator ListTherLoadFunctionIter;
/** @typedef Iterateur constant sur une std::list de ThermalLoad */
typedef ListTherLoadFunction::const_iterator ListTherLoadFunctionCIter;

#endif /* THERMALLOAD_H_ */
