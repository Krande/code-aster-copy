#ifndef THERMALLOAD_H_
#define THERMALLOAD_H_

/**
 * @file ThermalLoad.h
 * @brief Fichier entete de la classe ThermalLoad
 * @author Jean-Pierre Lefebvre
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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
#include "MemoryManager/JeveuxVector.h"
#include "Loads/ThermalLoadDescription.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ThermalLoadClass
 * @brief Classe definissant une charge thermique (issue d'AFFE_CHAR_THER)
 * @author Jean-Pierre Lefebvre
 */
template< class ConstantFieldOnCellsType>
class ThermalLoadClass : public DataStructure {
  private:
    /** @brief Vecteur Jeveux '.TYPE' */
    JeveuxVectorChar8 _type;
    /** @brief sd_char_chth '.CHTH' */
    ThermalLoadDescriptionPtr< ConstantFieldOnCellsType > _therLoadDesc;

  public:
    /**
     * @typedef ThermalLoadPtr
     * @brief Pointeur intelligent vers un ThermalLoad
     */
    typedef boost::shared_ptr< ThermalLoadClass > ThermalLoadPtr;

    /**
     * @brief Constructeur
     */
    ThermalLoadClass( void ) = delete;

    /**
     * @brief Constructeur
     */
    ThermalLoadClass( const ModelPtr &currentModel )
        : ThermalLoadClass( ResultNaming::getNewResultName(), currentModel ){};

    /**
     * @brief Constructeur
     */
    ThermalLoadClass( const std::string name, const ModelPtr &currentModel )
        : DataStructure( name, 8, "CHAR_THER" ),
        _therLoadDesc( boost::make_shared<
            ThermalLoadDescriptionClass< ConstantFieldOnCellsType > >(getName() + ".CHTH",
                                                                            currentModel )),
          _type( getName() + ".TYPE" ){};


    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const
    { return _therLoadDesc->getFiniteElementDescriptor(); };

    /**
     * @brief Get the model
     */
    ModelPtr getModel() const { return _therLoadDesc->getModel(); };

    /**
     * @brief Get the mesh
     */
    BaseMeshPtr getMesh() const { return _therLoadDesc->getMesh(); };
};

/**********************************************************/
/*  Explicit instantiation of template classes
/**********************************************************/

/** @typedef ThermalLoadRealClass Class d'une charge mécanique réelle */
typedef ThermalLoadClass< ConstantFieldOnCellsRealClass > ThermalLoadRealClass;
/** @typedef ThermalLoadFuncClass Class d'une charge mécanique de fonctions */
typedef ThermalLoadClass< ConstantFieldOnCellsChar24Class > ThermalLoadFunctionClass;

/** @typedef ThermalLoad  */
template< class ConstantFieldOnCellsType>
using ThermalLoadPtr =
    boost::shared_ptr< ThermalLoadClass< ConstantFieldOnCellsType > >;

typedef boost::shared_ptr< ThermalLoadRealClass > ThermalLoadRealPtr;
typedef boost::shared_ptr< ThermalLoadFunctionClass > ThermalLoadFunctionPtr;


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
