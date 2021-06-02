#ifndef MECHANICALLOAD_H_
#define MECHANICALLOAD_H_

/**
 * @file MechanicalLoad.h
 * @author Natacha Bereux
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
#include "DataFields/ListOfTables.h"
#include "DataStructures/DataStructure.h"
#include "Loads/MechanicalLoadDescription.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"


/**
 * @class MechanicalLoadClass
 * @brief Define a generic mechanical load
 * @author Nicolas Sellenet
 */
template< class ConstantFieldOnCellsType>
class MechanicalLoadClass : public DataStructure, public ListOfTablesClass {

  protected:
    /** @typedef Definition d'un pointeur intelligent sur un VirtualMeshEntity */
    typedef boost::shared_ptr< VirtualMeshEntity > MeshEntityPtr;

    /** @brief MeshEntity sur laquelle repose le "blocage" */
    MeshEntityPtr _meshEntity;
    /** @brief Vecteur Jeveux '.TYPE' */
    JeveuxVectorChar8 _type;
    /** @brief Vecteur Jeveux '.LISMA01' */
    JeveuxVectorLong _lisma01;
    /** @brief Vecteur Jeveux '.LISMA02' */
    JeveuxVectorLong _lisma02;
    /** @brief Vecteur Jeveux '.TRANS01' */
    JeveuxVectorReal _trans01;
    /** @brief Vecteur Jeveux '.TRANS02' */
    JeveuxVectorReal _trans02;
    /** @brief Vecteur Jeveux '.POIDS_MAILLE' */
    JeveuxVectorReal _poidsMaille;
    /** @brief sd_char_chme '.CHME' */
    MechanicalLoadDescriptionClassPtr< ConstantFieldOnCellsType > _mecaLoadDesc;


  public:
    /**
     * @typedef MechanicalLoadPtr
     * @brief Pointeur intelligent vers un MechanicalLoad
     */
    typedef boost::shared_ptr< MechanicalLoadClass > MechanicalLoadPtr;

    /**
     * @brief Constructor
     */
    MechanicalLoadClass( void ) = delete;

    /**
     * @brief Constructor
     */
    MechanicalLoadClass( const ModelPtr &currentModel )
        : MechanicalLoadClass( ResultNaming::getNewResultName(), currentModel ){};

    /**
     * @brief Constructor
     */
    MechanicalLoadClass( const std::string name, const ModelPtr &currentModel )
        : DataStructure( name, 8, "CHAR_MECA" ),
          ListOfTablesClass( name ),
            _mecaLoadDesc( boost::make_shared<
                MechanicalLoadDescriptionClass< ConstantFieldOnCellsType > >
                    (getName() + ".CHME", currentModel) ),
            _type( getName() + ".TYPE" ), _lisma01( getName() + ".LISMA01" ),
            _lisma02( getName() + ".LISMA02" ), _trans01( getName() + ".TRANS01" ),
            _trans02( getName() + ".TRANS02" ), _poidsMaille( getName() + ".POIDS_MAILLE" ){};

    /**
     * @brief Get the model
     */
    const MechanicalLoadDescriptionClassPtr< ConstantFieldOnCellsType >&
    getMechanicalLoadDescription() const { return _mecaLoadDesc; };

    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const
    { return _mecaLoadDesc->getFiniteElementDescriptor(); };

    /**
     * @brief Get the model
     */
    ModelPtr getModel() const {
        return _mecaLoadDesc->getModel();
    };

    /**
     * @brief Get the model
     */
    BaseMeshPtr getMesh() const {
        return _mecaLoadDesc->getMesh();
    };

    JeveuxVectorChar8 getType() const
    {
        return _type;
    }
};

/**********************************************************/
/*  Explicit instantiation of template classes
/**********************************************************/

/** @typedef MechanicalLoadRealClass Class d'une charge mécanique réelle */
typedef MechanicalLoadClass< ConstantFieldOnCellsRealClass > MechanicalLoadRealClass;
/** @typedef MechanicalLoadFuncClass Class d'une charge mécanique de fonctions */
typedef MechanicalLoadClass< ConstantFieldOnCellsChar24Class > MechanicalLoadFunctionClass;
/** @typedef MechanicalLoadComplexClass Class d'une charge mécanique de complexe */
typedef MechanicalLoadClass< ConstantFieldOnCellsComplexClass > MechanicalLoadComplexClass;

/** @typedef MechanicalLoad  */
template< class ConstantFieldOnCellsType>
using MechanicalLoadPtr =
    boost::shared_ptr< MechanicalLoadClass< ConstantFieldOnCellsType > >;

typedef boost::shared_ptr< MechanicalLoadRealClass > MechanicalLoadRealPtr;
typedef boost::shared_ptr< MechanicalLoadFunctionClass > MechanicalLoadFunctionPtr;
typedef boost::shared_ptr< MechanicalLoadComplexClass > MechanicalLoadComplexPtr;


/** @typedef std::list de MechanicalLoad */
typedef std::list< MechanicalLoadRealPtr > ListMecaLoadReal;
/** @typedef Iterateur sur une std::list de MechanicalLoad */
typedef ListMecaLoadReal::iterator ListMecaLoadRealIter;
/** @typedef Iterateur constant sur une std::list de MechanicalLoad */
typedef ListMecaLoadReal::const_iterator ListMecaLoadRealCIter;

/** @typedef std::list de MechanicalLoad */
typedef std::list< MechanicalLoadFunctionPtr > ListMecaLoadFunction;
/** @typedef Iterateur sur une std::list de MechanicalLoad */
typedef ListMecaLoadFunction::iterator ListMecaLoadFunctionIter;
/** @typedef Iterateur constant sur une std::list de MechanicalLoad */
typedef ListMecaLoadFunction::const_iterator ListMecaLoadFunctionCIter;

#endif /* MECHANICALLOAD_H_ */
