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
 * @class MechanicalLoad
 * @brief Define a generic mechanical load
 * @author Nicolas Sellenet
 */
template< class ConstantFieldOnCellsType>
class MechanicalLoad : public DataStructure, public ListOfTables {

  protected:
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
    MechanicalLoadDescriptionPtr< ConstantFieldOnCellsType > _mecaLoadDesc;

  public:
    /**
     * @typedef MechanicalLoadPtr
     * @brief Pointeur intelligent vers un MechanicalLoad
     */
    typedef boost::shared_ptr< MechanicalLoad > MechanicalLoadPtr;

    /**
     * @brief Constructor
     */
    MechanicalLoad( void ) = delete;

    /**
     * @brief Constructor
     */
    MechanicalLoad( const ModelPtr &currentModel )
        : MechanicalLoad( ResultNaming::getNewResultName(), currentModel ){};

    /**
     * @brief Constructor
     */
    MechanicalLoad( const std::string name, const ModelPtr &currentModel )
        : DataStructure( name, 8, "CHAR_MECA" ),
          ListOfTables( name ),
            _mecaLoadDesc( boost::make_shared<
                MechanicalLoadDescription< ConstantFieldOnCellsType > >
                    (getName() + ".CHME", currentModel) ),
            _type( getName() + ".TYPE" ), _lisma01( getName() + ".LISMA01" ),
            _lisma02( getName() + ".LISMA02" ), _trans01( getName() + ".TRANS01" ),
            _trans02( getName() + ".TRANS02" ), _poidsMaille( getName() + ".POIDS_MAILLE" ){};

    /**
     * @brief Get the model
     */
    const MechanicalLoadDescriptionPtr< ConstantFieldOnCellsType >&
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

    bool hasLoad(const std::string& load_name) const
    {
        return _mecaLoadDesc->hasLoad(load_name);
    }

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return true si la mise a jour s'est bien deroulee, false sinon
     */
    bool updateValuePointers() {
        bool retour = _mecaLoadDesc->updateValuePointers();
        retour = ( retour && _type->updateValuePointer() );
        retour = ( retour && _lisma01->updateValuePointer() );
        retour = ( retour && _lisma02->updateValuePointer() );
        retour = ( retour && _trans01->updateValuePointer() );
        retour = ( retour && _trans02->updateValuePointer() );
        retour = ( retour && _poidsMaille->updateValuePointer() );

        return retour;
    };
};

/**********************************************************/
/*  Explicit instantiation of template classes
/**********************************************************/

/** @typedef MechanicalLoadReal Class d'une charge mécanique réelle */
typedef MechanicalLoad< ConstantFieldOnCellsReal > MechanicalLoadReal;
/** @typedef MechanicalLoadFunc Class d'une charge mécanique de fonctions */
typedef MechanicalLoad< ConstantFieldOnCellsChar24 > MechanicalLoadFunction;
/** @typedef MechanicalLoadComplex Class d'une charge mécanique de complexe */
typedef MechanicalLoad< ConstantFieldOnCellsComplex > MechanicalLoadComplex;

/** @typedef MechanicalLoad  */
template< class ConstantFieldOnCellsType>
using MechanicalLoadPtr =
    boost::shared_ptr< MechanicalLoad< ConstantFieldOnCellsType > >;

typedef boost::shared_ptr< MechanicalLoadReal > MechanicalLoadRealPtr;
typedef boost::shared_ptr< MechanicalLoadFunction > MechanicalLoadFunctionPtr;
typedef boost::shared_ptr< MechanicalLoadComplex > MechanicalLoadComplexPtr;


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
