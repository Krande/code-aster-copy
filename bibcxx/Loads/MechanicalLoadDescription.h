#ifndef MECHANICALLOADDESCRIPTION_H_
#define MECHANICALLOADDESCRIPTION_H_

/**
 * @file MechanicalLoadDescription.h
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

#include <string>

#include "astercxx.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"


/**
 * @class MechanicalLoadDescription
 * @brief Encapsulation of sd_char_chme
 */
template< typename ConstantFieldOnCellsType>
class MechanicalLoadDescription : public DataStructure
{

  public:
    typedef boost::shared_ptr< ConstantFieldOnCellsType > ConstantFieldOnCellsTypePtr;

  private:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Vecteur Jeveux '.TEMPE.TEMP' */
    JeveuxVectorChar8 _temperatureField;
    /** @brief Vecteur Jeveux '.MODEL.NOMO' */
    JeveuxVectorChar8 _modelName;
    /** @brief Vecteur Jeveux '.VEASS' */
    JeveuxVectorChar8 _nameOfAssemblyVector;
    /** @brief Vecteur Jeveux '.VEISS' */
    JeveuxVectorChar8 _veiss;
    /** @brief Vecteur Jeveux '.EVOL.CHAR' */
    JeveuxVectorChar8 _evolChar;
    /** @brief Vecteur Jeveux '.LIGRE' */
    FiniteElementDescriptorPtr _FEDesc;
    /** @brief Carte '.CIMPO' */
    ConstantFieldOnCellsTypePtr _cimpo;
    /** @brief Carte '.CMULT' */
    ConstantFieldOnCellsRealPtr _cmult;
    /** @brief Carte '.DPGEN' */
    ConstantFieldOnCellsTypePtr _dpgen;
    /** @brief Carte '.EPSIN' */
    ConstantFieldOnCellsTypePtr _epsin;
    /** @brief Carte '.F1D2D' */
    ConstantFieldOnCellsTypePtr _f1d2d;
    /** @brief Carte '.F1D3D' */
    ConstantFieldOnCellsTypePtr _f1d3d;
    /** @brief Carte '.F2D3D' */
    ConstantFieldOnCellsTypePtr _f2d3d;
    /** @brief Carte '.FCO2D' */
    ConstantFieldOnCellsTypePtr _fco2d;
    /** @brief Carte '.FCO3D' */
    ConstantFieldOnCellsTypePtr _fco3d;
    /** @brief Carte '.FELEC' */
    ConstantFieldOnCellsTypePtr _felec;
    /** @brief Carte '.FL101' */
    ConstantFieldOnCellsTypePtr _fl101;
    /** @brief Carte '.FL102' */
    ConstantFieldOnCellsTypePtr _fl102;
    /** @brief Carte '.FORNO' */
    ConstantFieldOnCellsTypePtr _forno;
    /** @brief Carte '.IMPED' */
    ConstantFieldOnCellsTypePtr _imped;
    /** @brief Carte '.PESAN' */
    ConstantFieldOnCellsTypePtr _pesan;
    /** @brief Carte '.PRESS' */
    ConstantFieldOnCellsTypePtr _press;
    /** @brief Carte '.ROTAT' */
    ConstantFieldOnCellsTypePtr _rotat;
    /** @brief Carte '.SIGIN' */
    ConstantFieldOnCellsTypePtr _sigin;
    /** @brief Carte '.SIINT' */
    ConstantFieldOnCellsTypePtr _siint;
    /** @brief Carte '.VNOR' */
    ConstantFieldOnCellsTypePtr _vnor;
    /** @brief Carte '.ONDPL' */
    ConstantFieldOnCellsTypePtr _ondpl;
    /** @brief Carte '.ONDPR' */
    ConstantFieldOnCellsTypePtr _ondpr;

  public:

    MechanicalLoadDescription( void ) = delete;

    /**
     * @brief Constructor
     */
    MechanicalLoadDescription( const std::string name, const ModelPtr &currentModel )
        : DataStructure( name, 13, "CHAR_CHME" ),
            _model( currentModel ),
            _temperatureField( getName() + ".TEMPE.TEMP" ),
            _modelName( getName() + ".MODEL.NOMO" ),
            _nameOfAssemblyVector( getName() + ".VEASS" ),
            _veiss( getName() + ".VEISS" ),
            _evolChar( getName() + ".EVOL.CHAR" ),
            _FEDesc( boost::make_shared<FiniteElementDescriptor>(
                getName() + ".LIGRE", _model->getMesh() ) ),
            _cimpo(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".CIMPO", _FEDesc ) ),
            _cmult(boost::make_shared<ConstantFieldOnCellsReal>(
                                                                 getName() + ".CMULT", _FEDesc ) ),
            _dpgen(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".DPGEN", _FEDesc ) ),
            _epsin(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".EPSIN", _FEDesc ) ),
            _f1d2d(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".F1D2D", _FEDesc ) ),
            _f1d3d(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".F1D3D", _FEDesc ) ),
            _f2d3d(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".F2D3D", _FEDesc ) ),
            _fco2d(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".FCO2D", _FEDesc ) ),
            _fco3d(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".FCO3D", _FEDesc ) ),
            _felec(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".FELEC", _FEDesc ) ),
            _fl101(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".FL101", _FEDesc ) ),
            _fl102(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".FL102", _FEDesc ) ),
            _forno(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".FORNO", _FEDesc ) ),
            _imped(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".IMPED", _FEDesc ) ),
            _pesan(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".PESAN", _FEDesc ) ),
            _press(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".PRESS", _FEDesc ) ),
            _rotat(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".ROTAT", _FEDesc ) ),
            _sigin(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".SIGIN", _FEDesc ) ),
            _siint(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".SIINT", _FEDesc ) ),
            _vnor(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".VNOR", _FEDesc ) ),
            _ondpl(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".ONDPL", _FEDesc ) ),
            _ondpr(boost::make_shared<ConstantFieldOnCellsType>( getName() + ".ONDPR", _FEDesc ) )
            {};

    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    BaseMeshPtr getMesh() const { return _model->getMesh(); }

    ModelPtr getModel() const { return _model; }

    ConstantFieldOnCellsTypePtr getImpositionField() const { return _cimpo; }
    ConstantFieldOnCellsRealPtr getMultiplicativeField() const { return _cmult; }

    bool hasLoad( const std::string load_name) const
    {
        if( load_name == "IMPE_FACE")
        {
            _imped->updateValuePointers();
            return _imped->exists();
        }
        else
        {
            AS_ASSERT(false);
        }

        return false;
    }

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return true si la mise a jour s'est bien deroulee, false sinon
     */
    bool updateValuePointers() {
        bool retour = _temperatureField->updateValuePointer();
        retour = ( retour && _modelName->updateValuePointer() );
        retour = ( retour && _nameOfAssemblyVector->updateValuePointer() );
        retour = ( retour && _veiss->updateValuePointer() );
        retour = ( retour && _evolChar->updateValuePointer() );
        retour = ( retour && _cimpo->updateValuePointers() );
        retour = ( retour && _cmult->updateValuePointers() );
        retour = ( retour && _dpgen->updateValuePointers() );
        retour = ( retour && _epsin->updateValuePointers() );
        retour = ( retour && _f1d2d->updateValuePointers() );
        retour = ( retour && _f1d3d->updateValuePointers() );
        retour = ( retour && _f2d3d->updateValuePointers() );
        retour = ( retour && _fco2d->updateValuePointers() );
        retour = ( retour && _fco3d->updateValuePointers() );
        retour = ( retour && _felec->updateValuePointers() );
        retour = ( retour && _fl101->updateValuePointers() );
        retour = ( retour && _fl102->updateValuePointers() );
        retour = ( retour && _forno->updateValuePointers() );
        retour = ( retour && _imped->updateValuePointers() );
        retour = ( retour && _pesan->updateValuePointers() );
        retour = ( retour && _press->updateValuePointers() );
        retour = ( retour && _sigin->updateValuePointers() );
        retour = ( retour && _rotat->updateValuePointers() );
        retour = ( retour && _siint->updateValuePointers() );
        retour = ( retour && _vnor->updateValuePointers() );
        retour = ( retour && _ondpl->updateValuePointers() );
        retour = ( retour && _ondpr->updateValuePointers() );
        return retour;
    };
};

/**********************************************************/
/*  Explicit instantiation of template classes
/**********************************************************/

/** @typedef MechanicalLoadDescriptionReal Class d'une charge mécanique réelle */
typedef MechanicalLoadDescription< ConstantFieldOnCellsReal >
    MechanicalLoadDescriptionReal;
/** @typedef MechanicalLoadDescriptionFunc Class d'une charge mécanique de fonctions */
typedef MechanicalLoadDescription< ConstantFieldOnCellsChar24 >
MechanicalLoadDescriptionFunction;
/** @typedef MechanicalLoadDescriptionComplex Class d'une charge mécanique de complexe */
typedef MechanicalLoadDescription< ConstantFieldOnCellsComplex >
        MechanicalLoadDescriptionComplex;


template< typename ConstantFieldOnCellsType >
using MechanicalLoadDescriptionPtr =
    boost::shared_ptr< MechanicalLoadDescription< ConstantFieldOnCellsType > >;



#endif /* MECHANICALLOADDESCRIPTION_H_ */
