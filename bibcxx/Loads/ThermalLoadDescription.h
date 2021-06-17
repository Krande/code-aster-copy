#ifndef THERMALLOADDESCRIPTION_H_
#define THERMALLOADDESCRIPTION_H_

/**
 * @file ThermalLoad.h
 * @brief Fichier entete de la classe ThermalLoad
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
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ThermalLoadDescription
 * @brief Classe definissant une charge thermique sd_char_chth
 * @author Jean-Pierre Lefebvre
 */
template< typename ConstantFieldOnCellsType>
class ThermalLoadDescription : public DataStructure {

  public:
    typedef boost::shared_ptr< ConstantFieldOnCellsType > ConstantFieldOnCellsTypePtr;

  private:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Vecteur Jeveux '.MODEL.NOMO' */
    JeveuxVectorChar8 _modelName;
    /** @brief Vecteur Jeveux '.CONVE.VALE' */
    JeveuxVectorChar8 _convection;
    /** @brief Vecteur Jeveux '.LIGRE' */
    FiniteElementDescriptorPtr _FEDesc;
    /** @brief Carte '.CIMPO' */
    ConstantFieldOnCellsTypePtr _cimpo;
    /** @brief Carte '.CMULT' */
    ConstantFieldOnCellsRealPtr _cmult;
    /** @brief Carte '.COEFH' */
    ConstantFieldOnCellsTypePtr _coefh;
    /** @brief Carte '.FLUNL' */
    ConstantFieldOnCellsTypePtr _flunl;
    /** @brief Carte '.FLURE' */
    ConstantFieldOnCellsTypePtr _flure;
    /** @brief Carte '.GRAIN' */
    ConstantFieldOnCellsTypePtr _grain;
    /** @brief Carte '.HECHP' */
    ConstantFieldOnCellsTypePtr _hechp;
    /** @brief Carte '.SOURE' */
    ConstantFieldOnCellsTypePtr _soure;
    /** @brief Carte '.T_EXT' */
    ConstantFieldOnCellsTypePtr _tExt;

  public:

    ThermalLoadDescription( void ) = delete;

        /** @brief Constructeur */
    ThermalLoadDescription( const std::string &name, const ModelPtr &currentModel ):
        DataStructure( name, 13, "CHAR_CHTH" ),
        _model( currentModel ),
        _modelName( getName() + ".MODEL.NOMO" ), _convection( getName() + ".CONVE.VALE" ),
        _FEDesc( boost::make_shared< FiniteElementDescriptor >( getName() + ".LIGRE",
                                                                    _model->getMesh() ) ),
        _cimpo( boost::make_shared< ConstantFieldOnCellsType >( getName() + ".CIMPO", _FEDesc ) ),
        _cmult( boost::make_shared< ConstantFieldOnCellsReal >(
                                                                getName() + ".CMULT", _FEDesc ) ),
        _coefh( boost::make_shared< ConstantFieldOnCellsType >( getName() + ".COEFH", _FEDesc ) ),
        _flunl( boost::make_shared< ConstantFieldOnCellsType >( getName() + ".FLUNL", _FEDesc ) ),
        _flure( boost::make_shared< ConstantFieldOnCellsType >( getName() + ".FLURE", _FEDesc ) ),
        _grain( boost::make_shared< ConstantFieldOnCellsType >( getName() + ".GRAIN", _FEDesc ) ),
        _hechp( boost::make_shared< ConstantFieldOnCellsType >( getName() + ".HECHP", _FEDesc ) ),
        _soure( boost::make_shared< ConstantFieldOnCellsType >( getName() + ".SOURE", _FEDesc ) ),
        _tExt( boost::make_shared< ConstantFieldOnCellsType >( getName() + ".T_EXT", _FEDesc ) ){};

    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    /**
     * @brief Get the model
     */
    ModelPtr getModel() const { return _model; };

    /**
     * @brief Get the mesh
     */
    BaseMeshPtr getMesh() const { return _model->getMesh(); };
};


/**********************************************************/
/*  Explicit instantiation of template classes
/**********************************************************/

/** @typedef ThermalLoadDescriptionReal Class d'une charge mécanique réelle */
typedef ThermalLoadDescription< ConstantFieldOnCellsReal >
    ThermalLoadDescriptionReal;
/** @typedef ThermalLoadDescriptionFunc Class d'une charge mécanique de fonctions */
typedef ThermalLoadDescription< ConstantFieldOnCellsChar24 >
    ThermalLoadDescriptionFunction;


template< typename ConstantFieldOnCellsType >
using ThermalLoadDescriptionPtr =
    boost::shared_ptr< ThermalLoadDescription< ConstantFieldOnCellsType > >;

#endif /* THERMALLOADDESCRIPTION_H_ */
