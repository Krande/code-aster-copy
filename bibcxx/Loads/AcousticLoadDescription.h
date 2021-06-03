/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#ifndef ACOUSTICLOADDESCRIPTION_H_
#define ACOUSTICLOADDESCRIPTION_H_

/**
 * @file AcousticLoad.h
 * @brief Fichier entete de la classe AcousticLoad
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, eiAcou version 3 of the License, or
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


template< typename ConstantFieldOnCellsType>
class AcousticLoadDescriptionClass : public DataStructure {

  public:
    typedef boost::shared_ptr< ConstantFieldOnCellsType > ConstantFieldOnCellsTypePtr;

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
    /** @brief Carte '.CIMPED' */
    ConstantFieldOnCellsTypePtr _impedanceValues;
    /** @brief Carte '.VITFA' */
    ConstantFieldOnCellsTypePtr _speedValues;


  public:
    /**
     * @typedef AcousticLoadPtr
     * @brief Pointeur intelligent vers un AcousticLoad
     */
    typedef boost::shared_ptr< AcousticLoadDescriptionClass > AcousticLoadDescriptionPtr;

    /**
     * @brief Constructeur
     */
    AcousticLoadDescriptionClass( void ) = delete;

    /**
     * @brief Constructeur
     */
    AcousticLoadDescriptionClass( const ModelPtr &model )
        : AcousticLoadDescriptionClass( ResultNaming::getNewResultName(), model ){};

    /**
     * @brief Constructeur
     */
    AcousticLoadDescriptionClass( const std::string name, const ModelPtr &model )
        : DataStructure( name, 14, "CHAR_CHAC" ), _model( model ),
          _FEDesc( boost::make_shared< FiniteElementDescriptorClass >( getName() + ".LIGRE",
                                                                model->getMesh() ) ),
          _modelName( JeveuxVectorChar8( getName() + ".MODEL.NOMO" ) ),
          _imposedValues( boost::make_shared< ConstantFieldOnCellsType >(
                        getName() + ".CIMPO", _FEDesc ) ),
          _multiplier( boost::make_shared< ConstantFieldOnCellsComplexClass >(
                        getName() + ".CMULT", _FEDesc ) ),
          _impedanceValues( boost::make_shared< ConstantFieldOnCellsType >(
                        getName() + ".IMPED", _FEDesc ) ),
          _speedValues( boost::make_shared< ConstantFieldOnCellsType >(
                        getName() + ".VITFA", _FEDesc ) ){};

    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    /**
     * @brief Get the model
     */
    ModelPtr getModel() { return _model; };

    /**
     * @brief Get the model
     */
    BaseMeshPtr getMesh() { return _model->getMesh(); };
};

/**********************************************************/
/*  Explicit instantiation of template classes
/**********************************************************/


/** @typedef AcousticLoadDescriptionFuncClass Class d'une charge mécanique de fonctions */
typedef AcousticLoadDescriptionClass< ConstantFieldOnCellsComplexClass >
   AcousticLoadDescriptionComplexClass;


template< typename ConstantFieldOnCellsType >
using AcousticLoadDescriptionPtr =
    boost::shared_ptr<AcousticLoadDescriptionClass< ConstantFieldOnCellsType > >;


#endif /* ACOUSTICLOADDESCRIPTION_H_ */
