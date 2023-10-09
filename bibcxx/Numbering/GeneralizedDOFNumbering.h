#ifndef GENERALIZEDDOFNUMBERING_H_
#define GENERALIZEDDOFNUMBERING_H_

/**
 * @file GeneralizedDOFNumbering.h
 * @brief Fichier entete de la classe GeneralizedDOFNumbering
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

#include "DataStructures/DataStructure.h"
#include "LinearAlgebra/MatrixStorage.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/GeneralizedModel.h"
#include "Numbering/GeneralizedEquationNumbering.h"
#include "Results/ForwardGeneralizedModeResult.h"
#include "Results/ForwardModeResult.h"
#include "Supervis/ResultNaming.h"

/**
 * @class GeneralizedDOFNumbering
 * @brief Cette classe correspond a un sd_nume_ddl_gene
 * @author Nicolas Sellenet
 */
class GeneralizedDOFNumbering : public DataStructure {
  private:
    /** @brief Objet Jeveux '.BASE' */
    JeveuxVectorReal _base;
    /** @brief Objet Jeveux '.NOMS' */
    JeveuxVectorChar8 _noms;
    /** @brief Objet Jeveux '.TAIL' */
    JeveuxVectorLong _tail;
    /** @brief Objet Jeveux '.SMOS' */
    MorseStoragePtr _smos;
    /** @brief Objet Jeveux '.SLCS' */
    LigneDeCielPtr _slcs;
    /** @brief GeneralizedEquationNumbering */
    GeneralizedEquationNumberingPtr _nume;
    /** @brief GeneralizedModel */
    GeneralizedModelPtr _model;
    /** @brief modal basis */
    ForwardModeResultPtr _basis1;
    /** @brief modal basis */
    ForwardGeneralizedModeResultPtr _basis2;

  public:
    /**
     * @typedef GeneralizedDOFNumberingPtr
     * @brief Pointeur intelligent vers un GeneralizedDOFNumbering
     */
    typedef std::shared_ptr< GeneralizedDOFNumbering > GeneralizedDOFNumberingPtr;

    /**
     * @brief Constructeur
     */
    GeneralizedDOFNumbering() : GeneralizedDOFNumbering( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Constructeur
     */
    GeneralizedDOFNumbering( const std::string name )
        : DataStructure( name, 14, "NUME_DDL_GENE" ),
          _base( JeveuxVectorReal( getName() + ".ELIM.BASE" ) ),
          _noms( JeveuxVectorChar8( getName() + ".ELIM.NOMS" ) ),
          _tail( JeveuxVectorLong( getName() + ".ELIM.TAIL" ) ),
          _smos( new MorseStorage( getName() + ".SMOS" ) ),
          _slcs( new LigneDeCiel( getName() + ".SLCS" ) ),
          _nume( new GeneralizedEquationNumbering( getName() + ".NUME" ) ),
          _model( nullptr ),
          _basis1( nullptr ),
          _basis2( nullptr ) {};

    /**
     * @brief Get the GeneralizedModel
     */
    GeneralizedModelPtr getGeneralizedModel() const { return _model; };

    GeneralizedEquationNumberingPtr getDescription() const { return _nume; };

    /**
     * @brief Get modal basis
     */
    GeneralizedModeResultPtr getModalBasisFromGeneralizedModeResult() {
        if ( _basis2.isSet() )
            return _basis2.getPointer();
        return GeneralizedModeResultPtr( nullptr );
    };

    /**
     * @brief Get modal basis
     */
    ModeResultPtr getModalBasisFromModeResult() {
        if ( _basis1.isSet() )
            return _basis1.getPointer();
        return ModeResultPtr( nullptr );
    };

    /**
     * @brief Set the GeneralizedModel
     */
    bool setGeneralizedModel( const GeneralizedModelPtr &model ) {
        _model = model;
        return true;
    };

    /**
     * @brief Set modal basis
     */
    bool setModalBasis( const GeneralizedModeResultPtr &mecaModeC ) {
        if ( mecaModeC != nullptr ) {
            _basis2 = mecaModeC;
            _basis1 = nullptr;
            return true;
        }
        return false;
    };

    /**
     * @brief Set modal basis
     */
    bool setModalBasis( const ModeResultPtr &mecaModeC ) {
        if ( mecaModeC != nullptr ) {
            _basis1 = mecaModeC;
            _basis2 = nullptr;
            return true;
        }
        return false;
    };

    /**
     * @brief Get Morse storage
     */
    MorseStoragePtr getMorseStorage() const { return _smos; };
};

/**
 * @typedef GeneralizedDOFNumberingPtr
 * @brief Pointeur intelligent vers un GeneralizedDOFNumbering
 */
typedef std::shared_ptr< GeneralizedDOFNumbering > GeneralizedDOFNumberingPtr;

#endif /* GENERALIZEDDOFNUMBERING_H_ */
