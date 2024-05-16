#ifndef DYNAMICMACROELEMENT_H_
#define DYNAMICMACROELEMENT_H_

/**
 * @file DynamicMacroElement.h
 * @brief Fichier entete de la classe DynamicMacroElement
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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
#include "LinearAlgebra/AssemblyMatrix.h"
#include "LinearAlgebra/GeneralizedAssemblyMatrix.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Numbering/DOFNumbering.h"
#include "Results/ModeResult.h"
#include "Supervis/ResultNaming.h"

/**
 * @class DynamicMacroElement
 * @brief Cette classe correspond a un MACR_ELEM_DYNA
 * @author Nicolas Sellenet
 */
class DynamicMacroElement : public DataStructure {
  private:
    /** @brief Objet NUME_DDL */
    DOFNumberingPtr _numeDdl;
    /** @brief Objet Jeveux '.DESM' */
    JeveuxVectorLong _desm;
    /** @brief Objet Jeveux '.REFM' */
    JeveuxVectorChar8 _refm;
    /** @brief Objet Jeveux '.CONX' */
    JeveuxVectorLong _conx;
    /** @brief Objet Jeveux '.LINO' */
    JeveuxVectorLong _lino;
    /** @brief Objet Jeveux '.MAEL_DESC' */
    JeveuxVectorLong _maelDesc;
    /** @brief Objet Jeveux '.MAEL_REFE' */
    JeveuxVectorChar24 _maelRefe;
    /** @brief Objet Jeveux '.LICH' */
    JeveuxCollectionChar8 _lich;
    /** @brief Objet Jeveux '.LICA' */
    JeveuxCollectionReal _lica;
    /** @brief Objet Jeveux '.MAEL_RAID_DESC' */
    JeveuxVectorLong _maelRaidDesc;
    /** @brief Objet Jeveux '.MAEL_RAID_REFE' */
    JeveuxVectorChar24 _maelRaidRefe;
    /** @brief Objet Jeveux '.MAEL_RAID_VALE' */
    JeveuxCollectionReal _maelRaidVale;
    /** @brief Objet Jeveux '.MAEL_MASS_DESC' */
    JeveuxVectorLong _maelMassDesc;
    /** @brief Objet Jeveux '.MAEL_MASS_REFE' */
    JeveuxVectorChar24 _maelMassRefe;
    /** @brief Objet Jeveux '.MAEL_MASS_VALE' */
    JeveuxCollectionReal _maelMassVale;
    /** @brief Objet Jeveux '.MAEL_AMOR_DESC' */
    JeveuxVectorLong _maelAmorDesc;
    /** @brief Objet Jeveux '.MAEL_AMOR_REFE' */
    JeveuxVectorChar24 _maelAmorRefe;
    /** @brief Objet Jeveux '.MAEL_AMOR_VALE' */
    JeveuxCollectionReal _maelAmorVale;
    /** @brief Objet Jeveux '.MAEL_INER_REFE' */
    JeveuxVectorChar24 _maelInerRefe;
    /** @brief Objet Jeveux '.MAEL_INER_VALE' */
    JeveuxVectorReal _maelInterVale;
    /** @brief Mode Meca sur lequel repose le macro emement */
    ModeResultPtr _mechanicalMode;
    /** @brief double rigidity matrix */
    AssemblyMatrixDisplacementRealPtr _rigidityDMatrix;
    /** @brief complex rigidity matrix */
    AssemblyMatrixDisplacementComplexPtr _rigidityCMatrix;
    /** @brief mass matrix */
    AssemblyMatrixDisplacementRealPtr _massMatrix;
    /** @brief damping matrix */
    AssemblyMatrixDisplacementRealPtr _dampingMatrix;
    /** @brief MATR_IMPE */
    GeneralizedAssemblyMatrixComplexPtr _impeMatrix;
    /** @brief MATR_IMPE_RIGI */
    GeneralizedAssemblyMatrixComplexPtr _impeRigiMatrix;
    /** @brief MATR_IMPE_MASS */
    GeneralizedAssemblyMatrixComplexPtr _impeMassMatrix;
    /** @brief MATR_IMPE_AMOR */
    GeneralizedAssemblyMatrixComplexPtr _impeAmorMatrix;

  public:
    /**
     * @typedef DynamicMacroElementPtr
     * @brief Pointeur intelligent vers un DynamicMacroElement
     */
    typedef std::shared_ptr< DynamicMacroElement > DynamicMacroElementPtr;

    /**
     * @brief Constructeur
     */
    DynamicMacroElement() : DynamicMacroElement( ResultNaming::getNewResultName() ) {};

    DynamicMacroElement( const std::string name )
        : DataStructure( name, 8, "MACR_ELEM_DYNA" ),
          _numeDdl( new DOFNumbering( getName() + "      .NUME" ) ),
          _desm( JeveuxVectorLong( getName() + ".DESM" ) ),
          _refm( JeveuxVectorChar8( getName() + ".REFM" ) ),
          _conx( JeveuxVectorLong( getName() + ".CONX" ) ),
          _lino( JeveuxVectorLong( getName() + ".LINO" ) ),
          _maelDesc( JeveuxVectorLong( getName() + ".MAEL_DESC" ) ),
          _maelRefe( JeveuxVectorChar24( getName() + ".MAEL_REFE" ) ),
          _lich( JeveuxCollectionChar8( getName() + ".LICH" ) ),
          _lica( JeveuxCollectionReal( getName() + ".LICA" ) ),
          _maelRaidDesc( JeveuxVectorLong( getName() + ".MAEL_RAID_DESC" ) ),
          _maelRaidRefe( JeveuxVectorChar24( getName() + ".MAEL_RAID_REFE" ) ),
          _maelRaidVale( JeveuxCollectionReal( getName() + ".MAEL_RAID_VALE" ) ),
          _maelMassDesc( JeveuxVectorLong( getName() + ".MAEL_MASS_DESC" ) ),
          _maelMassRefe( JeveuxVectorChar24( getName() + ".MAEL_MASS_REFE" ) ),
          _maelMassVale( JeveuxCollectionReal( getName() + ".MAEL_MASS_VALE" ) ),
          _maelAmorDesc( JeveuxVectorLong( getName() + ".MAEL_AMOR_DESC" ) ),
          _maelAmorRefe( JeveuxVectorChar24( getName() + ".MAEL_AMOR_REFE" ) ),
          _maelAmorVale( JeveuxCollectionReal( getName() + ".MAEL_AMOR_VALE" ) ),
          _maelInerRefe( JeveuxVectorChar24( getName() + ".MAEL_INER_REFE" ) ),
          _maelInterVale( JeveuxVectorReal( getName() + ".MAEL_INER_VALE" ) ),
          _mechanicalMode( nullptr ),
          _rigidityDMatrix( nullptr ),
          _rigidityCMatrix( nullptr ),
          _massMatrix( nullptr ),
          _dampingMatrix( nullptr ),
          _impeMatrix( nullptr ),
          _impeRigiMatrix( nullptr ),
          _impeMassMatrix( nullptr ),
          _impeAmorMatrix( nullptr ) {};

    /**
     * @brief Get damping matrix
     */
    AssemblyMatrixDisplacementRealPtr getDampingMatrix() { return _dampingMatrix; };

    DOFNumberingPtr getDOFNumbering() const { return _numeDdl; };

    /**
     * @brief Get impedance matrix
     */
    GeneralizedAssemblyMatrixComplexPtr getImpedanceDampingMatrix() { return _impeAmorMatrix; };

    /**
     * @brief Get impedance matrix
     */
    GeneralizedAssemblyMatrixComplexPtr getImpedanceMatrix() { return _impeMatrix; };

    /**
     * @brief Get impedance matrix
     */
    GeneralizedAssemblyMatrixComplexPtr getImpedanceMassMatrix() { return _impeMassMatrix; };

    /**
     * @brief Get impedance matrix
     */
    GeneralizedAssemblyMatrixComplexPtr getImpedanceStiffnessMatrix() { return _impeRigiMatrix; };

    /**
     * @brief Get mass matrix
     */
    AssemblyMatrixDisplacementRealPtr getMassMatrix() { return _massMatrix; };

    /**
     * @brief Get size of object .LINO
     */
    int getNumberOfNodes() { return _lino->size(); };

    /**
     * @brief Get rigidity matrix
     */
    AssemblyMatrixDisplacementComplexPtr getStiffnessMatrixComplex() { return _rigidityCMatrix; };

    /**
     * @brief Get rigidity matrix
     */
    AssemblyMatrixDisplacementRealPtr getStiffnessMatrixReal() { return _rigidityDMatrix; };

    /**
     * @brief Get mechanical mode
     */
    ModeResultPtr getMechanicalMode() { return _mechanicalMode; };

    /**
     * @brief Set damping matrix
     * @param matrix AssemblyMatrixDisplacementRealPtr object
     */
    bool setDampingMatrix( const AssemblyMatrixDisplacementRealPtr &matrix ) {
        _dampingMatrix = matrix;
        return true;
    };

    /**
     * @brief Set impedance matrix
     * @param matrix AssemblyMatrixDisplacementRealPtr object
     */
    bool setImpedanceDampingMatrix( const GeneralizedAssemblyMatrixComplexPtr &matrix ) {
        _impeAmorMatrix = matrix;
        return true;
    };

    /**
     * @brief Set impedance matrix
     * @param matrix AssemblyMatrixDisplacementRealPtr object
     */
    bool setImpedanceMatrix( const GeneralizedAssemblyMatrixComplexPtr &matrix ) {
        _impeMatrix = matrix;
        return true;
    };

    /**
     * @brief Set impedance matrix
     * @param matrix AssemblyMatrixDisplacementRealPtr object
     */
    bool setImpedanceMassMatrix( const GeneralizedAssemblyMatrixComplexPtr &matrix ) {
        _impeMassMatrix = matrix;
        return true;
    };

    /**
     * @brief Set impedance matrix
     * @param matrix AssemblyMatrixDisplacementRealPtr object
     */
    bool setImpedanceStiffnessMatrix( const GeneralizedAssemblyMatrixComplexPtr &matrix ) {
        _impeRigiMatrix = matrix;
        return true;
    };

    /**
     * @brief Set mass matrix
     * @param matrix AssemblyMatrixDisplacementRealPtr object
     */
    bool setMassMatrix( const AssemblyMatrixDisplacementRealPtr &matrix ) {
        _massMatrix = matrix;
        return true;
    };

    /**
     * @brief Set the ModeResult
     * @param mechanicalMode ModeResultPtr object
     */
    bool setMechanicalMode( ModeResultPtr &mechanicalMode ) {
        _mechanicalMode = mechanicalMode;
        return true;
    };

    /**
     * @brief Set rigidity matrix
     * @param matrix AssemblyMatrixDisplacementComplexPtr object
     */
    bool setStiffnessMatrix( const AssemblyMatrixDisplacementComplexPtr &matrix ) {
        _rigidityCMatrix = matrix;
        _rigidityDMatrix = nullptr;
        return true;
    };

    /**
     * @brief Set rigidity matrix
     * @param matrix AssemblyMatrixDisplacementRealPtr object
     */
    bool setStiffnessMatrix( const AssemblyMatrixDisplacementRealPtr &matrix ) {
        _rigidityDMatrix = matrix;
        _rigidityCMatrix = nullptr;
        return true;
    };

    JeveuxCollectionReal getGeneralizedStiffnessMatrix() { return _maelRaidVale; }

    JeveuxCollectionReal getGeneralizedMassMatrix() { return _maelMassVale; }

    JeveuxCollectionReal getGeneralizedDampingMatrix() { return _maelAmorVale; }
};

/**
 * @typedef DynamicMacroElementPtr
 * @brief Pointeur intelligent vers un DynamicMacroElement
 */
typedef std::shared_ptr< DynamicMacroElement > DynamicMacroElementPtr;

#endif /* DYNAMICMACROELEMENT_H_ */
