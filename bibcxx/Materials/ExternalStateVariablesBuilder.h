#ifndef CALCULATIONEXTERNALVARIABLE_H_
#define CALCULATIONEXTERNALVARIABLE_H_

/**
 * @file BaseExternalStateVariables.h
 * @brief Fichier entete de la classe ExternalStateVariables
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Materials/CodedMaterial.h"
#include "Materials/MaterialField.h"
#include "Modeling/Model.h"
#include "Numbering/DOFNumbering.h"

/**
 * @class ExternalStateVariable
 * @brief Compute external state variables (AFFE_VARC)
 */
class ExternalStateVariablesBuilder : public DataStructure {
  private:
    /** @brief Model */
    ModelPtr _model;

    /** @brief Field of material parameters */
    MaterialFieldPtr _materialField;

    /** @brief Elementary characteristics */
    ElementaryCharacteristicsPtr _elemChara;

    /** @brief Coded material */
    CodedMaterialPtr _codedMaterial;

    /** @brief Nodal field of references values for external state variables */
    FieldOnCellsRealPtr _varcRefe;

    /** @brief Nodal field of values for external state variables */
    FieldOnCellsRealPtr _varcVale;

    /** @brief Field of time */
    ConstantFieldOnCellsRealPtr _timeField;

    /** @brief Time */
    ASTERDOUBLE _time;

    /** @brief Flags to detect external state variables */
    bool _pTot;
    bool _hydr;
    bool _sech;
    bool _temp;

  public:
    /** @typedef ExternalStateVariablesPtr */
    typedef boost::shared_ptr< ExternalStateVariablesBuilder > ExternalStateVariablesBuilderPtr;

    /** @brief Constructor */
    ExternalStateVariablesBuilder();

    /** @brief Constructor */
    ExternalStateVariablesBuilder( const ModelPtr &currModel,
                                   const MaterialFieldPtr &currMaterialField,
                                   const ElementaryCharacteristicsPtr &currElemChara,
                                   const CodedMaterialPtr &currCodedMaterial );

    /** @brief Destructor */
    ~ExternalStateVariablesBuilder(){};

    /**
     * @brief Compute external state variables at a given time
     * @param time Time
     */
    void build( const ASTERDOUBLE &time );

    /**
     * @brief Compute nodal field for external state variables
     * @param dofNUM DOF numbering
     * @return Nodal field for external state variables
     */
    FieldOnNodesRealPtr computeExternalStateVariablesLoad( const BaseDOFNumberingPtr &dofNUM );

    /** @brief Detect some external state variables */
    bool hasExternalStateVariables() { return _pTot || _hydr || _sech || _temp; };

    /**
     * @brief Get nodal field for external state variables
     * @return Nodal field for external state variables
     */
    FieldOnCellsRealPtr getExternalStateVariablesField() { return _varcVale; }
};

/** @typedef ExternalStateVariablesPtr*/
typedef boost::shared_ptr< ExternalStateVariablesBuilder > ExternalStateVariablesBuilderPtr;

#endif /* CALCULATIONEXTERNALVARIABLE_H_ */
