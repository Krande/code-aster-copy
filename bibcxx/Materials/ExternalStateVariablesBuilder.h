#ifndef CALCULATIONEXTERNALVARIABLE_H_
#define CALCULATIONEXTERNALVARIABLE_H_

/**
 * @file BaseExternalStateVariables.h
 * @brief Fichier entete de la classe ExternalStateVariables
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"
#include "Materials/MaterialField.h"
#include "Materials/CodedMaterial.h"
#include "Modeling/Model.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Numbering/DOFNumbering.h"

/**
 * @class ExternalStateVariable
 * @brief Calculation Input Variables
 * @author Nicolas Sellenet
 */
class ExternalStateVariablesBuilder : public DataStructure {
  private:
    ModelPtr _model;
    MaterialFieldPtr _mater;
    CodedMaterialPtr _codMater;
    ElementaryCharacteristicsPtr _elemCara;
    FieldOnCellsRealPtr _varRef;
    FieldOnCellsRealPtr _varInst;
    ConstantFieldOnCellsRealPtr _timeValue;
    ASTERDOUBLE _currentTime;
    bool _pTot;
    bool _hydr;
    bool _sech;
    bool _temp;

  public:
    /**
     * @typedef ExternalStateVariablesPtr
     * @brief Pointeur intelligent vers un ExternalStateVariables
     */
    typedef boost::shared_ptr< ExternalStateVariablesBuilder > ExternalStateVariablesBuilderPtr;

    /**
     * @brief Constructeur
     */
    ExternalStateVariablesBuilder( );

    /**
     * @brief Constructeur
     */
    ExternalStateVariablesBuilder( const ModelPtr &model, const MaterialFieldPtr &mater,
                                       const ElementaryCharacteristicsPtr &cara,
                                       const CodedMaterialPtr &codMater );

    /**
     * @brief Destructeur
     */
    ~ExternalStateVariablesBuilder() { return; };

    /**
     * @brief Compute Input Variables at a given time
     */
    void build( const ASTERDOUBLE &time );

    /**
     * @brief Compute Loads after computing of input variables
     */
    FieldOnNodesRealPtr computeExternalStateVariablesLoad( const BaseDOFNumberingPtr &dofNUM );

    bool hasExternalStateVariables() { return _pTot || _hydr || _sech || _temp; };
};

/**
 * @typedef ExternalStateVariablesPtr
 * @brief Pointeur intelligent vers un ExternalStateVariable
 */
typedef boost::shared_ptr< ExternalStateVariablesBuilder > ExternalStateVariablesBuilderPtr;

#endif /* CALCULATIONEXTERNALVARIABLE_H_ */
