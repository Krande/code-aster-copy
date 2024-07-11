#ifndef HHOMODEL_H_
#define HHOMODEL_H_

/**
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

#include "astercxx.h"

#include "DataFields/DataField.h"

/** @brief Forward declaration of FieldOnNodes */
template < class ValueType >
class FieldOnNodes;
typedef FieldOnNodes< ASTERDOUBLE > FieldOnNodesReal;
typedef std::shared_ptr< FieldOnNodesReal > FieldOnNodesRealPtr;

/** @brief Forward declaration of FieldOnCells */
template < class ValueType >
class FieldOnCells;
typedef FieldOnCells< ASTERDOUBLE > FieldOnCellsReal;
typedef std::shared_ptr< FieldOnCellsReal > FieldOnCellsRealPtr;

/**
 * @class HHOModel
 * @brief Datastructure for HHOModel
 */
class HHOModel {
  private:
    /** Fields for HHO */

    FieldOnNodesRealPtr _basis;

    FieldOnCellsRealPtr _op_stab;
    FieldOnCellsRealPtr _op_grad;

  public:
    HHOModel( const std::string &modelName );

    FieldOnNodesRealPtr getBasis() const { return _basis; };

    FieldOnCellsRealPtr getGradient() const { return _op_grad; };

    FieldOnCellsRealPtr getStabilization() const { return _op_stab; };
};

typedef std::shared_ptr< HHOModel > HHOModelPtr;

#endif /* HHOMODEL_H_ */
