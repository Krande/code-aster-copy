#ifndef XFEMMODEL_H_
#define XFEMMODEL_H_

/**
 * @file XfemModel.h
 * @brief Header for class XfemModel
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

#include "astercxx.h"

#include "DataFields/DataField.h"

/** @brief Forward declaration of FieldOnNodes */
template < class ValueType >
class FieldOnNodes;
typedef FieldOnNodes< ASTERINTEGER > FieldOnNodesLong;
typedef std::shared_ptr< FieldOnNodesLong > FieldOnNodesLongPtr;

/** @brief Forward declaration of FieldOnCells */
template < class ValueType >
class FieldOnCells;
typedef FieldOnCells< ASTERDOUBLE > FieldOnCellsReal;
typedef std::shared_ptr< FieldOnCellsReal > FieldOnCellsRealPtr;
typedef FieldOnCells< ASTERINTEGER > FieldOnCellsLong;
typedef std::shared_ptr< FieldOnCellsLong > FieldOnCellsLongPtr;

/**
 * @class XfemModel
 * @brief Datastructure for XfemModel (MODI_MODELE_XFEM)
 */
class XfemModel {
  private:
    /** Fields for XFEM */
    std::map< std::string, DataFieldPtr > _listfields;

    class SubElementTopology {
        const std::string _name;
        const std::string getName() const { return _name; };

      public:
        FieldOnCellsRealPtr pin, pai, pmi;
        FieldOnCellsLongPtr cns, hea, lon;
        SubElementTopology( const std::string );
    };
    class FacetTopology {
        const std::string _name;
        FieldOnCellsLongPtr _heaviside;
        const std::string getName() const { return _name; };

      public:
        FieldOnCellsRealPtr intersection_edge, intersection_pt, intersection_pt2, base;
        FieldOnCellsLongPtr connectivity, length;
        FacetTopology( const std::string );
    };
    class NodalTopology {
        const std::string _name;
        const std::string getName() const { return _name; };

      public:
        FieldOnCellsLongPtr hno, hfa, hse;
        NodalTopology( const std::string );
    };
    SubElementTopology _topose;
    FacetTopology _topofac;
    NodalTopology _topono;
    FieldOnCellsRealPtr _normal_levelset, _tangent_levelset, _local_basis;
    FieldOnCellsLongPtr _nodal_status, _crack_nodes, _crack_conn, _heaviside, _cracked_cells;
    FieldOnNodesLongPtr _xfem_nodes;
    JeveuxVectorLong _contact, _crack_number;
    JeveuxVectorChar8 _crack_names, _pre_cond, _thermic;

  public:
    XfemModel( const std::string modelName );

    DataFieldPtr getField( const std::string fieldType ) const {
        return _listfields.at( fieldType );
    };

    ASTERINTEGER getContact() const;
};

typedef std::shared_ptr< XfemModel > XfemModelPtr;

#endif /* XFEMMODEL_H_ */
