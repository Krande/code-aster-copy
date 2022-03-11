#ifndef FIELDBUILDER_H_
#define FIELDBUILDER_H_

/**
 * @file FieldBuilder.h
 * @brief Header of class FieldBuilder
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

#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Numbering/DOFNumbering.h"

/**
 * @class FieldBuilder
 * @brief This class builds FieldOnNodes and FieldOnCells with respect of
 *        FieldOnNodesDescription and FiniteElementDescriptor
 * @author Nicolas Sellenet
 */
class FieldBuilder {
  private:
    std::map< std::string, FieldOnNodesDescriptionPtr > _mapProfChno;
    std::map< std::string, FiniteElementDescriptorPtr > _mapLigrel;

  public:
    /**
     * @brief Constructor
     */
    FieldBuilder(){};

    /**
     * @brief Add a existing FieldOnNodesDescription in FieldBuilder
     */
    void addFieldOnNodesDescription( const FieldOnNodesDescriptionPtr &fond ) {
        AS_ASSERT( fond );

        _mapProfChno[trim( fond->getName() )] = fond;
    };

    /**
     * @brief Add a existing FiniteElementDescriptor in FieldBuilder
     */
    void addFiniteElementDescriptor( const FiniteElementDescriptorPtr &fed ) {
        AS_ASSERT( fed );

        _mapLigrel[trim( fed->getName() )] = fed;
    };

    /**
     * @brief Build a FieldOnCells with a FiniteElementDescriptor
     */
    template < typename ValueType >
    std::shared_ptr< FieldOnCells< ValueType > > buildFieldOnCells( const std::string &name,
                                                                      const BaseMeshPtr mesh ) {
        typedef FiniteElementDescriptor FEDDesc;
        typedef FiniteElementDescriptorPtr FEDDescP;

        std::shared_ptr< FieldOnCells< ValueType > > field =
            std::make_shared< FieldOnCells< ValueType > >( name );
        field->updateValuePointers();

        const std::string ligrel = trim( ( *( *field )._reference )[0].toString() );

        auto curIter = _mapLigrel.find( ligrel );
        FEDDescP curDesc;
        if ( curIter != _mapLigrel.end() ) {
            curDesc = curIter->second;
        } else {
            curDesc = std::make_shared< FEDDesc >( ligrel, mesh );
            _mapLigrel[ligrel] = curDesc;
        }
        field->setDescription( curDesc );

        return field;
    };

    /**
     * @brief Build a ConstantFieldOnCells with a FiniteElementDescriptor
     */
    template < typename ValueType >
    std::shared_ptr< ConstantFieldOnCells< ValueType > >
    buildConstantFieldOnCells( const std::string &name, const BaseMeshPtr mesh ) {

        std::shared_ptr< ConstantFieldOnCells< ValueType > > field =
            std::make_shared< ConstantFieldOnCells< ValueType > >( name, mesh );
        field->updateValuePointers();

        return field;
    };

    /**
     * @brief Build a FieldOnNodes with a FieldOnNodesDescription
     */
    template < typename ValueType >
    std::shared_ptr< FieldOnNodes< ValueType > > buildFieldOnNodes( std::string name ) {
        typedef FieldOnNodesDescription FONDesc;
        typedef FieldOnNodesDescriptionPtr FONDescP;

        std::shared_ptr< FieldOnNodes< ValueType > > field =
            std::make_shared< FieldOnNodes< ValueType > >( name );
        field->updateValuePointers();

        const std::string profchno = trim( ( *( *field )._reference )[1].toString() );

        auto curIter = _mapProfChno.find( profchno );
        FONDescP curDesc;
        if ( curIter != _mapProfChno.end() )
            curDesc = curIter->second;
        else {
            curDesc = std::make_shared< FONDesc >( profchno );
            _mapProfChno[profchno] = curDesc;
        }
        field->setDescription( curDesc );

        return field;
    };
};

#endif /* FIELDBUILDER_H_ */
