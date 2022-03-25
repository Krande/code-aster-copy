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

    // I use them to debug easily mutiple creation
    // I don't use map directly to avoid to keep in memory unnecessary objects
    static std::set< std::string > _setProfChno;
    static std::set< std::string > _setLigrel;

    /**
     * @brief Add a existing FiniteElementDescriptor in FieldBuilder
     */
    FiniteElementDescriptorPtr newFiniteElementDescriptor( const std::string &name,
                                                           const BaseMeshPtr mesh ) {
#ifdef ASTER_DEBUG_CXX
        if ( _setLigrel.count( trim( name ) ) > 0 ) {
            raiseAsterError( "LIGREL already exists: " + name );
        }
#endif

        auto curDesc = std::make_shared< FiniteElementDescriptor >( name, mesh );

        addFiniteElementDescriptor( curDesc );

        return curDesc;
    };

    /**
     * @brief Add a existing FieldOnNodesDescription in FieldBuilder
     */
    FieldOnNodesDescriptionPtr newFieldOnNodesDescription( const std::string &name ) {
#ifdef ASTER_DEBUG_CXX
        if ( _setProfChno.count( trim( name ) ) > 0 ) {
            raiseAsterError( "PROF_CHNO already exists: " + name );
        }
#endif

        auto curDesc = std::make_shared< FieldOnNodesDescription >( name );
        addFieldOnNodesDescription( curDesc );

        return curDesc;
    };

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
        _setProfChno.insert( trim( fond->getName() ) );
    };

    /**
     * @brief Add a existing FiniteElementDescriptor in FieldBuilder
     */
    void addFiniteElementDescriptor( const FiniteElementDescriptorPtr &fed ) {
        AS_ASSERT( fed );

        _mapLigrel[trim( fed->getName() )] = fed;
        _setLigrel.insert( trim( fed->getName() ) );
    };

    /**
     * @brief Build a FieldOnCells with a FiniteElementDescriptor
     */
    template < typename ValueType >
    std::shared_ptr< FieldOnCells< ValueType > > buildFieldOnCells( const std::string &name,
                                                                    const BaseMeshPtr mesh ) {
        std::shared_ptr< FieldOnCells< ValueType > > field =
            std::make_shared< FieldOnCells< ValueType > >( name );
        field->updateValuePointers();

        const std::string ligrel = trim( ( *( *field )._reference )[0].toString() );

        if ( !ligrel.empty() ) {
            auto curIter = _mapLigrel.find( ligrel );
            FiniteElementDescriptorPtr curDesc;
            if ( curIter != _mapLigrel.end() ) {
                curDesc = curIter->second;
            } else {
                curDesc = newFiniteElementDescriptor( ligrel, mesh );
            }

            field->setDescription( curDesc );
        }
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
        std::shared_ptr< FieldOnNodes< ValueType > > field =
            std::make_shared< FieldOnNodes< ValueType > >( name );
        field->updateValuePointers();

        const std::string profchno = trim( ( *( *field )._reference )[1].toString() );
        if ( !profchno.empty() ) {

            auto curIter = _mapProfChno.find( profchno );
            FieldOnNodesDescriptionPtr curDesc;
            if ( curIter != _mapProfChno.end() )
                curDesc = curIter->second;
            else {
                curDesc = newFieldOnNodesDescription( profchno );
            }
            field->setDescription( curDesc );
        }

        return field;
    };

    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptors() const {
        std::vector< FiniteElementDescriptorPtr > ret;

        for ( auto &[name, fed] : _mapLigrel )
            ret.push_back( fed );

        return ret;
    };

    std::vector< FieldOnNodesDescriptionPtr > getFieldOnNodesDescriptions() const {
        std::vector< FieldOnNodesDescriptionPtr > ret;

        for ( auto &[name, fnd] : _mapProfChno )
            ret.push_back( fnd );

        return ret;
    };
};

#endif /* FIELDBUILDER_H_ */
