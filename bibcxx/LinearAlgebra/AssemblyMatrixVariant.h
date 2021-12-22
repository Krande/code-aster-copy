#ifndef ASSEMBLYMATRIXVARIANT_H_
#define ASSEMBLYMATRIXVARIANT_H_

/**
 * @file AssemblyMatrixVariant.h
 * @brief Fichier entete de la classe AssemblyMatrix
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

#include "LinearAlgebra/AssemblyMatrix.h"

class AssemblyMatrixVariant {
  public:
    typedef boost::variant< AssemblyMatrixDisplacementRealPtr, AssemblyMatrixDisplacementComplexPtr,
                            AssemblyMatrixTemperatureRealPtr, AssemblyMatrixPressureRealPtr >
        MatrixVariant;

  private:
    class GetName : public boost::static_visitor< std::string > {
      public:
        template < typename T > std::string operator()( const T &operand ) const {
            return operand->getName();
        };
    };

    class GetFacto : public boost::static_visitor< bool > {
      public:
        template < typename T > bool operator()( const T &operand ) const {
            return operand->isFactorized();
        };
    };

    class IsFacto : public boost::static_visitor< void > {
      public:
        bool _facto;
        IsFacto( const bool &facto ) : _facto( facto ){};

        template < typename T > void operator()( const T &operand ) const {
            operand->isFactorized( _facto );
        };
    };

    class DelFacto : public boost::static_visitor< bool > {
      public:
        template < typename T > bool operator()( const T &operand ) const {
            return operand->deleteFactorizedMatrix();
        };
    };

    class SetSolverName : public boost::static_visitor< void > {
      public:
        std::string _name;
        SetSolverName( const std::string &name ) : _name( name ){};

        template < typename T > void operator()( const T &operand ) const {
            operand->setSolverName( _name );
        };
    };

    /* a matrix */
    MatrixVariant _matrix;
    /* has a matrix */
    bool _hasMatrix;

  public:
    typedef boost::shared_ptr< AssemblyMatrixVariant > AssemblyMatrixVariantPtr;

    AssemblyMatrixVariant() : _hasMatrix( false ){};

    bool hasMatrix() const { return _hasMatrix; };

    template < typename T >
    bool holds_alternative() const { // return std::holds_alternative< T >( _matrix ); when c++17
        return boost::get< T >( &_matrix ) != nullptr;
    };

    void operator=( const AssemblyMatrixDisplacementRealPtr &mat ) {
        _matrix = mat;
        _hasMatrix = true;
    }
    void operator=( const AssemblyMatrixDisplacementComplexPtr &mat ) {
        _matrix = mat;
        _hasMatrix = true;
    }
    void operator=( const AssemblyMatrixTemperatureRealPtr &mat ) {
        _matrix = mat;
        _hasMatrix = true;
    }
    void operator=( const AssemblyMatrixPressureRealPtr &mat ) {
        _matrix = mat;
        _hasMatrix = true;
    }

    void setMatrix( const AssemblyMatrixDisplacementRealPtr &mat ) {
        _matrix = mat;
        _hasMatrix = true;
    }

    void setMatrix( const AssemblyMatrixDisplacementComplexPtr &mat ) {
        _matrix = mat;
        _hasMatrix = true;
    }

    void setMatrix( const AssemblyMatrixTemperatureRealPtr &mat ) {
        _matrix = mat;
        _hasMatrix = true;
    }

    void setMatrix( const AssemblyMatrixPressureRealPtr &mat ) {
        _matrix = mat;
        _hasMatrix = true;
    }

    std::string getName() const { return boost::apply_visitor( GetName(), _matrix ); };

    bool isFactorized() const { return boost::apply_visitor( GetFacto(), _matrix ); };

    void isFactorized( const bool &facto ) { boost::apply_visitor( IsFacto( facto ), _matrix ); };

    bool deleteFactorizedMatrix() { return boost::apply_visitor( DelFacto(), _matrix ); };

    void setSolverName( const std::string &name ) {
        boost::apply_visitor( SetSolverName( name ), _matrix );
    }
};

typedef boost::shared_ptr< AssemblyMatrixVariant > AssemblyMatrixVariantPtr;

#endif /* ASSEMBLYMATRIXVARAINT_H_ */
