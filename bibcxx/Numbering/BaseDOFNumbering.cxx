/**
 * @file DOFNumbering.cxx
 * @brief Implementation de DOFNumbering
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

#include "astercxx.h"

#include "Numbering/BaseDOFNumbering.h"
#include "Supervis/ResultNaming.h"


BaseDOFNumbering::BaseDOFNumbering( const std::string &type, const JeveuxMemory memType )
    : DataStructure( ResultNaming::getNewResultName(), 14, type, memType ),
      _nameOfSolverDataStructure( JeveuxVectorChar24( getName() + ".NSLV" ) ),
      _globalNumbering( new GlobalEquationNumbering( getName() + ".NUME" ) ),
      _dofDescription( new FieldOnNodesDescription( getName() + ".NUME" ) ),
      _localNumbering( new LocalEquationNumbering( getName() + ".NUML" ) ),
      _model( ModelPtr( nullptr ) ), _listOfLoads( new ListOfLoads() ),
      _smos( new MorseStorage( getName() + ".SMOS" ) ),
      _slcs( new LigneDeCiel( getName() + ".SLCS" ) ),
      _mltf( new MultFrontGarbage( getName() + ".MLTF" ) ), _isEmpty( true ){};

BaseDOFNumbering::BaseDOFNumbering( const std::string name, const std::string &type,
                                              const JeveuxMemory memType )
    : DataStructure( name, 14, type, memType ),
      _nameOfSolverDataStructure( JeveuxVectorChar24( getName() + ".NSLV" ) ),
      _globalNumbering( new GlobalEquationNumbering( getName() + ".NUME" ) ),
      _dofDescription( new FieldOnNodesDescription( getName() + ".NUME" ) ),
      _localNumbering( new LocalEquationNumbering( getName() + ".NUML" ) ),
      _model( ModelPtr( nullptr ) ), _listOfLoads( new ListOfLoads() ),
      _smos( new MorseStorage( getName() + ".SMOS" ) ),
      _slcs( new LigneDeCiel( getName() + ".SLCS" ) ),
      _mltf( new MultFrontGarbage( getName() + ".MLTF" ) ), _isEmpty( true ){};

bool BaseDOFNumbering::computeNumbering() {
    if ( _model ) {
        if ( _model->isEmpty() )
            throw std::runtime_error( "Model is empty" );

        _listOfLoads->build(_model);

        const std::string base( "GG" );
        const std::string null( " " );
        CALLO_NUMERO_WRAP( getName(), base, null, null, _model->getName(),
                           _listOfLoads->getName() );

        const auto FEDescs = _listOfLoads->getFiniteElementDescriptors();
        this->addFiniteElementDescriptors(FEDescs);

    } else if ( _matrix.size() != 0 ) {
        _model = this->getModel();

        ASTERINTEGER nb_matr = _matrix.size();
        JeveuxVectorChar24 jvListOfMatr(ResultNaming::getNewResultName());
        bool retour = jvListOfMatr->allocate( nb_matr );

        int ind = 0;
        for ( const auto &mat : _matrix )
        {
            (*jvListOfMatr)[ind++] = ( boost::apply_visitor( ElementaryMatrixGetName(), mat ) );
            auto FEDescs = boost::apply_visitor( ElementaryMatrixGetFEDescrp(), mat );
            this->addFiniteElementDescriptors(FEDescs);
        }

        CALLO_NUME_DDL_MATR(getName(), jvListOfMatr->getName(), &nb_matr);
    } else
        throw std::runtime_error( "No matrix or model defined" );
    _isEmpty = false;

    return true;
};

std::string BaseDOFNumbering::getPhysicalQuantity() const {
    _globalNumbering->_informations->updateValuePointer();
    JeveuxChar24 physicalQuantity = ( *_globalNumbering->_informations )[1];
    return physicalQuantity.rstrip();
};

VectorLong BaseDOFNumbering::getDirichletBCDOFs( void ) const
{
    JeveuxVectorLong ccid( "&&NUME_CCID" );
    std::string base( "V" );

    // Il faudrait eventuellement rajouter une liste de charge en plus donné par le user
    CALLO_NUMCIMA(_listOfLoads->getName(), getName(), ccid->getName(), base);

    ccid->updateValuePointer();
    return ccid->toVector();
};
