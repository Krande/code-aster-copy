/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#ifndef BASEDOFNUMBERING_H_
#define BASEDOFNUMBERING_H_

/**
 * @file BaseDOFNumbering.h
 * @brief Fichier entete de la classe BaseDOFNumbering
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
 *   This file is part of code_aster.
 *
 *   code_aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   code_aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 */

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"
#include "boost/variant.hpp"

#include "DataStructures/DataStructure.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "LinearAlgebra/MatrixStorage.h"
#include "Loads/DirichletBC.h"
#include "Loads/ListOfLoads.h"
#include "Loads/MechanicalLoad.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Numbering/FieldOnNodesDescription.h"

/**
 * @class BaseDOFNumbering
 * @brief Class definissant un nume_ddl
 *        Cette classe est volontairement succinte car on n'en connait pas encore l'usage
 * @author Nicolas Sellenet
 */
class BaseDOFNumbering : public DataStructure {
  private:
    typedef boost::variant< ElementaryMatrixDisplacementRealPtr,
                            ElementaryMatrixDisplacementComplexPtr,
                            ElementaryMatrixTemperatureRealPtr, ElementaryMatrixPressureComplexPtr >
        MatrElem;

    class ElementaryMatrixGetModel : public boost::static_visitor< ModelPtr > {
      public:
        template < typename T > ModelPtr operator()( const T &operand ) const {
            return operand->getModel();
        };
    };

    class ElementaryMatrixGetName : public boost::static_visitor< std::string > {
      public:
        template < typename T > std::string operator()( const T &operand ) const {
            return operand->getName();
        };
    };

    class ElementaryMatrixGetFEDescrp : public boost::static_visitor<
                std::vector< FiniteElementDescriptorPtr > > {
      public:
        template < typename T > std::vector< FiniteElementDescriptorPtr >
            operator()( const T &operand ) const {
            return operand->getFiniteElementDescriptors();
        };
    };

  private:
    class MultFrontGarbage {
        /** @brief Objet Jeveux '.ADNT' */
        JeveuxVectorShort _adnt;
        /** @brief Objet Jeveux '.GLOB' */
        JeveuxVectorShort _glob;
        /** @brief Objet Jeveux '.LOCL' */
        JeveuxVectorShort _locl;
        /** @brief Objet Jeveux '.PNTI' */
        JeveuxVectorShort _pnti;
        /** @brief Objet Jeveux '.RENU' */
        JeveuxVectorChar8 _renu;
        /** @brief Objet Jeveux '.ADPI' */
        JeveuxVectorLong _adpi;
        /** @brief Objet Jeveux '.ADRE' */
        JeveuxVectorLong _adre;
        /** @brief Objet Jeveux '.ANCI' */
        JeveuxVectorLong _anci;
        /** @brief Objet Jeveux '.LFRN' */
        JeveuxVectorLong _debf;
        /** @brief Objet Jeveux '.DECA' */
        JeveuxVectorLong _deca;
        /** @brief Objet Jeveux '.DEFS' */
        JeveuxVectorLong _defs;
        /** @brief Objet Jeveux '.DESC' */
        JeveuxVectorLong _desc;
        /** @brief Objet Jeveux '.DIAG' */
        JeveuxVectorLong _diag;
        /** @brief Objet Jeveux '.FILS' */
        JeveuxVectorLong _fils;
        /** @brief Objet Jeveux '.FRER' */
        JeveuxVectorLong _frer;
        /** @brief Objet Jeveux '.LGBL' */
        JeveuxVectorLong _lgbl;
        /** @brief Objet Jeveux '.LGSN' */
        JeveuxVectorLong _lgsn;
        /** @brief Objet Jeveux '.NBAS' */
        JeveuxVectorLong _nbas;
        /** @brief Objet Jeveux '.NBLI' */
        JeveuxVectorLong _nbli;
        /** @brief Objet Jeveux '.NCBL' */
        JeveuxVectorLong _nbcl;
        /** @brief Objet Jeveux '.NOUV' */
        JeveuxVectorLong _nouv;
        /** @brief Objet Jeveux '.PARE' */
        JeveuxVectorLong _pare;
        /** @brief Objet Jeveux '.SEQU' */
        JeveuxVectorLong _sequ;
        /** @brief Objet Jeveux '.SUPN' */
        JeveuxVectorLong _supn;

        MultFrontGarbage( const std::string &DOFNumName )
            : _adnt( DOFNumName + ".ADNT" ), _glob( DOFNumName + ".GLOB" ),
              _locl( DOFNumName + ".LOCL" ), _pnti( DOFNumName + ".PNTI" ),
              _renu( DOFNumName + ".RENU" ), _adpi( DOFNumName + ".ADPI" ),
              _adre( DOFNumName + ".ADRE" ), _anci( DOFNumName + ".ANCI" ),
              _debf( DOFNumName + ".LFRN" ), _deca( DOFNumName + ".DECA" ),
              _defs( DOFNumName + ".DEFS" ), _desc( DOFNumName + ".DESC" ),
              _diag( DOFNumName + ".DIAG" ), _fils( DOFNumName + ".FILS" ),
              _frer( DOFNumName + ".FRER" ), _lgbl( DOFNumName + ".LGBL" ),
              _lgsn( DOFNumName + ".LGSN" ), _nbas( DOFNumName + ".NBAS" ),
              _nbli( DOFNumName + ".NBLI" ), _nbcl( DOFNumName + ".NCBL" ),
              _nouv( DOFNumName + ".NOUV" ), _pare( DOFNumName + ".PARE" ),
              _sequ( DOFNumName + ".SEQU" ), _supn( DOFNumName + ".SUPN" ){};
        friend class BaseDOFNumbering;
    };
    typedef boost::shared_ptr< MultFrontGarbage > MultFrontGarbagePtr;

    class GlobalEquationNumbering {
        /** @brief Objet Jeveux '.NEQU' */
        JeveuxVectorLong _numberOfEquations;
        /** @brief Objet Jeveux '.REFN' */
        JeveuxVectorChar24 _informations;
        /** @brief Objet Jeveux '.DELG' */
        JeveuxVectorLong _lagrangianInformations;

        GlobalEquationNumbering( const std::string &DOFNumName )
            : _numberOfEquations( DOFNumName + ".NEQU" ), _informations( DOFNumName + ".REFN" ),
              _lagrangianInformations( DOFNumName + ".DELG" ){};

        public:
            /**
             * @brief Returns a vector of information of the Lagrange multipliers
             */
            const JeveuxVectorLong getLagrangianInformations() const {
                return _lagrangianInformations;
            }
            /**
             * @brief Returns a vector of information on the numer of equations
             */
            const JeveuxVectorLong getNumberOfEquations() const {
                return _numberOfEquations;
            }

        friend class BaseDOFNumbering;
    };
    typedef boost::shared_ptr< GlobalEquationNumbering > GlobalEquationNumberingPtr;

    class LocalEquationNumbering {
        /** @brief Objet Jeveux '.NEQU' */
        JeveuxVectorLong _numberOfEquations;
        /** @brief Objet Jeveux '.DELG' */
        JeveuxVectorLong _lagrangianInformations;
        /** @brief Objet Jeveux '.PRNO' */
        JeveuxCollectionLong _componentsOnNodes;
        /** @brief Objet Jeveux '.NUEQ' */
        JeveuxVectorLong _indexationVector;
        /** @brief Objet Jeveux '.NULG' */
        JeveuxVectorLong _globalToLocal;
        /** @brief Objet Jeveux '.NUGL' */
        JeveuxVectorLong _LocalToGlobal;

        LocalEquationNumbering( const std::string &DOFNumName )
            : _numberOfEquations( DOFNumName + ".NEQU" ),
              _lagrangianInformations( DOFNumName + ".DELG" ),
              _componentsOnNodes( DOFNumName + ".PRNO" ), _indexationVector( DOFNumName + ".NUEQ" ),
              _globalToLocal( DOFNumName + ".NULG" ), _LocalToGlobal( DOFNumName + ".NUGL" ){};

        public:
            /**
             * @brief Returns the vector of local to global numbering
             */
            const JeveuxVectorLong getLocalToGlobal() const {
                return _LocalToGlobal;
            }
            /**
             * @brief Returns the vector of global to local numbering
             */
            const JeveuxVectorLong getGlobalToLocal() const {
                return _globalToLocal;
            }
        friend class BaseDOFNumbering;
    };
    typedef boost::shared_ptr< LocalEquationNumbering > LocalEquationNumberingPtr;

    // !!! Classe succinte car on ne sait pas comment elle sera utiliser !!!
    /** @brief Objet Jeveux '.NSLV' */
    JeveuxVectorChar24 _nameOfSolverDataStructure;
    /** @brief Objet '.NUME' */
    GlobalEquationNumberingPtr _globalNumbering;
    /** @brief Objet prof_chno */
    FieldOnNodesDescriptionPtr _dofDescription;
    /** @brief Objet '.NUML' */
    LocalEquationNumberingPtr _localNumbering;
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Matrices elementaires */
    std::vector< MatrElem > _matrix;
    /** @brief Chargements */
    ListOfLoadsPtr _listOfLoads;
    /** @brief Objet Jeveux '.SMOS' */
    MorseStoragePtr _smos;
    /** @brief Objet Jeveux '.SLCS' */
    LigneDeCielPtr _slcs;
    /** @brief Objet Jeveux '.MLTF' */
    MultFrontGarbagePtr _mltf;
    /** @brief Booleen permettant de preciser sur la sd est vide */
    bool _isEmpty;

    /** @brief Vectors of FiniteElementDescriptor */
    std::vector< FiniteElementDescriptorPtr > _FEDVector;
    std::set< std::string > _FEDNames;

  protected:
    /**
     * @brief Constructeur
     */
    BaseDOFNumbering( const std::string &type, const JeveuxMemory memType = Permanent );

    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le BaseDOFNumbering d'une sd_resu)
     */
    BaseDOFNumbering( const std::string name, const std::string &type,
                           const JeveuxMemory memType = Permanent );

  public:
    /**
     * @typedef BaseDOFNumberingPtr
     * @brief Pointeur intelligent vers un BaseDOFNumbering
     */
    typedef boost::shared_ptr< BaseDOFNumbering > BaseDOFNumberingPtr;

    /**
     * @brief Add a FiniteElementDescriptor to elementary matrix
     * @param FiniteElementDescriptorPtr FiniteElementDescriptor
     */
    bool addFiniteElementDescriptor( const FiniteElementDescriptorPtr &curFED ) {
        const auto name = trim( curFED->getName() );
        if ( _FEDNames.find( name ) == _FEDNames.end() ) {
            _FEDVector.push_back( curFED );
            _FEDNames.insert( name );
            return true;
        }
        return false;
    };

    /**
     * @brief Add a FiniteElementDescriptor to elementary matrix
     * @param FiniteElementDescriptorPtr FiniteElementDescriptor
     */
    bool addFiniteElementDescriptors( const std::vector<FiniteElementDescriptorPtr> &curFEDs ) {
        for( auto& curFED : curFEDs)
        {
            const bool ret = this->addFiniteElementDescriptor( curFED );
            if( !ret )
                return false;
        }

        return true;
    };

    /**
     * @brief Function d'ajout d'un chargement
     * @param Args... Liste d'arguments template
     */
    template < typename... Args > void addLoad( const Args &... a ) {
        _listOfLoads->addLoad( a... );
    };

    /**
     * @brief Returns the GlobalEquationNumberingPtr
     */
    const GlobalEquationNumberingPtr getGlobalNumbering() const {
        return _globalNumbering;
    }

    /**
     * @brief Returns the LocalEquationNumberingPtr
     */
    const LocalEquationNumberingPtr  getLocalNumbering() const {
        return _localNumbering;
    }

    /**
     * @brief Build the Numbering of DOFs
     */
    bool computeNumbering();

    /**
     * @brief Get Physical Quantity
     */
    std::string getPhysicalQuantity() const;

    /**
     * @brief Are Lagrange Multipliers used for BC or MPC
     */
    virtual bool useLagrangeMultipliers() const { throw std::runtime_error( "Not allowed" ); };

    /**
     * @brief Are Single Lagrange Multipliers used for BC or MPC
     */
    virtual bool useSingleLagrangeMultipliers() const {
        throw std::runtime_error( "Not allowed" ); };

    /**
     * @brief Get The Component Associated To A Given Row
     */
    virtual std::string getComponentAssociatedToRow(const ASTERINTEGER row) const {
        throw std::runtime_error( "Not allowed" ); };

    /**
     * @brief Get The Components Associated To A Given Node
     */
    virtual VectorString getComponentsAssociatedToNode(const ASTERINTEGER node) const {
        throw std::runtime_error( "Not allowed" ); };

    /**
     * @brief Get The Node Id Associated To A Given Row
     */
    virtual ASTERINTEGER getNodeAssociatedToRow(const ASTERINTEGER row) const {
        throw std::runtime_error( "Not allowed" ); };

    /**
     * @brief Get The total number of Dofs
     */
    virtual ASTERINTEGER getNumberOfDofs() const { throw std::runtime_error( "Not allowed" ); };

    /**
     * @brief get the Row index Associated To the Component of a Node
     */
    virtual ASTERINTEGER getRowAssociatedToNodeComponent(const ASTERINTEGER node,
                                                         const std::string comp) const {
        throw std::runtime_error( "Not allowed" );
    }

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    virtual VectorLong getRowsAssociatedToPhysicalDofs() const {
        throw std::runtime_error( "Not allowed" );
    }

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    virtual VectorLong getRowsAssociatedToLagrangeMultipliers() const {
        throw std::runtime_error( "Not allowed" );
    }

    /**
     * @brief Get Assigned Components
     */
    virtual VectorString getComponents() const {
        throw std::runtime_error( "Not allowed" );
    }

    /**
     * @brief Get FieldOnNodesDescription
     */
    FieldOnNodesDescriptionPtr getDescription() const { return _dofDescription; };

    /**
     * @brief Get all FiniteElementDescriptors
     * @return vector of all FiniteElementDescriptors
     */
    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptors() { return _FEDVector; };

    /**
     * @brief Get model
     */
    ModelPtr getModel() const {
        if ( _model != nullptr )
            return _model;
        else {
            if ( _matrix.size() != 0 )
                return boost::apply_visitor( ElementaryMatrixGetModel(), _matrix[0] );
        }
        return ModelPtr( nullptr );
    };

    /**
     * @brief Get mesh
     * @return Internal mesh
     */
    BaseMeshPtr getMesh() const {
        const auto model = this->getModel();
        if ( model != nullptr ) {
            return model->getMesh();
        }
        return nullptr;
    };

    /**
     * @brief Methode permettant de savoir si la numerotation est vide
     * @return true si la numerotation est vide
     */
    bool isEmpty() { return _isEmpty; };

    /**
     * @brief Methode permettant de savoir si l'objet est parallel
     * @return false
     */
    virtual bool isParallel() { return false; };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    virtual void setElementaryMatrix( const ElementaryMatrixDisplacementRealPtr &currentMatrix )

    {
        if ( _model )
            throw std::runtime_error(
                "It is not allowed to defined Model and ElementaryMatrix together" );
        _matrix.push_back( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    virtual void setElementaryMatrix( const ElementaryMatrixDisplacementComplexPtr &currentMatrix )

    {
        if ( _model )
            throw std::runtime_error(
                "It is not allowed to defined Model and ElementaryMatrix together" );
        _matrix.push_back( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    virtual void setElementaryMatrix( const ElementaryMatrixTemperatureRealPtr &currentMatrix )

    {
        if ( _model )
            throw std::runtime_error(
                "It is not allowed to defined Model and ElementaryMatrix together" );
        _matrix.push_back( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir les matrices elementaires
     * @param currentMatrix objet ElementaryMatrix
     */
    virtual void setElementaryMatrix( const ElementaryMatrixPressureComplexPtr &currentMatrix )

    {
        if ( _model )
            throw std::runtime_error(
                "It is not allowed to defined Model and ElementaryMatrix together" );
        _matrix.push_back( currentMatrix );
    };

    /**
     * @brief Methode permettant de definir la liste de charge
     * @param currentList Liste charge
     */
    void setListOfLoads( const ListOfLoadsPtr &currentList ) { _listOfLoads = currentList; };

    ListOfLoadsPtr getListOfLoads( void ) const { return _listOfLoads; };

    VectorLong getDirichletBCDOFs (void ) const;

    /**
     * @brief Methode permettant de definir le modele
     * @param currentModel Modele de la numerotation
     */
    virtual void setModel( const ModelPtr &currentModel ) {
        if ( _matrix.size() != 0 )
            throw std::runtime_error(
                "It is not allowed to defined Model and ElementaryMatrix together" );
        _model = currentModel;
        this->addFiniteElementDescriptor(_model->getFiniteElementDescriptor());
    };

    bool hasDirichletBC() const
    {
        return _listOfLoads->hasDirichletBC();
    }
};


/**
 * @typedef BaseDOFNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un BaseDOFNumbering
 * @author Nicolas Sellenet
 */
typedef boost::shared_ptr< BaseDOFNumbering > BaseDOFNumberingPtr;

#endif /* BASEDOFNUMBERING_H_ */
