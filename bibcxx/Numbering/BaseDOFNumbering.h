/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org             */
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

/**
 * @file BaseDOFNumbering.h
 * @brief Fichier entete de la classe BaseDOFNumbering
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

#pragma once

// Forward declaration
template < typename ValueType, PhysicalQuantityEnum PhysicalQuantity >
class ElementaryMatrix;

class DOFNumbering;
class ParallelDOFNumbering;

using ElementaryMatrixDisplacementRealPtr =
    std::shared_ptr< ElementaryMatrix< ASTERDOUBLE, Displacement > >;
using ElementaryMatrixDisplacementComplexPtr =
    std::shared_ptr< ElementaryMatrix< ASTERCOMPLEX, Displacement > >;
using ElementaryMatrixTemperatureRealPtr =
    std::shared_ptr< ElementaryMatrix< ASTERDOUBLE, Temperature > >;
using ElementaryMatrixPressureComplexPtr =
    std::shared_ptr< ElementaryMatrix< ASTERCOMPLEX, Pressure > >;

/**
 * @class BaseDOFNumbering
 * @brief Class definissant un nume_ddl
 * @author Nicolas Sellenet
 */
class BaseDOFNumbering : public DataStructure {
  public:
    typedef std::variant< ElementaryMatrixDisplacementRealPtr,
                          ElementaryMatrixDisplacementComplexPtr,
                          ElementaryMatrixTemperatureRealPtr, ElementaryMatrixPressureComplexPtr >
        MatrElem;

  private:
    class ElementaryMatrixGetModel {
      public:
        template < typename T >
        ModelPtr operator()( const T &operand ) const {
            return operand->getModel();
        };
    };

    class ElementaryMatrixGetName {
      public:
        template < typename T >
        std::string operator()( const T &operand ) const {
            return operand->getName();
        };
    };

    class ElementaryMatrixGetFEDescrp {
      public:
        template < typename T >
        std::vector< FiniteElementDescriptorPtr > operator()( const T &operand ) const {
            return operand->getFiniteElementDescriptors();
        };
    };

    class ElementaryMatrixGetMesh {
      public:
        template < typename T >
        BaseMeshPtr operator()( const T &operand ) const {
            return operand->getMesh();
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
            : _adnt( DOFNumName + ".ADNT" ),
              _glob( DOFNumName + ".GLOB" ),
              _locl( DOFNumName + ".LOCL" ),
              _pnti( DOFNumName + ".PNTI" ),
              _renu( DOFNumName + ".RENU" ),
              _adpi( DOFNumName + ".ADPI" ),
              _adre( DOFNumName + ".ADRE" ),
              _anci( DOFNumName + ".ANCI" ),
              _debf( DOFNumName + ".LFRN" ),
              _deca( DOFNumName + ".DECA" ),
              _defs( DOFNumName + ".DEFS" ),
              _desc( DOFNumName + ".DESC" ),
              _diag( DOFNumName + ".DIAG" ),
              _fils( DOFNumName + ".FILS" ),
              _frer( DOFNumName + ".FRER" ),
              _lgbl( DOFNumName + ".LGBL" ),
              _lgsn( DOFNumName + ".LGSN" ),
              _nbas( DOFNumName + ".NBAS" ),
              _nbli( DOFNumName + ".NBLI" ),
              _nbcl( DOFNumName + ".NCBL" ),
              _nouv( DOFNumName + ".NOUV" ),
              _pare( DOFNumName + ".PARE" ),
              _sequ( DOFNumName + ".SEQU" ),
              _supn( DOFNumName + ".SUPN" ){};
        friend class BaseDOFNumbering;
    };
    typedef std::shared_ptr< MultFrontGarbage > MultFrontGarbagePtr;

  protected:
    class GlobalEquationNumbering : public DataStructure {
      protected:
        /** @brief Objet Jeveux '.NEQU' */
        JeveuxVectorLong _numberOfEquations;
        /** @brief Objet Jeveux '.REFN' */
        JeveuxVectorChar24 _informations;
        /** @brief Objet Jeveux '.DELG' */
        JeveuxVectorLong _lagrangianInformations;

        GlobalEquationNumbering( const std::string &baseName )
            : DataStructure( baseName + ".NUME", 19, "NUME_EQUA" ),
              _numberOfEquations( getName() + ".NEQU" ),
              _informations( getName() + ".REFN" ),
              _lagrangianInformations( getName() + ".DELG" ){};

      public:
        /**
         * @brief Returns a vector of information of the Lagrange multipliers
         */
        const JeveuxVectorLong getLagrangianInformations() const { return _lagrangianInformations; }

        /**
         * @brief Returns a vector of information on the numer of equations
         */
        const JeveuxVectorLong getNumberOfEquations() const { return _numberOfEquations; }

        /**
         * @brief Returns the vector of local to global numbering
         */
        virtual const JeveuxVectorLong getLocalToGlobal() const {
            throw std::runtime_error( "Vector LocalToGlobal doesn't exist in sequential" );
            return JeveuxVectorLong( "RIEN" );
        };

        /**
         * @brief Returns the vector of the rank owning the local dof number
         */
        virtual const JeveuxVectorLong getLocalToRank() const {
            throw std::runtime_error( "Vector LocalToRank doesn't exist in sequential" );
            return JeveuxVectorLong( "RIEN" );
        };

        friend class BaseDOFNumbering;
        friend class DOFNumbering;
        friend class ParallelDOFNumbering;
    };
    typedef std::shared_ptr< GlobalEquationNumbering > GlobalEquationNumberingPtr;

    class LocalEquationNumbering : public DataStructure {
      protected:
        /** @brief Objet Jeveux '.NEQU' */
        JeveuxVectorLong _numberOfEquations;
        /** @brief Objet Jeveux '.DELG' */
        JeveuxVectorLong _lagrangianInformations;
        /** @brief Objet Jeveux '.PRNO' */
        JeveuxCollectionLong _componentsOnNodes;
        /** @brief Objet Jeveux '.NUEQ' */
        JeveuxVectorLong _indexationVector;
        /** @brief Objet Jeveux '.NULG' */
        JeveuxVectorLong _localToGlobal;
        /** @brief Objet Jeveux '.NUGL' */
        JeveuxVectorLong _globalToLocal;
        /** @brief Objet Jeveux '.PDDL' */
        JeveuxVectorLong _localToRank;

        LocalEquationNumbering( const std::string &baseName )
            : DataStructure( baseName + ".NUML", 19, "NUML_EQUA" ),
              _numberOfEquations( getName() + ".NEQU" ),
              _lagrangianInformations( getName() + ".DELG" ),
              _componentsOnNodes( getName() + ".PRNO" ),
              _indexationVector( getName() + ".NUEQ" ),
              _localToGlobal( getName() + ".NULG" ),
              _globalToLocal( getName() + ".NUGL" ),
              _localToRank( getName() + ".PDDL" ){};

      public:
        /**
         * @brief Returns the vector of local to global numbering
         */
        const JeveuxVectorLong getLocalToGlobal() const { return _localToGlobal; }
        /**
         * @brief Returns the vector of global to local numbering
         */
        const JeveuxVectorLong getGlobalToLocal() const { return _globalToLocal; }
        /**
         * @brief Returns the vector of the rank owning the local dof number
         */
        const JeveuxVectorLong getLocalToRank() const { return _localToRank; }

        friend class BaseDOFNumbering;
        friend class DOFNumbering;
        friend class ParallelDOFNumbering;
    };
    typedef std::shared_ptr< LocalEquationNumbering > LocalEquationNumberingPtr;

  private:
    /** @brief Objet Jeveux '.NSLV' */
    JeveuxVectorChar24 _nameOfSolverDataStructure;
    /** @brief Objet prof_chno */
    BaseMeshPtr _mesh;
    /** @brief Objet prof_chno */
    FieldOnNodesDescriptionPtr _dofDescription;
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
     * @param name nom souhaité de la sd (utile pour le BaseDOFNumbering d'une sd_resu)
     */
    BaseDOFNumbering( const std::string name, const std::string &type,
                      const FieldOnNodesDescriptionPtr fdof, const MeshPtr mesh );

    BaseDOFNumbering( const std::string name, const std::string &type );

  public:
    /**
     * @typedef BaseDOFNumberingPtr
     * @brief Pointeur intelligent vers un BaseDOFNumbering
     */
    typedef std::shared_ptr< BaseDOFNumbering > BaseDOFNumberingPtr;

    /**
     * @brief Add a FiniteElementDescriptor to elementary matrix
     * @param FiniteElementDescriptorPtr FiniteElementDescriptor
     */
    bool addFiniteElementDescriptor( const FiniteElementDescriptorPtr &curFED );

    /**
     * @brief Add a FiniteElementDescriptor to elementary matrix
     * @param FiniteElementDescriptorPtr FiniteElementDescriptor
     */
    bool addFiniteElementDescriptors( const std::vector< FiniteElementDescriptorPtr > &curFEDs );

    void setEmpty( const bool &empty ) { _isEmpty = empty; };

    /**
     * @brief Returns the GlobalEquationNumberingPtr
     */
    virtual GlobalEquationNumberingPtr getGlobalNumbering() const = 0;

    /**
     * @brief Build the Numbering of DOFs
     */
    virtual bool computeNumbering( const ModelPtr model, const ListOfLoadsPtr listOfLoads );

    /**
     * @brief Build the Numbering of DOFs
     */
    virtual bool computeNumbering( const std::vector< MatrElem > matrix );

    /**
     * @brief renumbering of DOFs
     */
    virtual bool computeRenumbering( const ModelPtr model, const ListOfLoadsPtr listOfLoads );

    /**
     * @brief Build the Numbering of DOFs
     */
    virtual bool computeNumberingWithLocalMode( const std::string &localMode );

    /**
     * @brief Get Physical Quantity
     */
    virtual std::string getPhysicalQuantity() const = 0;

    /**
     * @brief Are Lagrange Multipliers used for BC or MPC
     */
    virtual bool useLagrangeMultipliers() const { AS_ABORT( "Not allowed" ); };

    /**
     * @brief Are Single Lagrange Multipliers used for BC or MPC
     */
    virtual bool useSingleLagrangeMultipliers() const { AS_ABORT( "Not allowed" ); };

    /**
     * @brief Get The Component Associated To A Given Row
     */
    virtual std::string getComponentAssociatedToRow( const ASTERINTEGER row ) const {
        AS_ABORT( "Not allowed" );
    };

    /**
     * @brief Get The Components Associated To A Given Node
     */
    virtual VectorString getComponentsAssociatedToNode( const ASTERINTEGER node ) const {
        AS_ABORT( "Not allowed" );
    };

    /**
     * @brief Get The Node Id Associated To A Given Row
     */
    virtual ASTERINTEGER getNodeAssociatedToRow( const ASTERINTEGER row ) const {
        AS_ABORT( "Not allowed" );
    };

    /**
     * @brief Return true if a physical dof is Associated To A Given Row
     */
    virtual bool isRowAssociatedToPhysical( const ASTERINTEGER row ) const {
        AS_ABORT( "Not allowed" );
    };

    /**
     * @brief Get The total number of Dofs
     */
    virtual ASTERINTEGER getNumberOfDofs() const { AS_ABORT( "Not allowed" ); };

    /**
     * @brief Get The total number of Dofs
     */
    virtual ASTERINTEGER getNumberOfDofs( const bool local ) const { AS_ABORT( "Not allowed" ); };

    /**
     * @brief get the Row index Associated To the Component of a Node
     */
    virtual ASTERINTEGER getRowAssociatedToNodeComponent( const ASTERINTEGER node,
                                                          const std::string comp ) const {
        AS_ABORT( "Not allowed" );
    }

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    virtual VectorLong getRowsAssociatedToPhysicalDofs() const { AS_ABORT( "Not allowed" ); }

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    virtual VectorLong getRowsAssociatedToLagrangeMultipliers() const { AS_ABORT( "Not allowed" ); }

    /**
     * @brief Get Assigned Components
     */
    virtual VectorString getComponents() const { AS_ABORT( "Not allowed" ); }

    /**
     * @brief Get FieldOnNodesDescription
     */
    FieldOnNodesDescriptionPtr getDescription() const { return _dofDescription; };

    /**
     * @brief Get FieldOnNodesDescription
     */
    void setDescription( const FieldOnNodesDescriptionPtr dofd ) { _dofDescription = dofd; };

    /**
     * @brief Get all FiniteElementDescriptors
     * @return vector of all FiniteElementDescriptors
     */
    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptors() { return _FEDVector; };

    /**
     * @brief Get mesh
     * @return Internal mesh
     */
    BaseMeshPtr getMesh() const;

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
};

/**
 * @typedef BaseDOFNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un BaseDOFNumbering
 * @author Nicolas Sellenet
 */
typedef std::shared_ptr< BaseDOFNumbering > BaseDOFNumberingPtr;
