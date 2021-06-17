#ifndef EXTERNALVARIABLEDEFINITION_H_
#define EXTERNALVARIABLEDEFINITION_H_

/**
 * @file ExternalStateVariable.h
 * @brief Fichier entete de la classe ExternalStateVariable
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

#include "DataFields/DataField.h"
#include "Functions/Formula.h"
#include "Functions/Function.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/Mesh.h"
#include "Meshes/MeshEntities.h"
#include "Meshes/ParallelMesh.h"
#include "Meshes/Skeleton.h"
#include "Results/TransientResult.h"
#include "astercxx.h"

class MaterialFieldBuilder;

/**
 * @class EvolutionParameter
 * @brief Class EvolutionParameter to be used when ExternalStateVariable is time dependant
 * @author Nicolas Sellenet
 */
class EvolutionParameter {
  private:
    TransientResultPtr _evol;
    std::string _nomCham;
    std::string _prolGauche;
    std::string _prolDroite;
    FunctionPtr _foncInst;
    FormulaPtr _formInst;

  public:
    EvolutionParameter( const TransientResultPtr &evol )
        : _evol( evol ), _nomCham( "" ), _prolGauche( "EXCLU" ), _prolDroite( "EXCLU" ),
          _foncInst( nullptr ), _formInst( nullptr ){};

    std::string getFieldName() { return _nomCham; };

    std::string getLeftExtension() { return _prolGauche; };

    std::string getRightExtension() { return _prolDroite; };

    TransientResultPtr getTransientResult() { return _evol; };

    FormulaPtr getTimeFormula() { return _formInst; };

    FunctionPtr getTimeFunction() { return _foncInst; };

    void prohibitLeftExtension() { _prolGauche = "EXCLU"; };

    void prohibitRightExtension() { _prolDroite = "EXCLU"; };

    void setConstantLeftExtension() { _prolGauche = "CONSTANT"; };

    void setConstantRightExtension() { _prolDroite = "CONSTANT"; };

    void setFieldName( const std::string &name ) { _nomCham = name; };

    void setLinearLeftExtension() { _prolGauche = "LINEAIRE"; };

    void setLinearRightExtension() { _prolDroite = "LINEAIRE"; };

    void setTimeFunction( const FormulaPtr &func ) {
        _foncInst = nullptr;
        _formInst = func;
    };

    void setTimeFunction( const FunctionPtr &func ) {
        _formInst = nullptr;
        _foncInst = func;
    };
};

typedef boost::shared_ptr< EvolutionParameter > EvolutionParameterPtr;

struct TemperatureExternalStateVariablesTraits {
    constexpr static const char *name = "TEMP";
};

struct GeometryExternalStateVariablesTraits {
    constexpr static const char *name = "GEOM";
};

struct CorrosionExternalStateVariablesTraits {
    constexpr static const char *name = "CORR";
};

struct IrreversibleDeformationExternalStateVariablesTraits {
    constexpr static const char *name = "EPSA";
};

struct ConcreteHydratationExternalStateVariablesTraits {
    constexpr static const char *name = "HYDR";
};

struct IrradiationExternalStateVariablesTraits {
    constexpr static const char *name = "IRRA";
};

struct SteelPhasesExternalStateVariablesTraits {
    constexpr static const char *name = "M_ACIER";
};

struct ZircaloyPhasesExternalStateVariablesTraits {
    constexpr static const char *name = "M_ZIRC";
};

struct Neutral1ExternalStateVariablesTraits {
    constexpr static const char *name = "NEUT1";
};

struct Neutral2ExternalStateVariablesTraits {
    constexpr static const char *name = "NEUT2";
};

struct Neutral3ExternalStateVariablesTraits {
    constexpr static const char *name = "NEUT3";
};

struct ConcreteDryingExternalStateVariablesTraits {
    constexpr static const char *name = "SECH";
};

struct TotalFluidPressureExternalStateVariablesTraits {
    constexpr static const char *name = "PTOT";
};

struct VolumetricDeformationExternalStateVariablesTraits {
    constexpr static const char *name = "DIVU";
};

/**
 * @class BaseExternalStateVariable
 * @brief Input Variable Definition
 * @author Nicolas Sellenet
 */
class BaseExternalStateVariable {
  private:
    BaseMeshPtr _mesh;
    MeshEntityPtr _localization;
    ASTERDOUBLE _refValue;
    bool _refValueSet;
    DataFieldPtr _chamGd;
    EvolutionParameterPtr _evolParam;

  public:
    typedef boost::shared_ptr< BaseExternalStateVariable > BaseExternalStateVariablePtr;

    /**
     * @brief Constructeur
     */
    BaseExternalStateVariable( const BaseMeshPtr &mesh )
        : _mesh( mesh ), _localization( new AllMeshEntities() ), _refValue( 0. ),
          _refValueSet( false ){};

    /**
     * @brief Constructeur
     */
    BaseExternalStateVariable( const BaseMeshPtr &mesh, const std::string &nameOfGroup )
        : _mesh( mesh ), _localization( new GroupOfCells( nameOfGroup ) ), _refValue( 0. ),
          _refValueSet( false ) {
        if ( !_mesh->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );
    };

    /**
     * @brief Destructeur
     */
    ~BaseExternalStateVariable(){};

    /**
     * @brief Function to know if a reference value exists in input variable
     */
    bool hasReferenceValue() const { return _refValueSet; };

    /**
     * @brief Get the field of values of the input variable
     */
    EvolutionParameterPtr getEvolutionParameter() const { return _evolParam; };

    /**
     * @brief Get the field of values of the input variable
     */
    DataFieldPtr getValue() const { return _chamGd; };

    /**
     * @brief Get the reference value of input variable
     */
    ASTERDOUBLE getReferenceValue() const {
        if ( !_refValueSet )
            throw std::runtime_error( "Reference value not set" );
        return _refValue;
    };

    /**
     * @brief Get the name of the variable
     */
    virtual std::string getVariableName() const {
        throw std::runtime_error( "Not allowed" );
        return std::string( "NOTHING" );
    };

    /**
     * @brief Function to set the evolution parameter of the input variable
     */
    void setEvolutionParameter( const EvolutionParameterPtr &evol ) {
        _chamGd = nullptr;
        _evolParam = evol;
    };

    /**
     * @brief Function to set the field of values of the input variable
     */
    void setValue( const DataFieldPtr &field ) {
        _evolParam = nullptr;
        _chamGd = field;
    };

    /**
     * @brief Function to set the reference value of input variable
     */
    void setReferenceValue( const ASTERDOUBLE &value ) {
        _refValue = value;
        _refValueSet = true;
    };
};


/**
 * @class ExternalStateVariable
 * @brief Input Variable Definition
 * @author Nicolas Sellenet
 */
template < typename Parameter >
class ExternalStateVariable : public BaseExternalStateVariable {
  private:
    constexpr static const char *_varcName = Parameter::name;

  public:
    /**
     * @brief Constructeur
     */
    ExternalStateVariable( const BaseMeshPtr &mesh ) : BaseExternalStateVariable( mesh ){};

    /**
     * @brief Constructeur
     */
    ExternalStateVariable( const BaseMeshPtr &mesh, const std::string &nameOfGroup )
        : BaseExternalStateVariable( mesh, nameOfGroup ){};

    /**
     * @brief Destructeur
     */
    ~ExternalStateVariable(){};

    /**
     * @brief Get the name of the variable
     */
    std::string getVariableName() const { return std::string( _varcName ); };
};

typedef ExternalStateVariable< TemperatureExternalStateVariablesTraits >
    TemperatureExternalStateVariable;
typedef ExternalStateVariable< GeometryExternalStateVariablesTraits > GeometryExternalStateVariable;
typedef ExternalStateVariable< CorrosionExternalStateVariablesTraits >
    CorrosionExternalStateVariable;
typedef ExternalStateVariable< IrreversibleDeformationExternalStateVariablesTraits >
    IrreversibleDeformationExternalStateVariable;
typedef ExternalStateVariable< ConcreteHydratationExternalStateVariablesTraits >
    ConcreteHydratationExternalStateVariable;
typedef ExternalStateVariable< IrradiationExternalStateVariablesTraits >
    IrradiationExternalStateVariable;
typedef ExternalStateVariable< SteelPhasesExternalStateVariablesTraits >
    SteelPhasesExternalStateVariable;
typedef ExternalStateVariable< ZircaloyPhasesExternalStateVariablesTraits >
    ZircaloyPhasesExternalStateVariable;
typedef ExternalStateVariable< Neutral1ExternalStateVariablesTraits > Neutral1ExternalStateVariable;
typedef ExternalStateVariable< Neutral2ExternalStateVariablesTraits > Neutral2ExternalStateVariable;
typedef ExternalStateVariable< Neutral3ExternalStateVariablesTraits > Neutral3ExternalStateVariable;
typedef ExternalStateVariable< ConcreteDryingExternalStateVariablesTraits >
    ConcreteDryingExternalStateVariable;
typedef ExternalStateVariable< TotalFluidPressureExternalStateVariablesTraits >
    TotalFluidPressureExternalStateVariable;
typedef ExternalStateVariable< VolumetricDeformationExternalStateVariablesTraits >
    VolumetricDeformationExternalStateVariable;

typedef boost::shared_ptr< BaseExternalStateVariable > BaseExternalStateVariablePtr;
typedef boost::shared_ptr< TemperatureExternalStateVariable > TemperatureExternalStateVariablePtr;
typedef boost::shared_ptr< GeometryExternalStateVariable > GeometryExternalStateVariablePtr;
typedef boost::shared_ptr< CorrosionExternalStateVariable > CorrosionExternalStateVariablePtr;
typedef boost::shared_ptr< IrreversibleDeformationExternalStateVariable >
    IrreversibleDeformationExternalStateVariablePtr;
typedef boost::shared_ptr< ConcreteHydratationExternalStateVariable >
    ConcreteHydratationExternalStateVariablePtr;
typedef boost::shared_ptr< IrradiationExternalStateVariable > IrradiationExternalStateVariablePtr;
typedef boost::shared_ptr< SteelPhasesExternalStateVariable > SteelPhasesExternalStateVariablePtr;
typedef boost::shared_ptr< ZircaloyPhasesExternalStateVariable >
    ZircaloyPhasesExternalStateVariablePtr;
typedef boost::shared_ptr< Neutral1ExternalStateVariable > Neutral1ExternalStateVariablePtr;
typedef boost::shared_ptr< Neutral2ExternalStateVariable > Neutral2ExternalStateVariablePtr;
typedef boost::shared_ptr< Neutral3ExternalStateVariable > Neutral3ExternalStateVariablePtr;

typedef boost::shared_ptr< ConcreteDryingExternalStateVariable >
    ConcreteDryingExternalStateVariablePtr;
typedef boost::shared_ptr< TotalFluidPressureExternalStateVariable >
    TotalFluidPressureExternalStateVariablePtr;
typedef boost::shared_ptr< VolumetricDeformationExternalStateVariable >
    VolumetricDeformationExternalStateVariablePtr;


#endif /* EXTERNALVARIABLEDEFINITION_H_ */
