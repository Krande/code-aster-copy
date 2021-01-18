#ifndef KINEMATICSLOAD_H_
#define KINEMATICSLOAD_H_

/**
 * @file DirichletBC.h
 * @brief Fichier entete de la classe DirichletBC
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
#include <list>
#include <stdexcept>
#include <string>

#include "astercxx.h"
#include "Functions/Function.h"
#include "Loads/UnitaryLoad.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"

/**
 * @class DirichletBCClass
 * @brief Classe definissant une charge cinematique (issue d'AFFE_CHAR_CINE)
 * @author Nicolas Sellenet
 */
class DirichletBCClass : public DataStructure {
  protected:
    /** @typedef Pointeur intelligent sur un VirtualMeshEntity */
    typedef boost::shared_ptr< VirtualMeshEntity > MeshEntityPtr;

    /** @typedef std::list de RealLoadDisplacement */
    typedef std::list< RealLoadDisplacement > ListRealDisp;
    /** @typedef Iterateur sur ListRealDisp */
    typedef ListRealDisp::iterator ListRealDispIter;

    /** @typedef std::list de RealLoadTemperature */
    typedef std::list< RealLoadTemperature > ListRealTemp;
    typedef std::list< FunctionLoadTemperature > ListFunctionTemp;
    /** @typedef Iterateur sur ListRealTemp */
    typedef ListRealTemp::iterator ListRealTempIter;

    /** @brief Modele */
    ModelPtr _model;
    /** @brief Listes des valeurs imposees DEPL_R et TEMP_R */
    ListRealDisp _listOfRealImposedDisplacement;
    ListRealTemp _listOfRealImposedTemperature;
    ListFunctionTemp _listOfFunctionImposedTemperature;
    JeveuxVectorLong _intParam;
    JeveuxVectorChar8 _charParam;
    JeveuxVectorReal _doubleParam;
    /** @brief La SD est-elle vide ? */
    bool _isEmpty;

    DirichletBCClass( void ) = delete;

    /**
     * @brief Constructeur
     */
    DirichletBCClass( const std::string &type, const ModelPtr& model );

    /**
     * @brief Constructeur
     */
    DirichletBCClass( const std::string &name, const std::string &type, const ModelPtr& model );

  public:
    /**
     * @typedef DirichletBCPtr
     * @brief Pointeur intelligent vers un DirichletBC
     */
    typedef boost::shared_ptr< DirichletBCClass > DirichletBCPtr;

    /**
     * @brief Construction de la charge (appel a OP0101)
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool build() ;

    ModelPtr getModel() const
    {
        return _model;
    }

  private:
    /**
     * @brief Definition du modele
     * @param currentModel objet Model sur lequel la charge reposera
     */
    bool setModel( const ModelPtr &currentModel ) {
        if ( currentModel->isEmpty() )
            throw std::runtime_error( "Model is empty" );
        _model = currentModel;
        return true;
    };
};

/**
 * @class MechanicalDirichletBCClass
 * @brief Classe definissant une charge cinematique (issue d'AFFE_CHAR_CINE)
 * @author Nicolas Sellenet
 */
class MechanicalDirichletBCClass : public DirichletBCClass {
  public:
    /**
     * @brief Constructeur
     */
    MechanicalDirichletBCClass( void ) = delete;


    MechanicalDirichletBCClass( const ModelPtr& model)
        : DirichletBCClass( "_MECA", model ){};

    /**
     * @brief Constructeur
     */
    MechanicalDirichletBCClass( const std::string name, const ModelPtr& model )
        : DirichletBCClass( name, "_MECA", model ){};

    /**
     * @typedef MechanicalDirichletBCPtr
     * @brief Pointeur intelligent vers un MechanicalDirichletBC
     */
    typedef boost::shared_ptr< MechanicalDirichletBCClass > MechanicalDirichletBCPtr;

    /**
     * @brief Ajout d'une valeur mecanique imposee sur un groupe de mailles
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addImposedMechanicalDOFOnCells(
        const PhysicalQuantityComponent &coordinate, const double &value,
        const std::string &nameOfGroup ) {
        if ( !_model->getMesh()->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

        MeshEntityPtr meshEnt( new GroupOfCells( nameOfGroup ) );
        RealLoadDisplacement resu( meshEnt, coordinate, value );
        _listOfRealImposedDisplacement.push_back( resu );
        return true;
    };

    /**
     * @brief Ajout d'une valeur mecanique imposee sur un groupe de mailles
     * @param namesOfGroup Noms des groupes sur lequels imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addImposedMechanicalDOFOnCells(
        const PhysicalQuantityComponent &coordinate, const double &value,
        const VectorString &namesOfGroup ) {
        for ( const auto &nameOfGroup : namesOfGroup )
            addImposedMechanicalDOFOnCells( coordinate, value, nameOfGroup );
        return true;
    };

    /**
     * @brief Ajout d'une valeur mecanique imposee sur un groupe de noeuds
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool
    addImposedMechanicalDOFOnNodes( const PhysicalQuantityComponent &coordinate,
                                    const double &value,
                                    const std::string &nameOfGroup ) {
        if ( !_model->getMesh()->hasGroupOfNodes( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

        MeshEntityPtr meshEnt( new GroupOfNodes( nameOfGroup ) );
        RealLoadDisplacement resu( meshEnt, coordinate, value );
        _listOfRealImposedDisplacement.push_back( resu );
        return true;
    };

    /**
     * @brief Ajout d'une valeur mecanique imposee sur un groupe de noeuds
     * @param namesOfGroup Noms des groupe sur lequels imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addImposedMechanicalDOFOnNodes(
        const PhysicalQuantityComponent &coordinate, const double &value,
        const VectorString &namesOfGroup ) {
        for ( const auto &nameOfGroup : namesOfGroup )
            addImposedMechanicalDOFOnNodes( coordinate, value, nameOfGroup );
        return true;
    };
};

/**
 * @class ThermalDirichletBCClass
 * @brief Classe definissant une charge cinematique (issue d'AFFE_CHAR_CINE)
 * @author Nicolas Sellenet
 */
class ThermalDirichletBCClass : public DirichletBCClass {
  public:
    /**
     * @brief Constructeur
     */
    ThermalDirichletBCClass( void ) = delete;

    /**
     * @brief Constructeur
     */
    ThermalDirichletBCClass(const ModelPtr& model)
        : DirichletBCClass( "_THER", model ){};

    /**
     * @brief Constructeur
     */
    ThermalDirichletBCClass( const std::string name, const ModelPtr& model )
        : DirichletBCClass( name, "_THER", model ){};

    /**
     * @typedef ThermalDirichletBCPtr
     * @brief Pointeur intelligent vers un ThermalDirichletBC
     */
    typedef boost::shared_ptr< ThermalDirichletBCClass > ThermalDirichletBCPtr;

    /**
     * @brief Ajout d'une valeur thermique imposee sur un groupe de mailles
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool
    addImposedThermalDOFOnCells( const PhysicalQuantityComponent &coordinate,
                                    const double &value,
                                    const std::string &nameOfGroup ) {
        if ( !_model->getMesh()->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

        MeshEntityPtr meshEnt( new GroupOfCells( nameOfGroup ) );
        RealLoadTemperature resu( meshEnt, coordinate, value );
        _listOfRealImposedTemperature.push_back( resu );
        return true;
    };

    /**
     * @brief Ajout d'une valeur thermique imposee sur un groupe de mailles
     * @param namesOfGroup Noms des groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addImposedThermalDOFOnCells(
        const PhysicalQuantityComponent &coordinate, const double &value,
        const VectorString &namesOfGroup ) {
        for ( const auto &nameOfGroup : namesOfGroup )
            addImposedThermalDOFOnCells( coordinate, value, nameOfGroup );
        return true;
    };

    /**
     * @brief Ajout d'une valeur thermique imposee sur un groupe de noeuds
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addImposedThermalDOFOnNodes( const PhysicalQuantityComponent &coordinate,
                                      const double &value,
                                      const std::string &nameOfGroup ) {
        if ( !_model->getMesh()->hasGroupOfNodes( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

        MeshEntityPtr meshEnt( new GroupOfNodes( nameOfGroup ) );
        RealLoadTemperature resu( meshEnt, coordinate, value );
        _listOfRealImposedTemperature.push_back( resu );
        return true;
    };

    /**
     * @brief Ajout d'une valeur thermique imposee sur un groupe de noeuds
     * @param namesOfGroup Noms des groupes sur lequels imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addImposedThermalDOFOnNodes(
        const PhysicalQuantityComponent &coordinate, const double &value,
        const VectorString &namesOfGroup ) {
        for ( const auto &nameOfGroup : namesOfGroup )
            addImposedThermalDOFOnNodes( coordinate, value, nameOfGroup );
        return true;
    };

    /**
     * @brief Ajout d'une fonction thermique imposee sur un groupe de noeuds
     * @param nameOfGroup Nom du groupe sur lequel imposer la fonction
     * @param FunctionPtr function imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addImposedThermalDOFOnNodes( const PhysicalQuantityComponent &coordinate,
                                      const FunctionPtr &function,
                                      const std::string &nameOfGroup ) {
        if ( !_model->getMesh()->hasGroupOfNodes( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

        MeshEntityPtr meshEnt( new GroupOfNodes( nameOfGroup ) );
        FunctionLoadTemperature resu( meshEnt, coordinate, function );
        _listOfFunctionImposedTemperature.push_back( resu );
        return true;
    };

    /**
     * @brief Ajout d'une fonction thermique imposee sur un groupe de noeuds
     * @param namesOfGroup Noms des groupes sur lequels imposer la valeur
     * @param FunctionPtr function imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addImposedThermalDOFOnNodes(
        const PhysicalQuantityComponent &coordinate, const FunctionPtr &function,
        const VectorString &namesOfGroup ) {
        for ( const auto &nameOfGroup : namesOfGroup )
            addImposedThermalDOFOnNodes( coordinate, function, nameOfGroup );
        return true;
    };
};

/**
 * @class AcousticDirichletBCClass
 * @brief Classe definissant une charge cinematique (issue d'AFFE_CHAR_CINE)
 * @author Nicolas Sellenet
 */
class AcousticDirichletBCClass : public DirichletBCClass {
  public:

  /**
     * @brief Constructeur
     */
    AcousticDirichletBCClass( void ) = delete;

    /**
     * @brief Constructeur
     */
    AcousticDirichletBCClass(const ModelPtr& model)
        : DirichletBCClass( "_ACOU", model ){};

    /**
     * @brief Constructeur
     */
    AcousticDirichletBCClass( const std::string name, const ModelPtr& model )
        : DirichletBCClass( name, "_ACOU", model ){};

    /**
     * @typedef AcousticDirichletBCPtr
     * @brief Pointeur intelligent vers un AcousticDirichletBC
     */
    typedef boost::shared_ptr< AcousticDirichletBCClass > AcousticDirichletBCPtr;

    /**
     * @brief Ajout d'une valeur acoustique imposee sur un groupe de mailles
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addImposedAcousticDOFOnCells( const std::string &nameOfGroup,
                                          const double &value ) {
        throw std::runtime_error( "Not yet implemented" );
    };

    /**
     * @brief Ajout d'une valeur acoustique imposee sur un groupe de noeuds
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addImposedAcousticDOFOnNodes( const std::string &nameOfGroup,
                                       double value ) {
        throw std::runtime_error( "Not yet implemented" );
    };
};

/**
 * @typedef DirichletBC
 * @brief Pointeur intelligent vers un DirichletBCClass
 */
typedef boost::shared_ptr< DirichletBCClass > DirichletBCPtr;

/**
 * @typedef MechanicalDirichletBCPtr
 * @brief Pointeur intelligent vers un MechanicalDirichletBC
 */
typedef boost::shared_ptr< MechanicalDirichletBCClass > MechanicalDirichletBCPtr;

/**
 * @typedef ThermalDirichletBCPtr
 * @brief Pointeur intelligent vers un ThermalDirichletBC
 */
typedef boost::shared_ptr< ThermalDirichletBCClass > ThermalDirichletBCPtr;

/**
 * @typedef AcousticDirichletBCPtr
 * @brief Pointeur intelligent vers un AcousticDirichletBC
 */
typedef boost::shared_ptr< AcousticDirichletBCClass > AcousticDirichletBCPtr;

/** @typedef std::list de DirichletBC */
typedef std::list< DirichletBCPtr > ListKineLoad;
/** @typedef Iterateur sur une std::list de DirichletBC */
typedef ListKineLoad::iterator ListKineLoadIter;
/** @typedef Iterateur constant sur une std::list de DirichletBC */
typedef ListKineLoad::const_iterator ListKineLoadCIter;

#endif /* KINEMATICSLOAD_H_ */
