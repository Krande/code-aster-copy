#ifndef DIRICHLETBCSLOAD_H_
#define DIRICHLETBCSLOAD_H_

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

    /** @typedef std::list de DisplacementLoadReal */
    typedef std::list< DisplacementLoadReal > ListDispReal;
    /** @typedef Iterateur sur ListDispReal */
    typedef ListDispReal::iterator ListDispRealIter;

    /** @typedef std::list de TemperatureLoadReal */
    typedef std::list< TemperatureLoadReal > ListDispTemp;
    typedef std::list< TemperatureLoadFunction > ListFunctionTemp;
    /** @typedef Iterateur sur ListDispTemp */
    typedef ListDispTemp::iterator ListDispTempIter;

    /** @brief Modele */
    ModelPtr _model;
    /** @brief Listes des valeurs imposees DEPL_R et TEMP_R */
    ListDispReal _listOfRealImposedDisplacement;
    ListDispTemp _listOfRealImposedTemperature;
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

    virtual int getPhysics( void ) const {
        throw std::runtime_error( "Not allowed" );
    };

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

    virtual int getPhysics( void ) const {
        return Physics::Mechanics;
    };

    /**
     * @brief Ajout d'une valeur mecanique imposee sur un groupe de mailles
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addBCOnCells(
        const PhysicalQuantityComponent &coordinate, const ASTERDOUBLE &value,
        const std::string &nameOfGroup ) {
        if ( !_model->getMesh()->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

        MeshEntityPtr meshEnt( new GroupOfCells( nameOfGroup ) );
        DisplacementLoadReal resu( meshEnt, coordinate, value );
        _listOfRealImposedDisplacement.push_back( resu );
        return true;
    };

    /**
     * @brief Ajout d'une valeur mecanique imposee sur un groupe de mailles
     * @param namesOfGroup Noms des groupes sur lequels imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addBCOnCells(
        const PhysicalQuantityComponent &coordinate, const ASTERDOUBLE &value,
        const VectorString &namesOfGroup ) {
        for ( const auto &nameOfGroup : namesOfGroup )
            addBCOnCells( coordinate, value, nameOfGroup );
        return true;
    };

    /**
     * @brief Ajout d'une valeur mecanique imposee sur un groupe de noeuds
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool
    addBCOnNodes( const PhysicalQuantityComponent &coordinate,
                                    const ASTERDOUBLE &value,
                                    const std::string &nameOfGroup ) {
        if ( !_model->getMesh()->hasGroupOfNodes( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

        MeshEntityPtr meshEnt( new GroupOfNodes( nameOfGroup ) );
        DisplacementLoadReal resu( meshEnt, coordinate, value );
        _listOfRealImposedDisplacement.push_back( resu );
        return true;
    };

    /**
     * @brief Ajout d'une valeur mecanique imposee sur un groupe de noeuds
     * @param namesOfGroup Noms des groupe sur lequels imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addBCOnNodes(
        const PhysicalQuantityComponent &coordinate, const ASTERDOUBLE &value,
        const VectorString &namesOfGroup ) {
        for ( const auto &nameOfGroup : namesOfGroup )
            addBCOnNodes( coordinate, value, nameOfGroup );
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

    virtual int getPhysics( void ) const {
        return Physics::Thermal;
    };

    /**
     * @brief Ajout d'une valeur thermique imposee sur un groupe de mailles
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool
    addBCOnCells( const PhysicalQuantityComponent &coordinate,
                                    const ASTERDOUBLE &value,
                                    const std::string &nameOfGroup ) {
        if ( !_model->getMesh()->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

        MeshEntityPtr meshEnt( new GroupOfCells( nameOfGroup ) );
        TemperatureLoadReal resu( meshEnt, coordinate, value );
        _listOfRealImposedTemperature.push_back( resu );
        return true;
    };

    /**
     * @brief Ajout d'une valeur thermique imposee sur un groupe de mailles
     * @param namesOfGroup Noms des groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addBCOnCells(
        const PhysicalQuantityComponent &coordinate, const ASTERDOUBLE &value,
        const VectorString &namesOfGroup ) {
        for ( const auto &nameOfGroup : namesOfGroup )
            addBCOnCells( coordinate, value, nameOfGroup );
        return true;
    };

    /**
     * @brief Ajout d'une valeur thermique imposee sur un groupe de noeuds
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addBCOnNodes( const PhysicalQuantityComponent &coordinate,
                                      const ASTERDOUBLE &value,
                                      const std::string &nameOfGroup ) {
        if ( !_model->getMesh()->hasGroupOfNodes( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

        MeshEntityPtr meshEnt( new GroupOfNodes( nameOfGroup ) );
        TemperatureLoadReal resu( meshEnt, coordinate, value );
        _listOfRealImposedTemperature.push_back( resu );
        return true;
    };

    /**
     * @brief Ajout d'une valeur thermique imposee sur un groupe de noeuds
     * @param namesOfGroup Noms des groupes sur lequels imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addBCOnNodes(
        const PhysicalQuantityComponent &coordinate, const ASTERDOUBLE &value,
        const VectorString &namesOfGroup ) {
        for ( const auto &nameOfGroup : namesOfGroup )
            addBCOnNodes( coordinate, value, nameOfGroup );
        return true;
    };

    /**
     * @brief Ajout d'une fonction thermique imposee sur un groupe de noeuds
     * @param nameOfGroup Nom du groupe sur lequel imposer la fonction
     * @param FunctionPtr function imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addBCOnNodes( const PhysicalQuantityComponent &coordinate,
                                      const FunctionPtr &function,
                                      const std::string &nameOfGroup ) {
        if ( !_model->getMesh()->hasGroupOfNodes( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );

        MeshEntityPtr meshEnt( new GroupOfNodes( nameOfGroup ) );
        TemperatureLoadFunction resu( meshEnt, coordinate, function );
        _listOfFunctionImposedTemperature.push_back( resu );
        return true;
    };

    /**
     * @brief Ajout d'une fonction thermique imposee sur un groupe de noeuds
     * @param namesOfGroup Noms des groupes sur lequels imposer la valeur
     * @param FunctionPtr function imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addBCOnNodes(
        const PhysicalQuantityComponent &coordinate, const FunctionPtr &function,
        const VectorString &namesOfGroup ) {
        for ( const auto &nameOfGroup : namesOfGroup )
            addBCOnNodes( coordinate, function, nameOfGroup );
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

    virtual int getPhysics( void ) const {
        return Physics::Acoustic;
    };

    /**
     * @brief Ajout d'une valeur acoustique imposee sur un groupe de mailles
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addBCOnCells( const std::string &nameOfGroup,
                                          const ASTERDOUBLE &value ) {
        throw std::runtime_error( "Not yet implemented" );
    };

    /**
     * @brief Ajout d'une valeur acoustique imposee sur un groupe de noeuds
     * @param nameOfGroup Nom du groupe sur lequel imposer la valeur
     * @param value Valeur imposee
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool addBCOnNodes( const std::string &nameOfGroup,
                                       ASTERDOUBLE value ) {
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
typedef std::list< DirichletBCPtr > ListDiriBC;
/** @typedef Iterateur sur une std::list de DirichletBC */
typedef ListDiriBC::iterator ListDiriBCIter;
/** @typedef Iterateur constant sur une std::list de DirichletBC */
typedef ListDiriBC::const_iterator ListDiriBCCIter;

#endif /* DIRICHLETBCSLOAD_H_ */
