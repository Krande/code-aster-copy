#ifndef RESULTS_H_
#define RESULTS_H_

/**
 * @file Result.h
 * @brief Fichier entete de la classe Result
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

#include "DataStructures/DataStructure.h"
#include "Meshes/Mesh.h"
#include "Modeling/Model.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/NamesMap.h"
#include "DataFields/FieldOnNodes.h"
#include "DataFields/FieldOnCells.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/ParallelDOFNumbering.h"
#include "Supervis/ResultNaming.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Loads/ListOfLoads.h"
#include "DataFields/FieldBuilder.h"
#include "DataFields/ListOfTables.h"

#include "Python.h"

/**
 * @class Result
 * @brief Cette classe correspond a la sd_resultat de Code_Aster, elle stocke des champs
 * @author Nicolas Sellenet
 */
class Result : public DataStructure, public ListOfTables {
  protected:
    typedef std::map< ASTERINTEGER, FieldOnNodesRealPtr > MapOfFieldOnNodesReal;
    typedef std::map< ASTERINTEGER, FieldOnCellsRealPtr > MapOfFieldOnCellsReal;
    typedef std::map< ASTERINTEGER, ConstantFieldOnCellsRealPtr > MapOfConstantFieldOnCellsReal;

    typedef std::map< ASTERINTEGER, FieldOnNodesComplexPtr > MapOfFieldOnNodesComplex;
    typedef std::map< ASTERINTEGER, FieldOnCellsComplexPtr > MapOfFieldOnCellsComplex;

    typedef std::map< ASTERINTEGER, FieldOnCellsLongPtr > MapOfFieldOnCellsLong;

    typedef std::map< ASTERINTEGER, ConstantFieldOnCellsChar16Ptr > MapOfConstantFieldOnCellsChar16;

    /** @typedef std::map d'une chaine et des pointers vers toutes les DataStructure */
    typedef std::map< std::string, MapOfFieldOnNodesReal > mapStrMoFNR;
    typedef std::map< std::string, MapOfFieldOnCellsReal > mapStrMoFCR;
    typedef std::map< std::string, MapOfConstantFieldOnCellsReal > mapStrMoCFCR;

    typedef std::map< std::string, MapOfFieldOnNodesComplex > mapStrMoFNC;
    typedef std::map< std::string, MapOfFieldOnCellsComplex > mapStrMoFCC;

    typedef std::map< std::string, MapOfFieldOnCellsLong > mapStrMoFCI;

    typedef std::map< std::string, MapOfConstantFieldOnCellsChar16 > mapStrMoCFCK16;

    /** @typedef std::map du rang et des pointers vers ElementaryCharacteristicsPtr */
    typedef std::map< int, ElementaryCharacteristicsPtr > mapRankCaraElem;
    /** @typedef std::map du rang et des pointers vers ListOfLoadsPtr */
    typedef std::map< int, ListOfLoadsPtr > mapRankLoads;
    /** @typedef std::map du rang et des pointers vers MaterialFieldPtr */
    typedef std::map< int, MaterialFieldPtr > mapRankMaterial;
    /** @typedef std::map du rang et des pointers vers ModelPtr */
    typedef std::map< int, ModelPtr > mapRankModel;

    /** @brief Pointeur de nom Jeveux '.DESC' */
    NamesMapChar16 _symbolicNamesOfFields;
    /** @brief Collection '.TACH' */
    JeveuxCollectionChar24 _namesOfFields;
    /** @brief Pointeur de nom Jeveux '.NOVA' */
    NamesMapChar16 _accessVariables;
    /** @brief Collection '.TAVA' */
    JeveuxCollectionChar8 _calculationParameter;
    /** @brief Vecteur Jeveux '.ORDR' */
    JeveuxVectorLong _serialNumber;
    /** @brief Vecteur Jeveux '.RSPI' */
    JeveuxVectorLong _rspi;
    /** @brief Vecteur Jeveux '.RSPR' */
    JeveuxVectorReal _rspr;
    /** @brief Vecteur Jeveux '.RSP8' */
    JeveuxVectorChar8 _rsp8;
    /** @brief Vecteur Jeveux '.RS16' */
    JeveuxVectorChar16 _rs16;
    /** @brief Vecteur Jeveux '.RS24' */
    JeveuxVectorChar24 _rs24;

    /** @brief Liste des champs aux noeuds */
    mapStrMoFNR _dictOfMapOfFieldOnNodesReal;
    mapStrMoFNC _dictOfMapOfFieldOnNodesComplex;
    /** @brief Liste des champs aux éléments */
    mapStrMoFCR _dictOfMapOfFieldOnCellsReal;
    mapStrMoFCC _dictOfMapOfFieldOnCellsComplex;
    mapStrMoFCI _dictOfMapOfFieldOnCellsLong;
    /** @brief Liste des cartes K16 */
    mapStrMoCFCR _dictOfMapOfConstantFieldOnCellsReal;
    mapStrMoCFCK16 _dictOfMapOfConstantFieldOnCellsChar16;
    /** @brief Liste des NUME_DDL */
    std::vector< BaseDOFNumberingPtr > _listOfDOFNum;
    /** @brief List of ElementaryCharacteristicsPtr */
    mapRankCaraElem _mapElemCara;
    /** @brief List of ListOfLoadsPtr */
    mapRankLoads _mapLoads;
    /** @brief List of MaterialFieldPtr */
    mapRankMaterial _mapMaterial;
    /** @brief List of ModelPtr */
    mapRankModel _mapModel;

    /** @brief Maillage sur lequel repose la resultat */
    BaseMeshPtr _mesh;
    /** @brief Object to correctly manage fields and field descriptions */
    FieldBuilder _fieldBuidler;

    /**
     * @brief Get a name for field (wrap to rsexch.F90)
     * @param name Symbolic name of the field
     * @param rank Rank
     */
    std::pair< ASTERINTEGER, std::string > _getNewFieldName( const std::string &name,
                                                             const ASTERINTEGER &rank ) const;

    void _checkMesh( const BaseMeshPtr mesh) const;

  public:
    /**
     * @typedef ResultPtr
     * @brief Pointeur intelligent vers un Result
     */
    typedef boost::shared_ptr< Result > ResultPtr;

    /**
     * @brief Constructeur
     */
    Result( const std::string &resuTyp )
        : Result( ResultNaming::getNewResultName(), resuTyp ){};

    /**
     * @brief Constructeur
     */
    Result( const std::string &name, const std::string &resuTyp )
        : DataStructure( name, 19, resuTyp ),
          ListOfTables( name ),
          _symbolicNamesOfFields( NamesMapChar16( getName() + ".DESC" ) ),
          _namesOfFields( JeveuxCollectionChar24( getName() + ".TACH" ) ),
          _accessVariables( NamesMapChar16( getName() + ".NOVA" ) ),
          _calculationParameter( JeveuxCollectionChar8( getName() + ".TAVA" ) ),
          _serialNumber( JeveuxVectorLong( getName() + ".ORDR" ) ),
          _rspi( JeveuxVectorLong( getName() + ".RSPI" ) ),
          _rspr( JeveuxVectorReal( getName() + ".RSPR" ) ),
          _rsp8( JeveuxVectorChar8( getName() + ".RSP8" ) ),
          _rs16( JeveuxVectorChar16( getName() + ".RS16" ) ),
          _rs24( JeveuxVectorChar24( getName() + ".RS24" ) ),
          _mesh( nullptr ),
          _fieldBuidler( FieldBuilder() ){};

    /**
     * @brief Allouer une sd_resultat
     * @param nbRanks nombre de numéro d'ordre
     * @return true si l'allocation s'est bien passée
     */
    bool allocate( int nbRanks ) ;

    /**
     * @brief Add elementary characteristics to container
     * @param rank
     */
    void addElementaryCharacteristics( const ElementaryCharacteristicsPtr &,
                                       int rank ) ;

    /**
     * @brief Add a existing FieldOnNodesDescription in _fieldBuidler
     */
    void addFieldOnNodesDescription( const FieldOnNodesDescriptionPtr &fond )
    {
        _fieldBuidler.addFieldOnNodesDescription( fond );
    };

    /**
     * @brief Add elementary characteristics to container
     * @param rank
     */
    void addListOfLoads( const ListOfLoadsPtr &, int rank ) ;

    /**
     * @brief Add material definition
     * @param rank
     */
    void addMaterialField( const MaterialFieldPtr &, int rank ) ;

    /**
     * @brief Add model
     * @param rank
     */
    void addModel( const ModelPtr &, int rank ) ;

    /**
     * @brief Set model
     */
    void setMesh( const BaseMeshPtr &mesh );

    /**
     * @brief Add time value for one rank
     * @param rank
     */
    void addTimeValue( ASTERDOUBLE, int rank );

    /**
     * @brief Append a elementary characteristics on all rank of Result
     * @param ElementaryCharacteristicsPtr
     */
    void appendElementaryCharacteristicsOnAllRanks( const ElementaryCharacteristicsPtr& );

    /**
     * @brief Append a material on all rank of Result
     * @param MaterialFieldPtr
     */
    void appendMaterialFieldOnAllRanks( const MaterialFieldPtr & );

    /**
     * @brief Append a model on all rank of Result
     * @param ModelPtr
     */
    void appendModelOnAllRanks( const ModelPtr & );

    /**
     * @brief Obtenir un DOFNumbering à remplir
     * @return DOFNumbering à remplir
     */
    BaseDOFNumberingPtr getEmptyDOFNumbering();

/**
 * @brief Obtenir un DOFNumbering à remplir
 * @return DOFNumbering à remplir
 */
#ifdef ASTER_HAVE_MPI
    BaseDOFNumberingPtr getEmptyParallelDOFNumbering();
#endif /* ASTER_HAVE_MPI */

    /**
     * @brief Obtenir le dernier DOFNumbering
     * @return Dernier DOFNumbering
     */
    BaseDOFNumberingPtr getLastDOFNumbering() const {
        return _listOfDOFNum[_listOfDOFNum.size() - 1];
    };

    /**
     * @brief Add elementary characteristics to container
     * @param rank
     */
    ListOfLoadsPtr getListOfLoads( int rank ) ;

    /**
     * @brief Get elementary characteristics
     */
    ElementaryCharacteristicsPtr
    getElementaryCharacteristics() ;

    /**
     * @brief Get elementary characteristics
     */
    std::vector< ElementaryCharacteristicsPtr >
    getAllElementaryCharacteristics() const;

    bool hasElementaryCharacteristics() const;

    /**
     * @brief Get elementary characteristics
     * @param rank
     */
    ElementaryCharacteristicsPtr
    getElementaryCharacteristics( int rank ) ;

    /**
     * @brief Get elementary characteristics
     * @param rank
     */
    bool hasElementaryCharacteristics( ASTERINTEGER rank ) const;

    /**
     * @brief Get material
     */
    std::vector< MaterialFieldPtr > getMaterialFields() const;

    /**
     * @brief Get material
     */
    MaterialFieldPtr getMaterialField() ;

    /**
     * @brief Get material
     * @param rank
     */
    bool hasMaterialField( const ASTERINTEGER &rank ) const;

    /**
     * @brief Get material
     * @param rank
     */
    MaterialFieldPtr getMaterialField( int rank ) ;

    /**
     * @brief Get mesh
     */
    BaseMeshPtr getMesh();

    /**
     * @brief check for multiple models
     */
    bool hasMultipleModel() ;

    /**
     * @brief check for multiple models
     */
    bool hasModel( const ASTERINTEGER &rank ) const;

    /**
     * @brief Get models
     */
    std::vector< ModelPtr > getModels() const;

    /**
     * @brief Get model
     */
    ModelPtr getModel() ;

    /**
     * @brief Get model
     * @param rank
     */
    ModelPtr getModel( int rank ) ;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    FieldOnCellsRealPtr getFieldOnCellsReal( const std::string name,
                                             const ASTERINTEGER rank ) const;

    FieldOnCellsComplexPtr getFieldOnCellsComplex( const std::string name,
                                                   const ASTERINTEGER rank ) const;

    FieldOnCellsLongPtr getFieldOnCellsLong( const std::string name,
                                             const ASTERINTEGER rank ) const;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    ConstantFieldOnCellsChar16Ptr getConstantFieldOnCellsChar16( const std::string name,
                                                                 const ASTERINTEGER rank ) const;

    ConstantFieldOnCellsRealPtr getConstantFieldOnCellsReal( const std::string name,
                                                             const ASTERINTEGER rank ) const;

    /**
     * @brief Ajouter un champ par éléments réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    bool setField( const FieldOnCellsRealPtr field, const std::string &name,
                   const ASTERINTEGER rank );

    bool setField( const FieldOnCellsComplexPtr field, const std::string &name,
                   const ASTERINTEGER rank );

    bool setField( const FieldOnCellsLongPtr field, const std::string &name,
                   const ASTERINTEGER rank );

    /**
     * @brief Ajouter un champ par éléments réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    bool setField( const ConstantFieldOnCellsRealPtr field, const std::string &name,
                   const ASTERINTEGER rank );

    bool setField( const ConstantFieldOnCellsChar16Ptr field, const std::string &name,
                   const ASTERINTEGER rank );

    /**
    * @brief Get dict of access variables and their values
    * @return PyObject
    */
    PyObject *getAccessParameters() const;

    /**
     * @brief Get the list of fields on nodes
     * @return std::vector< string >
     */
    VectorString getFieldsOnNodesRealNames() const;

    VectorString getFieldsOnNodesComplexNames() const;

    /**
     * @brief Get the list of fields on elements
     * @return std::vector< string >
     */
    VectorString getFieldsOnCellsRealNames() const;

    VectorString getFieldsOnCellsComplexNames() const;

    VectorString getFieldsOnCellsLongNames() const;

    /**
     * @brief Get the list of constant fields on cells
     * @return std::vector< string >
     */
    VectorString getConstantFieldsOnCellsRealNames() const;

    VectorString getConstantFieldsOnCellsChar16Names() const;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnNodesRealPtr pointant vers le champ
     */
    FieldOnNodesRealPtr getFieldOnNodesReal( const std::string name,
                                             const ASTERINTEGER rank ) const;

    FieldOnNodesComplexPtr getFieldOnNodesComplex( const std::string name,
                                                   const ASTERINTEGER rank ) const;

    /**
     * @brief Ajouter un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnNodesRealPtr pointant vers le champ
     */
    bool setField( const FieldOnNodesRealPtr field, const std::string &name,
                   const ASTERINTEGER rank );

    bool setField( const FieldOnNodesComplexPtr field, const std::string &name,
                   const ASTERINTEGER rank );

    /**
     * @brief Impression de la sd au format MED
     * @param fileName Nom du fichier MED à imprimer
     * @return true
     * @todo revoir la gestion des mot-clés par défaut (ex : TOUT_ORDRE)
     * @todo revoir la gestion des unités logiques (notamment si fort.20 existe déjà)
     */
    bool printMedFile( const std::string fileName, std::string medName ) const ;

    bool printMedFile( const std::string fileName ) const
    { return printMedFile(fileName, std::string());} ;


    /**
    * @brief Get the number of steps stored in the Result
    * @return nbRanks
    */
    int getNumberOfRanks() const;

    /**
    * @brief Get the number of steps stored in the Result
    * @return nbRanks
    */
    VectorLong getRanks() const;

    /**
     * @brief Get all the fields stored in the Result
     * @return VectorString
     */
    VectorString getFieldsNames() const;

    /**
    * @brief Print all the fields stored in the Result
    * @return nbRanks
    */
    void listFields() const;

    /**
    * @brief Print informations about the Result content
    */
    void printInfo() const;

    /**
     * @brief Construire une sd_resultat à partir d'objet produit dans le Fortran
     * @return true si l'allocation s'est bien passée
     * @todo revoir l'agrandissement de dictOfVectorOfFieldsNodes et dictOfVectorOfFieldsCells
     */
    bool build() ;

    /**
     * @brief Update the  Result's size
     */
    bool resize( ASTERINTEGER nbRanks );
};

/**
 * @typedef ResultPtr
 * @brief Pointeur intelligent vers un Result
 */
typedef boost::shared_ptr< Result > ResultPtr;

#endif /* RESULTS_H_ */
