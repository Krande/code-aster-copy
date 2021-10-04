#ifndef RESULTS_H_
#define RESULTS_H_

/**
 * @file Result.h
 * @brief Fichier entete de la classe Result
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
    typedef std::vector< FieldOnNodesRealPtr > VectorOfFieldOnNodesReal;
    typedef std::vector< FieldOnCellsRealPtr > VectorOfFieldOnCellsReal;
    typedef std::vector< ConstantFieldOnCellsChar16Ptr > VectorOfConstantFieldOnCellsChar16;


    /** @typedef std::map d'une chaine et des pointers vers toutes les DataStructure */
    typedef std::map< std::string, VectorOfFieldOnNodesReal > mapStrVOFN;
    /** @typedef Iterateur sur le std::map */
    typedef mapStrVOFN::iterator mapStrVOFNIterator;
    /** @typedef Valeur contenue dans mapStrVOFN */
    typedef mapStrVOFN::value_type mapStrVOFNValue;

    /** @typedef std::map d'une chaine et des pointers vers toutes les DataStructure */
    typedef std::map< std::string, VectorOfFieldOnCellsReal > mapStrVOFE;
    /** @typedef Iterateur sur le std::map */
    typedef mapStrVOFE::iterator mapStrVOFEIterator;
    /** @typedef Valeur contenue dans mapStrVOFE */
    typedef mapStrVOFE::value_type mapStrVOFEValue;

        /** @typedef std::map d'une chaine et des pointers vers toutes les DataStructure */
    typedef std::map< std::string, VectorOfConstantFieldOnCellsChar16 > mapStrVOCF16;
    /** @typedef Iterateur sur le std::map */
    typedef mapStrVOCF16::iterator mapStrVOCF16Iterator;
    /** @typedef Valeur contenue dans mapStrVOFE */
    typedef mapStrVOCF16::value_type mapStrVOCF16Value;

    /** @typedef std::map du rang et des pointers vers ElementaryCharacteristicsPtr */
    typedef std::map< ASTERINTEGER, ElementaryCharacteristicsPtr > mapRankCaraElem;
    /** @typedef std::map du rang et des pointers vers ListOfLoadsPtr */
    typedef std::map< ASTERINTEGER, ListOfLoadsPtr > mapRankLoads;
    /** @typedef std::map du rang et des pointers vers MaterialFieldPtr */
    typedef std::map< ASTERINTEGER, MaterialFieldPtr > mapRankMaterial;
    /** @typedef std::map du rang et des pointers vers ModelPtr */
    typedef std::map< ASTERINTEGER, ModelPtr > mapRankModel;

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
    /** @brief Nombre de numéros d'ordre */
    ASTERINTEGER _nbRanks;
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
    mapStrVOFN _dictOfVectorOfFieldOnNodesReal;
    /** @brief Liste des champs aux éléments */
    mapStrVOFE _dictOfVectorOfFieldOnCellsReal;
    /** @brief Liste des cartes K16 */
    mapStrVOCF16 _dictOfVectorOfConstantFieldOnCellsChar16;
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
    std::pair< ASTERINTEGER, std::string> _getNewFieldName( const std::string& name,
                                                            const ASTERINTEGER& rank ) const;

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
          _serialNumber( JeveuxVectorLong( getName() + ".ORDR" ) ), _nbRanks( 0 ),
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
    bool allocate( ASTERINTEGER nbRanks ) ;

    /**
     * @brief Add elementary characteristics to container
     * @param rank
     */
    void setElementaryCharacteristics( const ElementaryCharacteristicsPtr &,
                                       ASTERINTEGER rank ) ;

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
    void setListOfLoads( const ListOfLoadsPtr &, ASTERINTEGER rank ) ;

    /**
     * @brief Add material definition
     * @param rank
     */
    void setMaterialField( const MaterialFieldPtr &, ASTERINTEGER rank ) ;

    /**
     * @brief Add model
     * @param rank
     */
    void setModel( const ModelPtr &, ASTERINTEGER rank ) ;

    /**
     * @brief Set model
     */
    void setMesh( const BaseMeshPtr &mesh ) { _mesh = mesh; };

    /**
     * @brief Add time value for one rank
     * @param rank
     */
    void setTimeValue( ASTERDOUBLE, ASTERINTEGER rank );

    ASTERDOUBLE getTimeValue( ASTERINTEGER rank );

    /**
     * @brief Append a elementary characteristics on all rank of Result
     * @param ElementaryCharacteristicsPtr
     */
    void setElementaryCharacteristics( const ElementaryCharacteristicsPtr& );

    /**
     * @brief Append a material on all rank of Result
     * @param MaterialFieldPtr
     */
    void setMaterialField( const MaterialFieldPtr & );

    /**
     * @brief Append a model on all rank of Result
     * @param ModelPtr
     */
    void setModel( const ModelPtr & );

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
     * @brief Obtenir un champ aux noeuds réel vide à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnNodesRealPtr pointant vers le champ
     */
    FieldOnNodesRealPtr getEmptyFieldOnNodesReal( const std::string name,
                                                      const ASTERINTEGER rank ) ;

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
    ListOfLoadsPtr getListOfLoads( ASTERINTEGER rank ) const;

    /**
     * @brief Get elementary characteristics
     */
    ElementaryCharacteristicsPtr
    getElementaryCharacteristics() const;

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
    getElementaryCharacteristics( ASTERINTEGER rank ) const;

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
    MaterialFieldPtr getMaterialField() const;

    /**
     * @brief Get material
     * @param rank
     */
    MaterialFieldPtr getMaterialField( ASTERINTEGER rank ) const;

    /**
     * @brief Get mesh
     */
    BaseMeshPtr getMesh() const;

    /**
     * @brief check for multiple models
     */
    bool hasMultipleModel() const;

    /**
     * @brief Get models
     */
    std::vector< ModelPtr > getModels() const;

    /**
     * @brief Get model
     */
    ModelPtr getModel() const;

    /**
     * @brief Get model
     * @param rank
     */
    ModelPtr getModel( ASTERINTEGER rank ) const;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    FieldOnCellsRealPtr getFieldOnCellsReal( const std::string name, const ASTERINTEGER rank )
    const;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    ConstantFieldOnCellsChar16Ptr getConstantFieldOnCellsChar16(
        const std::string name, const ASTERINTEGER rank )
    const;

    /**
     * @brief Ajouter un champ par éléments réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    bool setField( const FieldOnCellsRealPtr field, const std::string& name,
        const ASTERINTEGER rank );

     /**
     * @brief Ajouter un champ par éléments réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    bool setField( const ConstantFieldOnCellsChar16Ptr field, const std::string& name,
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
    VectorString getFieldsOnNodesNames() const;

    /**
    * @brief Get the list of fields on elements
    * @return std::vector< string >
    */
    VectorString getFieldsOnCellsNames() const;

    /**
    * @brief Get the list of constant fields on cells
    * @return std::vector< string >
    */
    VectorString getConstantFieldsOnCellsNames() const;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnNodesRealPtr pointant vers le champ
     */
    FieldOnNodesRealPtr getFieldOnNodesReal( const std::string name, const ASTERINTEGER rank )
    const;

    /**
     * @brief Ajouter un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param rank numéro d'ordre
     * @return FieldOnNodesRealPtr pointant vers le champ
     */
    bool setField( const FieldOnNodesRealPtr field,
                          const std::string& name, const ASTERINTEGER rank );

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
    ASTERINTEGER getNumberOfRanks() const;

    /**
    * @brief Get the number of steps stored in the Result
    * @return nbRanks
    */
    VectorLong getRanks() const;

    /**
    * @brief Print all the fields stored in the Result
    * @return nbRanks
    */
    void printListOfFields() const;

    /**
    * @brief Print informations about the Result content
    */
    void printInfo() const;

    /**
     * @brief Construire une sd_resultat à partir d'objet produit dans le Fortran
     * @return true si l'allocation s'est bien passée
     * @todo revoir l'agrandissement de dictOfVectorOfFieldOnNodesReal et
     *  dictOfVectorOfFieldOnCellsReal
     */
    bool build() ;

    /**
    * @brief Update the  Result's size
    */
    ASTERBOOL resize(ASTERINTEGER nbRanks);

};

/**
 * @typedef ResultPtr
 * @brief Pointeur intelligent vers un Result
 */
typedef boost::shared_ptr< Result > ResultPtr;

#endif /* RESULTS_H_ */
