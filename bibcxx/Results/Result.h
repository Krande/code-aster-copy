#ifndef RESULTS_H_
#define RESULTS_H_

/**
 * @file Result.h
 * @brief Fichier entete de la classe Result
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "Python.h"
#include "astercxx.h"

#include "DataFields/FieldBuilder.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "DataFields/ListOfTables.h"
#include "DataStructures/DataStructure.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Loads/ListOfLoads.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "MemoryManager/NamesMap.h"
#include "Meshes/Mesh.h"
#include "Modeling/Model.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/ParallelDOFNumbering.h"
#include "Supervis/ResultNaming.h"

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
    typedef std::map< ASTERINTEGER, ElementaryCharacteristicsPtr > mapIndexCaraElem;
    /** @typedef std::map du rang et des pointers vers ListOfLoadsPtr */
    typedef std::map< ASTERINTEGER, ListOfLoadsPtr > mapIndexLoads;
    /** @typedef std::map du rang et des pointers vers MaterialFieldPtr */
    typedef std::map< ASTERINTEGER, MaterialFieldPtr > mapIndexMaterial;
    /** @typedef std::map du rang et des pointers vers ModelPtr */
    typedef std::map< ASTERINTEGER, ModelPtr > mapIndexModel;

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

    /** @brief List of ElementaryCharacteristicsPtr */
    mapIndexCaraElem _mapElemCara;
    /** @brief List of ListOfLoadsPtr */
    mapIndexLoads _mapLoads;
    /** @brief List of MaterialFieldPtr */
    mapIndexMaterial _mapMaterial;
    /** @brief List of ModelPtr */
    mapIndexModel _mapModel;

    /** @brief Maillage sur lequel repose la resultat */
    BaseMeshPtr _mesh;
    /** @brief Object to correctly manage fields and field descriptions */
    FieldBuilder _fieldBuidler;

    /**
     * @brief Get a name for field (wrap to rsexch.F90)
     * @param name Symbolic name of the field
     * @param index Index
     */
    std::pair< ASTERINTEGER, std::string > _getNewFieldName( const std::string &name,
                                                             const ASTERINTEGER &index ) const;

    template < typename T >
    void
    _setFieldBase( const std::string &name, const ASTERINTEGER &index, std::shared_ptr< T > field,
                   std::map< std::string, std::map< ASTERINTEGER, std::shared_ptr< T > > > &dict );

    void _checkMesh( const BaseMeshPtr mesh ) const;

    ASTERINTEGER _getInternalIndex( const ASTERINTEGER &index ) const;

  public:
    /**
     * @typedef ResultPtr
     * @brief Pointeur intelligent vers un Result
     */
    typedef std::shared_ptr< Result > ResultPtr;

    /**
     * @brief Constructeur
     */
    Result( const std::string &resuTyp ) : Result( ResultNaming::getNewResultName(), resuTyp ){};

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
     * @brief Add a FiniteElementDescriptor to elementary matrix
     * @param FiniteElementDescriptorPtr FiniteElementDescriptor
     */
    bool addFiniteElementDescriptor( const FiniteElementDescriptorPtr curFED );

    /**
     * @brief Allouer une sd_resultat
     * @param nbIndexes nombre de numéro d'ordre
     * @return true si l'allocation s'est bien passée
     */
    void allocate( ASTERINTEGER nbIndexes );

    /**
     * @brief Add elementary characteristics to container
     * @param index
     */
    void setElementaryCharacteristics( const ElementaryCharacteristicsPtr &, ASTERINTEGER index );

    /**
     * @brief Add a existing GlobalEquationNumbering in _fieldBuidler
     */
    void addGlobalEquationNumbering( const GlobalEquationNumberingPtr &fond ) {
        _fieldBuidler.addGlobalEquationNumbering( fond );
    };

    /**
     * @brief Set list of loads at index
     * @param ListOfLoadsPtr, index
     */
    void setListOfLoads( const ListOfLoadsPtr &, ASTERINTEGER index );

    /**
     * @brief Add material definition
     * @param index
     */
    void setMaterialField( const MaterialFieldPtr &, ASTERINTEGER index );

    /**
     * @brief Add model
     * @param index
     */
    void setModel( const ModelPtr &, ASTERINTEGER index );

    /**
     * @brief Set model
     */
    void setMesh( const BaseMeshPtr &mesh );

    /**
     * @brief Add time value for one index
     * @param index
     */
    void setTimeValue( ASTERDOUBLE value, ASTERINTEGER index ) {
        this->setParameterValue( "INST", value, index );
    };

    /**
     * @brief Add parameter value for one index
     */
    void setParameterValue( std::string name, ASTERDOUBLE value, ASTERINTEGER index );

    ASTERDOUBLE getTimeValue( ASTERINTEGER index );

    /**
     * @brief Append a elementary characteristics on all index of Result
     * @param ElementaryCharacteristicsPtr
     */
    void setElementaryCharacteristics( const ElementaryCharacteristicsPtr & );

    /**
     * @brief Append a material on all index of Result
     * @param MaterialFieldPtr
     */
    void setMaterialField( const MaterialFieldPtr & );

    /**
     * @brief Append a model on all index of Result
     * @param ModelPtr
     */
    void setModel( const ModelPtr & );

    /**
     * @brief Get list of loads at index
     * @param index
     */
    ListOfLoadsPtr getListOfLoads( ASTERINTEGER index ) const;

    bool hasListOfLoads( const ASTERINTEGER &index ) const;

    bool hasListOfLoads() const;

    /**
     * @brief Get elementary characteristics
     */
    ElementaryCharacteristicsPtr getElementaryCharacteristics() const;

    /**
     * @brief Get elementary characteristics
     */
    std::vector< ElementaryCharacteristicsPtr > getAllElementaryCharacteristics() const;

    /**
     * @brief Get elementary characteristics
     * @param index
     */
    ElementaryCharacteristicsPtr getElementaryCharacteristics( ASTERINTEGER index ) const;

    /**
     * @brief Get elementary characteristics
     * @param index
     */
    bool hasElementaryCharacteristics( ASTERINTEGER index ) const;

    bool hasElementaryCharacteristics() const;

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
     * @param index
     */
    bool hasMaterialField( const ASTERINTEGER &index ) const;

    /**
     * @brief Get material
     * @param index
     */
    MaterialFieldPtr getMaterialField( ASTERINTEGER index ) const;

    /**
     * @brief Get mesh
     */
    BaseMeshPtr getMesh() const;

    /**
     * @brief check for multiple models
     */
    bool hasMultipleModel() const;

    /**
     * @brief check for multiple models
     */
    bool hasModel( const ASTERINTEGER &index ) const;

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
     * @param index
     */
    ModelPtr getModel( ASTERINTEGER index ) const;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param index numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    FieldOnCellsRealPtr getFieldOnCellsReal( const std::string name,
                                             const ASTERINTEGER index ) const;

    FieldOnCellsComplexPtr getFieldOnCellsComplex( const std::string name,
                                                   const ASTERINTEGER index ) const;

    FieldOnCellsLongPtr getFieldOnCellsLong( const std::string name,
                                             const ASTERINTEGER index ) const;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param index numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    ConstantFieldOnCellsChar16Ptr getConstantFieldOnCellsChar16( const std::string name,
                                                                 const ASTERINTEGER index ) const;

    ConstantFieldOnCellsRealPtr getConstantFieldOnCellsReal( const std::string name,
                                                             const ASTERINTEGER index ) const;

    /**
     * @brief Ajouter un champ par éléments réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param index numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    void setField( const FieldOnCellsRealPtr field, const std::string &name,
                   const ASTERINTEGER index );

    void setField( const FieldOnCellsComplexPtr field, const std::string &name,
                   const ASTERINTEGER index );

    void setField( const FieldOnCellsLongPtr field, const std::string &name,
                   const ASTERINTEGER index );

    /**
     * @brief Ajouter un champ par éléments réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param index numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    void setField( const ConstantFieldOnCellsRealPtr field, const std::string &name,
                   const ASTERINTEGER index );

    void setField( const ConstantFieldOnCellsChar16Ptr field, const std::string &name,
                   const ASTERINTEGER index );

    /**
     * @brief Get dict of access variables and their values
     * @return py::dict
     */
    py::dict getAccessParameters() const;

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
     * @param index numéro d'ordre
     * @return FieldOnNodesRealPtr pointant vers le champ
     */
    FieldOnNodesRealPtr getFieldOnNodesReal( const std::string name,
                                             const ASTERINTEGER index ) const;

    FieldOnNodesComplexPtr getFieldOnNodesComplex( const std::string name,
                                                   const ASTERINTEGER index ) const;

    /**
     * @brief Ajouter un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param index numéro d'ordre
     * @return FieldOnNodesRealPtr pointant vers le champ
     */
    void setField( const FieldOnNodesRealPtr field, const std::string &name,
                   const ASTERINTEGER index );

    void setField( const FieldOnNodesComplexPtr field, const std::string &name,
                   const ASTERINTEGER index );

    /**
     * @brief Impression de la sd au format MED
     * @param fileName Nom du fichier MED à imprimer
     * @return true
     * @todo revoir la gestion des mot-clés par défaut (ex : TOUT_ORDRE)
     * @todo revoir la gestion des unités logiques (notamment si fort.20 existe déjà)
     */
    void printMedFile( const std::string fileName, std::string medName = std::string(),
                       bool local = true ) const;

    /**
     * @brief Get the number of steps stored in the Result
     * @return nbIndexes
     */
    ASTERINTEGER getNumberOfIndexes() const;

    /**
     * @brief Get the number of steps stored in the Result
     * @return nbIndexes
     */
    VectorLong getIndexes() const;

    /**
     * @brief Get all the fields stored in the Result
     * @return VectorString
     */
    VectorString getFieldsNames() const;

    /**
     * @brief Print all the fields stored in the Result
     * @return nbIndexes
     */
    void printListOfFields() const;

    /**
     * @brief Print informations about the Result content
     */
    void printInfo() const;

    /**
     * @brief Construire une sd_resultat à partir d'objet produit dans le Fortran
     * @return true si l'allocation s'est bien passée
     * @todo revoir l'agrandissement de dictOfMapOfFieldOnNodesReal et
     *  dictOfMapOfFieldOnCellsReal
     */
    virtual bool build( const std::vector< FiniteElementDescriptorPtr > feds =
                            std::vector< FiniteElementDescriptorPtr >(),
                        const std::vector< GlobalEquationNumberingPtr > fnds =
                            std::vector< GlobalEquationNumberingPtr >() );

    /**
     * @brief Update the  Result's size
     */
    void resize( ASTERINTEGER nbIndexes );

    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptors() const;

    std::vector< GlobalEquationNumberingPtr > getGlobalEquationNumberings() const;

    void clear( const ASTERINTEGER &index );

    void clear();

    bool exists() const;
};

/**
 * @typedef ResultPtr
 * @brief Pointeur intelligent vers un Result
 */
typedef std::shared_ptr< Result > ResultPtr;

#endif /* RESULTS_H_ */
