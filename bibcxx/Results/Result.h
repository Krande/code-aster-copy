#ifndef RESULTS_H_
#define RESULTS_H_

/**
 * @file Result.h
 * @brief Fichier entete de la classe Result
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "DataFields/FieldBuilder.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "LinearAlgebra/GeneralizedAssemblyVector.h"
#include "Loads/ListOfLoads.h"
#include "Materials/MaterialField.h"
#include "Meshes/Mesh.h"
#include "Modeling/Model.h"

#include <filesystem>

/**
 * @class Result
 * @brief Cette classe correspond a la sd_resultat de Code_Aster, elle stocke des champs
 * @author Nicolas Sellenet
 */
class Result : public DataStructure, public ListOfTables {
  protected:
    using MapOfFieldOnNodesReal = std::map< ASTERINTEGER, FieldOnNodesRealPtr >;
    using MapOfFieldOnCellsReal = std::map< ASTERINTEGER, FieldOnCellsRealPtr >;
    using MapOfConstantFieldOnCellsReal = std::map< ASTERINTEGER, ConstantFieldOnCellsRealPtr >;
    using MapOfFieldOnNodesComplex = std::map< ASTERINTEGER, FieldOnNodesComplexPtr >;
    using MapOfFieldOnCellsComplex = std::map< ASTERINTEGER, FieldOnCellsComplexPtr >;
    using MapOfFieldOnCellsLong = std::map< ASTERINTEGER, FieldOnCellsLongPtr >;
    using MapOfConstantFieldOnCellsChar16 = std::map< ASTERINTEGER, ConstantFieldOnCellsChar16Ptr >;
    using MapOfGeneralizedVectorReal = std::map< ASTERINTEGER, GeneralizedAssemblyVectorRealPtr >;
    using MapOfGeneralizedVectorComplex =
        std::map< ASTERINTEGER, GeneralizedAssemblyVectorComplexPtr >;

    using mapStrMoFNR = std::map< std::string, MapOfFieldOnNodesReal >;
    using mapStrMoFCR = std::map< std::string, MapOfFieldOnCellsReal >;
    using mapStrMoCFCR = std::map< std::string, MapOfConstantFieldOnCellsReal >;
    using mapStrMoFNC = std::map< std::string, MapOfFieldOnNodesComplex >;
    using mapStrMoFCC = std::map< std::string, MapOfFieldOnCellsComplex >;
    using mapStrMoFCI = std::map< std::string, MapOfFieldOnCellsLong >;
    using mapStrMoCFCK16 = std::map< std::string, MapOfConstantFieldOnCellsChar16 >;
    using mapStrMoFGVR = std::map< std::string, MapOfGeneralizedVectorReal >;
    using mapStrMoFGVC = std::map< std::string, MapOfGeneralizedVectorComplex >;

    using mapIndexCaraElem = std::map< ASTERINTEGER, ElementaryCharacteristicsPtr >;
    using mapIndexLoads = std::map< ASTERINTEGER, ListOfLoadsPtr >;
    using mapIndexMaterial = std::map< ASTERINTEGER, MaterialFieldPtr >;
    using mapIndexModel = std::map< ASTERINTEGER, ModelPtr >;

    using mapParameters = std::map< std::string, std::string >;

    //
    // Objects for Fortran/JEVEUX
    //
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

    //
    // Objects for C++
    //
    /** @brief Fields */
    mapStrMoFNR _dictOfMapOfFieldOnNodesReal;
    mapStrMoFNC _dictOfMapOfFieldOnNodesComplex;
    mapStrMoFCR _dictOfMapOfFieldOnCellsReal;
    mapStrMoFCC _dictOfMapOfFieldOnCellsComplex;
    mapStrMoFCI _dictOfMapOfFieldOnCellsLong;
    mapStrMoCFCR _dictOfMapOfConstantFieldOnCellsReal;
    mapStrMoCFCK16 _dictOfMapOfConstantFieldOnCellsChar16;
    mapStrMoFGVR _dictOfMapOfGeneralizedVectorReal;
    mapStrMoFGVC _dictOfMapOfGeneralizedVectorComplex;

    /** @brief Parameters */
    mapParameters _dictParameters;

    /** @brief List of ElementaryCharacteristicsPtr */
    mapIndexCaraElem _mapElemCara;
    /** @brief List of ListOfLoadsPtr */
    mapIndexLoads _mapLoads;
    /** @brief List of MaterialFieldPtr */
    mapIndexMaterial _mapMaterial;
    /** @brief List of ModelPtr */
    mapIndexModel _mapModel;

    /** @brief Mesh */
    BaseMeshPtr _mesh;
    /** @brief Object to correctly manage fields and field descriptions */
    FieldBuilder _fieldBuilder;

  protected:
    /**
     * @brief Set field in datastructure
     * @param symbName Symbolic name of the field
     * @param storageIndex Index
     * @param field Field to save
     * @param dict Index
     */
    template < typename T >
    void
    _setFieldBase( const std::string &symbName, const ASTERINTEGER &storageIndex,
                   std::shared_ptr< T > field,
                   std::map< std::string, std::map< ASTERINTEGER, std::shared_ptr< T > > > &dict );

    /**
     * @brief Check consistency of mesh
     * @param BaseMeshPtr mesh
     */
    void _checkMesh( const BaseMeshPtr &mesh ) const;

    /**
     * @brief Get internal index
     * @param storageIndex Index to store field
     */
    ASTERINTEGER _getInternalIndex( const ASTERINTEGER &storageIndex ) const;

    /**
     * @brief Get a JEVEUX name for field
     * @param indexSymbName Index of symbolic name of the field
     * @param storageIndex Index to store field
     */
    std::string _generateFieldName( const ASTERINTEGER &indexSymbName,
                                    const ASTERINTEGER &storageIndex ) const;

    /**
     * @brief Prepare list of parameters
     */
    void _listOfParameters( void );

  public:
    using ResultPtr = std::shared_ptr< Result >;

    /**
     * @brief Constructeur
     */
    Result( const std::string &resuTyp ) : Result( ResultNaming::getNewResultName(), resuTyp ) {};

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
          _fieldBuilder( FieldBuilder() ) {};

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
     * @brief Set elementary characteristics to container
     * @param storageIndex
     */
    void setElementaryCharacteristics( const ElementaryCharacteristicsPtr &elemCara,
                                       ASTERINTEGER storageIndex, bool exists_ok );
    void setElementaryCharacteristics( const ElementaryCharacteristicsPtr &elemCara,
                                       ASTERINTEGER storageIndex ) {
        return setElementaryCharacteristics( elemCara, storageIndex, false );
    };

    /**
     * @brief Add a existing EquationNumbering in _fieldBuilder
     */
    void addEquationNumbering( const EquationNumberingPtr &fond ) {
        _fieldBuilder.addEquationNumbering( fond );
    };

    /**
     * @brief Set list of loads at storageIndex
     * @param ListOfLoadsPtr, storageIndex
     */
    void setListOfLoads( const ListOfLoadsPtr &load, ASTERINTEGER storageIndex );

    /**
     * @brief Add material definition
     * @param storageIndex
     */
    void setMaterialField( const MaterialFieldPtr &mate, ASTERINTEGER storageIndex,
                           bool exists_ok );
    void setMaterialField( const MaterialFieldPtr &mate, ASTERINTEGER storageIndex ) {
        return setMaterialField( mate, storageIndex, false );
    }

    /**
     * @brief Add model
     * @param storageIndex
     */
    void setModel( const ModelPtr &model, ASTERINTEGER storageIndex, bool exists_ok );
    void setModel( const ModelPtr &model, ASTERINTEGER storageIndex ) {
        return setModel( model, storageIndex, false );
    };

    /**
     * @brief Set model
     */
    void setMesh( const BaseMeshPtr &mesh );

    /**
     * @brief Add time for one storageIndex
     * @param storageIndex
     */
    void setTime( const ASTERDOUBLE &time, ASTERINTEGER storageIndex );

    /**
     * @brief Add storage index
     */
    void addStorageIndex( const ASTERINTEGER storageIndex );

    /**
     * @brief Add parameter value for one storageIndex
     */
    void setParameterValue( std::string paraName, ASTERDOUBLE paraValue,
                            ASTERINTEGER storageIndex );
    void setParameterValue( std::string paraName, std::string paraValue,
                            ASTERINTEGER storageIndex );

    /**
     * @brief Get internal index of parameter
     */
    ASTERINTEGER getParameterIndex( std::string paraName, std::string paraValue );

    /**
     * @brief Create index
     */
    ASTERINTEGER
    createIndexFromParameter( const std::string &paraName, const std::string &paraValue );

    ASTERDOUBLE getTime( ASTERINTEGER storageIndex ) const;

    /**
     * @brief Append a elementary characteristics on all index of Result
     * @param ElementaryCharacteristicsPtr
     */
    void setElementaryCharacteristics( const ElementaryCharacteristicsPtr &cara,
                                       bool exists_ok = false );

    /**
     * @brief Append a material on all index of Result
     * @param MaterialFieldPtr
     */
    void setMaterialField( const MaterialFieldPtr &mate, bool exists_ok = false );

    /**
     * @brief Append a model on all index of Result
     * @param ModelPtr
     */
    void setModel( const ModelPtr &model, bool exists_ok = false );

    /**
     * @brief Get list of loads at index
     * @param index
     */
    ListOfLoadsPtr getListOfLoads( ASTERINTEGER storageIndex ) const;

    bool hasListOfLoads( const ASTERINTEGER &storageIndex ) const;

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
    ElementaryCharacteristicsPtr getElementaryCharacteristics( ASTERINTEGER storageIndex ) const;

    /**
     * @brief Get elementary characteristics
     * @param index
     */
    bool hasElementaryCharacteristics( ASTERINTEGER storageIndex ) const;

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
     * @param storageIndex
     */
    bool hasMaterialField( const ASTERINTEGER &storageIndex ) const;

    /**
     * @brief Get material
     * @param storageIndex
     */
    MaterialFieldPtr getMaterialField( ASTERINTEGER storageIndex ) const;

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
    bool hasModel( const ASTERINTEGER &storageIndex ) const;

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
     * @param storageIndex
     */
    ModelPtr getModel( ASTERINTEGER storageIndex ) const;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param storageIndex numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    FieldOnCellsRealPtr getFieldOnCellsReal( const std::string name,
                                             const ASTERINTEGER storageIndex,
                                             const bool updatePtr = true ) const;

    FieldOnCellsComplexPtr getFieldOnCellsComplex( const std::string name,
                                                   const ASTERINTEGER storageIndex,
                                                   const bool updatePtr = true ) const;

    FieldOnCellsLongPtr getFieldOnCellsLong( const std::string name,
                                             const ASTERINTEGER storageIndex,
                                             const bool updatePtr = true ) const;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param storageIndex numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    ConstantFieldOnCellsChar16Ptr getConstantFieldOnCellsChar16( const std::string name,
                                                                 const ASTERINTEGER storageIndex,
                                                                 const bool updatePtr ) const;

    ConstantFieldOnCellsRealPtr getConstantFieldOnCellsReal( const std::string name,
                                                             const ASTERINTEGER storageIndex,
                                                             const bool updatePtr ) const;

    /**
     * @brief Ajouter un champ par éléments réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param storageIndex numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    void setField( const FieldOnCellsRealPtr field, const std::string &name,
                   const ASTERINTEGER storageIndex );

    void setField( const FieldOnCellsComplexPtr field, const std::string &name,
                   const ASTERINTEGER storageIndex );

    void setField( const FieldOnCellsLongPtr field, const std::string &name,
                   const ASTERINTEGER storageIndex );

    /**
     * @brief Ajouter un champ par éléments réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param storageIndex numéro d'ordre
     * @return FieldOnCellsRealPtr pointant vers le champ
     */
    void setField( const ConstantFieldOnCellsRealPtr field, const std::string &name,
                   const ASTERINTEGER storageIndex );

    void setField( const ConstantFieldOnCellsChar16Ptr field, const std::string &name,
                   const ASTERINTEGER storageIndex );

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
     * @brief Get the list of generalized vectors
     * @return std::vector< string >
     */
    VectorString getGeneralizedVectorRealNames() const;

    VectorString getGeneralizedVectorComplexNames() const;

    /**
     * @brief Obtenir un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param storageIndex numéro d'ordre
     * @return FieldOnNodesRealPtr pointant vers le champ
     */
    FieldOnNodesRealPtr getFieldOnNodesReal( const std::string name,
                                             const ASTERINTEGER storageIndex,
                                             const bool updatePtr = true ) const;

    FieldOnNodesComplexPtr getFieldOnNodesComplex( const std::string name,
                                                   const ASTERINTEGER storageIndex,
                                                   const bool updatePtr = true ) const;

    /**
     * @brief Ajouter un champ aux noeuds réel à partir de son nom et de son numéro d'ordre
     * @param name nom Aster du champ
     * @param storageIndex numéro d'ordre
     * @return FieldOnNodesRealPtr pointant vers le champ
     */
    void setField( const FieldOnNodesRealPtr field, const std::string &name,
                   const ASTERINTEGER storageIndex );

    void setField( const FieldOnNodesComplexPtr field, const std::string &name,
                   const ASTERINTEGER storageIndex );

    /**
     * @brief Interpolation
     */

    FieldOnNodesRealPtr
    interpolateFieldOnNodesReal( const std::string name, const ASTERDOUBLE value,
                                 const std::string para = "INST", const std::string left = "EXCLU",
                                 const std::string right = "EXCLU",
                                 const std::string crit = "RELATIF", const ASTERDOUBLE prec = 1.e-6,
                                 const bool updatePtr = true ) const;

    FieldOnCellsRealPtr
    interpolateFieldOnCellsReal( const std::string name, const ASTERDOUBLE value,
                                 const std::string para = "INST", const std::string left = "EXCLU",
                                 const std::string right = "EXCLU",
                                 const std::string crit = "RELATIF", const ASTERDOUBLE prec = 1.e-6,
                                 const bool updatePtr = true ) const;

    /**
     * @brief Impression de la sd au format MED
     * @param fileName Nom du fichier MED à imprimer
     * @return true
     * @todo revoir la gestion des mot-clés par défaut (ex : TOUT_ORDRE)
     * @todo revoir la gestion des unités logiques (notamment si fort.20 existe déjà)
     */
    virtual void printMedFile( const std::filesystem::path &fileName,
                               std::string medName = std::string(), bool local = true,
                               bool internalVar = true ) const;

    /**
     * @brief Get the number of steps stored in the Result
     * @return nbIndexes
     */
    ASTERINTEGER getNumberOfIndexes() const;

    /**
     * @brief Get last storageIndex stored in the Result
     * @return lastIndex
     */
    ASTERINTEGER getLastIndex() const;

    /**
     * @brief Get last time value
     * @return last time value
     */
    ASTERDOUBLE getLastTime() const;

    /**
     * @brief Get firs storageIndex stored in the Result
     * @return firstIndex
     */
    ASTERINTEGER getFirstIndex() const;

    /**
     * @brief Get the steps stored in the Result
     * @return nbIndexes
     */
    VectorLong getIndexes() const;

    /**
     * @brief Get the steps stored in the Result for a
     *        specific Field
     * @return List of Indexes
     */
    VectorLong getIndexesForFieldName( const std::string &name ) const;

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
    virtual bool
    build( const std::vector< FiniteElementDescriptorPtr > feds =
               std::vector< FiniteElementDescriptorPtr >(),
           const std::vector< EquationNumberingPtr > fnds = std::vector< EquationNumberingPtr >() );

    /**
     * @brief Update the  Result's size
     */
    void resize( ASTERINTEGER nbIndexes );

    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptors() const;

    std::vector< EquationNumberingPtr > getEquationNumberings() const;

    void clear( const ASTERINTEGER &storageIndex );

    void clear();

    bool exists() const;
};

using ResultPtr = std::shared_ptr< Result >;

#endif /* RESULTS_H_ */
