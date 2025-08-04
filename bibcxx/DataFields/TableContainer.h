#ifndef TABLECONTAINER_H_
#define TABLECONTAINER_H_

/**
 * @file TableContainer.h
 * @brief Fichier entete de la classe TableContainer
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

#include "DataFields/ConstantFieldOnCells.h"
#include "DataFields/DataField.h"
#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "DataFields/Table.h"
#include "Functions/Function.h"
#include "Functions/Function2D.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "LinearAlgebra/ElementaryVector.h"
#include "LinearAlgebra/GeneralizedAssemblyMatrix.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/Mesh.h"
#include "Results/ModeResult.h"
#include "Supervis/ResultNaming.h"

#include <map>
#include <string>

#ifdef ASTER_HAVE_MPI
#include "Meshes/ParallelMesh.h"
#endif

/**
 * @typedef TableContainer
 * @brief Definition of TableContainer (table_container)
 */
class TableContainer : public Table {
  private:
    JeveuxVectorChar16 _objectName;
    JeveuxVectorChar16 _objectType;
    JeveuxVectorChar24 _objectType24;
    JeveuxVectorChar8 _dsName;
    JeveuxVectorChar24 _dsName24;

    std::map< std::string, GeneralizedAssemblyMatrixRealPtr > _mapGAMD;
    std::map< std::string, ElementaryMatrixDisplacementRealPtr > _mapEMDD;
    std::map< std::string, ElementaryVectorDisplacementRealPtr > _mapEVDD;
    std::map< std::string, DataFieldPtr > _mapGDF;
    std::map< std::string, FieldOnNodesRealPtr > _mapFOND;
    std::map< std::string, ConstantFieldOnCellsRealPtr > _mapPCFOMD;
    std::map< std::string, FieldOnCellsRealPtr > _mapFOED;
    std::map< std::string, ModeResultPtr > _mapMMC;
    std::map< std::string, MeshPtr > _mapMesh;
    std::map< std::string, TablePtr > _mapT;
    std::map< std::string, FunctionPtr > _mapF;
    std::map< std::string, FunctionComplexPtr > _mapFC;
    std::map< std::string, Function2DPtr > _mapS;
#ifdef ASTER_HAVE_MPI
    std::map< std::string, ParallelMeshPtr > _mapPMesh;
#endif

  public:
    /**
     * @typedef TableContainerPtr
     * @brief Definition of a smart pointer to a TableContainer
     */
    typedef std::shared_ptr< TableContainer > TableContainerPtr;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux noeuds
     */
    TableContainer( const std::string &name ) : Table( name, "TABLE_CONTAINER" ) {};

    /**
     * @brief Constructeur
     */
    TableContainer() : TableContainer( ResultNaming::getNewResultName() ) {};

    /**
     * @brief Add ElementaryMatrixDisplacementReal to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &name, ElementaryMatrixDisplacementRealPtr );

    /**
     * @brief Add ElementaryVectorDisplacementReal to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, ElementaryVectorDisplacementRealPtr );

    /**
     * @brief Add FieldOnCellsReal to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, FieldOnCellsRealPtr );

    /**
     * @brief Add FieldOnNodesReal to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, FieldOnNodesRealPtr );

    /**
     * @brief Add Mesh to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, MeshPtr );

#ifdef ASTER_HAVE_MPI
    /**
     * @brief Add ParallelMesh to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, ParallelMeshPtr );
#endif

    /**
     * @brief Add Function to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, FunctionPtr );

    /**
     * @brief Add FunctionComplex to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, FunctionComplexPtr );

    /**
     * @brief Add generalized assembly matrix to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, GeneralizedAssemblyMatrixRealPtr );

    /**
     * @brief Add DataField to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, DataFieldPtr );

    /**
     * @brief Add ModeResult to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, ModeResultPtr );

    /**
     * @brief Add ConstantFieldOnCellsReal to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, ConstantFieldOnCellsRealPtr );

    /**
     * @brief Add Function2D to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, Function2DPtr );

    /**
     * @brief Add Table to TableContainer
     * @param name key used to find object
     */
    void addObject( const std::string &, TablePtr );

    /**
     * @brief Get ElementaryMatrixDisplacementReal stored in TableContainer
     * @param name key used to find object
     */
    ElementaryMatrixDisplacementRealPtr
    getElementaryMatrixDisplacementReal( const std::string &name ) const;

    /**
     * @brief Get ElementaryVectorDisplacementReal stored in TableContainer
     * @param name key used to find object
     */
    ElementaryVectorDisplacementRealPtr
    getElementaryVectorDisplacementReal( const std::string & ) const;

    /**
     * @brief Get FieldOnCellsReal stored in TableContainer
     * @param name key used to find object
     */
    FieldOnCellsRealPtr getFieldOnCellsReal( const std::string & ) const;

    /**
     * @brief Get FieldOnNodesReal stored in TableContainer
     * @param name key used to find object
     */
    FieldOnNodesRealPtr getFieldOnNodesReal( const std::string & ) const;

    /**
     * @brief Get mesh stored in TableContainer
     * @param name key used to find object
     */
    MeshPtr getMesh( const std::string & ) const;

#ifdef ASTER_HAVE_MPI
    /**
     * @brief Get mesh stored in TableContainer
     * @param name key used to find object
     */
    ParallelMeshPtr getParallelMesh( const std::string & ) const;
#endif

    /**
     * @brief Get Function stored in TableContainer
     * @param name key used to find object
     */
    FunctionPtr getFunction( const std::string & ) const;

    /**
     * @brief Get FunctionComplex stored in TableContainer
     * @param name key used to find object
     */
    FunctionComplexPtr getFunctionComplex( const std::string & ) const;

    /**
     * @brief Get generalized assembly matrix stored in TableContainer
     * @param name key used to find object
     */
    GeneralizedAssemblyMatrixRealPtr getGeneralizedAssemblyMatrix( const std::string & ) const;

    /**
     * @brief Get DataField stored in TableContainer
     * @param name key used to find object
     */
    DataFieldPtr getDataField( const std::string & ) const;

    /**
     * @brief Get ModeResult stored in TableContainer
     * @param name key used to find object
     */
    ModeResultPtr getModeResult( const std::string & ) const;

    /**
     * @brief Get ConstantFieldOnCellsReal stored in TableContainer
     * @param name key used to find object
     */
    ConstantFieldOnCellsRealPtr getConstantFieldOnCellsReal( const std::string & ) const;

    /**
     * @brief Get Function2D stored in TableContainer
     * @param name key used to find object
     */
    Function2DPtr getFunction2D( const std::string & ) const;

    /**
     * @brief Get Table stored in TableContainer
     * @param name key used to find object
     */
    TablePtr getTable( const std::string & ) const;

    /**
     * @brief Update the table
     * @todo add the case of ConstantFieldOnCells
     */
    bool build();
};

#endif /* TABLECONTAINER_H_ */
