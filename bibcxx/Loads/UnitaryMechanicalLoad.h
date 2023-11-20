#ifndef UNITARYMECHANICALLOAD_H_
#define UNITARYMECHANICALLOAD_H_

/**
 * @file UnitaryMechanicalLoad.h
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

#include "astercxx.h"

#include "aster_fort_superv.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "Loads/MechanicalLoad.h"
#include "Loads/MechanicalLoadDescription.h"
#include "Loads/PhysicalQuantity.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/MeshEntities.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

/**
 * @enum LoadEnum
 * @brief Inventory of all mechanical loads available in Code_Aster
 */
enum LoadEnum {
    NodalForce,
    ForceOnEdge,
    ForceOnFace,
    LineicForce,
    InternalForce,
    ForceOnBeam,
    ForceOnShell,
    PressureOnPipe,
    ImposedDoF,
    DistributedPressure,
    NormalSpeedOnFace,
    WavePressureOnFace,
    THMFlux
};

/**
 * @class LoadTraits
 * @brief Traits class for a Load
 */
// This is the most general case (defined but intentionally not implemented)
// It will be specialized for each load listed in the inventory

template < LoadEnum Load >
struct LoadTraits;

/*************************************************************/
/*  Loads consisting of a force applied to some localization */
/*************************************************************/
/**
 * @def LoadTraits<NodalForce>
 * @brief Declare specialization for NodalForce
 */
template <>
struct LoadTraits< NodalForce > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = false;
    static bool const isAllowedOnGroupOfCells = false;
    static bool const isAllowedOnGroupOfNodes = true;
};

/**
 * @def LoadTraits<ForceOnFace>
 * @brief Declare specialization for ForceOnFace
 */
template <>
struct LoadTraits< ForceOnFace > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = false;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
};

/**
 * @def LoadTraits<ForceOnEdge>
 * @brief Declare specialization for ForceOnEdge
 */
template <>
struct LoadTraits< ForceOnEdge > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = false;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
};

/**
 * @def LoadTraits<LineicForce>
 * @brief Declare specialization for LineicForce
 */
template <>
struct LoadTraits< LineicForce > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = false;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
};

/**
 * @def LoadTraits<InternalForce>
 * @brief Declare specialization for InternalForce
 */
template <>
struct LoadTraits< InternalForce > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = false;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
};

/**
 * @def LoadTraits<ForceOnBeam>
 * @brief Declare specialization for ForceOnBeam
 */
template <>
struct LoadTraits< ForceOnBeam > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = true;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
    /** @todo mot clé supplémentaire TYPE_CHARGE=FORCE */
};

/**
 * @def LoadTraits<ForceOnShell>
 * @brief Declare specialization for ForceOnShell
 */
template <>
struct LoadTraits< ForceOnShell > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = true;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
};

/**
 * @def LoadTraits<PressureOnPipe>
 * @brief Declare specialization for PressureOnPipe
 */
template <>
struct LoadTraits< PressureOnPipe > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = true;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
};

/**
 * @def LoadTraits<ImposedDoF>
 * @brief Declare specialization for ImposedDoF
 */
template <>
struct LoadTraits< ImposedDoF > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = true;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = true;
};

/**
 * @def LoadTraits<DistributedPressure>
 * @brief Declare specialization for DistributedPressure
 */
template <>
struct LoadTraits< DistributedPressure > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = true;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
};

/**
 * @def LoadTraits<NormalSpeedOnFace>
 * @brief Declare specialization for NormalSpeedOnFace
 */
template <>
struct LoadTraits< NormalSpeedOnFace > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = false;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
};

/**
 * @def LoadTraits<WavePressureOnFace>
 * @brief Declare specialization for WavePressureOnFace
 */
template <>
struct LoadTraits< WavePressureOnFace > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = false;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
};

/**
 * @def LoadTraits<THMFlux>
 * @brief Declare specialization for THMFlux
 */
template <>
struct LoadTraits< THMFlux > {
    // Mot clé facteur pour AFFE_CHAR_MECA
    static const std::string factorKeyword;
    // Authorized MeshEntity
    static bool const isAllowedOnWholeMesh = true;
    static bool const isAllowedOnGroupOfCells = true;
    static bool const isAllowedOnGroupOfNodes = false;
};

/**
 * @class MechanicalLoad
 * @brief Define a mechanical load
 * @author Natacha Bereux
 */
template < class PhysicalQuantity, LoadEnum Load >
class UnitaryMechanicalLoadReal : public MechanicalLoadReal {
  public:
    /** @typedef Traits Define the Traits type */
    typedef LoadTraits< Load > Traits;
    /** @typedef PhysicalQuantity Define the underlying PhysicalQuantity */
    typedef PhysicalQuantity PhysicalQuantityType;

  private:
    /** @typedef Definition d'un pointeur intelligent sur une PhysicalQuantity */
    typedef std::shared_ptr< PhysicalQuantity > PhysicalQuantityPtr;

    /** @typedef PhysicalQuantity que l'on veut imposer*/
    PhysicalQuantityPtr _physicalQuantity;

    /** @typedef Definition d'un pointeur intelligent sur un VirtualMeshEntity */
    typedef std::shared_ptr< VirtualMeshEntity > MeshEntityPtr;

    /** @brief MeshEntity sur laquelle repose le "blocage" */
    MeshEntityPtr _meshEntity;

  public:
    /**
     * @typedef MechanicalLoadPtr
     * @brief Pointeur intelligent vers un MechanicalLoad
     */
    typedef std::shared_ptr< UnitaryMechanicalLoadReal< PhysicalQuantity, Load > >
        UnitaryMechanicalLoadRealPtr;

    /**
     * @brief Constructor
     */
    UnitaryMechanicalLoadReal( const ModelPtr &model ) : MechanicalLoadReal( model ) {};

    /**
     * @brief Constructor
     */
    UnitaryMechanicalLoadReal( const std::string name, const ModelPtr &model )
        : MechanicalLoadReal( name, model ) {};

    /**
     * @brief Set a physical quantity on a MeshEntity (group of nodes
     *        or group of cells)
     * @param physPtr shared pointer to a PhysicalQuantity
     * @param nameOfGroup name of the group of cells
     * @return bool success/failure index
     */
    bool setValue( PhysicalQuantityPtr physPtr, std::string nameOfGroup = "" ) {
        // Check that the pointer to the model is not empty
        if ( ( !_mecaLoadDesc->getModel() ) || _mecaLoadDesc->getModel()->isEmpty() )
            throw std::runtime_error( "Model is empty" );

        // Get the type of MeshEntity
        BaseMeshPtr currentMesh = _mecaLoadDesc->getMesh();
        // If the MeshEntity is not given, the quantity is set on the whole mesh
        if ( nameOfGroup.size() == 0 && Traits::isAllowedOnWholeMesh ) {
            _meshEntity = MeshEntityPtr( new AllMeshEntities() );
        }
        // nameOfGroup is the name of a group of cells and
        // LoadTraits authorizes to base the current load on such a group
        else if ( currentMesh->hasGroupOfCells( nameOfGroup ) && Traits::isAllowedOnGroupOfCells ) {
            _meshEntity = MeshEntityPtr( new GroupOfCells( nameOfGroup ) );
        }
        // nameOfGroup is the name of a group of nodes and LoadTraits authorizes
        // to base the current load on such a group
        else if ( currentMesh->hasGroupOfNodes( nameOfGroup ) && Traits::isAllowedOnGroupOfNodes ) {
            _meshEntity = MeshEntityPtr( new GroupOfNodes( nameOfGroup ) );
        } else
            throw std::runtime_error(
                nameOfGroup + " does not exist in the mesh " +
                "or it is not authorized as a localization of the current load " );

        // Copy the shared pointer of the Physical Quantity
        _physicalQuantity = physPtr;
        return true;
    };

    /**
     * @brief appel de op0007
     */
    bool build() {
        // std::cout << " build " << std::endl;
        CommandSyntax cmdSt( "AFFE_CHAR_MECA" );
        cmdSt.setResult( ResultNaming::getCurrentName(), "CHAR_MECA" );

        SyntaxMapContainer dict;
        if ( !_mecaLoadDesc->getModel() )
            throw std::runtime_error( "Model is undefined" );
        dict.container["MODELE"] = _mecaLoadDesc->getModel()->getName();
        ListSyntaxMapContainer listeLoad;
        SyntaxMapContainer dict2;
        /* On itere sur les composantes de la "PhysicalQuantity" */
        typename PhysicalQuantityType::MapOfCompAndVal comp_val = _physicalQuantity->getMap();
        for ( typename PhysicalQuantityType::MapIt curIter( comp_val.begin() );
              curIter != comp_val.end(); ++curIter ) {
            const auto &tmp = ComponentNames.find( curIter->first );
            dict2.container[tmp->second] = curIter->second;
        }
        /* Caractéristiques du MeshEntity */
        if ( _meshEntity->getType() == AllMeshEntitiesType ) {
            dict2.container["TOUT"] = "OUI";
        } else {
            if ( _meshEntity->getType() == GroupOfNodesType ) {
                dict2.container["GROUP_NO"] = _meshEntity->getName();
                //        std::cout << "GROUP_NO " <<  _meshEntity->getName() << std::endl;
            } else if ( _meshEntity->getType() == GroupOfCellsType )
                dict2.container["GROUP_MA"] = _meshEntity->getName();
        }
        listeLoad.push_back( dict2 );
        // mot-clé facteur
        std::string kw = Traits::factorKeyword;
        dict.container[kw] = listeLoad;
        cmdSt.define( dict );
        // cmdSt.debugPrint();
        ASTERINTEGER op = 7;
        CALL_EXECOP( &op );
        return update_tables();
    };
};

/**********************************************************
 *  Explicit instantiation of template classes
 **********************************************************/

/* Appliquer une force nodale sur une modélisation 3D */
/** @typedef NodalForceReal  */
template class UnitaryMechanicalLoadReal< ForceReal, NodalForce >;
typedef UnitaryMechanicalLoadReal< ForceReal, NodalForce > NodalForceReal;
typedef std::shared_ptr< NodalForceReal > NodalForceRealPtr;

/* Appliquer une force nodale sur des éléments de structure */
/** @typedef NodalStructuralForceReal  */
template class UnitaryMechanicalLoadReal< StructuralForceReal, NodalForce >;
typedef UnitaryMechanicalLoadReal< StructuralForceReal, NodalForce > NodalStructuralForceReal;
typedef std::shared_ptr< NodalStructuralForceReal > NodalStructuralForceRealPtr;

/** @typedef ForceOnFaceReal  */
template class UnitaryMechanicalLoadReal< ForceReal, ForceOnFace >;
typedef UnitaryMechanicalLoadReal< ForceReal, ForceOnFace > ForceOnFaceReal;
typedef std::shared_ptr< ForceOnFaceReal > ForceOnFaceRealPtr;

/* Appliquer une force sur une arête d'élément volumique */
/** @typedef ForceOnEdgeReal  */
template class UnitaryMechanicalLoadReal< ForceReal, ForceOnEdge >;
typedef UnitaryMechanicalLoadReal< ForceReal, ForceOnEdge > ForceOnEdgeReal;
typedef std::shared_ptr< ForceOnEdgeReal > ForceOnEdgeRealPtr;

/* Appliquer une force sur une arête d'élément de structure (coque/plaque) */
/** @typedef StructuralForceOnEdgeReal  */
template class UnitaryMechanicalLoadReal< StructuralForceReal, ForceOnEdge >;
typedef UnitaryMechanicalLoadReal< StructuralForceReal, ForceOnEdge > StructuralForceOnEdgeReal;
typedef std::shared_ptr< StructuralForceOnEdgeReal > StructuralForceOnEdgeRealPtr;

/** @typedef LineicForceReal  */
template class UnitaryMechanicalLoadReal< ForceReal, LineicForce >;
typedef UnitaryMechanicalLoadReal< ForceReal, LineicForce > LineicForceReal;
typedef std::shared_ptr< LineicForceReal > LineicForceRealPtr;

/** @typedef InternalForceReal  */
template class UnitaryMechanicalLoadReal< ForceReal, InternalForce >;
typedef UnitaryMechanicalLoadReal< ForceReal, InternalForce > InternalForceReal;
typedef std::shared_ptr< InternalForceReal > InternalForceRealPtr;

/* Appliquer une force (définie dans le repère global) à une poutre */
/** @typedef StructuralForceOnBeamReal  */
template class UnitaryMechanicalLoadReal< StructuralForceReal, ForceOnBeam >;
typedef UnitaryMechanicalLoadReal< StructuralForceReal, ForceOnBeam > StructuralForceOnBeamReal;
typedef std::shared_ptr< StructuralForceOnBeamReal > StructuralForceOnBeamRealPtr;

/* Appliquer une force (définie dans le repère local) à une poutre */
/** @typedef LocalForceOnBeamReal  */
template class UnitaryMechanicalLoadReal< LocalBeamForceReal, ForceOnBeam >;
typedef UnitaryMechanicalLoadReal< LocalBeamForceReal, ForceOnBeam > LocalForceOnBeamReal;
typedef std::shared_ptr< LocalForceOnBeamReal > LocalForceOnBeamRealPtr;

/* Appliquer une force (définie dans le repère global) à une coque/plaque */
/** @typedef StructuralForceOnShellReal  */
template class UnitaryMechanicalLoadReal< StructuralForceReal, ForceOnShell >;
typedef UnitaryMechanicalLoadReal< StructuralForceReal, ForceOnShell > StructuralForceOnShellReal;
typedef std::shared_ptr< StructuralForceOnShellReal > StructuralForceOnShellRealPtr;

/* Appliquer une force (définie dans le repère local) à une coque/plaque */
/** @typedef LocalForceOnShellReal  */
template class UnitaryMechanicalLoadReal< LocalShellForceReal, ForceOnShell >;
typedef UnitaryMechanicalLoadReal< LocalShellForceReal, ForceOnShell > LocalForceOnShellReal;
typedef std::shared_ptr< LocalForceOnShellReal > LocalForceOnShellRealPtr;

/* Appliquer une pression à une coque/plaque */
/** @typedef PressureOnShellReal  */
template class UnitaryMechanicalLoadReal< PressureReal, ForceOnShell >;
typedef UnitaryMechanicalLoadReal< PressureReal, ForceOnShell > PressureOnShellReal;
typedef std::shared_ptr< PressureOnShellReal > PressureOnShellRealPtr;

/* Appliquer une pression à un tuyau */
/** @typedef PressureOnPipeReal  */
template class UnitaryMechanicalLoadReal< PressureReal, PressureOnPipe >;
typedef UnitaryMechanicalLoadReal< PressureReal, PressureOnPipe > PressureOnPipeReal;
typedef std::shared_ptr< PressureOnPipeReal > PressureOnPipeRealPtr;

/* Imposer un déplacement sur des noeuds */
/** @typedef ImposedDisplacementReal  */
template class UnitaryMechanicalLoadReal< DisplacementReal, ImposedDoF >;
typedef UnitaryMechanicalLoadReal< DisplacementReal, ImposedDoF > ImposedDisplacementReal;
typedef std::shared_ptr< ImposedDisplacementReal > ImposedDisplacementRealPtr;

/* Imposer une pression sur des noeuds */
/** @typedef ImposedPressureReal  */
template class UnitaryMechanicalLoadReal< PressureReal, ImposedDoF >;
typedef UnitaryMechanicalLoadReal< PressureReal, ImposedDoF > ImposedPressureReal;
typedef std::shared_ptr< ImposedPressureReal > ImposedPressureRealPtr;

/** @typedef DistributedPressureReal  */
template class UnitaryMechanicalLoadReal< PressureReal, DistributedPressure >;
typedef UnitaryMechanicalLoadReal< PressureReal, DistributedPressure > DistributedPressureReal;
typedef std::shared_ptr< DistributedPressureReal > DistributedPressureRealPtr;

/** @typedef NormalSpeedOnFaceReal  */
template class UnitaryMechanicalLoadReal< NormalSpeedReal, NormalSpeedOnFace >;
typedef UnitaryMechanicalLoadReal< NormalSpeedReal, NormalSpeedOnFace > NormalSpeedOnFaceReal;
typedef std::shared_ptr< NormalSpeedOnFaceReal > NormalSpeedOnFaceRealPtr;

/** @typedef WavePressureOnFaceReal  */
template class UnitaryMechanicalLoadReal< PressureReal, WavePressureOnFace >;
typedef UnitaryMechanicalLoadReal< PressureReal, WavePressureOnFace > WavePressureOnFaceReal;
typedef std::shared_ptr< WavePressureOnFaceReal > WavePressureOnFaceRealPtr;

/** @typedef DistributedHeatFluxReal  */
template class UnitaryMechanicalLoadReal< HeatFluxReal, THMFlux >;
typedef UnitaryMechanicalLoadReal< HeatFluxReal, THMFlux > DistributedHeatFluxReal;
typedef std::shared_ptr< DistributedHeatFluxReal > DistributedHeatFluxRealPtr;

/** @typedef DistributedHydraulicFluxReal  */
template class UnitaryMechanicalLoadReal< HydraulicFluxReal, THMFlux >;
typedef UnitaryMechanicalLoadReal< HydraulicFluxReal, THMFlux > DistributedHydraulicFluxReal;
typedef std::shared_ptr< DistributedHydraulicFluxReal > DistributedHydraulicFluxRealPtr;

#endif /* UNITARYMECHANICALLOAD_H_ */
