#ifndef PRESTRESSINGCABLE_H_
#define PRESTRESSINGCABLE_H_

/**
 * @file PrestressingCable.h
 * @brief Fichier entete de la classe PrestressingCable
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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
#include "DataFields/Table.h"
#include "DataStructures/DataStructure.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "Loads/ListOfLinearRelations.h"
#include "Materials/MaterialField.h"
#include "Meshes/Mesh.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

#include <stdexcept>

/**
 * @class PrestressingCable
 * @brief Produit une sd identique a celle produite par DEFI_CABLE_BP
 * @author Nicolas Sellenet
 * @todo ajouter un test pour valider des qu'on aura les macros
 */
class PrestressingCable : public DataStructure {
  private:
    ModelPtr _model;
    ElementaryCharacteristicsPtr _cara;
    MaterialFieldPtr _mater;
    MeshPtr _mesh;
    /** @brief Carte '.CHME.SIGIN' */
    ConstantFieldOnCellsRealPtr _sigin;
    /** @brief Table 'CABLEBP' */
    TablePtr _cableBP;
    /** @brief Table 'CABLEGL' */
    TablePtr _cableGL;
    /** @brief Table '.LIRELA' */
    ListOfLinearRelationsRealPtr _lirela;
    /** @brief Booleen indiquant si la sd a deja ete construite */
    bool _isBuilt;

  public:
    /**
     * @typedef PrestressingCable
     * @brief Pointeur intelligent vers un PrestressingCable
     */
    typedef std::shared_ptr< PrestressingCable > PrestressingCablePtr;

    /**
     * @brief Constructeur
     */
    PrestressingCable( const ModelPtr &, const MaterialFieldPtr &,
                       const ElementaryCharacteristicsPtr & );

    PrestressingCable( const std::string name, const ModelPtr &, const MaterialFieldPtr &,
                       const ElementaryCharacteristicsPtr & );

    // Since no constructor allows to have null or empty objects,
    // it is not necessary to check if they exist.
    ModelPtr getModel() const { return _model; };

    MaterialFieldPtr getMaterialField() const { return _mater; };

    ElementaryCharacteristicsPtr getElementaryCharacteristics() const { return _cara; };

    /**
     * @brief Methode permettant de savoir si l'objet a ete construit
     * @return true si le modele a deja ete construit
     */
    bool isBuilt() { return _isBuilt; };

    bool build();
};

/**
 * @typedef PrestressingCable
 * @brief Pointeur intelligent vers un PrestressingCable
 */
typedef std::shared_ptr< PrestressingCable > PrestressingCablePtr;

#endif /* PRESTRESSINGCABLE_H_ */
