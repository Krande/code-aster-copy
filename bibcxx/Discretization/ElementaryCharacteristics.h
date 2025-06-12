#ifndef ELEMENTARYCHARACTERISTICS_H_
#define ELEMENTARYCHARACTERISTICS_H_

/**
 * @file ElementaryCharacteristics.h
 * @brief Fichier entete de la classe ElementaryCharacteristics
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

#include "astercxx.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "DataFields/FieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

/**
 * @class ElementaryCharacteristics
 * @brief Class for elementary characteristics (from AFFE_CARA_ELEM)
 */
class ElementaryCharacteristics : public DataStructure {
  private:
    /** @brief Model */
    ModelPtr _model;
    /** @brief Mesh */
    BaseMeshPtr _mesh;
    /** @brief Objet Jeveux '.CARORIEN' for local basis */
    ConstantFieldOnCellsRealPtr _CARORIEN;
    /** @brief Objet Jeveux '.CARDISCK' for rigidity parameters for DIS_* elements */
    ConstantFieldOnCellsRealPtr _CARDISCK;
    /** @brief Objet Jeveux '.CARDISCM' for mass parameters for DIS_* elements */
    ConstantFieldOnCellsRealPtr _CARDISCM;
    /** @brief Objet Jeveux '.CARDISCA' for damping parameters for DIS_* elements */
    ConstantFieldOnCellsRealPtr _CARDISCA;
    /** @brief Objet Jeveux '.CARGENPO' for section properties for beam elements */
    ConstantFieldOnCellsRealPtr _CARGENPO;
    /** @brief Objet Jeveux '.CARGEOPO' for geometric properties for beam elements */
    ConstantFieldOnCellsRealPtr _CARGEOPO;
    /** @brief Objet Jeveux '.CARCOQUE' for properties of shell elements */
    ConstantFieldOnCellsRealPtr _CARCOQUE;
    /** @brief Objet Jeveux '.CARARCPO' for flexibility coefficients */
    ConstantFieldOnCellsRealPtr _CARARCPO;
    /** @brief Objet Jeveux '.CARCABLE' for properties of cable elements */
    ConstantFieldOnCellsRealPtr _CARCABLE;
    /** @brief Objet Jeveux '.CARGENBA' for properties of bar elements */
    ConstantFieldOnCellsRealPtr _CARGENBA;
    /** @brief Objet Jeveux '.CARMASSI' for orientation of material parameters */
    ConstantFieldOnCellsRealPtr _CARMASSI;
    /** @brief Objet Jeveux '.CARPOUFL' for properties of fluid beam elements */
    ConstantFieldOnCellsRealPtr _CARPOUFL;
    /** @brief Objet Jeveux '.CANBSP' for number of subpoints */
    FieldOnCellsLongPtr _CANBSP;
    /** @brief Objet Jeveux '.CAFIBR' for fibers */
    FieldOnCellsRealPtr _CAFIBR;
    /** @brief Objet Jeveux '.CARDINFO' for general parameters for DIS_* elements */
    ConstantFieldOnCellsRealPtr _CARDINFO;
    /** @brief V K8 '.MODELE' */
    JeveuxVectorChar8 _model_name;
    /** @brief Objet Jeveux '.CVENTCXF' */
    ConstantFieldOnCellsChar8Ptr _lineic;
    /** @brief Objet Jeveux '.CARDINFO' */
    ConstantFieldOnCellsRealPtr _infos;
    /** @brief Flag for empty datastructure */
    bool _isEmpty;

  public:
    /** @typedef ElementaryCharacteristicsPtr */
    typedef std::shared_ptr< ElementaryCharacteristics > ElementaryCharacteristicsPtr;

    /** @brief Constructor with a name */
    ElementaryCharacteristics( const std::string name, const ModelPtr &model );

    /** @brief Constructor with automatic name */
    ElementaryCharacteristics( const ModelPtr &model )
        : ElementaryCharacteristics( ResultNaming::getNewResultName(), model ) {};

    /** @brief Destructor */
    ~ElementaryCharacteristics() {};

    /** @brief Get the model */
    ModelPtr getModel() const;

    /** @brief Get the mesh */
    BaseMeshPtr getMesh() const;

    /**
     * @brief Detect state of datastructure
     * @return true if empty datastructure
     */
    bool isEmpty() const { return _isEmpty; };

    /** @brief Get fields */
    ConstantFieldOnCellsRealPtr getLocalBasis() { return _CARORIEN; };
    ConstantFieldOnCellsRealPtr getDiscreteRigidity() { return _CARDISCK; };
    ConstantFieldOnCellsRealPtr getDiscreteMass() { return _CARDISCM; };
    ConstantFieldOnCellsRealPtr getDiscreteDamping() { return _CARDISCA; };
    ConstantFieldOnCellsRealPtr getBeamGeometry() { return _CARGEOPO; };
    ConstantFieldOnCellsRealPtr getBeamSection() { return _CARGENPO; };
    ConstantFieldOnCellsRealPtr getShellParameters() { return _CARCOQUE; };
    ConstantFieldOnCellsRealPtr getFlexibilityCoefficients() { return _CARARCPO; };
    ConstantFieldOnCellsRealPtr getCableParameters() { return _CARCABLE; };
    ConstantFieldOnCellsRealPtr getBarParameters() { return _CARGENBA; };
    ConstantFieldOnCellsRealPtr getMaterialBase() { return _CARMASSI; };
    ConstantFieldOnCellsRealPtr getFluidBeamParameters() { return _CARPOUFL; };
    FieldOnCellsLongPtr getNumberOfSubpoints() { return _CANBSP; };
    FieldOnCellsRealPtr getFibers() { return _CAFIBR; };
    ConstantFieldOnCellsRealPtr getDiscreteParameters() { return _CARDINFO; };

    bool containsFieldOnCells() const {
        if ( _CANBSP->exists() || _CAFIBR->exists() ) {
            return true;
        }
        return false;
    }
};

/** @typedef ElementaryCharacteristicsPtr */
typedef std::shared_ptr< ElementaryCharacteristics > ElementaryCharacteristicsPtr;

#endif /* ELEMENTARYCHARACTERISTICS_H_ */
