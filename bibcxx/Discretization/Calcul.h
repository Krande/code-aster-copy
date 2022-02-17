#ifndef CALCUL_H_
#define CALCUL_H_

/**
 * @file Calcul.h
 * @brief Header of Calcul class
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

#include "astercxx.h"

#include "DataFields/ElementaryTerm.h"
#include "DataStructures/DataStructure.h"
#include "Discretization/ElementaryCharacteristics.h"
#include "MemoryManager/JeveuxObject.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Modeling/XfemModel.h"

/**
 * @class Calcul
 * @brief This class allow computation of elementary terms (CALCUL.F90)
 */
class Calcul {
  private:
    typedef std::map< std::string, DataFieldPtr > listFields;
    typedef std::map< std::string, ElementaryTermRealPtr > listElemTerms;
    typedef std::map< std::string, bool > listExists;

    /** @brief Option to compute */
    std::string _option;

    /** @brief Flag to complete elementary terms on all CPU */
    bool _completeField;

    /** @brief Flag for stop if option doesn't exist for one element */
    bool _stopCompute;

    /** @brief Where to compute elementary terms */
    bool _computeOnModel;
    BaseMeshPtr _mesh;
    ModelPtr _model;
    FiniteElementDescriptorPtr _FEDesc;

    /** @brief Input fields */
    listFields _inputFields;
    listElemTerms _inputElemTerms;

    /** @brief Output fields */
    listFields _outputFields;
    listExists _outputFieldsExist;

    /** @brief Output elementary terms */
    listElemTerms _outputElemTerms;
    listExists _outputElemTermsExist;

  private:
    /** @brief Create C++ objects for output fields */
    void postCompute();

  public:
    /** @typedef CalculPtr */
    typedef boost::shared_ptr< Calcul > CalculPtr;

    /** @brief Constructor */
    Calcul( const std::string &option );

    /** @brief Destructor */
    ~Calcul(){};

    /** @brief Set option to compute */
    void setOption( const std::string &option ) { _option = option; };

    /** @brief Set flag for complete field on all CPUs */
    void setCompleteField( const bool flag ) { _completeField = flag; };

    /** @brief Set flag to stop if option doesn't exist for one element */
    void setStopCompute( const bool flag ) { _stopCompute = flag; };

    /** @brief Set computation on finite element descriptor */
    void setFiniteElementDescriptor( const FiniteElementDescriptorPtr &FEDesc );

    /** @brief Set computation on model */
    void setModel( const ModelPtr &model );

    /** @brief Add input field */
    void addInputField( const std::string parameterName, const DataFieldPtr field );

    /** @brief Add output field */
    void addOutputField( const std::string parameterName, const DataFieldPtr field );

    /** @brief Add input elementary term */
    void addInputElementaryTerm( const std::string parameterName, const DataStructurePtr field );

    /** @brief Add output elementary term */
    void addOutputElementaryTerm( const std::string parameterName,
                                  const ElementaryTermRealPtr field );

    /** @brief Clear all input fields */
    void clearInputFields() { _inputFields.clear(); };

    /** @brief Clear all input elementary term */
    void clearInputElemTerms() { _inputElemTerms.clear(); };

    /** @brief Clear all output fields */
    void clearOutputFields() {
        _outputFields.clear();
        _outputFieldsExist.clear();
    };

    /** @brief Clear all output elementary term */
    void clearOutputElemTerms() {
        _outputElemTerms.clear();
        _outputElemTermsExist.clear();
    };

    /** @brief Clear all outputs  */
    void clearOutputs() {
        clearOutputFields();
        clearOutputElemTerms();
    };

    /** @brief Clear all inputs  */
    void clearInputs() {
        clearInputFields();
        clearInputElemTerms();
    };

    /** @brief Is output if is elementary term */
    bool hasOutputElementaryTerm( const std::string parameterName ) const;

    /** @brief Get output if is elementary term */
    ElementaryTermRealPtr getOutputElementaryTerm( const std::string parameterName ) const;

    /** @brief Compute option */
    void compute();

    /** @brief Add input fields for elementary characteristics */
    void addElementaryCharacteristicsField( const ElementaryCharacteristicsPtr &elemChara );

    /** @brief Create and add input field for Fourier */
    void addFourierModeField( const ASTERINTEGER nh );

    /** @brief Create and add input field for current time */
    void addTimeField( const ASTERDOUBLE time );

    /** @brief Create and add input fields for XFEM */
    void addXFEMField( const XfemModelPtr &xfemModel );
};

/**  @typedef CalculPtr */
using CalculPtr = boost::shared_ptr< Calcul >;

#endif /* CALCUL_H_ */
