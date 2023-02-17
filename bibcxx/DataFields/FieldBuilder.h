#ifndef FIELDBUILDER_H_
#define FIELDBUILDER_H_

/**
 * @file FieldBuilder.h
 * @brief Header of class FieldBuilder
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

#include "astercxx.h"

#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/GeneralizedGlobalEquationNumbering.h"

/**
 * @class FieldBuilder
 * @brief This class builds FieldOnNodes and FieldOnCells with respect of
 *        GlobalEquationNumbering and FiniteElementDescriptor
 * @author Nicolas Sellenet
 */
class FieldBuilder {
private:
  std::map<std::string, GlobalEquationNumberingPtr> _mapGlobNume;
  std::map<std::string, FiniteElementDescriptorPtr> _mapLigrel;

  // I use them to debug easily mutiple creation
  // I don't use map directly to avoid to keep in memory unnecessary objects
  static std::set<std::string> _setGlobNume;
  static std::set<std::string> _setLigrel;

  /**
   * @brief Add a existing FiniteElementDescriptor in FieldBuilder
   */
  FiniteElementDescriptorPtr
  newFiniteElementDescriptor(const std::string &name, const BaseMeshPtr mesh) {
    if (_setLigrel.count(trim(name)) > 0) {
      raiseAsterError("LIGREL already exists: " + name);
    }

    auto curDesc = std::make_shared<FiniteElementDescriptor>(name, mesh);

    addFiniteElementDescriptor(curDesc);

    return curDesc;
  };

  /**
   * @brief Add a existing GlobalEquationNumbering in FieldBuilder
   */
  GlobalEquationNumberingPtr
  newGlobalEquationNumbering(const std::string &name) {
    if (_setGlobNume.count(trim(name)) > 0) {
      raiseAsterError("NUME_EQUA already exists: " + name);
    }

    auto curDesc = std::make_shared<GlobalEquationNumbering>(name);
    addGlobalEquationNumbering(curDesc);

    return curDesc;
  };

  /**
   * @brief Add a existing generalizedGlobalEquationNumbering in FieldBuilder
   */
  GeneralizedGlobalEquationNumberingPtr
  newGeneralizedGlobalEquationNumbering(const std::string &name) {
    AS_ABORT(name);
    if (_setGlobNume.count(trim(name)) > 0) {
      raiseAsterError("PROF_GENE already exists: " + name);
    }

    auto curDesc = std::make_shared<GeneralizedGlobalEquationNumbering>(name);
    // addGlobalEquationNumbering(curDesc);

    return curDesc;
  };

public:
  /**
   * @brief Constructor
   */
  FieldBuilder(){};

  /**
   * @brief Add a existing GlobalEquationNumbering in FieldBuilder
   */
  void addGlobalEquationNumbering(const GlobalEquationNumberingPtr &fond) {
    AS_ASSERT(fond);

    _mapGlobNume[trim(fond->getName())] = fond;
    _setGlobNume.insert(trim(fond->getName()));
  };

  /**
   * @brief Add a existing FiniteElementDescriptor in FieldBuilder
   */
  void addFiniteElementDescriptor(const FiniteElementDescriptorPtr &fed) {
    AS_ASSERT(fed);

    _mapLigrel[trim(fed->getName())] = fed;
    _setLigrel.insert(trim(fed->getName()));
  };

  void clear() {
    _mapGlobNume.clear();
    _mapLigrel.clear();
  };

  /**
   * @brief Build a FieldOnCells with a FiniteElementDescriptor
   */
  template <typename ValueType>
  std::shared_ptr<FieldOnCells<ValueType>>
  buildFieldOnCells(const std::string &name, const BaseMeshPtr mesh) {
    std::shared_ptr<FieldOnCells<ValueType>> field =
        std::make_shared<FieldOnCells<ValueType>>(name);
    field->updateValuePointers();

    const std::string ligrel = trim((*(*field)._reference)[0].toString());

    if (!ligrel.empty()) {
      auto curIter = _mapLigrel.find(ligrel);
      FiniteElementDescriptorPtr curDesc;
      if (curIter != _mapLigrel.end()) {
        curDesc = curIter->second;
      } else {
        curDesc = newFiniteElementDescriptor(ligrel, mesh);
      }

      field->setDescription(curDesc);
    }
    return field;
  };

  /**
   * @brief Build a ConstantFieldOnCells with a FiniteElementDescriptor
   */
  template <typename ValueType>
  std::shared_ptr<ConstantFieldOnCells<ValueType>>
  buildConstantFieldOnCells(const std::string &name, const BaseMeshPtr mesh) {

    std::shared_ptr<ConstantFieldOnCells<ValueType>> field =
        std::make_shared<ConstantFieldOnCells<ValueType>>(name, mesh);
    field->updateValuePointers();

    return field;
  };

  /**
   * @brief Build a FieldOnNodes with a GlobalEquationNumbering
   */
  template <typename ValueType>
  std::shared_ptr<FieldOnNodes<ValueType>> buildFieldOnNodes(std::string name) {
    std::shared_ptr<FieldOnNodes<ValueType>> field =
        std::make_shared<FieldOnNodes<ValueType>>(name);
    field->updateValuePointers();

    const std::string globNume = trim((*(*field)._reference)[1].toString());
    AS_ASSERT(!globNume.empty());

    auto curIter = _mapGlobNume.find(globNume);
    GlobalEquationNumberingPtr curDesc;
    if (curIter != _mapGlobNume.end())
      curDesc = curIter->second;
    else {
      // .REFE de taille 2 pour les VGEN, voir vpstor.F90
      if (field->_reference->size() == 2) {
        auto curDesc2 = newGeneralizedGlobalEquationNumbering(globNume);
      } else {
        curDesc = newGlobalEquationNumbering(globNume);
      }
    }
    field->setDescription(curDesc);

    return field;
  };

  std::vector<FiniteElementDescriptorPtr> getFiniteElementDescriptors() const {
    std::vector<FiniteElementDescriptorPtr> ret;

    for (auto &[name, fed] : _mapLigrel)
      ret.push_back(fed);

    return ret;
  };

  std::vector<GlobalEquationNumberingPtr> getGlobalEquationNumberings() const {
    std::vector<GlobalEquationNumberingPtr> ret;

    for (auto &[name, fnd] : _mapGlobNume)
      ret.push_back(fnd);

    return ret;
  };
};

#endif /* FIELDBUILDER_H_ */
