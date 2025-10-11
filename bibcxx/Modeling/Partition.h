#ifndef PARTITION_H_
#define PARTITION_H_

/**
 * @file Partition.h
 * @brief Fichier entete de la classe Partition
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

#include "aster_fort_utils.h"

#include "DataStructures/DataStructure.h"

/**
 * @enum ModelSplitingMethod
 * @brief Types of partition for model
 */
enum ModelSplitingMethod { Centralized, SubDomain, GroupOfCellsSplit };
const int nbModelSplitingMethod = 3;

/**
 * @var ModelSplitingMethodNames
 * @brief Keyword for types of partition
 */
extern const char *const ModelSplitingMethodNames[nbModelSplitingMethod];

/**
 * @enum GraphPartitioner
 * @brief Graph partitioner
 */
enum GraphPartitioner { ScotchPartitioner, MetisPartitioner };
const int nbGraphPartitioner = 2;

/**
 * @var GraphPartitionerNames
 * @brief Keyword for graph partitioner
 */
extern const char *const GraphPartitionerNames[nbGraphPartitioner];

/**
 * @class Partition
 * @brief Datastructure for partition
 */

class Partition : public DataStructure {
    JeveuxVectorLong _prti;
    JeveuxVectorChar24 _prtk;
    JeveuxVectorLong _nupr;
    JeveuxVectorLong _fdim;
    JeveuxVectorLong _feta;
    JeveuxVectorChar8 _fref;

  public:
    Partition( const std::string );
    const std::string getMethod() const;
};

/**
 * @typedef Partition
 */
using PartitionPtr = std::shared_ptr< Partition >;

#endif /* PARTITION_H_ */
