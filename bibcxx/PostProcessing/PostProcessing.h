
/**
 * @file PostProcessing.h
 * @brief Header of class PostProcessing
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

#pragma once

#include "astercxx.h"

#include "DataFields/FieldOnCells.h"
#include "DataFields/FieldOnNodes.h"
#include "Results/Result.h"
#include "Studies/PhysicalProblem.h"

/**
 * @class PostProcessing
 * @brief Post-processing tools
 */
class PostProcessing {
  private:
    /** @brief Physical problem */
    PhysicalProblemPtr _phys_problem;

  public:
    /** @typedef PostProcessingPtr */
    using PostProcessingPtr = std::shared_ptr< PostProcessing >;

    /** @brief Default constructor disabled */
    PostProcessing( void ) = delete;

    /**
     * @brief Constructor
     * @param PhysicalProblemPtr study
     */
    PostProcessing( const PhysicalProblemPtr &currPhysProblem )
        : _phys_problem( currPhysProblem ) {};

    /** @brief Destructor */
    ~PostProcessing() {};

    /**
     * @brief Compute hydration
     * @return field for hydration
     */
    FieldOnCellsRealPtr computeHydration( const FieldOnNodesRealPtr temp_prev,
                                          const FieldOnNodesRealPtr temp_curr,
                                          const ASTERDOUBLE time_prev, const ASTERDOUBLE time_curr,
                                          const FieldOnCellsRealPtr hydr_prev ) const;

    /**
     * @brief Compute annealing
     * @return Internal state variables (VARI_ELGA)
     */
    FieldOnCellsRealPtr
    computeAnnealing( const FieldOnCellsRealPtr internVar, const ASTERDOUBLE &time_prev,
                      const ASTERDOUBLE &time_curr,
                      const FieldOnCellsRealPtr &externVarPrev = nullptr,
                      const FieldOnCellsRealPtr &externVarCurr = nullptr ) const;

    /** @brief Compute max value of EFGE_ELNO or EGRU_ELNO based on the maximum of the
     *         equivalent moment for piping studies*/
    FieldOnCellsRealPtr computeMaxResultantForPipe( const ResultPtr &resu,
                                                    const std::string &field_name ) const;
};

using PostProcessingPtr = std::shared_ptr< PostProcessing >;
