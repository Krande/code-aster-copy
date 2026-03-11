#ifndef RESULTMANAGER_H_
#define RESULTMANAGER_H_

/**
 * @file ResultManager.h
 * @brief Fichier entete de la classe ResultManager
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
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

#include "Results/Result.h"

/**
 * @class ResultManager
 * @brief Cette classe correspond a
 * @author Nicolas Sellenet
 */
class ResultManager {
  protected:
    using ResultManagerPtr = std::shared_ptr< ResultManager >;
    ResultPtr _currentResult;

    ResultManager() : _currentResult( nullptr ) {};

  public:
    ResultManager( ResultManager &other ) = delete;
    /**
     * Singletons should not be assignable.
     */
    void operator=( const ResultManager & ) = delete;

    void addFieldToResult( const std::string &nomSymb, const std::string &name,
                           ASTERINTEGER storageIndex );

    static ResultManagerPtr getInstance() {
        static ResultManagerPtr instance( new ResultManager() );
        return instance;
    };

    ResultPtr getCurrentResult() const { return _currentResult; };

    void releaseCurrentResult() { _currentResult = nullptr; };

    void setCurrentResult( const ResultPtr &result ) { _currentResult = result; };
};

using ResultManagerPtr = std::shared_ptr< ResultManager >;

#endif /* RESULTMANAGER_H_ */
