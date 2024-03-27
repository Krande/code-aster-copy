#ifndef SYNTAXSAVER_H_
#define SYNTAXSAVER_H_

/**
 * @file SyntaxSaver.h
 * @brief Fichier entete de la classe SyntaxSaver
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

#include "Supervis/CommandSyntax.h"
#include "Utilities/SyntaxDictionary.h"

/**
 * @class SyntaxSaver
 * @brief Object to save a syntax
 * @author Nicolas Sellenet
 */
class SyntaxSaver {
  private:
    ASTERINTEGER _op = 0;
    std::string _commandName = "";
    py::dict _keywords;

    void dictCopy( ListSyntaxMapContainer &, const py::dict & );

    void listCopy( ListSyntaxMapContainer &, const py::list & );

  public:
    SyntaxSaver( const std::string &commandName, const ASTERINTEGER &op, py::dict syntax );

    const std::string &commandName() const { return _commandName; };

    const py::dict &keywords() const { return _keywords; };

    const ASTERINTEGER &operatorNumber() const { return _op; };
};

/**
 * @typedef SyntaxSaverPtr
 * @brief Pointeur intelligent vers un SyntaxSaver
 */
typedef std::shared_ptr< SyntaxSaver > SyntaxSaverPtr;

#endif /* SYNTAXSAVER_H_ */
