#ifndef NonLinearResult_H_
#define NonLinearResult_H_

/**
 * @file NonLinearResult.h
 * @brief Fichier entete de la classe NonLinearResult
 * @author Natacha Béreux
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

/* person_in_charge: natacha.bereux at edf.fr */

#include "astercxx.h"

#include "Contact/Contact.h"
#include "Results/TransientResult.h"

#include <filesystem>

/**
 * @class NonLinearResult
 * @brief Cette classe correspond a un evol_noli, elle hérite de Result
          et stocke des champs
 * @author Natacha Béreux
 */
class NonLinearResult : public TransientResult {
  private:
    /** @typedef std::map du rang et des pointers vers ContactPtr */
    typedef std::map< ASTERINTEGER, ContactPtr > mapRankContact;

    /** @brief List of ContactPtr */
    mapRankContact _mapContact;

    static JeveuxVectorReal _mata;
    static JeveuxVectorReal _matc;

  public:
    /**
     * @brief Constructeur
     */
    NonLinearResult() : TransientResult( "EVOL_NOLI" ) {};

    /**
     * @brief Constructeur
     */
    NonLinearResult( const std::string name ) : TransientResult( name, "EVOL_NOLI" ) {};

    void setContact( const ContactPtr contact );

    void setContact( const ContactPtr contact, const ASTERINTEGER &rank );

    static VectorReal getTangentMatrix( const std::string & );

    void printMedFile( const std::filesystem::path &fileName, std::string medName = std::string(),
                       bool local = true, bool internalVar = true ) const;
};

/**
 * @typedef NonLinearResultPtr
 * @brief Pointeur intelligent vers un NonLinearResult
 */
typedef std::shared_ptr< NonLinearResult > NonLinearResultPtr;

#endif /* NonLinearResult_H_ */
