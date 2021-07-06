#ifndef CONTACT_H_
#define CONTACT_H_

/**
 * @file Contact.h
 * @brief Fichier entete de la class Contact
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include "Utilities/CapyConvertibleValue.h"
#include <stdexcept>
#include <list>
#include <string>
#include "astercxx.h"
#include "Modeling/Model.h"
#include "MemoryManager/JeveuxVector.h"
#include "DataStructures/DataStructure.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/ResultNaming.h"

class Contact : public DataStructure {
  private:
    /** @brief La SD est-elle vide ? */
    bool _isEmpty;
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Ligel ".CHME.LIGRE" */
    FiniteElementDescriptorPtr _FEDesc;

  public:
    /**
    * @typedef ContactPt
    * @brief Pointeur intelligent vers un Contact
    */
    typedef boost::shared_ptr< Contact > ContactPtr;
    /**
     * @brief Constructeur
     */
    Contact() : Contact( ResultNaming::getNewResultName() ){};

    /**
     * @brief Constructeur
     */
    Contact( const std::string name )
        : DataStructure( name, 8, "CHAR_CONTACT" ), _model( ModelPtr() ),
            _isEmpty( true ){};

    /**
     * @brief Constructeur
     */
    Contact( const std::string name, const ModelPtr model )
        : DataStructure( name, 8, "CHAR_CONTACT" ), _model( model ),
            _FEDesc( boost::make_shared< FiniteElementDescriptor >( getName() + ".LIGRE",
                                                                _model->getMesh() ) ),
            _isEmpty( true ){};

    /**
     * @brief Constructeur
     */
    Contact( const ModelPtr model )
        : Contact( ResultNaming::getNewResultName(),  model ){};

    void setModel( const ModelPtr model)
    {
        _model = model;
    }

    ModelPtr getModel( ) const
    {
        return _model;
    }

    FiniteElementDescriptorPtr getFiniteElementDescriptor() const
    {
        return _FEDesc;
    }
};

/**
* @typedef ContactPt
* @brief Pointeur intelligent vers un Contact
*/
typedef boost::shared_ptr< Contact > ContactPtr;

#endif /* CONTACT_H_ */
