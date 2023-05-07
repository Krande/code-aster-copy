#ifndef TIMESTEPPER_H_
#define TIMESTEPPER_H_

/**
 * @file TimesList.h
 * @brief Fichier entete de la classe TimesList
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

/* person_in_charge: nicolas.sellenet at edf.fr */
#include "astercxx.h"

#include "aster_pybind.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Steppers/GenericStepper.h"

typedef VectorReal::const_iterator VectorRealCIter;

/**
 * @class TimesList
 * @brief Cette classe permet de definir une liste d'instants
 * @author Nicolas Sellenet
 */
class TimesList : public DataStructure, public GenericStepper {
  private:
    /** @brief Liste des instants '.LIST'*/
    JeveuxVectorReal _values;
    /** @brief '.LIST.DITR' */
    JeveuxVectorReal _ditr;
    /** @brief '.LIST.INFOR' */
    JeveuxVectorReal _infor;
    /** @brief '.ECHE.EVENR' */
    JeveuxVectorReal _evenr;
    /** @brief '.ECHE.EVENK' */
    JeveuxVectorChar16 _evenk;
    /** @brief '.ECHE.SUBDR' */
    JeveuxVectorReal _subdr;
    /** @brief '.ADAP.EVENR' */
    JeveuxVectorReal _aevenr;
    /** @brief '.ADAP.TPLUK' */
    JeveuxVectorChar16 _tpluk;
    /** @brief '.ADAP.TPLUR' */
    JeveuxVectorReal _tplur;

  public:
    /**
     * @typedef TimesListPtr
     * @brief Pointeur intelligent vers un TimesList
     */
    typedef std::shared_ptr< TimesList > TimesListPtr;

    /* stored here before removing this object */
    py::object pyStepper;

    /**
     * @brief Constructeur
     */
    TimesList( const std::string name )
        : DataStructure( name, 8, "LIST_INST" ),
          _values( getName() + ".LIST.DITR" ),
          _infor( getName() + ".LIST.INFOR" ),
          _evenr( getName() + ".ECHE.EVENR" ),
          _evenk( getName() + ".ECHE.EVENK" ),
          _subdr( getName() + ".ECHE.SUBDR" ),
          _aevenr( getName() + ".ADAP.EVENR" ),
          _tpluk( getName() + ".ADAP.TPLUK" ),
          _tplur( getName() + ".ADAP.TPLUR" ) {};

    /**
     * @brief Constructeur
     */
    TimesList() : TimesList( DataStructureNaming::getNewName( 8 ) ) {};

    /**
     * @brief Destructeur
     */
    ~TimesList() {};

    struct const_iterator {
        ASTERDOUBLE *position;
        int rank;

        inline const_iterator() : position( NULL ), rank( 1 ) {};

        inline const_iterator( int curRank, ASTERDOUBLE *memoryPosition )
            : position( memoryPosition ), rank( curRank ) {};

        inline const_iterator( const const_iterator &iter )
            : position( iter.position ), rank( iter.rank ) {};

        inline const_iterator &operator=( const const_iterator &testIter ) {
            position = testIter.position;
            rank = testIter.rank;
            return *this;
        };

        inline const_iterator &operator++() {
            ++position;
            ++rank;
            return *this;
        };

        inline bool operator==( const const_iterator &testIter ) const {
            if ( testIter.position != position )
                return false;
            return true;
        };

        inline bool operator!=( const const_iterator &testIter ) const {
            if ( testIter.position != position )
                return true;
            return false;
        };

        inline const ASTERDOUBLE &operator->() const { return *position; };

        inline const ASTERDOUBLE &operator*() const { return *position; };
    };

    /**
     * @brief
     * @return
     */
    const_iterator begin() const { return const_iterator( 1, &( *_values )[0] ); };

    /**
     * @brief
     * @return
     */
    const_iterator end() const {
        return const_iterator( _values->size() + 1,
                               (double *)( _values->getDataPtr() + _values->size() ) );
    };

    /**
     * @brief Fonction permettant de mettre a jour le stepper
     * @return true si tout s'est bien passé
     */
    bool operator=( const VectorReal &vecReal ) { return setValues( vecReal ); };

    /**
     * @brief Fonction permettant de fixer la liste de pas de temps
     * @param values Liste des valeurs
     */
    bool setValues( const VectorReal &values );

    /**
     * @brief Fonction permettant de connaître le nombre de pas de temps
     * @return nombre de pas de temps
     */
    ASTERINTEGER size() const { return _values->size(); };

    VectorReal getValues() { return _values->toVector(); };

    /**
     * @brief Fonction permettant de mettre a jour le stepper
     * @return true si tout s'est bien passé
     */
    bool build() const {
        _values->updateValuePointer();
        return true;
    };
};

/**
 * @typedef TimesListPtr
 * @brief Pointeur intelligent vers un TimesList
 */
typedef std::shared_ptr< TimesList > TimesListPtr;

#endif /* TIMESTEPPER_H_ */
