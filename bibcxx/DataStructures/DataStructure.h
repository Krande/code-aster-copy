#ifndef DATASTRUCTURE_H_
#define DATASTRUCTURE_H_

/**
 * @file DataStructure.h
 * @brief Fichier entete de la classe DataStructure
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
/* person_in_charge: nicolas.sellenet at edf.fr */

#ifdef __cplusplus

#include "astercxx.h"

#include "aster_pybind.h"

#include "DataStructures/DataStructureNaming.h"
#include "MemoryManager/JeveuxAllowedTypes.h"
#include "MemoryManager/JeveuxVector.h"

/**
 * @class DataStructure
 * @brief Classe mere des classes representant des sd_aster
 * @author Nicolas Sellenet
 * @todo rendre cette classe virtuelle pure ?
 */
class DataStructure {
  public:
    /** @typedef shared_ptr d'une DataStructure */
    using DataStructurePtr = std::shared_ptr< DataStructure >;

  private:
    /** @brief Nom de la sd */
    std::string _name;
    /** @brief User variable name */
    std::string _user_name;
    /** @brief Object that stores the DataStructure type for jeveux requests */
    JeveuxVectorChar24 _tco;
    /** @brief Vector which contains reference to other DataStructure */
    std::vector< DataStructurePtr > _depsVector;
    /** @brief Python attributes for sdj/sdveri */
    py::object _sdj;
    /** @brief Python attributes to support values caching */
    py::object _cache;

  protected:
    /** @brief Object that stores the title of the DataStructure */
    JeveuxVectorChar80 _title;

  public:
    /**
     * @brief Constructeur
     * @param name Name of the jeveux datastructure
     * @param nameLength Length of the jeveux basename
     * @param type code_aster type of the datastructure

     */
    DataStructure( const std::string name, const int nameLength, const std::string type = "" );

    /**
     * @brief Constructeur
     * @param type code_aster type of the datastructure
     * @param nameLength Length of the jeveux basename
     */
    DataStructure( const int nameLength, const std::string type );

    /**
     * @brief Copy Constructor
     * @param other DataStructure to be copied
     */
    DataStructure( const DataStructure & );

    /**
     * @brief Move constructor
     */
    DataStructure( DataStructure && );

    /**
     * @brief Move assignment
     */
    DataStructure &operator=( DataStructure && );

    /**
     * @brief Destructeur
     */
    ~DataStructure();

    inline ASTERINTEGER id() { return (ASTERINTEGER)this; };

    /**
     * @brief Function to add a datastructure as a dependency
     * @param ds datastructure to reference
     */
    void addDependency( const DataStructurePtr & );

    /**
     * @brief Function to remove a datastructure from dependencies
     * @param ds datastructure to be removed
     */
    void removeDependency( const DataStructurePtr & );

    void resetDependencies();

    std::vector< DataStructurePtr > getDependencies() const;

    /**
     * @brief Function membre debugPrint
     * @param logicalUnit Unite logique d'impression
     * @param synchro Synchronisation des impressions en MPI
     */
    virtual void debugPrint( int logicalUnit = 6, bool synchro = true ) const;

    /**
     * @brief Function membre getName
     * @return une chaine contenant le nom de la sd
     */
    const std::string &getName() const { return _name; };

    /**
     * @brief Property userName
     * @return the user variable name
     */
    const std::string &getUserName() const { return _user_name; };

    void setUserName( const std::string );

    void setTitle( const std::string );

    std::string getTitle();

    /**
     * @brief Function membre getType
     * @return le type de la sd
     */
    const std::string getType() const {
        if ( !_tco.exists() )
            return "not_found";
        _tco->updateValuePointer();
        return ( *_tco )[0].rstrip();
    };

    /**
     * @brief Getter for SDJ property
     */
    const py::object &getSDJ() const { return _sdj; }

    /**
     * @brief Setter for SDJ property
     */
    void setSDJ( py::object &sdj ) { _sdj = sdj; }

    /**
     * @brief Getter for the cache property
     */
    const py::object &getCache() const { return _cache; }

    /**
     * @brief Setter for the cache property
     */
    void setCache( py::object &cache ) { _cache = cache; }

    /**
     * @brief Virtual function to update DataStructure
     */
    virtual bool build() { return true; };

  protected:
    /**
     * @brief Methode servant a fixer a posteriori le type d'une sd
     * @param newType chaine contenant le nouveau type
     */
    void setType( const std::string newType );
};
using DataStructurePtr = std::shared_ptr< DataStructure >;

class DSWithCppPickling : public DataStructure {
  public:
    using DSWithCppPicklingPtr = std::shared_ptr< DSWithCppPickling >;

    using DataStructure::DataStructure; // Inherit all constructors from DataStructure
};

using DSWithCppPicklingPtr = std::shared_ptr< DSWithCppPickling >;

#endif

#endif /* DATASTRUCTURE_H_ */
