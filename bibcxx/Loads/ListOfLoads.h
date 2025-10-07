#ifndef LISTOFLOADS_H_
#define LISTOFLOADS_H_

/**
 * @file ListOfLoads.h
 * @brief Fichier entete de la classe ListOfLoads
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

#include "Functions/Formula.h"
#include "Functions/Function.h"
#include "Functions/Function2D.h"
#include "Loads/AcousticLoad.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Loads/ParallelMechanicalLoad.h"
#include "Loads/ParallelThermalLoad.h"
#include "Loads/ThermalLoad.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"

typedef std::vector< GenericFunctionPtr > ListOfLoadFunctions;

/**
 * @class ListOfLoad
 * @brief Classe definissant une liste de charge
 * @author Nicolas Sellenet
 */
class ListOfLoads : public DataStructure {
  private:
    /** @brief Chargements cinematiques */
    ListDiriBC _listOfDirichletBCs;
    /** @brief List of functions for DirichletBCs */
    ListOfLoadFunctions _listOfDiriFun;
    /** @brief List of type_charge for  for DirichletBCs */
    VectorString _listOfDiriTyp;
    /** @brief Chargements Mecaniques */
    ListMecaLoadReal _listOfMechanicalLoadsReal;
    /** @brief List of functions for MechanicalLoads */
    ListOfLoadFunctions _listOfMechaFuncReal;
    /** @brief List of type_charge for MechanicalLoads */
    VectorString _listOfMechaTyp;
    /** @brief Chargements Mecaniques */
    ListMecaLoadFunction _listOfMechanicalLoadsFunction;
    /** @brief List of functions for MechanicalLoads */
    ListOfLoadFunctions _listOfMechaFuncFunction;
    /** @brief List of functions for MechanicalLoadsFunction */
    VectorString _listOfMechaFuncTyp;
    /** @brief Chargements Mecaniques */
    ListMecaLoadComplex _listOfMechanicalLoadsComplex;
    /** @brief List of functions for MechanicalLoads */
    ListOfLoadFunctions _listOfMechaFuncComplex;
    /** @brief Chargements Thermique */
    ListTherLoadReal _listOfThermalLoadsReal;
    /** @brief List of functions for ThermalLoads */
    ListOfLoadFunctions _listOfTherFuncReal;
    /** @brief Chargements Thermique */
    ListTherLoadFunction _listOfThermalLoadsFunction;
    /** @brief List of functions for ThermalLoads */
    ListOfLoadFunctions _listOfTherFuncFunction;
    /** @brief Chargements Acoustique */
    ListAcouLoadComplex _listOfAcousticLoadsComplex;
    /** @brief List of functions for AcousticLoads */
    ListOfLoadFunctions _listOfAcouFuncComplex;
#ifdef ASTER_HAVE_MPI
    /** @brief Chargements Mecaniques paralleles */
    ListParaMecaLoadReal _listOfParallelMechanicalLoadsReal;
    /** @brief List of functions for ParallelMechanicalLoads */
    ListOfLoadFunctions _listOfParaMechaFuncReal;
    /** @brief List of type_charge for ParallelMechanicalLoads */
    VectorString _listOfParaMechaTyp;
    /** @brief Chargements Mecaniques paralleles */
    ListParaMecaLoadFunction _listOfParallelMechanicalLoadsFunction;
    /** @brief List of functions for ParallelMechanicalLoads */
    ListOfLoadFunctions _listOfParaMechaFuncFunction;
    /** @brief List of type_charge for ParallelMechanicalLoads */
    VectorString _listOfParaMechaFuncTyp;
    /** @brief Chargements thermique paralleles */
    ListParaTherLoadReal _listOfParallelThermalLoadsReal;
    /** @brief List of functions for ParallelThermalLoads */
    ListOfLoadFunctions _listOfParaTherFuncReal;
    /** @brief Chargements thermique paralleles */
    ListParaTherLoadFunction _listOfParallelThermalLoadsFunction;
    /** @brief List of functions for ParallelThermalLoads */
    ListOfLoadFunctions _listOfParaTherFuncFunction;
#endif /* ASTER_HAVE_MPI */
    /** @brief .INFC */
    JeveuxVectorLong _loadInformations;
    /** @brief .LCHA */
    JeveuxVectorChar24 _list;
    /** @brief .FCHA */
    JeveuxVectorChar24 _listOfFunctions;
    /** @brief Le chargement a-t-il été construit ? */
    bool _isBuilt;
    /** @brief Model */
    ModelPtr _model;
    /** @brief differential displacement */
    FieldOnNodesRealPtr _differentialDisplacement;

  public:
    /**
     * @brief Constructeur
     */
    ListOfLoads();

    /**
     * @brief Constructeur
     */
    ListOfLoads( const std::string &name );

    /**
     * @brief Constructeur
     */
    ListOfLoads( const ModelPtr model );

    /**
     * @brief Constructeur
     */
    ListOfLoads( const std::string &name, const ModelPtr model );

    void addLoad( const DirichletBCPtr &currentLoad ) { addLoad( currentLoad, "FIXE_CSTE" ); };

    void addLoad( const DirichletBCPtr &currentLoad, const std::string &type ) {
        addLoad( currentLoad, emptyRealFunction, type );
    };

    /**
     * @brief Function d'ajout d'une charge cinematique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const DirichletBCPtr &currentLoad, const FunctionPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfDirichletBCs.push_back( currentLoad );
            _listOfDiriFun.push_back( func );
            _listOfDiriTyp.push_back( type );
        }
    };

    /**
     * @brief Function d'ajout d'une charge cinematique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const DirichletBCPtr &currentLoad, const FormulaPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfDirichletBCs.push_back( currentLoad );
            _listOfDiriFun.push_back( func );
            _listOfDiriTyp.push_back( type );
        }
    };

    /**
     * @brief Function d'ajout d'une charge cinematique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const DirichletBCPtr &currentLoad, const Function2DPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfDirichletBCs.push_back( currentLoad );
            _listOfDiriFun.push_back( func );
            _listOfDiriTyp.push_back( type );
        }
    };

    void addLoad( const MechanicalLoadRealPtr &currentLoad ) {
        addLoad( currentLoad, "FIXE_CSTE" );
    };

    void addLoad( const MechanicalLoadRealPtr &currentLoad, const std::string &type ) {
        addLoad( currentLoad, emptyRealFunction, type );
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const MechanicalLoadRealPtr &currentLoad, const FunctionPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfMechanicalLoadsReal.push_back( currentLoad );
            _listOfMechaFuncReal.push_back( func );
            _listOfMechaTyp.push_back( type );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const MechanicalLoadRealPtr &currentLoad, const FormulaPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfMechanicalLoadsReal.push_back( currentLoad );
            _listOfMechaFuncReal.push_back( func );
            _listOfMechaTyp.push_back( type );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const MechanicalLoadRealPtr &currentLoad, const Function2DPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfMechanicalLoadsReal.push_back( currentLoad );
            _listOfMechaFuncReal.push_back( func );
            _listOfMechaTyp.push_back( type );
        }
    };

    void addLoad( const MechanicalLoadComplexPtr &currentLoad ) {
        addLoad( currentLoad, emptyRealFunction );
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const MechanicalLoadComplexPtr &currentLoad, const FunctionPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfMechanicalLoadsComplex.push_back( currentLoad );
            _listOfMechaFuncComplex.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const MechanicalLoadComplexPtr &currentLoad, const FormulaPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfMechanicalLoadsComplex.push_back( currentLoad );
            _listOfMechaFuncComplex.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const MechanicalLoadComplexPtr &currentLoad, const Function2DPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfMechanicalLoadsComplex.push_back( currentLoad );
            _listOfMechaFuncComplex.push_back( func );
        }
    };

    void addLoad( const MechanicalLoadFunctionPtr &currentLoad ) {
        addLoad( currentLoad, "FIXE_CSTE" );
    };

    void addLoad( const MechanicalLoadFunctionPtr &currentLoad, const std::string &type ) {
        addLoad( currentLoad, emptyRealFunction, type );
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const MechanicalLoadFunctionPtr &currentLoad, const FunctionPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfMechanicalLoadsFunction.push_back( currentLoad );
            _listOfMechaFuncFunction.push_back( func );
            _listOfMechaFuncTyp.push_back( type );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const MechanicalLoadFunctionPtr &currentLoad, const FormulaPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfMechanicalLoadsFunction.push_back( currentLoad );
            _listOfMechaFuncFunction.push_back( func );
            _listOfMechaFuncTyp.push_back( type );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const MechanicalLoadFunctionPtr &currentLoad, const Function2DPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfMechanicalLoadsFunction.push_back( currentLoad );
            _listOfMechaFuncFunction.push_back( func );
            _listOfMechaFuncTyp.push_back( type );
        }
    };

#ifdef ASTER_HAVE_MPI
    void addLoad( const ParallelMechanicalLoadRealPtr &currentLoad ) {
        addLoad( currentLoad, "FIXE_CSTE" );
    };

    void addLoad( const ParallelMechanicalLoadRealPtr &currentLoad, const std::string &type ) {
        addLoad( currentLoad, emptyRealFunction, type );
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const ParallelMechanicalLoadRealPtr &currentLoad, const FunctionPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfParallelMechanicalLoadsReal.push_back( currentLoad );
            _listOfParaMechaFuncReal.push_back( func );
            _listOfParaMechaTyp.push_back( type );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const ParallelMechanicalLoadRealPtr &currentLoad, const FormulaPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfParallelMechanicalLoadsReal.push_back( currentLoad );
            _listOfParaMechaFuncReal.push_back( func );
            _listOfParaMechaTyp.push_back( type );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const ParallelMechanicalLoadRealPtr &currentLoad, const Function2DPtr &func,
                  const std::string &type ) {
        _isBuilt = false;
        this->setModel( currentLoad->getModel() );
        _listOfParallelMechanicalLoadsReal.push_back( currentLoad );
        _listOfParaMechaFuncReal.push_back( func );
        _listOfParaMechaTyp.push_back( type );
    };

    void addLoad( const ParallelMechanicalLoadFunctionPtr &currentLoad, const std::string &type ) {
        addLoad( currentLoad, emptyRealFunction, type );
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const ParallelMechanicalLoadFunctionPtr &currentLoad, const FunctionPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfParallelMechanicalLoadsFunction.push_back( currentLoad );
            _listOfParaMechaFuncFunction.push_back( func );
            _listOfParaMechaFuncTyp.push_back( type );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const ParallelMechanicalLoadFunctionPtr &currentLoad, const FormulaPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfParallelMechanicalLoadsFunction.push_back( currentLoad );
            _listOfParaMechaFuncFunction.push_back( func );
            _listOfParaMechaFuncTyp.push_back( type );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const ParallelMechanicalLoadFunctionPtr &currentLoad, const Function2DPtr &func,
                  const std::string &type ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfParallelMechanicalLoadsFunction.push_back( currentLoad );
            _listOfParaMechaFuncFunction.push_back( func );
            _listOfParaMechaFuncTyp.push_back( type );
        }
    };
#endif /* ASTER_HAVE_MPI */

    void addLoad( const ThermalLoadRealPtr &currentLoad ) {
        addLoad( currentLoad, emptyRealFunction );
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const ThermalLoadRealPtr &currentLoad, const FunctionPtr &func ) {
        _isBuilt = false;
        this->setModel( currentLoad->getModel() );
        _listOfThermalLoadsReal.push_back( currentLoad );
        _listOfTherFuncReal.push_back( func );
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const ThermalLoadRealPtr &currentLoad, const FormulaPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfThermalLoadsReal.push_back( currentLoad );
            _listOfTherFuncReal.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const ThermalLoadRealPtr &currentLoad, const Function2DPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfThermalLoadsReal.push_back( currentLoad );
            _listOfTherFuncReal.push_back( func );
        }
    };

    void addLoad( const ThermalLoadFunctionPtr &currentLoad ) {
        addLoad( currentLoad, emptyRealFunction );
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const ThermalLoadFunctionPtr &currentLoad, const FunctionPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfThermalLoadsFunction.push_back( currentLoad );
            _listOfTherFuncFunction.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const ThermalLoadFunctionPtr &currentLoad, const FormulaPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfThermalLoadsFunction.push_back( currentLoad );
            _listOfTherFuncFunction.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const ThermalLoadFunctionPtr &currentLoad, const Function2DPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfThermalLoadsFunction.push_back( currentLoad );
            _listOfTherFuncFunction.push_back( func );
        }
    };

#ifdef ASTER_HAVE_MPI

    void addLoad( const ParallelThermalLoadRealPtr &currentLoad ) {
        addLoad( currentLoad, emptyRealFunction );
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const ParallelThermalLoadRealPtr &currentLoad, const FunctionPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfParallelThermalLoadsReal.push_back( currentLoad );
            _listOfParaTherFuncReal.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const ParallelThermalLoadRealPtr &currentLoad, const FormulaPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfParallelThermalLoadsReal.push_back( currentLoad );
            _listOfParaTherFuncReal.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const ParallelThermalLoadRealPtr &currentLoad, const Function2DPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfParallelThermalLoadsReal.push_back( currentLoad );
            _listOfParaTherFuncReal.push_back( func );
        }
    };

    void addLoad( const ParallelThermalLoadFunctionPtr &currentLoad ) {
        addLoad( currentLoad, emptyRealFunction );
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const ParallelThermalLoadFunctionPtr &currentLoad, const FunctionPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfParallelThermalLoadsFunction.push_back( currentLoad );
            _listOfParaTherFuncFunction.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const ParallelThermalLoadFunctionPtr &currentLoad, const FormulaPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfParallelThermalLoadsFunction.push_back( currentLoad );
            _listOfParaTherFuncFunction.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge thermique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const ParallelThermalLoadFunctionPtr &currentLoad, const Function2DPtr &func ) {
        _isBuilt = false;
        this->setModel( currentLoad->getModel() );
        _listOfParallelThermalLoadsFunction.push_back( currentLoad );
        _listOfParaTherFuncFunction.push_back( func );
    };

#endif /* ASTER_HAVE_MPI */

    void addLoad( const AcousticLoadComplexPtr &currentLoad ) {
        addLoad( currentLoad, emptyRealFunction );
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function
     */
    void addLoad( const AcousticLoadComplexPtr &currentLoad, const FunctionPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfAcousticLoadsComplex.push_back( currentLoad );
            _listOfAcouFuncComplex.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier formula
     */
    void addLoad( const AcousticLoadComplexPtr &currentLoad, const FormulaPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfAcousticLoadsComplex.push_back( currentLoad );
            _listOfAcouFuncComplex.push_back( func );
        }
    };

    /**
     * @brief Function d'ajout d'une charge mécanique
     * @param currentLoad charge a ajouter a la sd
     * @param func multiplier function2d
     */
    void addLoad( const AcousticLoadComplexPtr &currentLoad, const Function2DPtr &func ) {
        if ( currentLoad ) {
            _isBuilt = false;
            this->setModel( currentLoad->getModel() );
            _listOfAcousticLoadsComplex.push_back( currentLoad );
            _listOfAcouFuncComplex.push_back( func );
        }
    };

    /**
     * @brief Construction de la liste de charge
     * @return Booleen indiquant que tout s'est bien passe
     */
    bool build( ModelPtr model = nullptr, std::string command_name = std::string() );

    /**
     * @brief Function de récupération des informations des charges
     * @return _loadInformations
     */
    const JeveuxVectorLong &getInformationVector() const { return _loadInformations; };

    /**
     * @brief Function de récupération de la liste des fonctions multiplicatrices
     * @return _listOfFunctions
     */
    const JeveuxVectorChar24 &getListOfFunctions() const { return _listOfFunctions; };

    /**
     * @brief Function de récupération de la liste des charges cinématiques
     * @return _listOfDirichletBCs
     */
    const ListDiriBC &getDirichletBCs() const { return _listOfDirichletBCs; };

    /**
     * @brief Function de récupération de la liste des charges mécaniques
     * @return _listOfMechanicalLoads
     */
    const ListMecaLoadReal &getMechanicalLoadsReal() const { return _listOfMechanicalLoadsReal; };

    /**
     * @brief Function de récupération de la liste des charges mécaniques
     * @return _listOfMechanicalLoads
     */
    const ListMecaLoadComplex &getMechanicalLoadsComplex() const {
        return _listOfMechanicalLoadsComplex;
    };

    /**
     * @brief Function de récupération de la liste des charges mécaniques
     * @return _listOfMechanicalLoads
     */
    const ListMecaLoadFunction &getMechanicalLoadsFunction() const {
        return _listOfMechanicalLoadsFunction;
    };

#ifdef ASTER_HAVE_MPI
    /**
     * @brief Function de récupération de la liste des charges mécaniques
     * @return _listOfMechanicalLoads
     */
    const ListParaMecaLoadReal &getParallelMechanicalLoadsReal() const {
        return _listOfParallelMechanicalLoadsReal;
    };

    /**
     * @brief Function de récupération de la liste des charges mécaniques
     * @return _listOfMechanicalLoads
     */
    const ListParaMecaLoadFunction &getParallelMechanicalLoadsFunction() const {
        return _listOfParallelMechanicalLoadsFunction;
    };
#endif /* ASTER_HAVE_MPI */

    /**
     * @brief Function de récupération de la liste des charges thermiques
     * @return _listOfThermalLoadsReal
     */
    const ListTherLoadReal &getThermalLoadsReal() const { return _listOfThermalLoadsReal; };

    /**
     * @brief Function de récupération de la liste des charges thermiques
     * @return _listOfThermalLoadsFunction
     */
    const ListTherLoadFunction &getThermalLoadsFunction() const {
        return _listOfThermalLoadsFunction;
    };

#ifdef ASTER_HAVE_MPI
    /**
     * @brief Function de récupération de la liste des charges thermiques
     * @return _listOfThermalLoadsReal
     */
    const ListParaTherLoadReal &getParallelThermalLoadsReal() const {
        return _listOfParallelThermalLoadsReal;
    };

    /**
     * @brief Function de récupération de la liste des charges thermiques
     * @return _listOfThermalLoadsFunction
     */
    const ListParaTherLoadFunction &getParallelThermalLoadsFunction() const {
        return _listOfParallelThermalLoadsFunction;
    };
#endif /* ASTER_HAVE_MPI */

    /**
     * @brief Function de récupération de la liste des charges thermiques
     * @return _listOfThermalLoadsReal
     */
    const ListAcouLoadComplex &getAcousticLoadsComplex() const {
        return _listOfAcousticLoadsComplex;
    };

    /**
     * @brief Function de récupération de la liste des charges
     * @return _list
     */
    JeveuxVectorChar24 getLoadNames() const { return _list; };

    /**
     * @brief Get all FiniteElementDescriptors
     * @return vector of all FiniteElementDescriptors
     */
    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptors() const;

    int getPhysics( void ) const;

    /**
     * @brief Methode permettant de savoir si le chargement a été construit
     * @return true si construit
     */
    bool isBuilt() { return _isBuilt; };

    /**
     * @brief Nombre de charges
     * @return taille de _listOfMechanicalLoads + taille de _listOfDirichletBCs
     */
    int getNumberOfLoads() const {
        return _listOfMechanicalLoadsReal.size() + _listOfMechanicalLoadsFunction.size() +
               _listOfMechanicalLoadsComplex.size() +
#ifdef ASTER_HAVE_MPI
               _listOfParallelMechanicalLoadsReal.size() +
               _listOfParallelMechanicalLoadsFunction.size() +
               _listOfParallelThermalLoadsReal.size() + _listOfParallelThermalLoadsFunction.size() +
#endif /* ASTER_HAVE_MPI */
               _listOfThermalLoadsReal.size() + _listOfThermalLoadsFunction.size() +
               _listOfAcousticLoadsComplex.size() + _listOfDirichletBCs.size();
    };

    bool hasDirichletBC() const { return _listOfDirichletBCs.size() > 0; }

    bool hasExternalLoad() const {
        return ( _listOfMechanicalLoadsReal.size() + _listOfMechanicalLoadsFunction.size() +
                 _listOfMechanicalLoadsComplex.size() +
#ifdef ASTER_HAVE_MPI
                 _listOfParallelMechanicalLoadsReal.size() +
                 _listOfParallelMechanicalLoadsFunction.size() +
                 _listOfParallelThermalLoadsReal.size() +
                 _listOfParallelThermalLoadsFunction.size() +
#endif /* ASTER_HAVE_MPI */
                 _listOfThermalLoadsReal.size() + _listOfThermalLoadsFunction.size() +
                 _listOfAcousticLoadsComplex.size() ) > 0;
    };

    /**
     * @brief Check that all loads have the same model
     * @return True if all loads have the same model
     */
    bool checkModelConsistency( const ModelPtr &model ) const;

    bool setModel( const ModelPtr &model );

    ModelPtr getModel( void ) const { return _model; };

    BaseMeshPtr getMesh() const {
        if ( _model ) {
            return _model->getMesh();
        }

        return nullptr;
    }

    VectorString getListOfMechaTyp();
    VectorString getListOfMechaFuncTyp();
    VectorString getListOfDiriTyp();

    void setDifferentialDisplacement( const FieldOnNodesRealPtr );

    const FieldOnNodesRealPtr getDifferentialDisplacement();

    bool hasDifferential();
    bool hasDifferentialLoads();
    bool hasDifferentialDirichletBC();
    bool hasDifferentialDisplacement();
};

/**
 * @typedef ListOfLoad
 * @brief Pointeur intelligent vers un ListOfLoad
 */
typedef std::shared_ptr< ListOfLoads > ListOfLoadsPtr;

#endif /* LISTOFLOADS_H_ */
