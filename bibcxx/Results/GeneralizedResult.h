#ifndef GENERALIZEDRESULTSCONTAINER_H_
#define GENERALIZEDRESULTSCONTAINER_H_

/**
 * @file GeneralizedResult.h
 * @brief Fichier entete de la classe GeneralizedResult
 * @author Natacha Béreux
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

/* person_in_charge: natacha.bereux at edf.fr */

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modal/StaticMacroElement.h"
#include "Results/DynamicResultsIndexing.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/GeneralizedDOFNumbering.h"
#include "Supervis/ResultNaming.h"

/**
 * @class GeneralizedResult
 * @brief Cette classe correspond a la sd_dyna_gene de Code_Aster.
 * Un objet sd_dyna_gene est un concept produit par un opérateur
 * dynamique sur base généralisée.
 * @author Natacha Béreux
 */
template <class ValueType>
class GeneralizedResult: public DataStructure
{
private:
    /** @brief DynamicResultsIndexing */
    DynamicResultsIndexingPtr _index;
    /** @brief Vecteur Jeveux '.DESC' */
    JeveuxVectorLong          _desc;
    /** @brief Vecteur Jeveux '.DISC' */
    /* Valeur des instants/fréquences sauvegardées */
    JeveuxVectorReal        _abscissasOfSamples;
    /** @brief Vecteur Jeveux '.ORDR' */
    JeveuxVectorLong          _indicesOfSamples;
    /** @brief Vecteur Jeveux '.DEPL' */
    JeveuxVector<ValueType>   _displacement;
    /** @brief Vecteur Jeveux '.VITE' */
    JeveuxVector<ValueType>   _velocity;
    /** @brief Vecteur Jeveux '.ACCE' */
    JeveuxVector<ValueType>   _acceleration;
    /** @brief si résulte d'un proj_mesu_modal */
    ProjMesuPtr               _projM;
    /** @brief Generalized DOFNumbering */
    GeneralizedDOFNumberingPtr _genDOFNum;
    /** @brief DOFNumbering */
    DOFNumberingPtr _DOFNum;
public:
    /**
     * @brief Constructeur
     */
    GeneralizedResult( const std::string &name, const std::string &resuTyp ):
        DataStructure( name, 19, resuTyp ),
        _index( new DynamicResultsIndexing( getName(), resuTyp )),
        _desc( JeveuxVectorLong( getName() + ".DESC" ) ),
        _abscissasOfSamples( JeveuxVectorReal( getName() +".DISC"  ) ),
        _indicesOfSamples( JeveuxVectorLong ( getName() +".ORDR"  ) ),
        _displacement( JeveuxVector<ValueType>( getName() +".DEPL"  ) ),
        _velocity( JeveuxVector<ValueType>( getName() +".VITE"  ) ),
        _acceleration( JeveuxVector<ValueType>( getName() +".ACCE"  ) ),
        _projM( new ProjMesu( getName() + ".PROJM" ) ),
        _genDOFNum( nullptr ),
        _DOFNum( nullptr )
    {};

    GeneralizedResult( const std::string &resuTyp ):
        GeneralizedResult( ResultNaming::getNewResultName(), resuTyp )
    {};

    GeneralizedDOFNumberingPtr getGeneralizedDOFNumbering() const
    {
        return _genDOFNum;
    };

    bool setGeneralizedDOFNumbering( const GeneralizedDOFNumberingPtr& genDOFNum )
    {
        _genDOFNum = genDOFNum;

        return true;
    };

    DOFNumberingPtr getDOFNumbering() const
    {
        return _DOFNum;
    };

    bool setDOFNumbering( const DOFNumberingPtr& DOFNum )
    {
        _DOFNum = DOFNum;

        return true;
    };
};

/** @typedef Définition d'un résultat généralisé à valeurs réelles */
template class GeneralizedResult< ASTERDOUBLE >;
typedef GeneralizedResult< ASTERDOUBLE >
    GeneralizedResultReal;
typedef boost::shared_ptr< GeneralizedResultReal >
    GeneralizedResultRealPtr;

/** @typedef Définition d'un résultat généralisé à valeurs complexes */
template class GeneralizedResult< ASTERCOMPLEX >;
typedef GeneralizedResult< ASTERCOMPLEX >
    GeneralizedResultComplex;
typedef boost::shared_ptr< GeneralizedResultComplex >
    GeneralizedResultComplexPtr;

class NonLinearDescriptor
{
    private:
    /** @brief Vecteur Jeveux '.NL.TYPE' */
    JeveuxVectorLong           _type;
    /** @brief Vecteur Jeveux '.NL.VINT' */
    JeveuxVectorReal          _internalVar;
    /** @brief Vecteur Jeveux '.NL.VIND' */
    JeveuxVectorLong            _vIndi;
    /** @brief Vecteur Jeveux '.NL.INTI */
    JeveuxVectorChar24          _vInti;
    public:
    /**
     * @brief Constructeur
    */
    NonLinearDescriptor(const std::string &name):
    _type( JeveuxVectorLong( name +".NL.TYPE")),
    _internalVar( JeveuxVectorReal ( name + ".NL.VINT" )),
    _vIndi( JeveuxVectorLong( name + ".NL.VIND")),
    _vInti( JeveuxVectorChar24( name + ".NL.INTI"))
    {};
};

class TransientGeneralizedResult:
    public GeneralizedResultReal
{
private:
    /** @brief Vecteur Jeveux '.PTEM' */
    /*  valeur du pas de temps aux instants de calcul sauvegardés*/
    JeveuxVectorReal  _timeSteps;
    /** @brief Vecteur Jeveux '.FACC' */
    /* Nom et type des fonctions d’excitation de type accélération */
    JeveuxVectorChar8   _acceExcitFunction;
    /** @brief Vecteur Jeveux '.FVIT' */
    /* Nom et type des fonctions d’excitation de type vitesse */
    JeveuxVectorChar8   _veloExcitFunction;
    /** @brief Vecteur Jeveux '.FDEP' */
    /* Nom et type des fonctions d’excitation de type vitesse */
    JeveuxVectorChar8   _displExcitFunction;
    /** @brief Vecteur Jeveux '.IPSD' */
    JeveuxVectorLong    _ipsd;
    /** @brief Description des nonlinéarités (si mot-clé COMPORTEMENT) */
    NonLinearDescriptor _nonLinDesc;

public:
    /**
     * @typedef TransientGeneralizedResultPtr
     * @brief Pointeur intelligent vers un TransientGeneralizedResult
     */
    typedef boost::shared_ptr< TransientGeneralizedResult >
        TransientGeneralizedResultPtr;

    /**
     * @brief Constructeur
     */
    TransientGeneralizedResult():
        TransientGeneralizedResult( ResultNaming::getNewResultName() )
    {};

    /**
     * @brief Constructeur
     */
    TransientGeneralizedResult( const std::string &name ):
        GeneralizedResultReal( name, "TRAN_GENE" ),
        _nonLinDesc( name ),
        _timeSteps( JeveuxVectorReal( getName() +".PTEM"  ) ),
        _acceExcitFunction(  JeveuxVectorChar8( getName() +".FACC"  ) ),
        _veloExcitFunction(  JeveuxVectorChar8( getName() +".FVIT"  ) ),
        _displExcitFunction(  JeveuxVectorChar8( getName() +".FDEP"  ) ),
        _ipsd( JeveuxVectorLong( getName() + ".IPSD" ) )
    {};


};
typedef boost::shared_ptr< TransientGeneralizedResult >
    TransientGeneralizedResultPtr;

/**
 * @class HarmoGeneralizedResult
 * @brief Cette classe correspond aux concepts  harm_gene,
 * résultats de calcul dynamique harmonique sur base généralisée
 * @author Natacha Béreux
 */
class HarmoGeneralizedResult: public GeneralizedResultComplex
{
private:


public:
    /**
     * @typedef HarmoGeneralizedResultPtr
     * @brief Pointeur intelligent vers un HarmoGeneralizedResult
     */
    typedef boost::shared_ptr< HarmoGeneralizedResult >
        HarmoGeneralizedResultPtr;

    /**
     * @brief Constructeur
     */
    HarmoGeneralizedResult():
        HarmoGeneralizedResult( ResultNaming::getNewResultName() )
    {};

    /**
     * @brief Constructeur
     */
    HarmoGeneralizedResult( const std::string &name ):
        GeneralizedResultComplex( name, "HARM_GENE" )
    {};
};

typedef boost::shared_ptr< HarmoGeneralizedResult >
    HarmoGeneralizedResultPtr;

#endif /* GENERALIZEDRESULTSCONTAINER_H_ */
