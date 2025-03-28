/**
 * @file LibAsterGC.cxx
 * @brief Main of LibAsterGC
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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

extern "C" {
void dintels( double *cequi, double *ht, double *bw, double *enrobi, double *enrobs, double *scmaxi,
              double *scmaxs, double *ssmax, long *uc, long *ntot, double *dnsinf = nullptr,
              double *dnssup = nullptr, double *nrd = nullptr, double *mrd = nullptr,
              long *ndemi = nullptr );

void dintelu( long *typco, double *alphacc, double *ht, double *bw, double *enrobi, double *enrobs,
              double *facier, double *fbeton, double *gammas, double *gammac, long *clacier,
              double *eys, long *typdiag, long *uc, long *ntot, double *dnsinf = nullptr,
              double *dnssup = nullptr, double *nrd = nullptr, double *mrd = nullptr,
              long *ndemi = nullptr );
}

const std::tuple< VectorReal, VectorReal >
dintels_wrapper( double cequi, double ht, double bw, double enrobi, double enrobs, double scmaxi,
                 double scmaxs, double ssmax, long uc, double dnsinf, double dnssup ) {
    long ntot = -1;
    // get size of output vectors
    dintels( &cequi, &ht, &bw, &enrobi, &enrobs, &scmaxi, &scmaxs, &ssmax, &uc, &ntot );

    VectorReal vect_nrd( ntot, 0. );
    VectorReal vect_mrd( ntot, 0. );
    // compute and fill vectors
    dintels( &cequi, &ht, &bw, &enrobi, &enrobs, &scmaxi, &scmaxs, &ssmax, &uc, &ntot, &dnsinf,
             &dnssup, vect_nrd.data(), vect_mrd.data() );

    return std::make_tuple( vect_nrd, vect_mrd );
}

const std::tuple< VectorReal, VectorReal >
dintelu_wrapper( long typco, double alphacc, double ht, double bw, double enrobi, double enrobs,
                 double facier, double fbeton, double gammas, double gammac, long clacier,
                 double eys, long typdiag, long uc, double dnsinf, double dnssup ) {
    long ntot = -1;
    // get size of output vectors
    dintelu( &typco, &alphacc, &ht, &bw, &enrobi, &enrobs, &facier, &fbeton, &gammas, &gammac,
             &clacier, &eys, &typdiag, &uc, &ntot );

    VectorReal vect_nrd( ntot, 0. );
    VectorReal vect_mrd( ntot, 0. );
    // compute and fill vectors
    dintelu( &typco, &alphacc, &ht, &bw, &enrobi, &enrobs, &facier, &fbeton, &gammas, &gammac,
             &clacier, &eys, &typdiag, &uc, &ntot, &dnsinf, &dnssup, vect_nrd.data(),
             vect_mrd.data() );

    return std::make_tuple( vect_nrd, vect_mrd );
}

PYBIND11_MODULE( libAsterGC, mod ) {
    mod.doc() = "This module provides some utilities for reinforced concrete structures";

    mod.def( "dintels", &dintels_wrapper, R"(
Construction du diagramme d'interaction d'une section ferraillée

Vérification d'un ferraillage existant selon le critère : limitation des contraintes (ELS).

Args:
    cequi (float): coefficient d'équivalence acier/beton
    ht (float): hauteur de la section
    bw (float): largeur de la section
    enrobi (float): enrobage des armatures inférieures
    enrobs (float): enrobage des armatures supérieures
    scmaxi (float): contrainte de compression maxi du beton en fibre inf
    scmaxs (float): contrainte de compression maxi du beton en fibre sup
    ssmax (float): contrainte maxi de l'acier de flexion
    uc (int): unite des contraintes : 0 en Pa, 1 en MPa
    dnsinf (float): densité de l'acier inférieur
    dnssup (float): densité de l'acier supérieur

Returns:
    tuple (list[float], list[float]):
    vecteur des efforts normaux résistants (diag inter) et
    vecteur des moments résistants (diag inter).
    )",
             py::arg( "cequi" ), py::arg( "ht" ), py::arg( "bw" ), py::arg( "enrobi" ),
             py::arg( "enrobs" ), py::arg( "scmaxi" ), py::arg( "scmaxs" ), py::arg( "ssmax" ),
             py::arg( "uc" ), py::arg( "dnsinf" ), py::arg( "dnssup" ) );

    mod.def( "dintelu", &dintelu_wrapper, R"(
Construction du diagramme d'interaction d'une section ferraillée

Vérification d'un ferraillage existant selon le critère : limitation des déformations (ELU).

Args:
    typco (int): codification utilisée (1 = bael91, 2 = ec2)
    alphacc (float): coefficient de sécurité sur la résistance de calcul du béton en surpression
    ht (float): hauteur de la section
    bw (float): largeur de la section
    enrobi (float): enrobage des armatures inférieures
    enrobs (float): enrobage des armatures supérieures
    facier (float): limite d'élasticité des aciers (contrainte)
    fbeton (float): résistance en surpression du béton (contrainte)
    gammas (float): coefficient de sécurité sur la résistance de calcul des aciers
    gammac (float): coefficient de sécurité sur la résistance de calcul du béton
    clacier (int): classe de ductilité des aciers (utilisé pour ec2).
        0: acier peu ductile (classe a),
        1: acier moyennement ductile (classe b),
        3: acier fortement ductile (classe c)
    eys (float): module d'young de l'acier
    typdiag (int): type de diagramme utilisé pour l'acier.
        typdiag = 1 ("b1" ==> palier incliné),
        typdiag = 2 ("b2" ==> palier horizontal)
    uc (int): unité des contraintes : 0 en Pa, 1 en MPa
    dnsinf (float): densité de l'acier inférieur
    dnssup (float): densité de l'acier supérieur

Returns:
    tuple (list[float], list[float]):
    vecteur des efforts normaux résistants (diag inter) et
    vecteur des moments résistants (diag inter).
    )",
             py::arg( "typco" ), py::arg( "alphacc" ), py::arg( "ht" ), py::arg( "bw" ),
             py::arg( "enrobi" ), py::arg( "enrobs" ), py::arg( "facier" ), py::arg( "fbeton" ),
             py::arg( "gammas" ), py::arg( "gammac" ), py::arg( "clacier" ), py::arg( "eys" ),
             py::arg( "typdiag" ), py::arg( "uc" ), py::arg( "dnsinf" ), py::arg( "dnssup" )

    );
};
