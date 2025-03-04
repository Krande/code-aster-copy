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
void dintels( ASTERDOUBLE *cequi, ASTERDOUBLE *ht, ASTERDOUBLE *bw, ASTERDOUBLE *enrobi,
              ASTERDOUBLE *enrobs, ASTERDOUBLE *scmaxi, ASTERDOUBLE *scmaxs, ASTERDOUBLE *ssmax,
              ASTERINTEGER *uc, ASTERDOUBLE *dnsinf, ASTERDOUBLE *dnssup, ASTERINTEGER *ntot,
              ASTERDOUBLE *nrd, ASTERDOUBLE *mrd );
}

const std::tuple< VectorReal, VectorReal >
dintels_wrapper( ASTERDOUBLE cequi, ASTERDOUBLE ht, ASTERDOUBLE bw, ASTERDOUBLE enrobi,
                 ASTERDOUBLE enrobs, ASTERDOUBLE scmaxi, ASTERDOUBLE scmaxs, ASTERDOUBLE ssmax,
                 ASTERINTEGER uc, ASTERDOUBLE dnsinf, ASTERDOUBLE dnssup, ASTERINTEGER ntot ) {
    VectorReal vect_nrd( ntot, 0. );
    VectorReal vect_mrd( ntot, 0. );
    dintels( &cequi, &ht, &bw, &enrobi, &enrobs, &scmaxi, &scmaxs, &ssmax, &uc, &dnsinf, &dnssup,
             &ntot, vect_nrd.data(), vect_mrd.data() );
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
    uc (float): unite des contraintes : 0 en Pa, 1 en MPa
    dnsinf (float): densité de l'acier inférieur
    dnssup (float): densité de l'acier supérieur
    ntot (int): dimensions des vecteurs

Returns:
    tuple (list[float], list[float]):
    vecteur des efforts normaux résistants (diag inter) et
    vecteur des moments résistants (diag inter).
    )",
             py::arg( "cequi" ), py::arg( "ht" ), py::arg( "bw" ), py::arg( "enrobi" ),
             py::arg( "enrobs" ), py::arg( "scmaxi" ), py::arg( "scmaxs" ), py::arg( "ssmax" ),
             py::arg( "uc" ), py::arg( "dnsinf" ), py::arg( "dnssup" ), py::arg( "ntot" ) );
};
