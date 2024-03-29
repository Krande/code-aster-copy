! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine xmfrot(algofr, coeffr, coeffp, ddlm, ddls, &
                  ffc, ffp, idepd, idepm, indco, &
                  jac, lact, mmat, mu, nd, &
                  ndim, nfh, nfiss, nno, nnol, &
                  nnos, nvit, pla, seuil, &
                  singu, fk, tau1, tau2)
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
! IN ALGOFR : ALGO FROTTEMENT (1:LAG, 2:PENA, 0:RIEN)
! IN CFACE  : CONNECTIVITE FACETTES DE CONTACT
! IN COEFFR : COEF AUGMENTATION FROT
! IN COEFFP : COEF PENALISATION FROT
! IN DDLM   : NOMBRE DE DDLS A CHAQUE NOEUD MILIEU
! IN DDLS   : NOMBRE DE DDLS A CHAQUE NOEUD SOMMET
! IN FFC    : FONCTIONS DE FORME DE CONTACT
! IN FFP    : FONCTIONS DE FORME ELEMENT PARENT
! IN IDEPD  : ADRESSE INCREMENT DEPLACEMENT COURANT
! IN IDEPM  : ADRESSE DEPLACEMENT INSTANT -
! IN IFA    : NUMERO FACETTE DE CONTACT
! IN INDCO  : ETAT DE CONTACT POINT DE GAUSS
! IN IPGF   : NUMERO POINT DE GAUSS DE CONTACT
! IN IVFF   : ADRESSE FONCTION DE FORME EL PARENT
! IN JAC    : PRODUIT JACOBIEN*POIDS
! IN LACT   : DDL DE LAGRANGE ACTIF OU NON
! OUT MMAT  : MATRICE ELEMENTAIRE DE CONTACT
! IN MU     : COEFFICIENT DE COULOMB
! IN ND     : NORMALE A LA SURFACE DE CONTACT AU PG
! IN NDIM   : DIMENSION DU MODELE
! IN NFH    : NOMBRE DE DDL HEAVISIDE
! IN NFISS  : NOMBRE DE FISSURES
! IN NNO    : NOMBRE DE NOEUDS TOTAL ELEMENT PARENT
! IN NNOL   : NOMBRE DE NOEUDS EL PARENT PORTEURS DE DDL LAGRANGE
! IN NNOS   : NOMBRE DE NOEUDS SOMMET ELEMENT PARENT
! IN NOEUD  : FORMULATION AUX NOEUDS
! IN NVIT   : ARETE VITALE OU NON
! IN PLA    : PLACE DES DDLS DE LAGRANGE
! IN SINGU  : ELEMENT ENRICHI CTIP OU ON
! IN TAU1   : 1ERE TANGENTE SURFACE DE CONTACT
! IN TAU2   : 2EME TANGENTE (3D)
#include "asterfort/assert.h"
#include "asterfort/xmmab3.h"
#include "asterfort/xmmab4.h"
#include "asterfort/xmmab5.h"
#include "asterfort/xmmab6.h"
#include "asterfort/xmmbp3.h"
#include "asterfort/xmmbp5.h"
#include "asterfort/xmmsa1.h"
    integer :: algofr, ddlm, ddls
    integer :: idepd, idepm, indco
    integer :: lact(8), ndim
    integer :: nfh, nfiss, nno, nnol, nnos, nvit
    integer :: pla(27), singu
    real(kind=8) :: coeffp, coeffr, ik(3, 3), jac, knp(3, 3)
    real(kind=8) :: mmat(216, 216), mu, nd(3), p(3, 3), ffc(8), ffp(27)
    real(kind=8) :: ptknp(3, 3), seuil, tau1(3), tau2(3)
    real(kind=8) :: fk(27, 3, 3)
    aster_logical :: adher
!
    if (mu .eq. 0.d0 .or. seuil .eq. 0.d0) indco = 0
    if (nfiss .gt. 1) indco = 0
!
    if (indco .eq. 0) then
        if (nvit .ne. 0) then
!
! --- CALCUL DE LA MATRICE F - CAS SANS CONTACT
!
            call xmmab6(ndim, nnol, pla, ffc, jac, &
                        tau1, tau2, lact, mmat)
!
        end if
    else if (indco .eq. 1) then
!
! --- CALCUL DES INCREMENTS - DÉPLACEMENTS&
! --- SEMI-MULTIPLICATEUR DE FROTTEMENT
!
        call xmmsa1(algofr, ndim, nno, nnos, nnol, &
                    pla, ffc, ffp, idepd, idepm, &
                    nfh, nd, tau1, tau2, singu, &
                    fk, lact, ddls, ddlm, coeffr, &
                    coeffp, p, adher, knp, ptknp, &
                    ik)
!
! --- CALCUL DE B, BT
!
        if (algofr .eq. 1) then
            call xmmab3(ndim, nno, nnos, nnol, pla, &
                        ffc, ffp, jac, knp, nfh, &
                        seuil, tau1, tau2, mu, singu, &
                        fk, lact, ddls, ddlm, mmat)
!
! --- CALCUL DE B_U
!
            call xmmab4(ndim, nno, nnos, ffp, jac, &
                        ptknp, nfh, seuil, mu, singu, &
                        fk, coeffr, ddls, ddlm, mmat)
!
! --- CALCUL DE F (NULLE EN ADHERENCE)
!
            if (.not. adher) then
                call xmmab5(ndim, nnol, pla, ffc, jac, &
                            coeffr, seuil, tau1, tau2, mu, &
                            ik, lact, mmat)
            end if
        else if (algofr .eq. 2) then
            call xmmbp3(ndim, nno, nnos, nnol, pla, &
                        ffc, ffp, jac, knp, nfh, &
                        seuil, tau1, tau2, mu, singu, &
                        fk, lact, ddls, ddlm, mmat)
!
! --- CALCUL DE B_U
!
            call xmmab4(ndim, nno, nnos, ffp, jac, &
                        ptknp, nfh, seuil, mu, singu, &
                        fk, coeffp, ddls, ddlm, mmat)
!
! --- CALCUL DE F - CAS GLISSANT OU PENALISATION (NON SEULE)
!
            call xmmbp5(ndim, nnol, pla, ffc, jac, &
                        coeffp, seuil, tau1, tau2, mu, &
                        lact, mmat)
        else
            ASSERT(algofr .eq. 0)
        end if
!
    else
        ASSERT(indco .eq. 0 .or. indco .eq. 1)
    end if
end subroutine
