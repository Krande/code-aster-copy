! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
!
subroutine te0500(option, nomte)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/thmGetElemPara.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/jevech.h"
#include "asterfort/lxlgut.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/assert.h"
!
    character(len=16) :: option, nomte
!    - FONCTION REALISEE:  CALCUL DE L'ESTIMATEUR D'ERREUR TEMPORELLE
!      SUR UN ELEMENT ISOPARAMETRIQUE POUR LES MODELISATIONS HM SATUREES
!
!      --> OPTION 'ERRE_TEMPS_THM'
!
! IN OPTION : NOM DE L'OPTION
! IN NOMTE  : NOM DU TYPE D'ELEMENT
!   -------------------------------------------------------------------
!     SUBROUTINES APPELLEES :
!       MESSAGE              : UTMESS,UTMESK.
!       CHAMPS LOCAUX        : JEVECH,TECACH,TECAEL.
!       ENVIMA               : R8MIEM.
!       MATERIAUX/CHARGES    : RCVALB,RCCOMA.
!       DEDIEES A TE0500     : CAETHM
!     FONCTIONS INTRINSEQUES : SQRT.
!   -------------------------------------------------------------------
!     ASTER INFORMATIONS :
!       25/02/08 (SM) : CREATION POUR CALCUL INDICATEUR D'ERREUR
!                       TEMPORELLE EN INSTATIONNAIRE .
!----------------------------------------------------------------------
!
!
! DECLARATION VARIABLES LOCALES
!
    integer(kind=8) :: nbre1, nbrr1
    parameter(nbre1=1, nbrr1=2)
!
    integer(kind=8) :: nbre2
    parameter(nbre2=2)
!
    integer(kind=8) :: nbre3, nbrr3
    parameter(nbre3=1, nbrr3=2)
!
    integer(kind=8) :: ndim, nno
!
  integer(kind=8) :: ipi, kpi, iaux, npg, igeom, jv_gano, imate, ierre, igrdca, iret, isigap, isigam
    integer(kind=8) :: itab(7), nbcmp
    integer(kind=8) :: dimdep, dimdef, dimcon
    integer(kind=8) :: jv_poids, jv_func, jv_dfunc, jv_poids2, jv_func2, jv_dfunc2
    integer(kind=8) :: nddl_meca, npi, nddl_p1, nddl_p2, nddl_2nd, nnos, nddls, nddlm
    integer(kind=8) :: mecani(5), press1(7), press2(7), tempe(5), second(5), dimuel
!
    real(kind=8) :: poids2
    real(kind=8) :: ovfl
    real(kind=8) :: valre1(nbre1), valrr1(nbrr1), valre2(nbre2), valre3(nbre3)
    real(kind=8) :: valrr3(nbrr3), valrr2(nbrr1), valrr4(nbrr3)
    real(kind=8) :: longc, presc, myoung
    real(kind=8) :: valpar(1), time, raux, rholiq, viscli, permin
    real(kind=8) :: fluhpx, fluhmx, fluhpy, fluhmy
    real(kind=8) :: tertps
!
    aster_logical :: l_axi
!
    integer(kind=8) :: codme1(nbre1), codmr1(nbrr1), codme2(nbre2), codmr2(nbrr1), codme3(nbre3), kpg, spt
    integer(kind=8) :: codmr3(nbrr3), codmr4(nbrr3)
    character(len=3) :: inte_type
    character(len=4) :: nompar(1)
    character(len=8) :: type_elem(2), valk
    character(len=8) :: nomre1(nbre1), nomrr1(nbrr1)
    character(len=8) :: nomrr3(nbrr3), nomrr4(nbrr3), nomre2(nbre2), nomre3(nbre3), fami, poum
    character(len=8) :: nomrr2(nbrr1)
!
    data nomre1/'PERM_IN'/
    data nomrr1/'PERMIN_L', 'PERMIN_N'/
    data nomrr2/'PERMIN_L', 'PERMIN_T'/
    data nomre2/'RHO', 'VISC'/
    data nomre3/'E'/
    data nomrr3/'E_L', 'E_N'/
    data nomrr4/'E_L', 'E_T'/
    type(THM_DS) :: ds_thm
!
! ------------------------------------------------------------------
!
    ovfl = r8miem()
!
! - Get all parameters for current element
!
    call thmGetElemPara(ds_thm, l_axi, &
                        type_elem, inte_type, ndim, &
                        mecani, press1, press2, tempe, second, &
                        dimdep, dimdef, dimcon, dimuel, &
                        nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                        nno, nnos, &
                        npi, npg, &
                        jv_poids, jv_func, jv_dfunc, &
                        jv_poids2, jv_func2, jv_dfunc2, &
                        jv_gano)

! =====================================================================
! 2. RECUPERATION DES PARAMETRES TEMPORELS
! =====================================================================
    call tecach('ONO', 'PINSTR', 'L', iret, iad=itab(1))
    if (iret .eq. 0) then
        time = zr(itab(1))
    else
        call utmess('F', 'INDICATEUR_11')
    end if
! =====================================================================
! 3. INITIALISATIONS/RECUPERATION DE LA GEOMETRIE ET DES CHAMPS LOCAUX
! =====================================================================
!
! 3.1. GEOMETRIE (IGEOM)
!
    call jevech('PGEOMER', 'L', igeom)
!
! 3.2. MATERIAU (IMATE)
!
    call jevech('PMATERC', 'L', imate)
!
! 3.3 CONTRAINTES ( T- ET T+ )
!
    call jevech('PCONTGM', 'L', isigam)
    call tecach('ONO', 'PCONTGP', 'L', iret, nval=3, itab=itab)
!
    isigap = itab(1)
    nbcmp = itab(2)/npi
!
! 3.4  RECHERCHE DES VALEURS NECESSAIRES AU CALCUL DE L'INDICATEUR
!
! --- A. GRANDEURS CARACTERISTIQUES
!
    call jevech('PGRDCA', 'L', igrdca)
    longc = zr(igrdca)
    presc = zr(igrdca+1)
!
! --- B. PERMEABILITE INTRINSEQUE DU MILIEU
!
! => PERMIN SI ISOTROPE
! => PERMIN_X,PERMIN_Y ET PERMIN_Z SINON
!
    nompar(1) = 'INST'
    valpar(1) = time
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'THM_DIFFU', 1, nompar, [valpar], &
                nbre1, nomre1, valre1, codme1, 0, nan='OUI')
!
    if (codme1(1) .eq. 0) then
        permin = valre1(1)
    else if (codme1(1) .eq. 1) then
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'THM_DIFFU', 1, nompar, [valpar], &
                    nbrr1, nomrr1, valrr1, codmr1, 0, nan='OUI')
        if ((codmr1(1) .eq. 0) .and. (codmr1(2) .eq. 0)) then
            permin = sqrt(valrr1(1)**2+valrr1(2)**2+valrr1(1)**2)
        else
            call rcvalb(fami, kpg, spt, poum, zi(imate), &
                        ' ', 'THM_DIFFU', 1, nompar, [valpar], &
                        nbrr1, nomrr2, valrr2, codmr2, 0, nan='OUI')
            if ((codmr2(1) .eq. 0) .and. (codmr2(2) .eq. 0)) then
                permin = sqrt(valrr2(1)**2+valrr2(2)**2)
            else
                ASSERT(.false.)
            end if
        end if
    else
        call utmess('F', 'ELEMENTS4_78', sk=nomre1(1))
    end if
!
! --- C. MASSE VOLUMIQUE DU LIQUIDE
!        VISCOSITE DYNAMIQUE DU LIQUIDE
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'THM_LIQU', 1, nompar, [valpar], &
                nbre2, nomre2, valre2, codme2, 1)
!
    if (codme2(1) .eq. 0 .and. codme2(2) .eq. 0) then
        rholiq = valre2(1)
        viscli = valre2(2)
    else
        call utmess('F', 'ELEMENTS4_69', sk=nomre2(1)//nomre2(2))
    end if
!
! --- D. MODULE DE YOUNG
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'ELAS', 1, nompar, [valpar], &
                nbre3, nomre3, valre3, codme3, 0)
!
    if (codme3(1) .eq. 0) then
        myoung = valre3(1)
    else if (codme3(1) .eq. 1) then
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'ELAS_ISTR', 1, nompar, [valpar], &
                    nbrr3, nomrr3, valrr3, codmr3, 0)
        if ((codmr3(1) .eq. 0) .and. (codmr3(2) .eq. 0)) then
            myoung = sqrt(valrr3(1)**2+valrr3(2)**2)
        else
            call rcvalb(fami, kpg, spt, poum, zi(imate), &
                        ' ', 'ELAS_ORTH', 1, nompar, [valpar], &
                        nbrr3, nomrr4, valrr4, codmr4, 0)
            if ((codmr4(1) .eq. 0) .and. (codmr4(2) .eq. 0)) then
                myoung = sqrt(valrr4(1)**2+valrr4(2)**2)
            else
                ASSERT(.false.)
            end if
        end if
    else
        call utmess('F', 'ELEMENTS4_71', sk=nomre3(1))
    end if
!
! 3.5 CALCUL DU COEFFICIENT D'ADIMENSIONNEMENT
!
    if (abs(longc) .gt. ovfl) then
        raux = (myoung*viscli)/(longc**ndim)
!
        if (rholiq .gt. ovfl) then
            raux = raux/(rholiq**2)
!
            if (presc .gt. ovfl) then
                raux = raux/(presc**2)
!
                if (permin .gt. ovfl) then
                    raux = raux/permin
!
                else
                    call utmess('F', 'INDICATEUR_20')
                end if
            else
                valk = 'pression'
                call utmess('F', 'INDICATEUR_21', sk=valk)
            end if
        else
            call utmess('F', 'INDICATEUR_22')
        end if
    else
        valk = 'longueur'
        call utmess('F', 'INDICATEUR_21', sk=valk)
    end if
!
! =====================================================================
! 4. CALCUL DE L'INDICATEUR TEMPOREL
! =====================================================================
!
! 4.1. INITIALISATION
!
    tertps = 0.d0
!
! --- BOUCLE SUR LES POINTS DE GAUSS
!
    do ipi = 1, npg
!
        kpi = ipi
!
        if (ndim .eq. 2) then
! =====================================================================
! => EN DIMENSION 2
! =====================================================================
!
! 4.2. ON RECUPERE LES POIDS D'INTEGRATION AUX POINTS DE GAUSS
!
            call dfdm2d(nnos, kpi, jv_poids2, jv_dfunc2, zr(igeom), &
                        poids2)
!
            iaux = nbcmp*(kpi-1)
!
            fluhpx = zr(isigap+iaux+8+5)
            fluhmx = zr(isigam+iaux+8+5)
!
            fluhpy = zr(isigap+iaux+9+5)
            fluhmy = zr(isigam+iaux+9+5)
!
            tertps = tertps+raux*poids2*((fluhpx-fluhmx)**2+(fluhpy-fluhmy)**2)
!
        else
            iaux = lxlgut(option)
            call utmess('F', 'INDICATEUR_92', sk=option(1:iaux))
        end if
!
    end do
!
! --- FIN BOUCLE SUR LES POINTS DE GAUSS
!
! =====================================================================
! 5. STOCKAGE DE L'ERREUR
! =====================================================================
!
    call jevech('PERREUR', 'E', ierre)
!
    zr(ierre) = tertps
!
!
end subroutine
