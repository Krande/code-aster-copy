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
subroutine cgcrtb(table, option, ndim, typfis, nxpara, &
                  lmoda, nbpara, linopa, litypa)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cgajpa.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
    integer(kind=8) :: nbpara, nxpara, ndim
    aster_logical :: lmoda
    character(len=*) :: litypa(nxpara), linopa(nxpara)
    character(len=8) :: table, typfis
    character(len=16) :: option
!
! person_in_charge: samuel.geniaut at edf.fr
!
!     SOUS-ROUTINE DE L'OPERATEUR CALC_G
!
!     BUT : CREATION DE LA TABLE ISSUE DE CALC_G
!           ET AFFECTATION DES PARAMETRES
!
! ----------------------------------------------
!  IN :
!     TABLE : NOM DE LA TABLE
!     OPTION : OPTION DE CALCUL
!     NDIM    : DIMENSION DU CALCUL : 2 OU 3
!     TYPFIS : TYPE D'OBJET POUR DECRIRE LE FOND DE FISSURE
!              'FONDFISS' OU 'FISSURE' OU 'THETA'
!     NXPARA : NOMBRE MAXI DE PARAMETRES DE LA TABLE
!     LMODA  : .TRUE.  SI TYPE SD RESULTAT = MODE_MECA
!              .FALSE. SINON
!
!  OUT :
!     NBPARA : NOMBRE DE PARAMETRES
!     linopa : NOMS DES PARAMETRES
!     litypa : TYPES DES PARAMETRES
! ----------------------------------------------
!
    integer(kind=8) :: i
!
    aster_logical :: debug
!
    nbpara = 0
    debug = .false.
!
!   --------------------
!   EXCLUSION DES OPTIONS A SUPPRIMER G_BILI(_GLOB) et G_MAX(_GLOB) et CALC_K_MAX
!   ---------------------
    if ((option .ne. 'G_BILI') .and. (option .ne. 'G_BILI_GLOB') .and. (option .ne. 'G_MAX') &
        .and. (option .ne. 'G_MAX_GLOB') .and. (option .ne. 'CALC_K_MAX')) then
!-------------------------------------------------------------------------
!   --------------------
!   1. IDENTIFICATEURS
!   --------------------
!   --------------------
!   1.1 FOND DE FISSURE
!   ---------------------
        if (typfis .ne. 'THETA') then
            call cgajpa('NUME_FOND', 'I', nbpara, linopa, litypa, &
                        nxpara)
        end if
!   --------------------
!   1.2 TEMPOREL/CHARGEMENT
!   ---------------------
        if (lmoda) then
            call cgajpa('NUME_MODE', 'I', nbpara, linopa, litypa, &
                        nxpara)
        else
            call cgajpa('NUME_ORDRE', 'I', nbpara, linopa, litypa, &
                        nxpara)
            call cgajpa('INST', 'R', nbpara, linopa, litypa, &
                        nxpara)
        end if
!   --------------------
!   1.3 POINT DU FOND DE FISSURE
!   ---------------------
        if (typfis .ne. 'THETA') then
            if (typfis .eq. 'FONDFISS') then
                call cgajpa('NOEUD', 'K8', nbpara, linopa, litypa, &
                            nxpara)
            end if
!            if (typfis.ne.'THETA') then
            call cgajpa('COOR_X', 'R', nbpara, linopa, litypa, &
                        nxpara)
            call cgajpa('COOR_Y', 'R', nbpara, linopa, litypa, &
                        nxpara)
            if (ndim .eq. 3) then
                call cgajpa('COOR_Z', 'R', nbpara, linopa, litypa, &
                            nxpara)
                call cgajpa('NUM_PT', 'I', nbpara, linopa, litypa, &
                            nxpara)
                call cgajpa('ABSC_CURV', 'R', nbpara, linopa, litypa, &
                            nxpara)
            end if
!            endif
        end if
!   --------------------
!   2. OPTIONS DE CALCUL
!   --------------------
!   --------------------
!   2.1 G COMMUN A TOUTES LES OPTIONS
!   ---------------------
        call cgajpa('G', 'R', nbpara, linopa, litypa, &
                    nxpara)
!   --------------------
!   2.2 CALC_K_G
!   ---------------------
        if (option .eq. 'CALC_K_G') then
            call cgajpa('K1', 'R', nbpara, linopa, litypa, &
                        nxpara)
            call cgajpa('K2', 'R', nbpara, linopa, litypa, &
                        nxpara)
            call cgajpa('G_IRWIN', 'R', nbpara, linopa, litypa, &
                        nxpara)
            if (ndim .eq. 3) then
                call cgajpa('K3', 'R', nbpara, linopa, litypa, &
                            nxpara)
            end if
        end if
!-------------------------------------------------------------------------
!   --------------------
!   2.3 OPTIONS A SUPPRIMER (G_BILI(_GLOB) et G_MAX(_GLOB) et CALC_K_MAX)
!   ---------------------
    else if (option .eq. 'CALC_K_MAX') then
        nbpara = 14
        linopa(1) = 'NUME_FOND'
        litypa(1) = 'I'
        linopa(2) = 'NUME_ORDRE'
        litypa(2) = 'I'
        linopa(3) = 'INST'
        litypa(3) = 'R'
!
        if (typfis .eq. 'FONDFISS') then
            linopa(4) = 'NOEUD'
            litypa(4) = 'K8'
        end if
        linopa(5) = 'COOR_X'
        litypa(5) = 'R'
        linopa(6) = 'COOR_Y'
        litypa(6) = 'R'
        linopa(7) = 'COOR_Z'
        litypa(7) = 'R'
        linopa(8) = 'NUM_PT'
        litypa(8) = 'I'
        linopa(9) = 'ABSC_CURV'
        litypa(9) = 'R'
        linopa(10) = 'K1'
        litypa(10) = 'R'
        linopa(11) = 'K2'
        litypa(11) = 'R'
        linopa(12) = 'K3'
        litypa(12) = 'R'
        linopa(13) = 'G'
        litypa(13) = 'R'
        linopa(14) = 'G_IRWIN'
        litypa(14) = 'R'
    else if ((option .eq. 'G_BILI') .or. (option .eq. 'G_MAX')) then
        nbpara = 6
        linopa(1) = 'INST'
        litypa(1) = 'R'
        linopa(2) = 'NUME_CMP_I'
        litypa(2) = 'I'
        linopa(3) = 'NUME_CMP_J'
        litypa(3) = 'I'
        linopa(4) = 'NOEUD'
        litypa(4) = 'K8'
        linopa(5) = 'ABSC_CURV'
        litypa(5) = 'R'
        linopa(6) = 'G_BILI_LOCAL'
        litypa(6) = 'R'
    else if ((option .eq. 'G_BILI_GLOB') .or. (option .eq. 'G_MAX_GLOB')) then
        nbpara = 3
        linopa(1) = 'NUME_CMP_I'
        litypa(1) = 'I'
        linopa(2) = 'NUME_CMP_J'
        litypa(2) = 'I'
        linopa(3) = 'G_BILIN'
        litypa(3) = 'R'
    end if
!
!   --------------------
!   3. CREATION DE LA TABLE
!   --------------------
    call tbcrsd(table, 'G')
    call tbajpa(table, nbpara, linopa, litypa)
!
!   --------
!   4. DEBUG
!   --------
    if (debug) then
        write (6, *) 'OPTION = ', option
        write (6, *) 'NOMBRE DE PARAMETRES DE LA TABLE = ', nbpara
        write (6, *) 'NO_PARA, NOM_PARA, TYP_PARA'
        do i = 1, nbpara
            write (6, *) i, linopa(i), litypa(i)
        end do
    end if
!
end subroutine
