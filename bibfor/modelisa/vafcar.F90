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

subroutine vafcar(tpgz, imclf, nmobjz, nutyel, ntyele, car, ncar, ivr, kioc, ier)
!
!
! --------------------------------------------------------------------------------------------------
!
!       VERIFICATION DE LA BONNE AFFECTATION DES DONNEES :
!         CARAC POUTRE      >  ELEMENT POUTRE
!         CARAC DISCRET     >  ELEMENT DISCRET DE TYPE L OU N
!         CARAC COQUE       >  ELEMENT COQUE
!         CARAC ORIENTATION >  ELEMENT DISCRET OU POUTRE
!         CARAC DEFI_ARC    >  ELEMENT POUTRE COURBE
!         CARAC CABLE       >  ELEMENT CABLE
!         CARAC BARRE       >  ELEMENT BARRE
!         CARAC MASSIF      >  ELEMENT THERMIQUE
!         CARAC GRILLE      >  ELEMENT GRILLE
!         CARAC MEMBRANE    >  ELEMENT MEMBRANE
!
! --------------------------------------------------------------------------------------------------
!
    use cara_elem_parameter_module
    implicit none
    integer(kind=8) :: ntyele(*), ivr(*), ncar, nutyel, ier, imclf
    character(len=6) :: kioc
    character(len=*) :: tpgz, nmobjz, car(*)
!
#include "asterfort/utmess.h"
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: ii, ll0, ll1
    integer(kind=8) :: npd, npf
    character(len=8) :: tpg, nomobj, carz
    character(len=17) :: tpge
    character(len=24) :: valk(4)
! --------------------------------------------------------------------------------------------------
    tpg = tpgz
    nomobj = nmobjz
!
!   Vérification de l'affectation de la maille par un élément
    tpge = tpg//' '//nomobj
    if (nutyel .eq. 0) then
        if (ivr(1) .eq. 1) then
            valk(1) = kioc
            valk(2) = ACE_MCLEF(imclf)
            valk(3) = tpge
            call utmess('A', 'MODELISA7_63', nk=3, valk=valk)
            ier = ier+1
        end if
        goto 999
    end if
!
!   Vérification du bon type de l'élément
    npd = 0; npf = 0
    if ((imclf .eq. ACE_POUTRE) .or. (imclf .eq. ACE_POUTRE_FLUI)) then
        npd = 1
        npf = ACE_NB_POUTRE
    else if ((imclf .eq. ACE_DISCRET) .or. (imclf .eq. ACE_DISCRET_2D) .or. &
             (imclf .eq. ACE_RIGI_PARASOL) .or. &
             (imclf .eq. ACE_MASS_AJOU) .or. &
             (imclf .eq. ACE_MASS_REP)) then
        npd = ACE_NB_POUTRE+1
        npf = ACE_NB_POUTRE+ACE_NB_DISCRET
    else if (imclf .eq. ACE_ORIENTATION) then
        npd = 1
        npf = ACE_NB_POUTRE+ACE_NB_DISCRET
    else if (imclf .eq. ACE_COQUE) then
        npd = ACE_NB_POUTRE+ACE_NB_DISCRET+1
        npf = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE
    else if (imclf .eq. ACE_CABLE) then
        npd = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+1
        npf = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ACE_NB_CABLE
    else if (imclf .eq. ACE_BARRE) then
        npd = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ACE_NB_CABLE+1
        npf = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ACE_NB_CABLE+ACE_NB_BARRE
    else if (imclf .eq. ACE_MASSIF) then
        npd = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ACE_NB_CABLE+ACE_NB_BARRE+1
        npf = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ACE_NB_CABLE+ACE_NB_BARRE+ &
              ACE_NB_MASSIF
    else if (imclf .eq. ACE_GRILLE) then
        npd = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ACE_NB_CABLE+ACE_NB_BARRE+ &
              ACE_NB_MASSIF+1
        npf = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ACE_NB_CABLE+ACE_NB_BARRE+ &
              ACE_NB_MASSIF+ACE_NB_GRILLE
    else if (imclf .eq. ACE_MEMBRANE) then
        npd = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ACE_NB_CABLE+ACE_NB_BARRE+ &
              ACE_NB_MASSIF+ACE_NB_GRILLE+1
        npf = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ACE_NB_CABLE+ACE_NB_BARRE+ &
              ACE_NB_MASSIF+ACE_NB_GRILLE+ACE_NB_MEMBRANE
    else
        valk(1) = ACE_MCLEF(imclf)
        call utmess('A', 'MODELISA9_11', sk=valk(1))
    end if
!
    do ii = npd, npf
        if (nutyel .eq. ntyele(ii)) goto 20
    end do
    if (ivr(1) .eq. 1) then
        valk(1) = kioc
        valk(2) = ACE_MCLEF(imclf)
        valk(3) = tpge
        call utmess('A', 'MODELISA7_64', nk=3, valk=valk)
        ier = ier+1
    end if
    goto 999
20  continue
!
! --- CAS PARTICULIER DES ELEMENTS DISCRETS
    if ((imclf .eq. ACE_DISCRET) .or. (imclf .eq. ACE_DISCRET_2D)) then
        ll0 = 0; ll1 = 0
        do ii = 1, ncar
            carz = car(ii)
            if (carz(3:4) .eq. 'T_') then
                ll0 = 1
                ll1 = 5
            else if (carz(3:4) .eq. 'TR') then
                ll0 = 3
                ll1 = 7
            end if
            if (((carz(5:5) .eq. 'N' .or. carz(6:6) .eq. 'N' .or. &
                  carz(7:7) .eq. 'N' .or. carz(8:8) .eq. 'N') &
                 .and. (nutyel .ne. ntyele(ACE_NB_POUTRE+ll0) .and. &
                        nutyel .ne. ntyele(ACE_NB_POUTRE+ll1))) .or. &
                ((carz(5:5) .eq. 'L' .or. carz(6:6) .eq. 'L' .or. carz(7:7) .eq. 'L' .or. &
                  carz(8:8) .eq. 'L') &
                 .and. (nutyel .ne. ntyele(ACE_NB_POUTRE+1+ll0) .and. &
                        nutyel .ne. ntyele(ACE_NB_POUTRE+1+ll1))) .or. &
                ((carz(5:5) .eq. 'D' .or. carz(6:6) .eq. 'D') &
                 .and. (nutyel .ne. ntyele(ACE_NB_POUTRE+ll0) .and. &
                        nutyel .ne. ntyele(ACE_NB_POUTRE+ll1) .and. &
                        nutyel .ne. ntyele(ACE_NB_POUTRE+1+ll0) .and. &
                        nutyel .ne. ntyele(ACE_NB_POUTRE+1+ll1))) .or. &
                ((carz(1:5) .eq. 'M_T_D' .or. carz(1:6) .eq. 'M_TR_D') &
                 .and. (nutyel .eq. ntyele(ACE_NB_POUTRE+1+ll1) .and. &
                        nutyel .ne. ntyele(ACE_NB_POUTRE+ll1)))) then
                if (ivr(1) .eq. 1) then
                    valk(1) = kioc
                    valk(2) = ACE_MCLEF(imclf)
                    valk(3) = tpge
                    valk(4) = carz
                    call utmess('A', 'MODELISA7_65', nk=4, valk=valk)
                    ier = ier+1
                end if
            end if
        end do
    end if
!
999 continue
end subroutine
