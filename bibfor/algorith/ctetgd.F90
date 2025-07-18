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

subroutine ctetgd(basmod, numd, numg, nbsec, teta, &
                  nbtet)
    implicit none
!
!***********************************************************************
!    P. RICHARD     DATE 11/03/91
!-----------------------------------------------------------------------
!  BUT:     < CALCUL DE LA MATRICE TETA GAUCHE-DROITE >
!
!   CALCUL DE LA MATRICE TETA PERMETTANT DE PASSER DES  DDL DE
!  L'INTERFACE DROITE A CEUX DE L'INTERFACE GAUCHE
!
!-----------------------------------------------------------------------
!
! BASMOD   /I/: NOM UTLISATEUR DE LA BASE MODALE
! NUMD     /I/: NUMERO DE L'INTERFACE DROITE
! NUMG     /I/: NUMERO DE L'INTERFACE GAUCHE
! NBSEC    /I/: NOMBRE DE SECTEURS
! TETA     /O/: MATRICE CARREE DE CHANGEMENT DE REPERE RECHERCHE
! NBTET    /I/: DIMENSION DELA MATRICE TETA
!
!
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/amppr.h"
#include "asterfort/bmnodi.h"
#include "asterfort/dismoi.h"
#include "asterfort/intet0.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid(1), icomp, iloci, ilocj, inod
    integer(kind=8) :: inog, j, k, lldesc, llnod, llnog
    integer(kind=8) :: nbcmp, nbcpmx, nbdcou, nbddr, nbdga, nbec
    integer(kind=8) :: nbnod, nbnog, nbnot, nbsec, nbtet, noer, numd
    integer(kind=8) :: numg
    real(kind=8) :: angle, pi, x
!-----------------------------------------------------------------------
    parameter(nbcpmx=300)
    character(len=24) :: valk(2)
    character(len=8) :: basmod, mailla, nomnoe, tyd, intf, kbid
    real(kind=8) :: xd(10), xg(10), xtd(10), xtg(10), tet0(10, 10)
    real(kind=8) :: teta(nbtet, nbtet)
    aster_logical :: nook
    integer(kind=8) :: idecd(nbcpmx), idecg(nbcpmx)
    integer(kind=8) :: vali(2), jnocmp
!
!-----------------------------------------------------------------------
!
    data nook/.false./
!
!-----------------------------------------------------------------------
!
    call jemarq()
    pi = r8pi()
!
!-----------------RECUPERATION DES CONCEPTS AMONT-----------------------
!
!
!
    call dismoi('REF_INTD_PREM', basmod, 'RESU_DYNA', repk=intf)
    call dismoi('NOM_MAILLA', intf, 'INTERF_DYNA', repk=mailla)
!
!----------------RECUPERATION DU NOMBRE D'ENTIERS CODES-----------------
!
    call dismoi('NB_CMP_MAX', intf, 'INTERF_DYNA', repi=nbcmp)
    call dismoi('NB_EC', intf, 'INTERF_DYNA', repi=nbec)
    if (nbec .gt. 10) then
        call utmess('F', 'MODELISA_94')
    end if
!
!
!
!-------------------REQUETTE DESCRIPTEUR DES DEFORMEES STATIQUES--------
!
    call jeveuo(intf//'.IDC_DEFO', 'L', lldesc)
    call jelira(intf//'.IDC_DEFO', 'LONMAX', nbnot)
!**************************************************************
    nbnot = nbnot/(2+nbec)
!      NBNOT=NBNOT/3
!**************************************************************
!
!-----------REQUETTE SUR DEFINITION INTERFACES DROITE ET GAUCHE---------
!
    call jeveuo(jexnum(intf//'.IDC_LINO', numd), 'L', llnod)
    call jeveuo(jexnum(intf//'.IDC_LINO', numg), 'L', llnog)
!
!
!--------------RECUPERATION NOMBRE DE NOEUDS AUX INTERFACES-------------
!
    call jelira(jexnum(intf//'.IDC_LINO', numd), 'LONMAX', nbnod)
!
    call jelira(jexnum(intf//'.IDC_LINO', numg), 'LONMAX', nbnog)
!
    if (nbnod .ne. nbnog) then
        vali(1) = nbnod
        vali(2) = nbnog
        call utmess('F', 'ALGORITH14_99', ni=2, vali=vali)
    end if
!
!
!--------------RECUPERATION NOMBRE DE DDL AUX INTERFACES----------------
!
    kbid = ' '
    call bmnodi(basmod, kbid, '          ', numd, 0, &
                ibid, nbddr)
    kbid = ' '
    call bmnodi(basmod, kbid, '          ', numg, 0, &
                ibid, nbdga)
    if (nbdga .ne. nbddr) then
        vali(1) = nbddr
        vali(2) = nbdga
        call utmess('F', 'ALGORITH15_1', ni=2, vali=vali)
    end if
!
!
    if (nbddr .ne. nbtet) then
        vali(1) = nbddr
        vali(2) = nbtet
        call utmess('F', 'ALGORITH15_2', ni=2, vali=vali)
    end if
!
!----------------------CALCUL DU TETA ELEMENTAIRE-----------------------
!
    angle = 2*pi/nbsec
    call intet0(angle, tet0, 3)
!
!
    nbdcou = 0
    call jeveuo(jexnom('&CATA.GD.NOMCMP', 'DEPL_R'), 'L', jnocmp)
    do i = 1, nbnod
        inod = zi(llnod+i-1)
!******************************************************************
!        ICODD=ZI(LLDESC+2*NBNOT+INOD-1)
        inog = zi(llnog+i-1)
!        ICODG=ZI(LLDESC+2*NBNOT+INOG-1)
        call isdeco(zi(lldesc+2*nbnot+(inod-1)*nbec+1-1), idecd, nbcmp)
        call isdeco(zi(lldesc+2*nbnot+(inog-1)*nbec+1-1), idecg, nbcmp)
!******************************************************************
        do j = 1, 10
            if (idecd(j) .eq. 1) then
                xd(j) = 1.d0
            else
                xd(j) = 0.d0
            end if
!
            if (idecg(j) .eq. 1) then
                xg(j) = 1.d0
            else
                xg(j) = 0.d0
            end if
        end do
!
!
        do j = 1, 10
            xtd(j) = 0.d0
            xtg(j) = 0.d0
            do k = 1, 10
                xtd(j) = xtd(j)+abs(tet0(j, k))*xd(k)
                xtg(j) = xtg(j)+abs(tet0(k, j))*xg(k)
            end do
        end do
!
!
!    VERIFICATION SUR COHERENCE DES DDL INTERFACES
!
        do j = 1, 6

            if (xtd(j) .gt. 0.d0 .and. xg(j) .eq. 0.d0) then
                noer = zi(lldesc+inog-1)
                nomnoe = int_to_char8(noer)
                tyd = zk8(jnocmp-1+j)
                call utmess('E', 'ALGORITH15_3')
                valk(1) = tyd
                valk(2) = nomnoe
                call utmess('E', 'ALGORITH15_4', nk=2, valk=valk)
                nook = .true.
            end if
            if (xtg(j) .gt. 0.d0 .and. xd(j) .eq. 0.d0) then
                noer = zi(lldesc+inod-1)
                nomnoe = int_to_char8(noer)
                tyd = zk8(jnocmp-1+j)
                call utmess('E', 'ALGORITH15_3')
                valk(1) = tyd
                valk(2) = nomnoe
                call utmess('E', 'ALGORITH15_6', nk=2, valk=valk)
                nook = .true.
            end if
        end do
!
        if (nook) then
            call utmess('F', 'ALGORITH15_7')
        end if

        do j = 7, nbcmp
            if (idecd(j) .eq. 1.d0) then
                noer = zi(lldesc+inod-1)
                nomnoe = int_to_char8(noer)
                tyd = zk8(jnocmp-1+j)
                valk(1) = tyd
                valk(2) = nomnoe
                call utmess('E', 'ALGORITH15_5', nk=2, valk=valk)
                nook = .true.
            end if
            if (idecg(j) .eq. 1.d0) then
                noer = zi(lldesc+inog-1)
                nomnoe = int_to_char8(noer)
                tyd = zk8(jnocmp-1+j)
                valk(1) = tyd
                valk(2) = nomnoe
                call utmess('E', 'ALGORITH15_5', nk=2, valk=valk)
                nook = .true.
            end if
        end do
!
        if (nook) then
            call utmess('F', 'ALGORITH15_9')
        end if
!
        iloci = 0
        icomp = 0
        do j = 1, 10
            if (idecg(j) .gt. 0) then
                iloci = iloci+1
                ilocj = 0
                icomp = icomp+1
                do k = 1, 10
                    if (idecd(k) .gt. 0) then
                        ilocj = ilocj+1
                        x = tet0(j, k)
                        call amppr(teta, nbddr, nbddr, [x], 1, &
                                   1, nbdcou+iloci, nbdcou+ilocj)
                    end if
                end do
            end if
        end do
!
        nbdcou = nbdcou+icomp
!
    end do
!
    call jedema()
end subroutine
