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

subroutine vtcop1(chin, chout, kstop, codret)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/vrrefe.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
    character(len=*) :: chin, chout
    character(len=1) :: kstop
    integer :: codret
! person_in_charge: nicolas.sellenet at edf.fr
!     APPELLE PAR LA ROUTINE CHAPEAU VTCOPY
!     RECOPIE LES VALEURS DU CHAM_NO CHIN DANS LE CHAM_NO CHOUT
!     CETTE ROUTINE PERMET DE CHANGER LA NUMEROTATION D'UN CHAM_NO
!     SI KSTOP.NE.'F' EN ENTREE ET QUE CODRET != 0 EN SORTIE, ALORS
!     DES COMPOSANTES DU CHAMP CHIN N'ONT PAS PU ETRE RECOPIEES ET
!     ONT ETE MISES A ZERO (ON ARRETE LA ROUTINE SI KSTOP .EQ. 'F')
!
!     PRECAUTIONS D'EMPLOI :
!     - LES CHAM_NOS DOIVENT EXISTER.
!     - LES DDLS DE "LAGRANGE" SONT MIS A ZERO DANS CHOUT.
!
!     ------------------------------------------------------------------
!
!
!
    integer :: iret, ieq1, ieq2, neq1, jvale1, jvale2
    integer :: neq2
    integer :: nnomx, ncpmx, nuno2, nucp2, nuno1, nucp1
    integer :: jcmpgd, ncmpmx, icmp
    character(len=1) :: typ1, typ2
    character(len=8) :: nomgd, mesh1, mesh2
    character(len=24) :: valk(4)
    character(len=19) :: ch1, ch2, pfchno
    integer, pointer :: trav1(:) => null()
    aster_logical, pointer :: trav2(:) => null()
    character(len=24), pointer :: refe1(:) => null()
    character(len=24), pointer :: refe2(:) => null()
    integer, pointer :: deeq1(:) => null()
    integer, pointer :: deeq2(:) => null()
    integer, pointer :: deeq(:) => null()
!
!     ------------------------------------------------------------------
!
    call jemarq()
    codret = 0
    ch1 = chin
    ch2 = chout
!
    call vrrefe(ch1, ch2, iret)
!
!
!     1. SI LES 2 CHAMPS ONT LES MEMES NUMEROTATIONS :
!     -------------------------------------------------
    if (iret .eq. 0) then
        call jelira(ch1//'.VALE', 'TYPE', cval=typ1)
        call jelira(ch1//'.VALE', 'LONMAX', neq1)
        call jeveuo(ch1//'.VALE', 'L', jvale1)
        call jelira(ch2//'.VALE', 'TYPE', cval=typ2)
        call jeveuo(ch2//'.VALE', 'E', jvale2)
        call dismoi('NUME_EQUA', ch2, 'CHAM_NO', repk=pfchno)
        call jeveuo(pfchno//'.DEEQ', 'L', vi=deeq)
        if (typ1 .eq. typ2) then
            if (typ1 .eq. 'R') then
                do ieq1 = 0, neq1-1
                    if (deeq(2*ieq1+2) .le. 0) then
                        zr(jvale2+ieq1) = 0.d0
                    else
                        zr(jvale2+ieq1) = zr(jvale1+ieq1)
                    end if
                end do
            else if (typ1 .eq. 'C') then
                do ieq1 = 0, neq1-1
                    if (deeq(2*ieq1+2) .le. 0) then
                        zc(jvale2+ieq1) = dcmplx(0.d0, 0.d0)
                    else
                        zc(jvale2+ieq1) = zc(jvale1+ieq1)
                    end if
                end do
            else
                valk(1) = ch1
                valk(2) = ch2
                valk(3) = typ1
                call utmess('F', 'ALGELINE3_93', nk=3, valk=valk)
            end if
        else
            if (typ1 .eq. 'R' .and. typ2 .eq. 'C') then
                do ieq1 = 0, neq1-1
                    if (deeq(2*ieq1+2) .le. 0) then
                        zc(jvale2+ieq1) = dcmplx(0.d0, 0.d0)
                    else
                        zc(jvale2+ieq1) = zr(jvale1+ieq1)
                    end if
                end do
            else
                valk(1) = ch1
                valk(2) = typ1
                valk(3) = ch2
                valk(4) = typ2
                call utmess('F', 'ALGELINE3_94', nk=4, valk=valk)
            end if
        end if
        goto 999
    end if
!
!
!     2. SI LES 2 CHAMPS N'ONT PAS LES MEMES NUMEROTATIONS :
!     ------------------------------------------------------
    call jelira(ch1//'.VALE', 'TYPE', cval=typ1)
    call jelira(ch1//'.VALE', 'LONMAX', neq1)
    call jeveuo(ch1//'.VALE', 'L', jvale1)
    call jelira(ch2//'.VALE', 'TYPE', cval=typ2)
    call jelira(ch2//'.VALE', 'LONMAX', neq2)
    call jeveuo(ch2//'.VALE', 'E', jvale2)
!
    call jeveuo(ch1//'.REFE', 'L', vk24=refe1)
    call jeveuo(ch2//'.REFE', 'L', vk24=refe2)
    call dismoi("NOM_MAILLA", ch1, "CHAM_NO", repk=mesh1)
    call dismoi("NOM_MAILLA", ch2, "CHAM_NO", repk=mesh2)
    if (mesh1 .ne. mesh2) then
        valk(1) = ch1
        valk(2) = ch2
        valk(3) = mesh1
        valk(4) = mesh2
        call utmess('F', 'CALCULEL2_1', nk=4, valk=valk)
    end if
    call jeveuo(refe1(2) (1:19)//'.DEEQ', 'L', vi=deeq1)
    call jeveuo(refe2(2) (1:19)//'.DEEQ', 'L', vi=deeq2)
!
!
!     2.1 ON CHERCHE LE NUMERO DE CMP LE PLUS GRAND ET
!     LE NUMERO DE NOEUD LE PLUS GRAND DANS CH2->.DEEQ2 :
!     ---------------------------------------------------------------
!
    nnomx = 0
    ncpmx = 0
    do ieq2 = 1, neq2
        nnomx = max(nnomx, deeq2(2*(ieq2-1)+1))
        ncpmx = max(ncpmx, deeq2(2*(ieq2-1)+2))
    end do
!
!
!     2.2 ON REMPLIT UN OBJET DE TRAVAIL :
!     ------------------------------------
    AS_ALLOCATE(vi=trav1, size=nnomx*ncpmx)
    AS_ALLOCATE(vl=trav2, size=neq2)
    do ieq2 = 1, neq2
        nuno2 = deeq2(2*(ieq2-1)+1)
        nucp2 = deeq2(2*(ieq2-1)+2)
        if (nucp2 .gt. 0) trav1((nuno2-1)*ncpmx+nucp2) = ieq2
        trav2(ieq2) = .false.
    end do
!
!
!     2.3 ON RECOPIE LES VALEURS DE CH1 DANS CH2 :
!     -------------------------------------------
    if (typ1 .eq. typ2) then
        if (typ1 .eq. 'R') then
            do ieq1 = 1, neq1
                nuno1 = deeq1(2*(ieq1-1)+1)
                nucp1 = deeq1(2*(ieq1-1)+2)
                if ((nucp1 .gt. 0) .and. (nuno1 .le. nnomx) .and. (nucp1 .le. ncpmx)) then
                    ieq2 = trav1((nuno1-1)*ncpmx+nucp1)
                    if (ieq2 .gt. 0) then
                        trav2(ieq2) = .true.
                        zr(jvale2-1+ieq2) = zr(jvale1-1+ieq1)
                    end if
                end if
            end do
        else if (typ1 .eq. 'C') then
            do ieq1 = 1, neq1
                nuno1 = deeq1(2*(ieq1-1)+1)
                nucp1 = deeq1(2*(ieq1-1)+2)
                if ((nucp1 .gt. 0) .and. (nuno1 .le. nnomx)) then
                    ieq2 = trav1((nuno1-1)*ncpmx+nucp1)
                    if (ieq2 .gt. 0) then
                        trav2(ieq2) = .true.
                        zc(jvale2-1+ieq2) = zc(jvale1-1+ieq1)
                    end if
                end if
            end do
        else
            valk(1) = ch1
            valk(2) = ch2
            valk(3) = typ1
            call utmess('F', 'ALGELINE3_93', nk=3, valk=valk)
        end if
!
    else if (typ1 .eq. 'R' .and. typ2 .eq. 'C') then
        do ieq1 = 1, neq1
            nuno1 = deeq1(2*(ieq1-1)+1)
            nucp1 = deeq1(2*(ieq1-1)+2)
            if ((nucp1 .gt. 0) .and. (nuno1 .le. nnomx)) then
                ieq2 = trav1((nuno1-1)*ncpmx+nucp1)
                if (ieq2 .gt. 0) then
                    trav2(ieq2) = .true.
                    zc(jvale2-1+ieq2) = zr(jvale1-1+ieq1)
                end if
            end if
        end do
!
    else
        valk(1) = ch1
        valk(2) = typ1
        valk(3) = ch2
        valk(4) = typ2
        call utmess('F', 'ALGELINE3_94', nk=4, valk=valk)
    end if
!
!     A CAUSE DE LA SOUS-STRUCTURATION STATIQUE, ON DOIT AJOUTER
!     UNE GLUTE POUR OUBLIER LA COMPOSANTE 'LAGR' POUR LA VERIF
    call dismoi('NOM_GD', ch2, 'CHAM_NO', repk=nomgd)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', jcmpgd)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgd), 'LONMAX', ncmpmx)
    icmp = -200
    icmp = indik8(zk8(jcmpgd), 'LAGR', 1, ncmpmx)
!
    do ieq2 = 1, neq2
        nuno2 = deeq2(2*(ieq2-1)+1)
        nucp2 = deeq2(2*(ieq2-1)+2)
!       NUCP2.NE.ICMP == GLUTE POUR LA SOUS-STRUCTURATION STATIQUE
        if (nucp2 .gt. 0 .and. nucp2 .ne. icmp .and. .not. trav2(ieq2)) then

            valk(1) = zk8(jcmpgd+nucp2-1)
            valk(2) = int_to_char8(nuno2)
            valk(3) = ch1
            codret = 1

            if (kstop .eq. 'F') then
                call utmess('F', 'ALGELINE7_20', nk=3, valk=valk)
            else
                call utmess('A', 'ALGELINE7_20', nk=3, valk=valk)
                exit
            end if
        end if
    end do
    AS_DEALLOCATE(vi=trav1)
    AS_DEALLOCATE(vl=trav2)
!
999 continue
    call jedema()
end subroutine
