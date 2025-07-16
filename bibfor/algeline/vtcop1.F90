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
subroutine vtcop1(fieldNodeInZ, fieldNodeOutZ, codret)
!
    implicit none
!
#include "asterc/indik8.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/vrrefe.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: fieldNodeInZ, fieldNodeOutZ
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
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
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, ieq1, ieq2, nbEquaIn, jvValeIn, jvValeOut
    integer(kind=8) :: nbEquaOut
    integer(kind=8) :: nnomx, ncpmx, nuno2, nucp2, nuno1, nucp1
    integer(kind=8) :: jcmpgd, ncmpmx, cmpLagr
    character(len=1) :: scalTypeIn, scalTypeOut
    character(len=8) :: nomgd, meshIn, meshOut
    character(len=19) :: fieldNodeIn, fieldNodeOut, pfchno
    integer(kind=8), pointer :: trav1(:) => null()
    aster_logical, pointer :: trav2(:) => null()
    character(len=24), pointer :: refeIn(:) => null()
    character(len=24), pointer :: refeOut(:) => null()
    integer(kind=8), pointer :: deeqIn(:) => null()
    integer(kind=8), pointer :: deeqOut(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    codret = 0
    fieldNodeIn = fieldNodeInZ
    fieldNodeOut = fieldNodeOutZ

! - Check REFE
    call vrrefe(fieldNodeIn, fieldNodeOut, iret)

! - Access to fields
    call jelira(fieldNodeIn//'.VALE', 'TYPE', cval=scalTypeIn)
    call jelira(fieldNodeIn//'.VALE', 'LONMAX', nbEquaIn)
    call jeveuo(fieldNodeIn//'.REFE', 'L', vk24=refeIn)
    call jeveuo(fieldNodeIn//'.VALE', 'L', jvValeIn)
    call jelira(fieldNodeOut//'.VALE', 'TYPE', cval=scalTypeOut)
    call jelira(fieldNodeOut//'.VALE', 'LONMAX', nbEquaOut)
    call jeveuo(fieldNodeOut//'.REFE', 'L', vk24=refeOut)
    call jeveuo(fieldNodeOut//'.VALE', 'E', jvValeOut)

! - Check meshes
    call dismoi("NOM_MAILLA", fieldNodeIn, "CHAM_NO", repk=meshIn)
    call dismoi("NOM_MAILLA", fieldNodeOut, "CHAM_NO", repk=meshOut)
    if (meshIn .ne. meshOut) then
        call utmess('F', 'FIELD0_2')
    end if

    if (iret .eq. 0) then
! ----- Same numbering: copy values (set zero on Lagrange multipliers for BC)
        call dismoi('NUME_EQUA', fieldNodeOut, 'CHAM_NO', repk=pfchno)
        call jeveuo(pfchno//'.DEEQ', 'L', vi=deeq)
        if (scalTypeIn .eq. scalTypeOut) then
            if (scalTypeIn .eq. 'R') then
                do ieq1 = 0, nbEquaIn-1
                    if (deeq(2*ieq1+2) .le. 0) then
                        zr(jvValeOut+ieq1) = 0.d0
                    else
                        zr(jvValeOut+ieq1) = zr(jvValeIn+ieq1)
                    end if
                end do
            else if (scalTypeIn .eq. 'C') then
                do ieq1 = 0, nbEquaIn-1
                    if (deeq(2*ieq1+2) .le. 0) then
                        zc(jvValeOut+ieq1) = dcmplx(0.d0, 0.d0)
                    else
                        zc(jvValeOut+ieq1) = zc(jvValeIn+ieq1)
                    end if
                end do
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            if (scalTypeIn .eq. 'R' .and. scalTypeOut .eq. 'C') then
                do ieq1 = 0, nbEquaIn-1
                    if (deeq(2*ieq1+2) .le. 0) then
                        zc(jvValeOut+ieq1) = dcmplx(0.d0, 0.d0)
                    else
                        zc(jvValeOut+ieq1) = zr(jvValeIn+ieq1)
                    end if
                end do
            else
                call utmess("F", "FIELD0_1")
            end if
        end if
    else
! ----- Different numbering: copy values (set zero on Lagrange multipliers for BC)
        call jeveuo(refeIn(2) (1:19)//'.DEEQ', 'L', vi=deeqIn)
        call jeveuo(refeOut(2) (1:19)//'.DEEQ', 'L', vi=deeqOut)
!
!
!     2.1 ON CHERCHE LE NUMERO DE CMP LE PLUS GRAND ET
!     LE NUMERO DE NOEUD LE PLUS GRAND DANS CH2->.DEEQ2 :
!     ---------------------------------------------------------------
!
        nnomx = 0
        ncpmx = 0
        do ieq2 = 1, nbEquaOut
            nnomx = max(nnomx, deeqOut(2*(ieq2-1)+1))
            ncpmx = max(ncpmx, deeqOut(2*(ieq2-1)+2))
        end do
!
!
!     2.2 ON REMPLIT UN OBJET DE TRAVAIL :
!     ------------------------------------
        AS_ALLOCATE(vi=trav1, size=nnomx*ncpmx)
        AS_ALLOCATE(vl=trav2, size=nbEquaOut)
        do ieq2 = 1, nbEquaOut
            nuno2 = deeqOut(2*(ieq2-1)+1)
            nucp2 = deeqOut(2*(ieq2-1)+2)
            if (nucp2 .gt. 0) trav1((nuno2-1)*ncpmx+nucp2) = ieq2
            trav2(ieq2) = .false.
        end do

!     2.3 ON RECOPIE LES VALEURS DE CH1 DANS CH2
        if (scalTypeIn .eq. scalTypeOut) then
            if (scalTypeIn .eq. 'R') then
                do ieq1 = 1, nbEquaIn
                    nuno1 = deeqIn(2*(ieq1-1)+1)
                    nucp1 = deeqIn(2*(ieq1-1)+2)
                    if ((nucp1 .gt. 0) .and. (nuno1 .le. nnomx) .and. (nucp1 .le. ncpmx)) then
                        ieq2 = trav1((nuno1-1)*ncpmx+nucp1)
                        if (ieq2 .gt. 0) then
                            trav2(ieq2) = .true.
                            zr(jvValeOut-1+ieq2) = zr(jvValeIn-1+ieq1)
                        end if
                    end if
                end do
            else if (scalTypeIn .eq. 'C') then
                do ieq1 = 1, nbEquaIn
                    nuno1 = deeqIn(2*(ieq1-1)+1)
                    nucp1 = deeqIn(2*(ieq1-1)+2)
                    if ((nucp1 .gt. 0) .and. (nuno1 .le. nnomx)) then
                        ieq2 = trav1((nuno1-1)*ncpmx+nucp1)
                        if (ieq2 .gt. 0) then
                            trav2(ieq2) = .true.
                            zc(jvValeOut-1+ieq2) = zc(jvValeIn-1+ieq1)
                        end if
                    end if
                end do
            else
                ASSERT(ASTER_FALSE)
            end if
!
        else if (scalTypeIn .eq. 'R' .and. scalTypeOut .eq. 'C') then
            do ieq1 = 1, nbEquaIn
                nuno1 = deeqIn(2*(ieq1-1)+1)
                nucp1 = deeqIn(2*(ieq1-1)+2)
                if ((nucp1 .gt. 0) .and. (nuno1 .le. nnomx)) then
                    ieq2 = trav1((nuno1-1)*ncpmx+nucp1)
                    if (ieq2 .gt. 0) then
                        trav2(ieq2) = .true.
                        zc(jvValeOut-1+ieq2) = zr(jvValeIn-1+ieq1)
                    end if
                end if
            end do
        else
            call utmess("F", "FIELD0_1")
        end if
!
!     A CAUSE DE LA SOUS-STRUCTURATION STATIQUE, ON DOIT AJOUTER
!     UNE GLUTE POUR OUBLIER LA COMPOSANTE 'LAGR' POUR LA VERIF
        call dismoi('NOM_GD', fieldNodeOut, 'CHAM_NO', repk=nomgd)
        call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', jcmpgd)
        call jelira(jexnom('&CATA.GD.NOMCMP', nomgd), 'LONMAX', ncmpmx)
        cmpLagr = -200
        cmpLagr = indik8(zk8(jcmpgd), 'LAGR', 1, ncmpmx)
!
        do ieq2 = 1, nbEquaOut
            nuno2 = deeqOut(2*(ieq2-1)+1)
            nucp2 = deeqOut(2*(ieq2-1)+2)
!       NUCP2.NE.ICMP == GLUTE POUR LA SOUS-STRUCTURATION STATIQUE
            if (nucp2 .gt. 0 .and. nucp2 .ne. cmpLagr .and. .not. trav2(ieq2)) then
                codret = 1
            end if
        end do
        AS_DEALLOCATE(vi=trav1)
        AS_DEALLOCATE(vl=trav2)
    end if
!
    call jedema()
end subroutine
