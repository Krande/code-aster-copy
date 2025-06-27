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
subroutine caflux(load, model, mesh, geomDime, valeType)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/char_affe_neum.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/tbexp2.h"
#include "asterfort/tbliva.h"
#include "asterfort/tecart.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: load, mesh, model
    integer(kind=8), intent(in) :: geomDime
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads 'FLUX_REP'
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
! In  mesh             : mesh
! In  geomDime         : dimension of space
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'FLUX_REP'
    integer(kind=8) :: ibid, nflux, jvalv1, jvalv2, iocc, n, n1, n2, n3
    integer(kind=8) :: n4, n5, n6, n7, n8, n11, n12, ngr, ncmp, ncmp1, ncmps(2)
    integer(kind=8) :: ncmp2, iret
    real(kind=8) :: r8b, aire, xlong
    complex(kind=8) :: c16b
    aster_logical :: icre1, icre2
    character(len=8) :: k8b, nomtab
    character(len=19) :: cart1, cart2, cartes(2)
    character(len=24) :: para, mongrm
    character(len=24) :: valk(2)
    character(len=8), pointer :: vncmp1(:) => null()
    character(len=8), pointer :: vncmp2(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    icre1 = .false.
    icre2 = .false.
    call getfac(keywordFact, nflux)
!
    do iocc = 1, nflux
        n5 = 0
        if (valeType .eq. 'REEL') then
            call getvr8(keywordFact, 'FLUN', iocc=iocc, nbval=0, nbret=n11)
            call getvr8(keywordFact, 'FLUN_INF', iocc=iocc, nbval=0, nbret=n2)
            call getvr8(keywordFact, 'FLUN_SUP', iocc=iocc, nbval=0, nbret=n3)
            call getvid(keywordFact, 'CARA_TORSION', iocc=iocc, nbval=0, nbret=n12)
            n1 = n11+n12
        else if (valeType .eq. 'FONC') then
            call getvid(keywordFact, 'FLUN', iocc=iocc, nbval=0, nbret=n1)
            call getvid(keywordFact, 'FLUN_INF', iocc=iocc, nbval=0, nbret=n2)
            call getvid(keywordFact, 'FLUN_SUP', iocc=iocc, nbval=0, nbret=n3)
            call getvid(keywordFact, 'FLUX_X', iocc=iocc, nbval=0, nbret=n6)
            call getvid(keywordFact, 'FLUX_Y', iocc=iocc, nbval=0, nbret=n7)
            call getvid(keywordFact, 'FLUX_Z', iocc=iocc, nbval=0, nbret=n8)
            n5 = n6+n7+n8
        else
            ASSERT(.false.)
        end if
        n4 = n1+n2+n3
        if ((n5 .ne. 0) .and. (n4 .ne. 0)) then
            if (valeType .eq. 'FONC') then
                call utmess('F', 'MODELISA2_64')
            end if
        end if
        if (n4 .ne. 0) icre1 = .true.
        if (n5 .ne. 0) icre2 = .true.
    end do
!
!     ALLOCATION EVENTUELLE DES CARTES CART1 ET CART2 :
!
    cart1 = load//'.CHTH.FLURE'
    cart2 = load//'.CHTH.FLUR2'
    if (valeType .eq. 'REEL') then
        if (icre1) call alcart('G', cart1, mesh, 'FLUN_R')
        if (icre2) call alcart('G', cart2, mesh, 'FLUX_R')
    else if (valeType .eq. 'FONC') then
        if (icre1) call alcart('G', cart1, mesh, 'FLUN_F')
        if (icre2) call alcart('G', cart2, mesh, 'FLUX_F')
    else
        ASSERT(.false.)
    end if
!
    if (icre1) then
        call jeveuo(cart1//'.NCMP', 'E', vk8=vncmp1)
        call jeveuo(cart1//'.VALV', 'E', jvalv1)
    end if
    if (icre2) then
        call jeveuo(cart2//'.NCMP', 'E', vk8=vncmp2)
        call jeveuo(cart2//'.VALV', 'E', jvalv2)
    end if
!
!      STOCKAGE DE FLUX NULS SUR TOUT LE MAILLAGE
!
    if (icre1) then
        ncmp = 3
        vncmp1(1) = 'FLUN'
        vncmp1(2) = 'FLUN_INF'
        vncmp1(3) = 'FLUN_SUP'
        if (valeType .eq. 'REEL') then
            zr(jvalv1-1+1) = 0.d0
            zr(jvalv1-1+2) = 0.d0
            zr(jvalv1-1+3) = 0.d0
        else
            zk8(jvalv1-1+1) = '&FOZERO'
            zk8(jvalv1-1+2) = '&FOZERO'
            zk8(jvalv1-1+3) = '&FOZERO'
        end if
        call nocart(cart1, 1, ncmp)
    end if
!
    if (icre2) then
        ncmp = 3
        vncmp2(1) = 'FLUX'
        vncmp2(2) = 'FLUY'
        vncmp2(3) = 'FLUZ'
        if (valeType .eq. 'REEL') then
            zr(jvalv2-1+1) = 0.d0
            zr(jvalv2-1+2) = 0.d0
            zr(jvalv2-1+3) = 0.d0
        else
            zk8(jvalv2-1+1) = '&FOZERO'
            zk8(jvalv2-1+2) = '&FOZERO'
            zk8(jvalv2-1+3) = '&FOZERO'
        end if
        call nocart(cart2, 1, ncmp)
    end if
!
!     STOCKAGE DANS LES CARTES
!
    do iocc = 1, nflux
        ncmp1 = 0
        ncmp2 = 0
!
        if (valeType .eq. 'REEL') then
!
            call getvid(keywordFact, 'CARA_TORSION', iocc=iocc, scal=nomtab, nbret=n)
            if (n .eq. 1) then
!              VERIFICATION DES PARAMETRES DE LA TABLE 'NOMTAB'
                call tbexp2(nomtab, 'AIRE')
                call tbexp2(nomtab, 'LONGUEUR')
                call tbexp2(nomtab, 'GROUP_MA')
!
                call getvem(mesh, 'GROUP_MA', keywordFact, 'GROUP_MA', iocc, &
                            1, mongrm, ngr)
                para = 'AIRE'
                call tbliva(nomtab, 1, 'GROUP_MA', [ibid], [r8b], &
                            [c16b], mongrm, k8b, [r8b], para, &
                            k8b, ibid, aire, c16b, k8b, &
                            iret)
                if (iret .eq. 1) then
                    valk(1) = para
                    valk(2) = nomtab
                    call utmess('F', 'MODELISA8_34', nk=2, valk=valk)
                else if (iret .eq. 2) then
                    valk(1) = para
                    call utmess('F', 'MODELISA8_35', sk=valk(1))
                else if (iret .eq. 3) then
                    valk(1) = mongrm
                    call utmess('F', 'MODELISA8_36', sk=valk(1))
                end if
                para = 'LONGUEUR'
                call tbliva(nomtab, 1, 'GROUP_MA', [ibid], [r8b], &
                            [c16b], mongrm, k8b, [r8b], para, &
                            k8b, ibid, xlong, c16b, k8b, &
                            iret)
                if (iret .eq. 1) then
                    valk(1) = para
                    valk(2) = nomtab
                    call utmess('F', 'MODELISA8_34', nk=2, valk=valk)
                else if (iret .eq. 2) then
                    valk(1) = para
                    call utmess('F', 'MODELISA8_35', sk=valk(1))
                else if (iret .eq. 3) then
                    valk(1) = mongrm
                    call utmess('F', 'MODELISA8_36', sk=valk(1))
                end if
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'FLUN'
                zr(jvalv1-1+ncmp1) = 2.0d0*aire/xlong
            end if
            call getvr8(keywordFact, 'FLUN', iocc=iocc, scal=r8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'FLUN'
                zr(jvalv1-1+ncmp1) = r8b
            end if
            call getvr8(keywordFact, 'FLUN_INF', iocc=iocc, scal=r8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'FLUN_INF'
                zr(jvalv1-1+ncmp1) = r8b
            end if
            call getvr8(keywordFact, 'FLUN_SUP', iocc=iocc, scal=r8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'FLUN_SUP'
                zr(jvalv1-1+ncmp1) = r8b
            end if
!
        else
            call getvid(keywordFact, 'FLUN', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'FLUN'
                zk8(jvalv1-1+ncmp1) = k8b
            end if
            call getvid(keywordFact, 'FLUN_INF', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'FLUN_INF'
                zk8(jvalv1-1+ncmp1) = k8b
            end if
            call getvid(keywordFact, 'FLUN_SUP', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp1 = ncmp1+1
                vncmp1(ncmp1) = 'FLUN_SUP'
                zk8(jvalv1-1+ncmp1) = k8b
            end if
!
            call getvid(keywordFact, 'FLUX_X', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp2 = ncmp2+1
                vncmp2(ncmp2) = 'FLUX'
                zk8(jvalv2-1+ncmp2) = k8b
            end if
            call getvid(keywordFact, 'FLUX_Y', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp2 = ncmp2+1
                vncmp2(ncmp2) = 'FLUY'
                zk8(jvalv2-1+ncmp2) = k8b
            end if
            call getvid(keywordFact, 'FLUX_Z', iocc=iocc, scal=k8b, nbret=n)
            if (n .eq. 1) then
                ncmp2 = ncmp2+1
                vncmp2(ncmp2) = 'FLUZ'
                zk8(jvalv2-1+ncmp2) = k8b
            end if
        end if
!
!
        cartes(1) = cart1
        cartes(2) = cart2
        ncmps(1) = ncmp1
        ncmps(2) = ncmp2
        call char_affe_neum(model, mesh, geomDime, keywordFact, iocc, 2, &
                            cartes, ncmps)
!
    end do
!
    if (icre1) call tecart(cart1)
    if (icre2) call tecart(cart2)
!
    call jedema()
end subroutine
