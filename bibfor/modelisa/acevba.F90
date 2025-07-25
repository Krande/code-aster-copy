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
subroutine acevba(nbocc, nlg, ier)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8maem.h"
#include "asterfort/acedat.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/codent.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbocc, nlg, ier
!     AFFE_CARA_ELEM
!     VERIFICATION DES MOTS CLES POUR L'ELEMENT BARRE
! ----------------------------------------------------------------------
! IN  : NBOCC  : NOMBRE D'OCCURENCE
! OUT : NLG    : NOMBRE TOTAL DE GROUPE DE MAILLE
! ----------------------------------------------------------------------
!     NSECBA : NOMBRE DE SECTIONS PAR BARRE
!     NTYPSE : NOMBRE DE TYPE DE SECTION
! ----------------------------------------------------------------------
    real(kind=8) :: tst
    character(len=8) :: k8b, kioc, ki, nomu
    character(len=24) :: valk(3)
    character(len=16) :: k16b, sec, concep, cmd
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ioc, irece, irech
    integer(kind=8) :: l, nbcar
    integer(kind=8) :: nbo, nbval, nc, ncar, ncara, ncmax, ndim
    integer(kind=8) :: ng, ns, nsec, nsecba, ntypse
    integer(kind=8) :: nv, nval
    character(len=8), pointer :: cara(:) => null()
    character(len=8), pointer :: carbar(:) => null()
    character(len=8), pointer :: expbar(:) => null()
    integer(kind=8), pointer :: ncp(:) => null()
    character(len=8), pointer :: tabbar(:) => null()
    integer(kind=8), pointer :: tab_para(:) => null()
    character(len=16), pointer :: typ_sect(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call getres(nomu, concep, cmd)
!
    AS_ALLOCATE(vi=tab_para, size=10)
    call acedat('BARRE', 0, tab_para, k16b, k8b, &
                k8b, k8b)
    nsecba = tab_para(1)
    ntypse = tab_para(1+1)
    nbo = tab_para(1+2)
    nbcar = tab_para(1+3)
    nbval = tab_para(1+4)
    AS_ALLOCATE(vi=ncp, size=ntypse)
    do i = 1, ntypse
        ncp(i) = tab_para(1+4+i)
    end do
    ndim = ncp(1+1)*ntypse
    AS_ALLOCATE(vk16=typ_sect, size=ntypse)
    AS_ALLOCATE(vk8=expbar, size=nbo)
    AS_ALLOCATE(vk8=tabbar, size=nbo)
    AS_ALLOCATE(vk8=carbar, size=ndim)
    call acedat('BARRE', 1, tab_para, typ_sect, expbar, &
                tabbar, carbar)
    AS_ALLOCATE(vk8=cara, size=nbcar)
    AS_ALLOCATE(vr=vale, size=nbval)
!
    tst = r8maem()
    nlg = 0
    do ioc = 1, nbocc
        call codent(ioc, 'G', kioc)
        call getvtx('BARRE', 'GROUP_MA', iocc=ioc, nbval=0, nbret=ng)
        call getvtx('BARRE', 'SECTION', iocc=ioc, nbval=0, nbret=ns)
        call getvtx('BARRE', 'SECTION', iocc=ioc, scal=sec, nbret=nsec)
        call getvtx('BARRE', 'CARA', iocc=ioc, nbval=0, nbret=nc)
        call getvtx('BARRE', 'CARA', iocc=ioc, nbval=nbcar, vect=cara, &
                    nbret=ncar)
        call getvr8('BARRE', 'VALE', iocc=ioc, nbval=0, nbret=nv)
        call getvr8('BARRE', 'VALE', iocc=ioc, nbval=nbval, vect=vale, &
                    nbret=nval)
!
! -- CARA
        if (ncar .gt. 0) then
            ncara = ncar
            do l = 1, ntypse
                if (sec .eq. typ_sect(l)) then
                    ncmax = ncp(l)*nsecba
                    call codent(ncmax, 'G', ki)
                    if (ncar .gt. ncmax .and. l .ne. 2) then
                        valk(1) = kioc
                        valk(2) = ki
                        valk(3) = typ_sect(l)
                        call utmess('E', 'MODELISA_44', nk=3, valk=valk)
                        ier = ier+1
                    end if
                    if (l .eq. 2) then
                        if (ncar .gt. 4) then
                            valk(1) = kioc
                            valk(2) = typ_sect(l)
                            call utmess('E', 'MODELISA_45', nk=2, valk=valk)
                            ier = ier+1
                        end if
                        irech = 0
                        irece = 0
                        do i = 1, ncar
                            if (cara(i) (1:2) .eq. 'H ') then
                                if (irech .eq. 2) then
                                    valk(1) = kioc
                                    valk(2) = typ_sect(l)
                                    call utmess('E', 'MODELISA_46', nk=2, valk=valk)
                                    ier = ier+1
                                end if
                                irech = 1
                            end if
                            if (cara(i) (1:2) .eq. 'HY' .or. cara(i) (1:2) .eq. 'HZ') then
                                if (irech .eq. 1) then
                                    valk(1) = kioc
                                    valk(2) = typ_sect(l)
                                    call utmess('E', 'MODELISA_47', nk=2, valk=valk)
                                    ier = ier+1
                                end if
                                irech = 2
                            end if
                            if (cara(i) (1:3) .eq. 'EP ') then
                                if (irece .eq. 1) then
                                    valk(1) = kioc
                                    valk(2) = typ_sect(l)
                                    call utmess('E', 'MODELISA_48', nk=2, valk=valk)
                                    ier = ier+1
                                end if
                                irece = 2
                            end if
                            if (cara(i) (1:3) .eq. 'EPX' .or. cara(i) (1:3) .eq. 'EPY') then
                                if (irece .eq. 2) then
                                    valk(1) = kioc
                                    valk(2) = typ_sect(l)
                                    call utmess('E', 'MODELISA_49', nk=2, valk=valk)
                                    ier = ier+1
                                end if
                                irece = 1
                            end if
                        end do
                    end if
                end if
            end do
        end if
!
! -- VALE
        if (nval .gt. 0) then
            if (nval .ne. ncara) then
                call codent(ncara, 'G', ki)
                valk(1) = kioc
                valk(2) = ki
                call utmess('E', 'MODELISA_50', nk=2, valk=valk)
                ier = ier+1
            else
                do i = 1, nval
                    call codent(i, 'G', ki)
                    if (vale(i) .eq. tst) then
                        valk(1) = kioc
                        valk(2) = typ_sect(l)
                        valk(3) = ki
                        call utmess('E', 'MODELISA_51', nk=3, valk=valk)
                        ier = ier+1
                    end if
                end do
            end if
        end if
!
        nlg = max(nlg, -ng)
!
    end do
!
    AS_DEALLOCATE(vi=tab_para)
    AS_DEALLOCATE(vi=ncp)
    AS_DEALLOCATE(vk16=typ_sect)
    AS_DEALLOCATE(vk8=expbar)
    AS_DEALLOCATE(vk8=tabbar)
    AS_DEALLOCATE(vk8=carbar)
    AS_DEALLOCATE(vk8=cara)
    AS_DEALLOCATE(vr=vale)
!
    call jedema()
end subroutine
