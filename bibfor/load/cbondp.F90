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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine cbondp(load, mesh, ndim, valeType)
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/char_read_val.h"
#include "asterfort/getelem.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/utmess.h"
#include "asterfort/vetyma.h"
!
    character(len=8), intent(in) :: load, mesh
    integer(kind=8), intent(in) :: ndim
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load ONDE_PLANE
!
! --------------------------------------------------------------------------------------------------
!
!
! In  mesh      : name of mesh
! In  load      : name of load
! In  ndim      : space dimension
! In  valeType  : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'ONDE_PLANE'
    character(len=24), parameter :: listCell = '&&CBONDP.LIST_ELEM'
    complex(kind=8) :: c16dummy
    real(kind=8) :: r8dummy
    character(len=8) :: k8dummy
    character(len=16) :: k16dummy
    real(kind=8) :: wave_dire(3), point_source(3), point_reflechi(3), wave_type_r
    character(len=8) :: signal, signde
    character(len=16) :: wave_type
    integer(kind=8) :: jvalv
    integer(kind=8) :: iocc, ndir, val_nb, nondp
    integer(kind=8) :: jvCell
    integer(kind=8) :: nbCell
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer(kind=8) :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    wave_dire = 0.d0
    point_source = 0.d0
    point_reflechi = 0.d0
!
    call getfac(keywordfact, nondp)
    if (nondp .eq. 0) goto 99
!
! - Initializations
!
    ASSERT(valeType .eq. 'FONC')
!
! - Creation and initialization to zero of <CARTE>
!
    call char_crea_cart('MECANIQUE', keywordfact, load, mesh, valeType, &
                        nbMap, map, nbCmp)
    ASSERT(nbMap .eq. 2)
!
! - Loop on factor keyword
!
    do iocc = 1, nondp

! ----- Read mesh affectation
        call getelem(mesh, keywordfact, iocc, 'A', listCell, nbCell)
        if (nbCell .ne. 0) then
            call jeveuo(listCell, 'L', jvCell)

! --------- Get wave function
            call char_read_val(keywordfact, iocc, 'FONC_SIGNAL', 'FONC', val_nb, &
                               r8dummy, signal, c16dummy, k16dummy)
            ASSERT(val_nb .eq. 1)
            call char_read_val(keywordfact, iocc, 'DEPL_IMPO', 'FONC', val_nb, &
                               r8dummy, signde, c16dummy, k16dummy)
            if (val_nb .ne. 1) then
                signde = '&FOZERO'
            end if

! --------- Affectation of values in <CARTE> - Wave function
            call jeveuo(map(1)//'.VALV', 'E', jvalv)
            zk8(jvalv-1+1) = signal
            zk8(jvalv-1+2) = signde
            call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCell, &
                        limanu=zi(jvCell))

! --------- Get direction
            wave_dire(1) = 0.d0
            wave_dire(2) = 0.d0
            wave_dire(3) = 0.d0
            call getvr8(keywordfact, 'DIRECTION', iocc=iocc, nbval=0, nbret=ndir)
            ndir = -ndir
            ASSERT(ndir .eq. 3)
            call getvr8(keywordfact, 'DIRECTION', iocc=iocc, nbval=ndir, vect=wave_dire)

! --------- Get wave type
            call char_read_val(keywordfact, iocc, 'TYPE_ONDE', 'TEXT', val_nb, &
                               r8dummy, k8dummy, c16dummy, wave_type)
            ASSERT(val_nb .eq. 1)
            if (ndim .eq. 3) then
                if (wave_type .eq. 'P ') then
                    wave_type_r = 0.d0
                else if (wave_type .eq. 'SV') then
                    wave_type_r = 1.d0
                else if (wave_type .eq. 'SH') then
                    wave_type_r = 2.d0
                else if (wave_type .eq. 'S ') then
                    call utmess('F', 'CHARGES2_61')
                else
                    ASSERT(ASTER_FALSE)
                end if
            else if (ndim .eq. 2) then
                if (wave_type .eq. 'P ') then
                    wave_type_r = 0.d0
                else if (wave_type .eq. 'S ') then
                    wave_type_r = 1.d0
                else if (wave_type .eq. 'SV' .or. wave_type .eq. 'SH') then
                    call utmess('A', 'CHARGES2_62')
                    wave_type_r = 1.d0
                else
                    ASSERT(ASTER_FALSE)
                end if
            else
                ASSERT(ASTER_FALSE)
            end if

! --------- Affectation of values in <CARTE> - Wave type and direction
            call jeveuo(map(2)//'.VALV', 'E', jvalv)
            zr(jvalv-1+1) = wave_dire(1)
            zr(jvalv-1+2) = wave_dire(2)
            zr(jvalv-1+3) = wave_dire(3)
            zr(jvalv-1+4) = wave_type_r

! --------- Get source point coordinates
            zr(jvalv-1+5) = r8vide()
            zr(jvalv-1+6) = r8vide()
            zr(jvalv-1+7) = r8vide()
            call getvr8(keywordfact, 'COOR_SOURCE', iocc=iocc, nbval=0, nbret=ndir)
            ndir = -ndir
            if (ndir .ne. 0) then
                call getvr8(keywordfact, 'COOR_SOURCE', iocc=iocc,&
                            &nbval=ndir, vect=point_source)
                zr(jvalv-1+5) = point_source(1)
                zr(jvalv-1+6) = point_source(2)
                if (ndim .eq. 3) then
                    zr(jvalv-1+7) = point_source(3)
                end if
            end if

! --------- Get reflection point coordinates
            zr(jvalv-1+8) = r8vide()
            zr(jvalv-1+9) = r8vide()
            zr(jvalv-1+10) = r8vide()
            call getvr8(keywordfact, 'COOR_REFLECHI', iocc=iocc, nbval=0, nbret=ndir)
            ndir = -ndir
            if (ndir .ne. 0) then
                call getvr8(keywordfact, 'COOR_REFLECHI', iocc=iocc,&
                            &nbval=ndir, vect=point_reflechi)
                zr(jvalv-1+8) = point_reflechi(1)
                zr(jvalv-1+9) = point_reflechi(2)
                if (ndim .eq. 3) then
                    zr(jvalv-1+10) = point_reflechi(3)
                end if
            end if

            call nocart(map(2), 3, nbCmp(2), mode='NUM', nma=nbCell, &
                        limanu=zi(jvCell))
        end if

! ----- Check elements
        call vetyma(mesh, ndim, keywordfact, listCell, nbCell)
!
        call jedetr(listCell)
!
    end do
!
99  continue
    call jedema()
end subroutine
