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
subroutine vedime(modelZ, loadNameJvZ, loadInfoJvZ, &
                  timeCurr, scalarType, vectElemZ, &
                  lCumul_, jvBase_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/load_list_info.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/vemare.h"
#include "asterfort/reajre.h"
!
    character(len=*), intent(in) :: modelZ
    character(len=*), intent(in) :: loadNameJvZ, loadInfoJvZ
    real(kind=8), intent(in) :: timeCurr
    character(len=1), intent(in) :: scalarType
    character(len=24), intent(inout) :: vectElemZ
    aster_logical, optional, intent(in) :: lCumul_
    character(len=1), optional, intent(in) :: jvBase_
!
! --------------------------------------------------------------------------------------------------
!
! Compute Dirichlet loads
!
! For Lagrange elements (AFFE_CHAR_MECA) - U(given)
!
! --------------------------------------------------------------------------------------------------
!
! In  model             : name of model
! In  loadNameJv        : name of object for list of loads name
! In  loadInfoJv        : name of object for list of loads info
! In  timeCurr          : current time
! In  scalarType        : type of coefficients (real or complex)
! IO  vectElem          : name of vectElem result
! In  jvBase            : JEVEUX base to create vector
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbin = 3, nbout = 1
    character(len=8) :: lpain(nbin), lpaout(nbout)
    character(len=24) :: lchin(nbin), lchout(nbout)
    character(len=8) :: loadName, newnom
    character(len=16) :: option
    character(len=24) :: vectElem, resuElem
    character(len=24), parameter :: chtime = '&&VEDIME.CH_INST_R'
    character(len=24) :: loadLigrel, chgeom
    integer(kind=8) :: loadNume, nbLoad, iLoad
    character(len=24), pointer :: listLoadName(:) => null()
    integer(kind=8), pointer :: listLoadInfo(:) => null()
    aster_logical :: noLoadInList, lCumul
    character(len=1) :: jvBase
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    newnom = '.0000000'
    resuElem = '&&VEDIME.???????'
    jvBase = 'V'
    if (present(jvBase_)) then
        jvBase = jvBase_
    end if
    lCumul = ASTER_FALSE
    if (present(lCumul_)) then
        lCumul = lCumul_
    end if
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "

! - Result name for vectElem
    vectElem = vectElemZ(1:19)
    if (vectElem .eq. ' ') then
        vectElem = '&&VEDIME'
    end if

! - Get loads
    call load_list_info(noLoadInList, nbLoad, listLoadName, listLoadInfo, &
                        loadNameJvZ, loadInfoJvZ)

! - Allocate result
    if (.not. lCumul) then
        call detrsd('VECT_ELEM', vectElem)
        call vemare(jvBase, vectElem, modelZ)
        call reajre(vectElem, ' ', jvBase)
    end if
    if (noLoadInList) then
        goto 99
    end if

! - Geometry field
    call megeom(modelZ, chgeom)

! - Time field
    call mecact('V', chtime, 'MODELE', modelZ, 'INST_R  ', &
                ncmp=1, nomcmp='INST', sr=timeCurr)

! - Input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PINSTR'
    lchin(2) = chtime

! - Output field
    if (scalarType .eq. 'R') then
        lpaout(1) = 'PVECTUR'
    else
        lpaout(1) = 'PVECTUC'
    end if

! - Computation
    do iLoad = 1, nbLoad
        loadName = listLoadName(iLoad) (1:8)
        loadNume = listLoadInfo(iLoad+1)
        if ((loadNume .gt. 0) .and. (loadNume .le. 4)) then
            loadLigrel = loadName//'.CHME.LIGRE'

! --------- Input field
            lchin(3) = loadName//'.CHME.CIMPO'
            if (loadNume .eq. 1) then
                if (scalarType .eq. 'R') then
                    option = 'MECA_DDLI_R'
                    lpain(3) = 'PDDLIMR'
                else
                    option = 'MECA_DDLI_C'
                    lpain(3) = 'PDDLIMC'
                end if
            else if (loadNume .eq. 2) then
                option = 'MECA_DDLI_F'
                lpain(3) = 'PDDLIMF'
            else if (loadNume .eq. 3) then
                option = 'MECA_DDLI_F'
                lpain(3) = 'PDDLIMF'
            else if (loadNume .eq. 4) then
                ASSERT(scalarType .eq. 'R')
                option = 'MECA_DDLI_R'
                lpain(3) = 'PDDLIMR'
            else
                ASSERT(ASTER_FALSE)
            end if

! --------- Generate new RESU_ELEM name
            call gcnco2(newnom)
            resuElem(10:16) = newnom(2:8)
            call corich('E', resuElem, ichin_=iLoad)
            lchout(1) = resuElem

! --------- Computation
            call calcul('S', option, loadLigrel, &
                        nbin, lchin, lpain, &
                        nbout, lchout, lpaout, &
                        jvBase, 'OUI')

! --------- Copying output field
            call reajre(vectElem, lchout(1), jvBase)
        end if
    end do
!
99  continue
!
    vectElemZ = vectElem(1:19)//'.RELR'
!
    call jedema()
end subroutine
