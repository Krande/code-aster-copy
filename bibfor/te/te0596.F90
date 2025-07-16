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
subroutine te0596(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/niinit.h"
#include "asterfort/nofnpd.h"
#include "asterfort/nufnlg.h"
#include "asterfort/nufnpd.h"
#include "asterfort/teattr.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
! FONCTION REALISEE:  CALCUL DE L'OPTION FORC_NODA POUR LES ELEMENTS
!                     INCOMPRESSIBLES A 2 CHAMPS UP
!                     EN 3D/D_PLAN/AXI
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ----------------------------------------------------------------------
!
    aster_logical :: mini
    integer(kind=8) :: ndim, nno1, nno2, nnos, npg, jgn, ntrou
    integer(kind=8) :: iw, ivf1, ivf2, idf1, idf2
    integer(kind=8) :: vu(3, 27), vg(27), vp(27), vpi(3, 27)
    integer(kind=8) :: igeom, jvSief, jvDisp, imate, ivectu
    integer(kind=8) :: ibid
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: lielrf(10), typmod(2), alias8
    character(len=24) :: valk
! ----------------------------------------------------------------------
!
! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elref2(nomte, 10, lielrf, ntrou)
    ASSERT(ntrou .ge. 2)
    call elrefe_info(elrefe=lielrf(2), fami='RIGI', ndim=ndim, nno=nno2, nnos=nnos, &
                     npg=npg, jpoids=iw, jvf=ivf2, jdfde=idf2, jgano=jgn)
    call elrefe_info(elrefe=lielrf(1), fami='RIGI', ndim=ndim, nno=nno1, nnos=nnos, &
                     npg=npg, jpoids=iw, jvf=ivf1, jdfde=idf1, jgano=jgn)
!
! - TYPE DE MODELISATION
    if (ndim .eq. 2 .and. lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS  '
    else if (ndim .eq. 2 .and. lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN  '
    else if (ndim .eq. 3) then
        typmod(1) = '3D'
    else
        call utmess('F', 'ELEMENTS_34', sk=nomte)
    end if
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PSIEFR', 'L', jvSief)
    call jevech('PDEPLAR', 'L', jvDisp)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PVECTUR', 'E', ivectu)
!
! - CALCUL DES FORCES INTERIEURES
    if (compor(DEFO) .eq. 'PETIT ') then
        if (lteatt('INCO', 'C2 ')) then
!
! - MINI ELEMENT ?
            call teattr('S', 'ALIAS8', alias8, ibid)
            if (alias8(6:8) .eq. 'TR3' .or. alias8(6:8) .eq. 'TE4') then
                mini = .true.
            else
                mini = .false.
            end if
! --------- Get index of dof
            call niinit(typmod, ndim, nno1, 0, &
                        nno2, 0, vu, vg, vp, &
                        vpi)
!
            call nufnpd(ndim, nno1, nno2, npg, iw, &
                        zr(ivf1), zr(ivf2), idf1, vu, vp, &
                        typmod, zi(imate), compor, zr(igeom), zr(jvSief), &
                        zr(jvDisp), mini, zr(ivectu))
        else if (lteatt('INCO', 'C2O')) then
! --------- Get index of dof
            call niinit(typmod, ndim, nno1, 0, &
                        nno2, nno2, vu, vg, vp, &
                        vpi)
!
            call nofnpd(ndim, nno1, nno2, nno2, npg, &
                        iw, zr(ivf1), zr(ivf2), zr(ivf2), idf1, &
                        vu, vp, vpi, typmod, zi(imate), &
                        compor, zr(igeom), nomte, zr(jvSief), zr(jvDisp), &
                        zr(ivectu))
        else
            valk = compor(DEFO)
            call utmess('F', 'MODELISA10_17', sk=valk)
        end if
    else if (compor(DEFO) .eq. 'GDEF_LOG') then
        if (lteatt('INCO', 'C2 ')) then
!
! - MINI ELEMENT ?
            call teattr('S', 'ALIAS8', alias8, ibid)
            if (alias8(6:8) .eq. 'TR3' .or. alias8(6:8) .eq. 'TE4') then
! - PAS ENCORE INTRODUIT
                valk = compor(DEFO)
                call utmess('F', 'MODELISA10_18', sk=valk)
            end if
! --------- Get index of dof
            call niinit(typmod, ndim, nno1, 0, &
                        nno2, 0, vu, vg, vp, &
                        vpi)
!
            call nufnlg(ndim, nno1, nno2, npg, iw, &
                        zr(ivf1), zr(ivf2), idf1, vu, vp, &
                        typmod, zi(imate), compor, zr(igeom), zr(jvSief), &
                        zr(jvDisp), zr(ivectu))
        else
            valk = compor(DEFO)
            call utmess('F', 'MODELISA10_17', sk=valk)
        end if
    else
        call utmess('F', 'ELEMENTS3_16', sk=compor(DEFO))
    end if
!
end subroutine
