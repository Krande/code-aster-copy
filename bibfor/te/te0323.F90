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

subroutine te0323(option, nomte)
!
! person_in_charge: jerome.laverne at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
!
#include "asterfort/ejfono.h"
#include "asterfort/ejfore.h"
#include "asterfort/ejinit.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/terefe.h"
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
!    OPTION FORC_NODA ET REFE_FORC_NODA POUR LES JOINTS QUADRA ET HYME
! ----------------------------------------------------------------------
!
    character(len=8) :: lielrf(10)
    aster_logical :: axi
    integer(kind=8) :: nno1, nno2, npg, ivf2, idf2, nnos, jgn, nddl
    integer(kind=8) :: iw, ivf1, idf1, igeom, ivectu, jvSief, ndim, ntrou
    integer(kind=8) :: iu(3, 16), ip(8)
    real(kind=8) :: sigref, fhyref
!
    call elref2(nomte, 2, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(1), fami='RIGI', ndim=ndim, nno=nno1, nnos=nnos, &
                     npg=npg, jpoids=iw, jvf=ivf1, jdfde=idf1, jgano=jgn)
    call elrefe_info(elrefe=lielrf(1), fami='RIGI', ndim=ndim, nno=nno2, nnos=nnos, &
                     npg=npg, jpoids=iw, jvf=ivf2, jdfde=idf2, jgano=jgn)
    ndim = ndim+1
    nddl = 2*ndim*nno1+nno2
    axi = lteatt('AXIS', 'OUI')
!
! - DECALAGE D'INDICE POUR LES ELEMENTS DE JOINT
    call ejinit(nomte, iu, ip)
!
    call jevech('PVECTUR', 'E', ivectu)
    call jevech('PGEOMER', 'L', igeom)
!
!      OPTIONS FORC_NODA ET REFE_FORC_NODA
!
    if (option .eq. 'FORC_NODA') then
!
        call jevech('PSIEFR', 'L', jvSief)
!
        call ejfono(ndim, nddl, axi, nno1, nno2, &
                    npg, iw, zr(iw), zr(ivf1), zr(ivf2), &
                    idf2, zr(idf2), zr(igeom), iu, ip, &
                    zr(jvSief), zr(ivectu))
!
    else if (option .eq. 'REFE_FORC_NODA') then
!
        call terefe('SIGM_REFE', 'THM_JOINT', sigref)
!
!      EN MECA PURE ON IMPOSE LA VALEUR DE FLUX DE REFERENCE A 1
        if (lteatt('TYPMOD2', 'EJ_HYME')) then
            call terefe('FLUX_HYD1_REFE', 'THM_JOINT', fhyref)
        else if (lteatt('TYPMOD2', 'ELEMJOIN')) then
            fhyref = 1.D0
        end if
!
        call ejfore(ndim, nddl, axi, nno1, nno2, &
                    npg, iw, zr(iw), zr(ivf1), zr(ivf2), &
                    idf2, zr(idf2), zr(igeom), iu, ip, &
                    sigref, fhyref, zr(ivectu))
    end if
!
!
end subroutine
