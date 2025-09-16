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
subroutine nmdecp(sddisc, iterNewt, iEvenActi, cutStepType, nbStep, &
                  deltac, ratio, optdec, ldcext, durdec, &
                  retdec)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/nmdcae.h"
#include "asterfort/utdidt.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: iterNewt, iEvenActi
    character(len=4), intent(out) :: cutStepType
    integer(kind=8), intent(out) :: nbStep
    real(kind=8), intent(out) :: deltac, ratio
    character(len=16), intent(out) :: optdec
    aster_logical, intent(out) :: ldcext
    real(kind=8), intent(out) :: durdec
    integer(kind=8), intent(out) :: retdec
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (GESTION DES EVENEMENTS - DECOUPE)
!
! PARAMETRES DE DECOUPE
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  iterNewt         : index of current Newton iteration
! In  iEvenActi     : index of active event
! OUT RATIO  : RATIO DU PREMIER PAS DE TEMPS
! OUT TYPDEC : TYPE DE DECOUPE
!              'SUBD' - SUBDIVISION PAR UN NOMBRE DE PAS DONNE
!              'DELT' - SUBDIVISION PAR UN INCREMENT DONNE
! OUT NBRPAS : NOMBRE DE PAS DE TEMPS
! OUT DELTAC : INCREMENT DE TEMPS CIBLE
! OUT OPTDEC : OPTION DE DECOUPE
!     'UNIFORME'   - DECOUPE REGULIERE ET UNIFORME
!     'PROGRESSIF' - DECOUPE EN DEUX ZONES, UN PAS LONG+ UNE SERIE
!                    DE PAS UNIFORMES
!     'DEGRESSIF'  - DECOUPE EN DEUX ZONES, UNE SERIE DE PAS
!                    UNIFORMES + UN PAS LONG
! OUT RETDEC : CODE RETOUR DECOUPE
!     0 - ECHEC DE LA DECOUPE
!     1 - ON A DECOUPE
!     2 - PAS DE DECOUPE
! OUT LDCEXT : .TRUE. SI ON DOIT CONTINUER LA DECOUPE
! OUT DURDEC : DUREEE DE LA DECOUPE APRES (SI LDCEXT =.TRUE.)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: subaut
!
! --------------------------------------------------------------------------------------------------
!
    retdec = 0
    optdec = ' '
    ratio = 0.d0
    nbStep = -1
    deltac = -1.d0
    cutStepType = ' '
    optdec = ' '
    ldcext = ASTER_FALSE
    durdec = -1.d0

! - TYPE DE DECOUPAGE AUTO
    call utdidt('L', sddisc, 'ECHE', 'SUBD_METHODE_AUTO', index_=iEvenActi, &
                valk_=subaut)

! - PARAMETRES SUIVANT DECOUPE
    if (subaut .eq. 'EXTRAPOLE') then
        call nmdcae(sddisc, iterNewt, cutStepType, nbStep, ratio, &
                    optdec, retdec)
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
