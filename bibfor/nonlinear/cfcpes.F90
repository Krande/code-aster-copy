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

subroutine cfcpes(resoco, jsecmb)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/compute_ineq_conditions_vector.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=24) :: resoco
    integer(kind=8) :: jsecmb
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (RESOLUTION - PENALISATION)
!
! CALCUL DU SECOND MEMBRE POUR LE CONTACT -E_N.[Ac]T.{JEU}
!
! ----------------------------------------------------------------------
!
!
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  JSECMB : ADRESSE VERS LE SECOND MEMBRE
!
!
!
!
    integer(kind=8) :: iliac
    character(len=19) :: mu
    integer(kind=8) :: jmu
    character(len=24) :: apcoef, apddl, appoin
    integer(kind=8) :: japcoe, japddl, japptr
    character(len=24) :: tacfin
    integer(kind=8) :: jtacf
    character(len=24) :: jeux
    integer(kind=8) :: jjeux
    integer(kind=8) :: ztacf
    integer(kind=8) :: nbliai, neq
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATION DES VARIABLES
!
    nbliai = cfdisd(resoco, 'NBLIAI')
    neq = cfdisd(resoco, 'NEQ')
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    appoin = resoco(1:14)//'.APPOIN'
    apddl = resoco(1:14)//'.APDDL'
    apcoef = resoco(1:14)//'.APCOEF'
    jeux = resoco(1:14)//'.JEUX'
    tacfin = resoco(1:14)//'.TACFIN'
    mu = resoco(1:14)//'.MU'
    call jeveuo(appoin, 'L', japptr)
    call jeveuo(apddl, 'L', japddl)
    call jeveuo(apcoef, 'L', japcoe)
    call jeveuo(jeux, 'L', jjeux)
    call jeveuo(tacfin, 'L', jtacf)
    call jeveuo(mu, 'E', jmu)
    ztacf = cfmmvd('ZTACF')

    call compute_ineq_conditions_vector(jsecmb, nbliai, neq, &
                                        japptr, japddl, japcoe, &
                                        jjeux, jtacf, jmu, &
                                        3, ztacf, iliac)
!
    call jedema()
!
end subroutine
