! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine nxcerr(sddisc)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmecrr.h"
#include "asterfort/nxdocn.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: sddisc
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear thermics - Discretization management
!
! Management of parameters - Create
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc          : datastructure for time discretization
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: newtKrylResi = 0.9d0
    real(kind=8) :: resiGlobRela, resiGlobMaxi
    real(kind=8), dimension(2) :: parcrr
    integer(kind=8), dimension(3) :: parcri
    integer(kind=8) :: nbIter, iterGlobMaxi
    character(len=24) :: sddiscIfcvName
    real(kind=8), pointer :: sddiscIfcv(:) => null()
    character(len=24) :: sddiscIfreName
    real(kind=8), pointer :: sddiscIfre(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - LECTURE DES PARAMETRES DU FICHIER DE COMMANDE
    call nxdocn(parcri, parcrr)
    iterGlobMaxi = parcri(3)
    resiGlobRela = parcrr(2)
    resiGlobMaxi = parcrr(1)
    nbIter = iterGlobMaxi

! - CREATION DU VECTEUR D'INFORMATIONS SUR LA CONVERGENCE
    sddiscIfcvName = sddisc(1:19)//'.IFCV'
    call wkvect(sddiscIfcvName, 'V V R', 10, vr=sddiscIfcv)

! - Set parameters
    call nmecrr(sddisc, 'MXITER', paraValeI_=iterGlobMaxi)
    call nmecrr(sddisc, 'NBITER', paraValeI_=nbIter)
    call nmecrr(sddisc, 'RESI_GLOB_RELA', paraValeR_=resiGlobRela)
    call nmecrr(sddisc, 'RESI_GLOB_MAXI', paraValeR_=resiGlobMaxi)
    call nmecrr(sddisc, 'INIT_NEWTON_KRYLOV', paraValeR_=newtKrylResi)
    call nmecrr(sddisc, 'ITER_NEWTON_KRYLOV', paraValeR_=newtKrylResi)

! - CREATION DU VECTEUR DE STOCKAGE DES RESIDUS
    nbIter = nbIter+1
    sddiscIfreName = sddisc(1:19)//'.IFRE'
    call wkvect(sddiscIfreName, 'V V R', 3*nbIter, vr=sddiscIfre)
!
    call jedema()
end subroutine
