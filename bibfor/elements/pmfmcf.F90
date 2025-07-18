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

subroutine pmfmcf(ip, nbgf, nbfib, nugf, sdcomp, &
                  crit, option, instam, instap, icdmat, &
                  nbvalc, defam, defap, varim, varimp, &
                  contm, defm, defp, epsm, modf, &
                  sigf, varip, codret)
!
! aslint: disable=W1504
! --------------------------------------------------------------------------------------------------
!
!       APPEL AU COMPORTEMENT DU GROUPE DE FIBRE
!
! --------------------------------------------------------------------------------------------------
!
!
    implicit none

#include "MultiFiber_type.h"
#include "asterfort/pmfcom.h"
!
    integer(kind=8) :: ip, nbgf, nbfib, nbvalc, nugf(*), icdmat, codret
    character(len=16) :: option
    character(len=24) :: sdcomp(*)
    real(kind=8) :: varim(*), varimp(*), varip(*), contm(*), defm(*), defp(*)
    real(kind=8) :: crit(*), instam, instap, defap(*), defam(*), epsm
    real(kind=8) :: sigf(*), modf(*)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ig, ngf, nbfig
    integer(kind=8) :: idcipv, idcipc, idecc, idecv, icp
    integer(kind=8) :: iposv, iposc
    integer(kind=8) :: codrep
!
! --------------------------------------------------------------------------------------------------
!
    codrep = 0
!
    idcipc = nbfib*(ip-1)
    idcipv = nbvalc*idcipc
    idecc = 1
    idecv = 1
    do ig = 1, nbgf
!       numéro du groupe de fibre
        ngf = nugf(ig)
        icp = 1+(ngf-1)*MULTI_FIBER_SIZEK
!       nombre de fibres de ce groupe
        read (sdcomp(icp-1+MULTI_FIBER_NBFI), '(I24)') nbfig
!       aiguillage suivant comportement :
!           module et contrainte sur chaque fibre
!           attention à la position du pointeur contrainte et variables internes
        iposv = idecv+idcipv
        iposc = idecc+idcipc
        call pmfcom(ip, idecc, option, sdcomp(icp), crit, &
                    nbfig, instam, instap, icdmat, nbvalc, &
                    defam, defap, varim(iposv), varimp(iposv), contm(iposc), &
                    defm(idecc), defp(idecc), epsm, modf(idecc), sigf(idecc), &
                    varip(iposv), codrep)
        if (codrep .ne. 0) then
            codret = codrep
!           code 3: on continue et on le renvoie a la fin. Autres codes: sortie immediate
            if (codrep .ne. 3) goto 900
        end if
        idecc = idecc+nbfig
        idecv = idecv+nbvalc*nbfig
    end do
!
900 continue
end subroutine
