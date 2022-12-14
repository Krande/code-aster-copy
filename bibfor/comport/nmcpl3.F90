! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine nmcpl3(compor, option, crit, deps, dsidep, &
                  ndim,   sigp,   vip,  cpl,  icp,    &
                  conv)
!     CONTRAINTES PLANES PAR LA METHODE DE BORST / CONDENSATION STATIQUE
!     POUR LES COMPORTEMENTS QUI N'INTEGRENT PAS LES CONTRAINTES PLANES
!     ATTENTION : POUR BIEN CONVERGER, IL FAUT REACTUALISER LA MATRICE
!     TANGENTE. DE PLUS, IL FAUT AJOUTER 4 VARIABLES INTERNES
!
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
!                               (3) = VALEUR TOLERANCE DE CONVERGENCE
!                                     (RESI_INTE_RELA == RESCREL)
!     DEPS    : INCREMENT DE DEFORMATION TOTALE :
!               DEPS(T) = DEPS(MECANIQUE(T)) + DEPS(DILATATION(T))
!     ICP     : NUMERO DE L'ITERATION
! VAR DSIDEP  : MATRICE TANGENTE CARREE
! VAR SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! VAR VIP     : LES 4 DERNIERES SONT RELATIVES A LA METHODE DE BORST
!
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16) :: option, compor(*)
    integer :: k, ndim, ncpmax, icp, cpl
    aster_logical :: conv, vecteu
    real(kind=8) :: vip(*), deps(*), crit(*), dsidep(6, 6), sigp(4), sigpeq
    real(kind=8) :: prec, signul, precr, ddezz
!
    vecteu = (option(1:9) .eq. 'FULL_MECA') .or. (option(1:9) .eq. 'RAPH_MECA')
    ncpmax = nint(crit(ITER_DEBORST_MAX))
!
    signul = crit(RESI_INTE_RELA)
    prec   = crit(RESI_DEBORST_MAX)
    conv   = .true.
!
    if (vecteu) then
!       DANS LE CAS D=1 ON NE FAIT RIEN CAR LES CONTRAINTES SONT NULLES
        if (compor(RELA_NAME) .eq. 'ENDO_ISOT_BETON') then
            if (vip(2) .gt. 1.5d0) goto 999
        endif
!
        if (prec .gt. 0.d0) then
            ! PRECISION RELATIVE
            sigpeq=0.d0
            do k = 1, 2*ndim
                sigpeq = sigpeq + sigp(k)**2
            enddo
            sigpeq = sqrt(sigpeq)
            if (sigpeq .lt. signul) then
                precr=prec
            else
                precr=prec*sigpeq
            endif
        else
            ! PRECISION ABSOLUE
            precr=abs(prec)
        endif
        conv = (icp.ge.ncpmax) .or. (abs(sigp(3)).lt.precr)
        !
        if (.not. conv) then
            if (cpl .eq. 2) then
                if (abs(dsidep(3,3)) .gt. precr) then
                    deps(3) = deps(3) - sigp(3)/dsidep(3,3)
                else
                    conv=.true.
                endif
            else if (cpl .eq. 1) then
                if (abs(dsidep(2,2)) .gt. precr) then
                    ddezz   = -(sigp(3) - dsidep(3,2)/dsidep(2,2)*sigp(2))/ &
                               (dsidep(3,3)- dsidep(3,2)*dsidep(2,3)/dsidep( 2,2))
                    deps(3) = deps(3) + ddezz
                    deps(2) = deps(2)-(sigp(2)+dsidep(2,3)*ddezz)/dsidep(2,2)
                else
                    conv=.true.
                endif
            endif
        endif
    endif
!
999 continue
end subroutine
