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
subroutine compt(nbpt, fn, offset, t, elapse, &
                 nbchoc, tchocm, tchmax, tchmin, nbrebo, &
                 trebom, tchoct, nbinst)
!        COMPTAGE DES CHOCS AMV
!
! IN  : NBPT   : NB DE POINTS DU SIGNAL
! IN  : FN     : TABLEAU DU SIGNAL
! IN  : T      : TABLEAU DU TEMPS
! IN  : OFFSET : VALEUR DU SEUIL DE DETECTION D UN CHOC
! IN  : ELAPSE : TEMPS MINIMUM POUR VRAI FIN DE CHOC
! OUT : NBCHOC : NB DE CHOC GLOBAUX ( CRITERE ELAPSE )
! OUT : NBREBO : NB DE REBONDS ( RETOUR AU SEUIL )
! OUT : TCHOCM : TEMPS DE CHOC GLOBAL MOYEN
! OUT : TREBOM : TEMPS DE REBOND MOYEN
! OUT : TCHOCT : TEMPS DE CHOC CUMULE
! IN  : NBINST : NB D'INSTANTS TOTAL DU CALCUL TRANSITOIRE
! ----------------------------------------------------------------------
!
    implicit none
    real(kind=8) :: fn(*), t(*)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ichoc, idebur, idebut, idech, ifin, ifinr
    integer(kind=8) :: irebo, j, jfin, nbchoc, nbinst, nbpas, nbpt
    integer(kind=8) :: nbrebo
    real(kind=8) :: dt, elapse, offset, tchmax, tchmin, tchoc, tchocm
    real(kind=8) :: tchoct, trebo, trebom, zero
!-----------------------------------------------------------------------
    zero = 0.d0
    nbchoc = 0
    nbrebo = 0
    tchocm = zero
    trebom = zero
    tchoct = zero
    tchmax = zero
    tchmin = 1.0d20
    irebo = 0
    ichoc = 0
    idebut = 1
    idebur = 1
    ifin = 1
    dt = t(4)-t(3)
    nbpas = max(1, nint(elapse/dt))
!
    do i = 1, nbpt
!
        if (abs(fn(i)) .le. offset) then
!
            if (irebo .eq. 1) then
                ifinr = i
                trebo = t(ifinr)-t(idebur)
                trebom = trebom+trebo
                nbrebo = nbrebo+1
            end if
!
            idech = 0
            jfin = min(i+nbpas, nbinst)
            if (jfin .gt. (i+1)) then
                do j = i+1, jfin
                    if (abs(fn(j)) .gt. offset) idech = 1
                end do
            end if
!
            if (idech .eq. 0 .and. ichoc .eq. 1) then
!
                ifin = i
                tchoc = t(ifin)-t(idebut)
                tchocm = tchocm+tchoc
!
                if (tchoc .gt. tchmax) tchmax = tchoc
!
                if (tchoc .lt. tchmin) tchmin = tchoc
!
                nbchoc = nbchoc+1
                ichoc = 0
!
            end if
!
            irebo = 0
!
        else
!
            if (ichoc .eq. 0) idebut = i
!
            if (irebo .eq. 0) idebur = i
            irebo = 1
            ichoc = 1
!
        end if
!
    end do
!
    tchoct = tchocm
    if (nbchoc .ne. 0) then
        tchocm = tchocm/nbchoc
    else
        tchocm = zero
    end if
!
    if (nbrebo .ne. 0) then
        trebom = trebom/nbrebo
    else
        trebom = zero
    end if
!
end subroutine
