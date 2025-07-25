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
subroutine fcweib(nrupt, cals, sk, sigw, nur, &
                  nt, nbres, indtp, nbtp, m, &
                  fc, dfc)
    implicit none
#include "asterf_types.h"
#include "asterc/r8maem.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nrupt, nur(*), nt(*), nbres, indtp(*), nbtp
    real(kind=8) :: sigw(*), m, fc, dfc, s1, s2, sk(*)
    aster_logical :: cals
!     AUTEUR : M. BONNAMY
!     ----------------------------------------------------------------
!
!     BUT: CALCUL DE RECALAGE DES PARAMETRES DE WEIBULL PAR LA
!          METHODE DU MAXIMUM DE VRAISSEMBLANCE
!
!     ----------------------------------------------------------------
!
!     NRUPT        /IN/:NOMBRE DE CONTRAINTES
!     CALS         /IN/:TRUE SI SIGMA_U EST FIXE
!     SK           /IN/:PARAMETRE SIGMA-U(K) DE WEIBULL
!     SIGW         /IN/:CONTRAINTES DE WEIBULL AUX INSTANTS DE RUPTURE
!     NUR          /IN/:NUMERO DE RESULTAT ASSOCIEE A
!                       LA CONTRAINTE SIGW(I)
!     NT           /IN/:DIMENSION DE LA SOUS-BASE CORRESPONDANT A LA
!                       TEMPERATURE T
!     NBRES        /IN/:NOMBRE DE BASES DE RESULTATS
!     INDTP        /IN/:INDICE DE TEMPERATURE POUR CHAQUE RESULTAT
!     NBTP         /IN/:NOMBRE DE TEMPERATURE DIFFERENTES
!
!     M            /OUT/:PARAMETRE M(K+1)DE WEIBULL
!     FC           /OUT/:FONCTION F(M) DE WEIBULL
!     DFC          /OUT/:DERIVEE DE LA FONCTION DF(M) DE WEIBULL
!
!
!     ----------------------------------------------------------------
!
    real(kind=8) :: swm, slw, slwm, sl2wm, sl2bwm, snt, maxr, maxm
    real(kind=8) :: valr
    integer(kind=8) :: i, itp, ir, vali
    character(len=4) :: valk
!
!     ----------------------------------------------------------------
!
    slw = 0.d0
    sl2bwm = 0.d0
    maxr = r8maem()
    maxm = log(maxr)/log(nrupt*sigw(nrupt))
    if (m .ge. maxm) then
        valr = maxm
        call utmess('F', 'UTILITAI8_22', sr=valr)
    end if
    if (m .le. 0.d0) then
        valr = m
        call utmess('F', 'UTILITAI8_23', sr=valr)
    end if
!
    do i = 1, nrupt
!
        if (sigw(i) .eq. 0.d0) then
            vali = i
            valk = 'SIGW'
            call utmess('F', 'UTILITAI8_25', si=vali, sk=valk)
        end if
        if (cals) then
            slw = slw+(log(sigw(i)/sk(1)))*(1.d0-(sigw(i)/sk(1))**m)
            sl2bwm = sl2bwm+( &
                     (sigw(i)/sk(1))**m)*(log(sigw(i)/sk(1))*log(sigw(i)/sk(1)))
        else
            slw = slw+log(sigw(i))
        end if
!
    end do
!
    s1 = 0.d0
    s2 = 0.d0
    do itp = 1, nbtp
!
        snt = 0.d0
!
        do ir = 1, nbres
!
            if (indtp(ir) .eq. itp) snt = snt+nt(ir)
!
        end do
!
        swm = 0.d0
        slwm = 0.d0
        sl2wm = 0.d0
        do i = 1, nrupt
!
            if (indtp(nur(i)) .eq. itp) then
                swm = swm+sigw(i)**m
                slwm = slwm+(sigw(i)**m)*(log(sigw(i)))
                sl2wm = sl2wm+(sigw(i)**m)*(log(sigw(i))*log(sigw(i)))
            end if
!
        end do
!
        s1 = s1+snt*slwm/swm
        s2 = s2+snt*((sl2wm/swm)*swm-(slwm/swm)*slwm)/swm
!
    end do
!
    if (cals) then
        fc = nrupt
        fc = fc/m+slw
        dfc = nrupt
        dfc = -dfc*(1.d0/(m*m))-sl2bwm
    else
        fc = nrupt
        fc = fc/m+slw-s1
        dfc = nrupt
        dfc = -dfc*(1.d0/(m*m))-s2
    end if
!
!
end subroutine
