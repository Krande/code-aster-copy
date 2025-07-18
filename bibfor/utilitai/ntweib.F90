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
subroutine ntweib(nrupt, cals, sk, sigw, nur, &
                  nt, nbres, x1, x2, xacc, &
                  rtsafe, impr, ifm, indtp, nbtp)
    implicit none
#include "asterf_types.h"
#include "asterfort/fcweib.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nrupt, nur(*), nt(*), nbres, indtp(*), nbtp, ifm
    real(kind=8) :: sigw(*), x1, x2, xacc, rtsafe, sk(*)
    aster_logical :: cals, impr
!     AUTEUR : M. BONNAMY
!     ----------------------------------------------------------------
!
!     BUT: CALCUL DE RECALAGE DES PARAMETRES DE WEIBULL PAR LA
!          METHODE DU MAXIMUM DE VRAISSEMBLANCE : METHODE DE NEWTON
!          COMBINEE AVEC METHODE DE LA BISSECTRICE
!
!     ----------------------------------------------------------------
!
!     NRUPT        /IN/:NOMBRE TOTAL DE CONTRAINTES DE WEIBULL
!     CALS         /IN/:TRUE SI SIGMA_U EST FIXE
!     SK           /IN/:PARAMETRE SIGMA-U(K) DE WEIBULL
!     SIGW         /IN/:CONTRAINTES DE WEIBULL AUX INSTANTS DE RUPTURE
!     NUR          /IN/:NUMERO DE RESULTAT ASSOCIEE A
!                       LA CONTRAINTE SIGW(I)
!     NT           /IN/:DIMENSION DE LA SOUS-BASE CORRESPONDANT A LA
!                       TEMPERATURE T
!     NBRES        /IN/:NOMBRE DE BASES DE RESULTATS
!     IMPR         /IN/:IMPRESSION DETAILLEE
!     INDTP        /IN/:INDICE DE TEMPERATURE POUR CHAQUE RESULTAT
!     NBTP         /IN/:NOMBRE DE TEMPERATURE DIFFERENTES
!     X1           /IN/:BORNE GAUCHE DE RECHERCHE
!     X2           /IN/:BORNE DROITE DE RECHERCHE
!     XACC         /IN/:PRECISION VOULUE
!
!     X,Y          /OUT/:VALEUR DES FONCTIONS LOG(SIGMAW)
!                        ET LOG(LOG(1/(1-PF)))
!     MKP          /OUT/:PARAMETRE M(K+1)DE WEIBULL
!     SKP          /OUT/:PARAMETRE SIGMA-U(K+1) DE WEIBULL
!     RTSAFE       /OUT/:PARAMETRE M(K+1)DE WEIBULL
!
!     ----------------------------------------------------------------
!
    real(kind=8) :: df, dx, dxold, f, fh, fl, temp, xh, xl, dfl, dfh
    real(kind=8) :: valr(4)
    integer(kind=8) :: maxit, j
    integer(kind=8) :: vali
    parameter(maxit=100)
!     ----------------------------------------------------------------
!
4   continue
    call fcweib(nrupt, cals, sk, sigw, nur, &
                nt, nbres, indtp, nbtp, x1, &
                fl, dfl)
    if (impr) write (ifm, *) 'F,DF,X1 SUR BORNE GAUCHE : ', fl, dfl, x1
5   continue
    call fcweib(nrupt, cals, sk, sigw, nur, &
                nt, nbres, indtp, nbtp, x2, &
                fh, dfh)
    if (impr) write (ifm, *) 'F,DF,X2 SUR BORNE DROITE : ', fh, dfh, x2
!
!     RECHERCHE DE LA BORNE DROITE SI F(X1) ET F(X2) DE MEME SIGNE
!
    if (((fl .gt. 0.d0 .and. fh .gt. 0.d0) .and. dfh .lt. 0.d0) .or. &
        ((fl .lt. 0.d0 .and. fh .lt. 0.d0) .and. dfh .gt. 0.d0)) then
        x2 = x2+0.9d0
        goto 5
    end if
    if (((fl .gt. 0.d0 .and. fh .gt. 0.d0) .and. dfl .gt. 0.d0) .or. &
        ((fl .lt. 0.d0 .and. fh .lt. 0.d0) .and. dfl .lt. 0.d0)) then
        x1 = x1-0.9d0
        goto 4
    end if
    if (fl .eq. 0.d0) then
        rtsafe = x1
        goto 999
    else if (fh .eq. 0.d0) then
        rtsafe = x2
        goto 999
    else if (fl .lt. 0.d0) then
        xl = x1
        xh = x2
    else
        xh = x1
        xl = x2
    end if
    rtsafe = 0.5d0*(x1+x2)
    dxold = abs(x2-x1)
    dx = dxold
    call fcweib(nrupt, cals, sk, sigw, nur, &
                nt, nbres, indtp, nbtp, rtsafe, &
                f, df)
    if (impr) write (ifm, *) 'F ET DF MILIEU INTERVALLE :', rtsafe, f, df
    do j = 1, maxit
        if (impr) write (ifm, *) '*** ITERATION DE NEWTON NO', j
        if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) .gt. 0.d0 .or. abs(2.d0*f) .gt. &
            abs(dxold*df)) then
            dxold = dx
            dx = 0.5d0*(xh-xl)
            rtsafe = xl+dx
            if (xl .eq. rtsafe) goto 999
            if (impr) write (ifm, *) 'INCREMENT - SOLUTION : ', dx, rtsafe
        else
            dxold = dx
            dx = f/df
            temp = rtsafe
            rtsafe = rtsafe-dx
            if (temp .eq. rtsafe) goto 999
            if (impr) write (ifm, *) 'INCREMENT - SOLUTION : ', dx, rtsafe
        end if
        if (abs(dx) .lt. xacc) goto 999
        call fcweib(nrupt, cals, sk, sigw, nur, &
                    nt, nbres, indtp, nbtp, rtsafe, &
                    f, df)
        if (impr) write (ifm, *) 'SOLUTION/F-DF : ', rtsafe, f, df
        if (f .lt. 0.d0) then
            xl = rtsafe
        else
            xh = rtsafe
        end if
    end do
    call utmess('F', 'UTILITAI2_53')
999 continue
    if (impr) then
        valr(1) = rtsafe
        valr(2) = f
        valr(3) = dx
        valr(4) = xacc
        vali = j
        call utmess('I', 'UTILITAI6_48', si=vali, nr=4, valr=valr)
    end if
!
end subroutine
