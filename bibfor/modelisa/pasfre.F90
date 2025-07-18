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
subroutine pasfre(disc, freq, pasf, dim, nbm, &
                  iv, imodi, freqi, freqf, nb)
    implicit none
!  CALCUL DE LA DISCRETISATION FREQUENTIELLE
!  -----------------------------------------
!  POUR CHAQUE FREQUENCE MODALE FREQ(I) PRISE EN COMPTE
!  DANS LA BANDE DE FREQUENCES, UNE ZONE CALCULEE
!  EN FONCTION DE FREQ(I-1) ET FREQ(I+1) EST DISCRETISEE
!  EN 1024 POINTS (NB = 1024 , EN PARAMETER DANS OP0146)
!  CHAQUE ZONE EST DIVISEE EN 4
!  2 ZONES A PAS FIN ET 2 ZONES A PAS LACHE DETERMINEES PAR :
!   FREQ_DEBUT = (FREQ(I)+FREQ(I-1))/2
!   F2         = FREQ(I) - DF
!   F3         = FREQ(I)
!   F4         = FREQ(I) + DF
!   FREQ_FIN   = (FREQ(I)+FREQ(I+1))/2
!       AVEC DF = 2*PI*FREQ(I)*AMOR(I)
!-----------------------------------------------------------------------
!     IN  : FREQ  : CARACT. MODALES DE LA BASE DE CONCEPT MELASFLU
!     OUT : PASF  : DISCRETISATION FREQUENTIELLE CALCULEE
!     IN  : DIM   : NOMBRE DE MODES PRIS EN COMPTE DANS LA BANDE DE FREQ
!     IN  : NBM   : NOMBRE DE MODES DSE LA BASE DE CONCEPT MELASFLU
!     IN  : IV    : INDICATEUR DE LA VITESSE ETUDIEE
!     IN  : IMODI : NUMERO D ORDRE DU PREMIER MODE PRIS EN COMPTE
!     IN  : FREQI : FREQUENCE INITIALE DE LA BANDE DE FREQ.
!     IN  : FREQF : FREQUENCE FINALE   DE LA BANDE DE FREQ.
!     IN  : NB    : NOMBRE DE POINTS PAR MODE POUR LA DISCR. FREQ.
!     IN  : DISC  : TABLEAU DE TRAVAIL
!
!
#include "jeveux.h"
#include "asterc/r8pi.h"
    integer(kind=8) :: dim, nb, nbm
    real(kind=8) :: freq(2, nbm, *), pasf(dim*nb), disc(2, *)
!-----------------------------------------------------------------------
!C
!-----------------------------------------------------------------------
    integer(kind=8) :: if, im, imodi, ip, iv, iz
    integer(kind=8) :: nb4, nbpf, nbz, numo
    real(kind=8) :: df, dff, fmax, fmin, freqf, freqi, pas
    real(kind=8) :: pi
!-----------------------------------------------------------------------
    nb4 = nb/4
    pi = r8pi()
    nbpf = nb*dim
    nbz = 4*dim
!
    disc(1, 1) = freqi
    disc(2, nbz) = freqf
!
    do im = 1, dim
        numo = imodi+(im-1)
!
        if (freq(2, numo, iv) .lt. 0.d0) then
            df = 20.d0*freq(1, numo, iv)*1.d-06
        else
            df = 2.d0*pi*freq(1, numo, iv)*freq(2, numo, iv)
        end if
!
        if (im .gt. 1) then
            ip = im-1
            disc(1, ip*4+1) = disc(2, ip*4)
        end if
!
        if (im .lt. dim) then
            disc(2, (im-1)*4+4) = (freq(1, numo, iv)+freq(1, numo+1, iv))/2.d0
        end if
!
        dff = freq(1, numo, iv)-disc(1, (im-1)*4+1)
        if (df .ge. dff) then
            disc(2, (im-1)*4+1) = freq(1, numo, iv)-(dff/2.d0)
            disc(1, (im-1)*4+2) = freq(1, numo, iv)-(dff/2.d0)
        else
            disc(2, (im-1)*4+1) = freq(1, numo, iv)-df
            disc(1, (im-1)*4+2) = freq(1, numo, iv)-df
        end if
!
        disc(2, (im-1)*4+2) = freq(1, numo, iv)
        disc(1, (im-1)*4+3) = freq(1, numo, iv)
!
        dff = disc(2, (im-1)*4+4)-freq(1, numo, iv)
        if (df .ge. dff) then
            disc(2, (im-1)*4+3) = freq(1, numo, iv)+(dff/2.d0)
            disc(1, (im-1)*4+4) = freq(1, numo, iv)+(dff/2.d0)
        else
            disc(2, (im-1)*4+3) = freq(1, numo, iv)+df
            disc(1, (im-1)*4+4) = freq(1, numo, iv)+df
        end if
    end do
!
    if = 1
    do iz = 1, nbz
        if (iz .lt. nbz) then
            fmin = disc(1, iz)
            fmax = disc(2, iz)
            pas = (fmax-fmin)/dble(nb4)
            pasf(if) = fmin
            do ip = 1, nb4-1
                if = if+1
                pasf(if) = fmin+pas*ip
            end do
            if = if+1
        else
            fmin = disc(1, iz)
            fmax = disc(2, iz)
            pas = (fmax-fmin)/dble(nb4-1)
            pasf(if) = fmin
            do ip = 1, nb4-2
                if = if+1
                pasf(if) = fmin+pas*ip
            end do
        end if
    end do
    pasf(nbpf) = fmax
end subroutine
