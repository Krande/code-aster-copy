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
subroutine cavini(ndim, nno, geom, vim, npg, &
                  lgpg, imate)
!
! CAVINI :
! CALCUL DES CONTRAINTES DE RUPTURE POUR MODELE ENDO_HETEROGENE
! VIM(3,GG) = CONTRAINTE D AMORCAGE AU PT DE GAUSS GG
! VIM(4,GG) = CONTRAINTE DE PROPAGATION AU PT DE GAUSS GG
!
    implicit none
#include "asterc/getran.h"
#include "asterc/r8pi.h"
#include "asterfort/casurf.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndim, nno, npg, lgpg, imate, zz, zzz, zzzz, nono, nitert, ntirmx
    real(kind=8) :: geom(1:ndim, 1:nno)
    real(kind=8) :: vim(1:lgpg, 1:npg), gr, pi
    real(kind=8) :: lc(1), mm, echp, ki, epai, ct1, ct2, randd, surff
    integer(kind=8) :: icodre(5)
    integer(kind=8) :: k2(1), kpg, spt
    character(len=16) :: nomres(5)
    character(len=8) :: fami, poum
    real(kind=8) :: valres(5), sa, sp, sc
!
!
    nitert = 0
    pi = r8pi()
567 continue
!
    nono = 0
! RMQ NICO : INITIALISATION DE LA CONTRAINTE D AMORCAGE
! SI NON PRECISEE
!
    call casurf(ndim, nno, geom, surff)
    nomres(1) = 'SY'
    nomres(2) = 'WEIBULL'
    nomres(3) = 'KI'
    nomres(4) = 'EPAI'
    nomres(5) = 'GR'
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, imate, &
                ' ', 'ENDO_HETEROGENE', 0, ' ', [0.d0], &
                5, nomres, valres, icodre, 1)
    call rcvalb(fami, kpg, spt, poum, imate, &
                ' ', 'NON_LOCAL', 0, ' ', [0.d0], &
                1, 'LONG_CARA', lc, k2, 1)
!  FACTEUR D ECHELLE
    echp = valres(1)
!  MODULE DE WEIBULL
    mm = valres(2)
! EPAISSEUR
!
    ki = valres(3)
    epai = valres(4)
! GRAINE
    gr = valres(5)
!
    if (gr .gt. 0.d0) then
!
! GRAINE NON NULLE : TIRAGE UNIQUE
! ON NE VERIFIE LE SEUIL QU'UNE SEULE FOIS
        ntirmx = 1
    else
        ntirmx = 25
    end if
!
    if (vim(3, 1) .lt. 0.0001d0) then
        call getran(randd)
        ct1 = 0.d0
        ct1 = 0.d0-log(1.d0-randd)
        sa = 0.d0
        sa = echp*((lc(1)**3.d0)**(1.d0/mm))/((surff*epai)**(1.d0/mm))*( &
             ct1**(1.d0/mm))
        do zz = 1, npg
            vim(3, zz) = sa
        end do
    end if
!  INITIALISATION DE LA CONTRAINTE DE PROPAGATION
! SI NON PRECISEE
    if (vim(4, 1) .lt. 0.0001d0) then
!
! TENACITE
!
        ct2 = 0.d0
        ct2 = 0.5736d0
        sp = 0.d0
        sp = ct2*((ki**2.d0/(pi*lc(1)))**(0.5d0))
        do zzz = 1, npg
            vim(4, zzz) = sp
        end do
    end if
!
!  VERIFICATION DE LA COHERENCE DES DEUX SEUILS
    sc = ((2.d0)**(0.5d0))*vim(4, 1)
    if (sc .gt. vim(3, 1)) then
        do zzzz = 1, npg
            vim(3, zzzz) = 0.d0
            vim(4, zzzz) = 0.d0
        end do
        nono = 1
        nitert = nitert+1
    end if
!
!
! ON N AUTORISE NTIRMX TIRAGES SINON L'UTILISATEUR DOIT
! REVOIR L'INITIALISATION DE SES SEUILS
    if ((nono .eq. 1) .and. (nitert .lt. ntirmx)) then
        goto 567
    else if ((nono .eq. 1) .and. (nitert .ge. ntirmx)) then
        call utmess('F', 'COMPOR2_14')
    end if
!
end subroutine
