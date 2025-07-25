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
function lcroy1()
    implicit none
    real(kind=8) :: lcroy1
!
! ********************************************************************
! *       INTEGRATION DE LA LOI DE ROUSSELIER LOCAL                  *
! *  CALCUL DE BORNES INF ET SUP DE LA FONCTION SEUIL(Y) QUAND S(0)>0*
! *  ET RESOLUTION DE SEUIL(Y)=0                                     *
! *  PAR UNE METHODE DE NEWTON AVEC BORNES CONTROLEES ET DICHOTOMIE  *
! *  - LA BORNE INFERIEURE EST TELLE QUE:                            *
! *    2*MU*EQTR-R(PM+DPINF)=3*MU*DPINF      => ON EN DEDUIT DPINF   *
! *    YINF*EXP(YINF)=K*FONC*DPINF/SIG1      => ON EN DEDUIT YINF    *
! *  - LA BORNE SUPERIEURE EST TELLE QUE:                            *
! *    YSUP*EXP(YSUP) =K*FONC*(2*MU*EQTR-S(YINF))/3*MU*SIG1          *
! ********************************************************************
!
! OUT  Y   : VALEUR DE Y TEL QUE SEUIL(Y)=0
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!  COMMON LOI DE COMPORTEMENT ROUSSELIER
!
#include "asterfort/lcrofg.h"
#include "asterfort/lcroty.h"
#include "asterfort/rcfonc.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: itemax, jprolp, jvalep, nbvalp
    real(kind=8) :: prec, young, nu, sigy, sig1, rousd, f0, fcr, acce
    real(kind=8) :: pm, rpm, fonc, fcd, dfcddj, dpmaxi, typoro
    common/lcrou/prec, young, nu, sigy, sig1, rousd, f0, fcr, acce,&
     &               pm, rpm, fonc, fcd, dfcddj, dpmaxi, typoro,&
     &               itemax, jprolp, jvalep, nbvalp
! ----------------------------------------------------------------------
!  COMMON GRANDES DEFORMATIONS CANO-LORENTZ
!
    integer(kind=8) :: ind1(6), ind2(6)
    real(kind=8) :: kr(6), rac2, rc(6)
    real(kind=8) :: lambda, mu, deuxmu, unk, troisk, cother
    real(kind=8) :: jm, dj, jp, djdf(3, 3)
    real(kind=8) :: etr(6), dvetr(6), eqetr, tretr, detrdf(6, 3, 3)
    real(kind=8) :: dtaude(6, 6)
!
    common/gdclc/&
     &          ind1, ind2, kr, rac2, rc,&
     &          lambda, mu, deuxmu, unk, troisk, cother,&
     &          jm, dj, jp, djdf,&
     &          etr, dvetr, eqetr, tretr, detrdf,&
     &          dtaude
!
! ----------------------------------------------------------------------
    integer(kind=8) :: iter
    real(kind=8) :: seuil, dseuil, s
    real(kind=8) :: y, dp, yinf, ysup, t, rp, pente, aire
!
! 1 - CALCUL DU MINORANT
!
    if (deuxmu*eqetr-rpm .le. 0) then
!
! LE TERME DE TRACE N'EST PAS NEGLIGEABLE
!
        yinf = 0
    else
!
! RESOLUTION DE L'EQUATION SANS LE TERME DE TRACE
! LCROTY RESOUD UNE EQUATION DU TYPE Y*EXP(Y)=CONSTANTE
!
        call rcfonc('E', 1, jprolp, jvalep, nbvalp, &
                    e=young, nu=nu, p=pm, rp=rp, rprim=pente, &
                    airerp=aire, sieleq=2.d0*mu*eqetr, dp=dp)
        yinf = lcroty(dp*unk*fonc/sig1, prec, itemax)
    end if
!
! 2 - CALCUL DU MAJORANT
! LCROTY RESOUD UNE EQUATION DU TYPE Y*EXP(Y)=CONSTANTE
!
    call lcrofg(yinf, dp, s, seuil, dseuil)
    if (seuil .lt. 0.d0) yinf = 0.d0
    call lcrofg(yinf, dp, s, seuil, dseuil)
!
    t = (deuxmu*eqetr-s)/(3.d0*mu)
    ysup = lcroty(t*unk*fonc/sig1, prec, itemax)
    call lcrofg(ysup, dp, s, seuil, dseuil)
!
! 3 - RESOLUTION PAR UNE METHODE DE NEWTON ENTRE LES BORNES
!
    y = ysup
    do iter = 1, itemax
        if (abs(seuil)/sigy .le. prec) goto 100
!
        y = y-seuil/dseuil
        if (y .le. yinf .or. y .ge. ysup) y = (yinf+ysup)/2
!
        call lcrofg(y, dp, s, seuil, dseuil)
!
        if (seuil .ge. 0) yinf = y
        if (seuil .le. 0) ysup = y
!
    end do
    call utmess('F', 'ALGORITH3_55')
!
!
100 continue
    lcroy1 = y
end function
