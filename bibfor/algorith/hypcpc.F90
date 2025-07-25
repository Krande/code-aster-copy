! --------------------------------------------------------------------
! Copyright (C) 2005 UCBL LYON1 - T. BARANGER     WWW.CODE-ASTER.ORG
! Copyright (C) 2007 - 2025 - EDF R&D - www.code-aster.org
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

subroutine hypcpc(c11, c22, c33, c12, k, &
                  c10, c01, c20, nitmax, epsi, &
                  sig, codret)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
    real(kind=8) :: c11, c22, c33, c12
    real(kind=8) :: k
    real(kind=8) :: c10, c01, c20
    integer(kind=8) :: nitmax
    real(kind=8) :: epsi
    real(kind=8) :: sig(6)
    integer(kind=8) :: codret
!
! ----------------------------------------------------------------------
!
! LOI DE COMPORTEMENT HYPERELASTIQUE DE SIGNORINI
!
! C_PLAN - CALCUL DES CONTRAINTES
!
! ----------------------------------------------------------------------
!
!
! IN  C11,C22,C33,C12: ELONGATIONS
! IN  C10,C01,C20:     CARACTERISTIQUES MATERIAUX
! IN  K      : MODULE DE COMPRESSIBILITE
! IN  NITMAX : NOMBRE MAXI D'ITERATIONS
! IN  EPSI   : CRITERE DE CONVERGENCE
! OUT SIG    : CONTRAINTES
! OUT CODRET : CODE RETOUR ERREUR INTEGRATION (1 SI PROBLEME, 0 SINON)
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: t1, t3, t5, t6, t8, t12, t13
    real(kind=8) :: t17, t20, t15, t24, t61
    real(kind=8) :: t19, t38, t7, t37, t58, t9, t41, t46
    real(kind=8) :: sn, dc33, dsn, test
    integer(kind=8) :: iter, facmul
!
!-----------------------------------------------------------------------
!
    facmul = 10
    iter = 1
10  continue
!
! --- DEBUT DE BOUCLE
!
    t1 = c11*c22
    t3 = c12**2
    t5 = t1*c33-t3*c33
    if (t5 .le. 0.d0) then
        codret = 1
        goto 99
    end if
    t6 = t5**(1.d0/3.d0)
    if (t6 .eq. 0.d0) then
        codret = 1
        goto 99
    end if
    t7 = 1.d0/t6
    t8 = c11+c22+c33
    t12 = t1-t3
    t15 = t7-t8/t6/t5*t12/3.d0
    t19 = t6**2
    if (t19 .eq. 0.d0) then
        codret = 1
        goto 99
    end if
    t38 = sqrt(t5)
    sn = 2.d0*c10*t15+2.d0*c01*( &
         (c11+c22)/t19-2.d0/3.d0*(t1+c11*c33+c22*c33-t3)/t19/t5*t12)+4.d0*c20*(t8*t7-3.d0)*t15+&
         &k*(t38-1.d0 &
         )/t38*t12
    t1 = c11*c22
    t3 = c12**2
    t5 = t1*c33-t3*c33
    t6 = t5**(1.d0/3.d0)
    t8 = 1.d0/t6/t5
    t9 = t1-t3
    t12 = c11+c22+c33
    t13 = t5**2
    t17 = t9**2
    t20 = -2.d0/3.d0*t8*t9+4.d0/9.d0*t12/t6/t13*t17
    t24 = t6**2
    t41 = 1.d0/t6
    t46 = (t41-t12*t8*t9/3.d0)**2.d0
    t58 = sqrt(t5)
    t61 = t58**2
    dsn = 2.d0*c10*t20+2.d0*c01*( &
          -4.d0/3.d0*(c11+c22)/t24/t5*t9+10.d0/9.d0*(t1+c11*c33+c22*c33-t3)/t24/t13*t17)+4.d0*c&
          &20*t46+4.d0*c20*(t41*t12-3.d0)*t20+k/t5*t17/2.d0-k*(t58-1.d0 &
          )/t61/t58*t17/2.d0
    if (dsn .eq. 0.d0) then
        codret = 1
        goto 99
    else
        dc33 = -sn/dsn
    end if
!
    c33 = c33+dc33
!
    iter = iter+1
    if (iter .lt. nitmax*facmul) then
        test = abs(sn)
        if (test .lt. epsi) then
            codret = 0
            goto 200
        else
            goto 10
        end if
    else
        codret = 1
        goto 99
    end if
!
! --- FIN DE BOUCLE
!
200 continue
!
    t1 = c11*c22
    t3 = c12**2
    t5 = t1*c33-t3*c33
    t6 = t5**(1.d0/3.d0)
!
    if ((t5 .le. 0.d0)) then
        codret = 1
        goto 99
    end if
!
    t7 = 1.d0/t6
    t8 = c11+c22+c33
    t12 = c22*c33
    t15 = t7-t8/t6/t5*t12/3.d0
    t19 = t6**2
    t37 = sqrt(t5)
    sig(1) = 2.d0*c10*t15+2.d0*c01*( &
             (c22+c33)/t19-2.d0/3.d0*(t1+c11*c33+t12-t3)/t19/t5*t12)+4.d0*c20*(t8*t7-3.d0)*t15+&
             & k*(t37-1.d0 &
             )/t37*c22*c33
!
    t1 = c11*c22
    t3 = c12**2
    t5 = t1*c33-t3*c33
    t6 = t5**(1.d0/3.d0)
    t7 = 1.d0/t6
    t8 = c11+c22+c33
    t12 = c11*c33
    t15 = t7-t8/t6/t5*t12/3.d0
    t19 = t6**2
    t37 = sqrt(t5)
!
    sig(2) = 2.d0*c10*t15+2.d0*c01*( &
             (c11+c33)/t19-2.d0/3.d0*(t1+c22*c33+t12-t3)/t19/t5*t12)+4.d0*c20*(t8*t7-3.d0)*t15+&
             & k*(t37-1.d0 &
             )/t37*c11*c33
!
    t1 = c11+c22+c33
    t3 = c11*c22
    t5 = c12**2
    t7 = t3*c33-t5*c33
    t8 = t7**(1.d0/3.d0)
    t12 = 1.d0/t8/t7*c12*c33
    t15 = t8**2
    t38 = sqrt(t7)
    sig(4) = 4.d0/3.d0*c10*t1*t12+2.d0*c01*( &
             -2.d0*c12/t15+4.d0/3.d0*(t3+c11*c33+c22*c33-t5)/t15/t7*c12*c33)+8.d0/3.d0*c20*(t1/&
             &t8-3.d0)*t1*t12-2.d0*k*(t38-1.d0 &
             )/t38*c12*c33
    sig(3) = 0.d0
    sig(5) = 0.d0
    sig(6) = 0.d0
!
99  continue
end subroutine
