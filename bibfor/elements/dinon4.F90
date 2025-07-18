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
subroutine dinon4(neq, ul, dul, utl, nno, &
                  nbcomp, varimo, raide, nbpar, param, &
                  okdire, varipl)
! ----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
#include "asterc/r8miem.h"
    integer(kind=8) :: neq, nbcomp, nno, nbpar
    real(kind=8) :: ul(neq), dul(neq), utl(neq)
    real(kind=8) :: varimo(nbcomp*1), varipl(nbcomp*1)
    real(kind=8) :: raide(nbcomp), param(6, nbpar)
    aster_logical :: okdire(6)
!
! ======================================================================
!
!     RELATION DE COMPORTEMENT "BI-LINEAIRE" (DISCRET NON LINEAIRE).
!
!    f < FPRE
!       df = KDEB*dU
!    f > FPRE
!       df = KFIN*dU
!
!        KDEB  : raideur au début
!        KFIN  : raideur à la fin
!        FPREC : effort de pre-tension
!
!======================================================================
!
! IN  :
!       NEQ    : NOMBRE DE DDL DE L'ELEMENT
!       UL     : DEPLACEMENT PRECEDENT REPERE LOCAL (DIM NEQ)
!       DUL    : INCREMENT DE DEPLACEMENT REPERE LOCAL (DIM NEQ)
!       UTL    : DEPLACEMENT COURANT REPERE LOCAL (DIM NEQ)
!       NNO    : NOMBRE DE NOEUDS
!       NBCOMP : NOMBRE DE COMPOSANTES
!       VARIMO : VARIABLES INTERNES A T- (1 PAR COMPOSANTES)
!       RAIDE  : RAIDEUR ELASTIQUE DES DISCRETS
!       NBPAR  : NOMBRE MAXIMAL DE PARAMETRE DE LA LOI
!       PARAM  : PARAMETRES DE LA LOI
!       OKDIRE : VRAI SI LE COMPORTEMENT AFFECTE CETTE DIRECTION
!
! OUT :
!       RAIDE  : RAIDEUR QUASI-TANGENTE AU COMPORTEMENT DES DISCRETS
!       VARIPL : VARIABLES INTERNES INTERNES A T+ (1 PAR COMPOSANTES)
!
!***************** DECLARATION DES VARIABLES LOCALES *******************
!
    integer(kind=8) :: ii
    real(kind=8) :: ulel, dulel, utlel, r8min
!
    integer(kind=8) :: ivari
    real(kind=8) :: kdeb, kfin, fpre, useuil, fplus, fmoins, depl
!
    real(kind=8) :: zero
    parameter(zero=0.0d0)
!
!************ FIN DES DECLARATIONS DES VARIABLES LOCALES ***************
    r8min = r8miem()
!
    do ii = 1, nbcomp
!        INDEX DES VARIABLES INTERNES
        ivari = ii
!        PAR DEFAUT LES VARIABLES N'EVOLUENT PAS
        varipl(ivari) = varimo(ivari)
!        SI LE COMPORTEMENT EST BI-LINEAIRE DANS CETTE DIRECTION
        if (.not. okdire(ii)) goto 20
        kdeb = param(ii, 1)
        kfin = param(ii, 2)
        fpre = param(ii, 3)
!
!          write(*,'(I2,3(2X,E12.5))') II,KDEB,KFIN,FPRE
        if (abs(kdeb) .le. r8min) goto 20
!
        if (nno .eq. 1) then
            dulel = dul(ii)
            ulel = ul(ii)
            utlel = utl(ii)
        else
            dulel = dul(ii+nbcomp)-dul(ii)
            ulel = ul(ii+nbcomp)-ul(ii)
            utlel = utl(ii+nbcomp)-utl(ii)
        end if
!        SEUIL EN DEPLACEMENT
        useuil = abs(fpre/kdeb)
        if (abs(dulel) .gt. r8min) then
!           A L'INSTANT MOINS
            depl = abs(ulel)
            if (depl .le. useuil) then
                fmoins = kdeb*depl
            else
                fmoins = fpre+kfin*(depl-useuil)
            end if
            if (ulel .lt. zero) fmoins = -fmoins
!           A L'INSTANT PLUS
            depl = abs(utlel)
            if (depl .le. useuil) then
                fplus = kdeb*depl
                varipl(ivari) = 1.0d0
            else
                fplus = fpre+kfin*(depl-useuil)
                varipl(ivari) = 2.0d0
            end if
            if (utlel .lt. zero) fplus = -fplus
            raide(ii) = abs(fplus-fmoins)/abs(dulel)
        else
!           CAS OU DUEL=0 ==> UTLEL=ULEL
            depl = abs(utlel)
            if (depl .le. useuil) then
                raide(ii) = kdeb
            else
                raide(ii) = kfin
            end if
        end if
20      continue
    end do
!
end subroutine
