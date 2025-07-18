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
subroutine ancrca(icabl, nbno, s, alpha, f0, &
                  delta, ea, frco, frli, sa, &
                  d, f)
    implicit none
!  DESCRIPTION : CALCUL DE LA TENSION LE LONG D'UN CABLE EN PRENANT EN
!  -----------   COMPTE LES PERTES PAR RECUL DE L'ANCRAGE
!                APPELANT : TENSK1, TENSK2
!
!  IN     : ICABL  : INTEGER , SCALAIRE
!                    NUMERO DU CABLE
!  IN     : NBNO   : INTEGER , SCALAIRE
!                    NOMBRE DE NOEUDS DU CABLE
!  IN     : S      : REAL*8 , VECTEUR DE DIMENSION NBNO
!                    CONTIENT LES VALEURS DE L'ABSCISSE CURVILIGNE
!                    LE LONG DU CABLE
!  IN     : ALPHA  : REAL*8 , VECTEUR DE DIMENSION NBNO
!                    CONTIENT LES VALEURS DE LA DEVIATION ANGULAIRE
!                    CUMULEE LE LONG DU CABLE
!  IN     : F0     : REAL*8 , SCALAIRE
!                    VALEUR DE LA TENSION APPLIQUEE A L'UN OU AUX DEUX
!                    ANCRAGES ACTIFS DU CABLE
!  IN     : DELTA  : REAL*8 , SCALAIRE
!                    VALEUR DU RECUL DE L'ANCRAGE
!  IN     : EA     : REAL*8 , SCALAIRE
!                    VALEUR DU MODULE D'YOUNG DE L'ACIER
!  IN     : FRCO   : REAL*8 , SCALAIRE
!                    VALEUR DU COEFFICIENT DE FROTTEMENT EN COURBE
!                    (CONTACT ENTRE LE CABLE ACIER ET LE MASSIF BETON)
!  IN     : FRLI   : REAL*8 , SCALAIRE
!                    VALEUR DU COEFFICIENT DE FROTTEMENT EN LIGNE
!                    (CONTACT ENTRE LE CABLE ACIER ET LE MASSIF BETON)
!  IN     : SA     : REAL*8 , SCALAIRE
!                    VALEUR DE L'AIRE DE LA SECTION DROITE DU CABLE
!  OUT    : D      : REAL*8 , SCALAIRE
!                    VALEUR DE LA LONGUEUR SUR LAQUELLE ON DOIT
!                    PRENDRE EN COMPTE LES PERTES DE TENSION
!                    PAR RECUL DE L'ANCRAGE
!  OUT    : F      : REAL*8 , VECTEUR DE DIMENSION NBNO
!                    CONTIENT LES VALEURS DE LA TENSION LE LONG DU CABLE
!                    APRES PRISE EN COMPTE DES PERTES PAR FROTTEMENT ET
!                    DES PERTES PAR RECUL DE L'ANCRAGE
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
! ARGUMENTS
! ---------
#include "asterfort/utmess.h"
#include "asterfort/wdefca.h"
    integer(kind=8) :: icabl, nbno
    real(kind=8) :: s(*), alpha(*), f0, delta, ea, frco, frli, sa, d, f(*)
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: iinf, ino, isup
    real(kind=8) :: alphad, alphai, alphas, df, ds, epsw, f2, pente, wcr, wdef
    real(kind=8) :: winf, wsup
    character(len=3) :: k3b
!
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    if (delta .eq. 0.0d0) then
        d = 0.0d0
        goto 999
    end if
!
    epsw = 1.0d-04
    wdef = ea*sa*delta
    wcr = wdefca(nbno, s, alpha, f0, frco, frli)
!
    if (wcr .lt. wdef) then
!
        write (k3b, '(I3)') icabl
        call utmess('A', 'MODELISA2_3', sk=k3b)
!
        d = s(nbno)
        if (wcr/wdef .lt. epsw) then
            df = wdef/d
            do ino = 1, nbno
                f(ino) = f(ino)-df
            end do
        else
            df = (wdef-wcr)/d
            f2 = f(nbno)
            f2 = f2*f2
            do ino = 1, nbno
                f(ino) = f2/f(ino)-df
            end do
        end if
!
    else if (wcr .eq. wdef) then
!
        write (k3b, '(I3)') icabl
        call utmess('A', 'MODELISA2_4', sk=k3b)
!
        d = s(nbno)
        f2 = f(nbno)
        f2 = f2*f2
        do ino = 1, nbno
            f(ino) = f2/f(ino)
        end do
!
    else
!
        iinf = 1
        isup = nbno
40      continue
        ino = (iinf+isup)/2
        if (wdefca(ino, s, alpha, f0, frco, frli) .gt. wdef) then
            isup = ino
        else
            iinf = ino
        end if
        if (isup-iinf .gt. 1) goto 40
!
        winf = wdefca(iinf, s, alpha, f0, frco, frli)
        wsup = wdefca(isup, s, alpha, f0, frco, frli)
        ds = s(isup)-s(iinf)
        pente = (wsup-winf)/ds
        d = s(iinf)+(wdef-winf)/pente
!
        alphai = alpha(iinf)
        alphas = alpha(isup)
        pente = (alphas-alphai)/ds
        alphad = alphai+pente*(d-s(iinf))
!
        f2 = f0*dble(exp(-frco*alphad-frli*d))
        f2 = f2*f2
        do ino = 1, iinf
            f(ino) = f2/f(ino)
        end do
!
    end if
!
999 continue
!
! --- FIN DE ANCRCA.
end subroutine
