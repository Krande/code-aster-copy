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
subroutine tensk2(icabl, nbno, s, alpha, f0, &
                  delta, ea, frco, frli, sa, &
                  f)
    implicit none
!  DESCRIPTION : CALCUL DE LA TENSION LE LONG D'UN CABLE EN PRENANT EN
!  -----------   COMPTE LES PERTES PAR FROTTEMENT ET LES PERTES PAR
!                RECUL DES ANCRAGES
!                CAS DE DEUX ANCRAGES ACTIFS
!                APPELANT : TENSCA
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
!                    VALEUR DE LA TENSION APPLIQUEE AUX DEUX ANCRAGES
!                    ACTIFS DU CABLE
!  IN     : DELTA  : REAL*8 , SCALAIRE
!                    VALEUR DU RECUL DES DEUX ANCRAGES
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
!  OUT    : F      : REAL*8 , VECTEUR DE DIMENSION NBNO
!                    CONTIENT LES VALEURS DE LA TENSION LE LONG DU CABLE
!                    APRES PRISE EN COMPTE DES PERTES PAR FROTTEMENT ET
!                    DES PERTES PAR RECUL DES DEUX ANCRAGES
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
!
! ARGUMENTS
! ---------
#include "jeveux.h"
#include "asterfort/ancrca.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    integer(kind=8) :: icabl, nbno
    real(kind=8) :: s(*), alpha(*), f0, delta, ea, frco, frli, sa, f(*)
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: ino
    real(kind=8) :: alphal, d1, d2, df, long, wcr
    real(kind=8), pointer :: absc2(:) => null()
    real(kind=8), pointer :: alpha2(:) => null()
    real(kind=8), pointer :: f1(:) => null()
    real(kind=8), pointer :: f2(:) => null()
    real(kind=8), pointer :: fmax(:) => null()
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call jemarq()
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   CREATION DES OBJETS DE TRAVAIL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    AS_ALLOCATE(vr=absc2, size=nbno)
    AS_ALLOCATE(vr=alpha2, size=nbno)
    AS_ALLOCATE(vr=f1, size=nbno)
    AS_ALLOCATE(vr=f2, size=nbno)
    AS_ALLOCATE(vr=fmax, size=nbno)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   CREATION DES DISCRETISATIONS DE L'ABSCISSE CURVILIGNE ET DE LA
!     DEVIATION ANGULAIRE CUMULEE CORRESPONDANT AU SENS DE PARCOURS
!     INVERSE LE LONG DU CABLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    long = s(nbno)
    absc2(1) = 0.0d0
    do ino = 2, nbno
        absc2(ino) = long-s(nbno-ino+1)
    end do
!
    alphal = alpha(nbno)
    alpha2(1) = 0.0d0
    do ino = 2, nbno
        alpha2(ino) = alphal-alpha(nbno-ino+1)
    end do
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 3   PRISE EN COMPTE DES PERTES DE TENSION PAR FROTTEMENT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 3.1 TENSION APPLIQUEE AU PREMIER ANCRAGE ACTIF
! ---
    do ino = 1, nbno
        f1(ino) = f0*dble(exp(-frco*alpha(ino)-frli*s(ino)))
    end do
!
! 3.2 TENSION APPLIQUEE AU SECOND ANCRAGE ACTIF
! ---
    do ino = 1, nbno
        f2(ino) = f0*dble(exp(-frco*alpha2(ino)-frli*absc2(ino)))
    end do
!
! 3.3 VALEUR MAXIMALE INDUITE PAR LES TENSIONS APPLIQUEES AUX DEUX
! --- ANCRAGES ACTIFS APRES PRISE EN COMPTE DES PERTES PAR FROTTEMENT
!
    do ino = 1, nbno
        fmax(ino) = dble(max(f1(ino), f2(1+nbno-ino)))
    end do
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 4   PRISE EN COMPTE DES PERTES DE TENSION PAR RECUL DES DEUX ANCRAGES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 4.1 TENSION APPLIQUEE AU PREMIER ANCRAGE ACTIF
! ---
    call ancrca(icabl, nbno, s, alpha, f0, &
                delta, ea, frco, frli, sa, &
                d1, f1)
!
! 4.2 TENSION APPLIQUEE AU SECOND ANCRAGE ACTIF
! ---
    call ancrca(icabl, nbno, absc2, alpha2, f0, &
                delta, ea, frco, frli, sa, &
                d2, f2)
!
! 4.3 VALEUR FINALE INDUITE PAR LES TENSIONS APPLIQUEES AUX DEUX
! --- ANCRAGES ACTIFS APRES PRISE EN COMPTE DES PERTES PAR RECUL
!     DES DEUX ANCRAGES
!
    if (d1+d2 .lt. long) then
        do ino = 1, nbno
            f(ino) = dble(max(f1(ino), f2(1+nbno-ino)))
        end do
    else
        do ino = 1, nbno
            f(ino) = dble(min(f1(ino), f2(1+nbno-ino)))
        end do
    end if
!
! 4.4 CORRECTION SI RECOUVREMENT DES LONGUEURS D'APPLICATION DES PERTES
! --- DE TENSION PAR RECUL DES DEUX ANCRAGES
!
    if (d1+d2 .gt. long) then
        wcr = 0.0d0
        do ino = 1, nbno-1
            wcr = wcr+( &
                  (fmax(ino)-f(ino))+(fmax(ino+1)-f(ino+1)))*(s(ino+1)-s(ino) &
                                                              )/2.0d0
        end do
        df = (ea*sa*2.0d0*delta-wcr)/long
        do ino = 1, nbno
            f(ino) = f(ino)-df
        end do
    end if
!
! --- MENAGE
    AS_DEALLOCATE(vr=absc2)
    AS_DEALLOCATE(vr=alpha2)
    AS_DEALLOCATE(vr=f1)
    AS_DEALLOCATE(vr=f2)
    AS_DEALLOCATE(vr=fmax)
!
    call jedema()
!
! --- FIN DE TENSK2.
end subroutine
