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
subroutine impvoi(texte, nbma, iaddvo, iadvoi)
    implicit none
#include "jeveux.h"
    character(len=*) :: texte
    integer(kind=8) :: nbma, iaddvo, iadvoi
    integer(kind=8) :: ima, iv, is
!-----------FONCTIONS  D ACCES A VGE -----------------------------------
!     IADDVO : ADRESSE JEVEUX DU TABLEAU DE POINTEURS DANS LA SD EL_VOIS
!     IADVOI : ADRESSE JEVEUX DE LA SD EL_VOIS
!
!     DES DONNEES DES VOISINS DE LA MAILLE IMA (0 SI MAILLE PAS ACTIVE)
#define zzadvo(ima) zi(iaddvo+ima-1)+iadvoi
!
!     NOMBBRE DE VOISINS DE IMA EME MAILLE
#define zznbvo(ima) zi(zzadvo(ima)-1+1)
!
!     POUR LA MAILLE IMA
!     POUR LE VOISIN IV
!     ADRESSE DES DONNEES
! FAUX ???      ZZADVE(IMA,IV) = ZI(ZZADVO(IMA)-1+1+IV)+IADVOI-1
#define zzadve(ima,iv) zi(zzadvo(ima)-1+1+iv)+zzadvo(ima)-1
!
!     POUR LA MAILLE IMA
!     POUR LE VOISIN IV
!     TYPE DE VOISINAGE :
!        3D PAR FACE    : F3 : 1
!        2D PAR FACE    : F2 : 2
!        3D PAR ARRETE  : A3 : 3
!        2D PAR ARRETE  : A2 : 4
!        1D PAR ARRETE  : A1 : 5
!        3D PAR SOMMET  : S3 : 6
!        2D PAR SOMMET  : S2 : 7
!        1D PAR SOMMET  : S1 : 8
!        0D PAR SOMMET  : S0 : 9
#define zztyvo(ima,iv) zi(zzadve(ima,iv)-1+1)
!
!     POUR LA MAILLE IMA
!     POUR LE VOISIN IV
!     NUMERO DE MAILLE
#define zzmavo(ima,iv) zi(zzadve(ima,iv)-1+2)
!
!     POUR LA MAILLE IMA
!     POUR LE VOISIN IV
!     NOMBRE DE NOEUDS DE MAILLE
#define zznbno(ima,iv) zi(zzadve(ima,iv)-1+3)
!
!
!     POUR LA MAILLE IMA
!     POUR LE VOISIN IV
!        NOMBRE DE SOMMETS COMMUNS
#define zznbsc(ima,iv) zi(zzadve(ima,iv)-1+4)
!
!
!     POUR LA MAILLE IMA
!     POUR LE VOISIN IV
!     POUR LE SOMMET COMMUN IS
!     NUMERO LOCAL DANS IMA
#define zzloc1(ima,iv,is) zi(zzadve(ima,iv)-1+4+1+2*(is-1))
!
!
!
!     POUR LA MAILLE IMA
!     POUR LE VOISIN IV
!    POUR LE SOMMET COMMUN IS
!     NUMERO LOCAL DANS IV
#define zzloc2(ima,iv,is) zi(zzadve(ima,iv)-1+4+1+2*(is-1)+1)
!-----------FIN FONCTIONS  D ACCES A VGE -------------------------------
!
!
    write (6, *)
    write (6, *) ' IMPRESSION OBJET VOISINAGE VGE '
    write (6, *) texte
    write (6, *)
!
    do ima = 1, nbma
        write (6, 9000) ima, zznbvo(ima)
        do iv = 1, zznbvo(ima)
            write (6, 9010) iv, zztyvo(ima, iv), zzmavo(ima, iv), zznbno(ima, iv), zznbsc(ima, iv)
            do is = 1, zznbsc(ima, iv)
                write (6, 9020) is, zzloc1(ima, iv, is), zzloc2(ima, iv, is)
            end do
        end do
    end do
    write (6, *) ' FIN IMPRESSION OBJET VOISINAGE VGE '
    write (6, *)
!
9000 format(' MAILLE ', i8, ' NB VOIS ', i2)
9010 format(' VOISIN ', i2, ' TYPE ', i2, ' MAILLE ', i8, ' NB NOEUDS ', i2,&
   &       ' NB SOMMETS COMMUN ', i2)
9020 format(' IS ', i2, ' NUMLOC ', i2, ' NUMLOC DANS VOISIN ', i2)
end subroutine
