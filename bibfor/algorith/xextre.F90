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
subroutine xextre(iptbor, vectn, nbfacb, jbas, jborl, &
                  jdirol, jnvdir)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "blas/ddot.h"
    integer(kind=8) :: iptbor(2), nbfacb
    integer(kind=8) :: jbas, jborl, jdirol, jnvdir
    real(kind=8) :: vectn(12)
!
!
!           CALCUL DES VECTEURS DE PROPAGATION AUX EXTREMITES DU FOND
!           DE FISSURE
!
!     ENTREE
!       IPTBOR   : VECTEUR CONTENANT LES INDICES DU OU DES POINTS DE
!                  BORD DE LA MAILLE
!       VECTN    : VECTEUR CONTENANT LES VECTEURS NORMAUX DES FACES DE
!                  BORD DE LA MAILLE
!       NBFACB   : NOMBRE DE FACES DE BORD DANS LA MAILLE
!       JBORL    : ADRESSE DU VECTEUR PERMETTANT DE SAVOIR SI LE VECTEUR
!                  DE DIRECTION DE PROPAGATION A DEJA ETE RECALCULE OU
!                  NON AUX POINTS EXTREMITES DE FONFIS (POUR SAVOIR SI
!                  ON DOIT REMPLACER LA VALEUR EXISTANTE OU LA LUI
!                  AJOUTER)
!       JDIROL   : ADRESSE DES VECTEURS DIRECTIONS DE PROPAGATION
!                  INITIAUX (CAD SANS MODIFICATION DES VECTEURS AUX
!                  POINTS EXTREMITES DE FONFIS)
!       JNVDIR   : ADRESSE DU VECTEUR CONTENANT 0 OU 1 AUX POINTS
!                  EXTREMITES DE FONFIS:
!                  0: LE PRODUIT SCALAIRE ENTRE LA NORMALE A LA FACE DE
!                     BORD ET LE VDIR INITIAL ESI INFERIEUR A 0
!                  1: LE PRODUIT SCALAIRE EST SUPERIEUR OU EGAL A 0
!     SORTIE
!       JBAS     : ADRESSE DU VECTEUR 'BASEFOND'
!
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: h, i, ind, k, nptbom, signe
    real(kind=8) :: maxi, norm, proj, sens, temp
    real(kind=8) :: normal(3), vdir(3), vdirol(3), vnor(3)
    aster_logical :: change, vecmax
    blas_int :: b_incx, b_incy, b_n
! ----------------------------------------------------------------------
    call jemarq()
!
    maxi = 0.d0
    change = .true.
    vecmax = .false.
    nptbom = 1
    if (iptbor(2) .ne. 0) nptbom = 2
!
!     BOUCLE SUR LE NOMBRE DE POINTS DU FOND DE LA MAILLE
!     QUI SONT SUR UNE FACE DE BORD
!     (CAS GENERAL NPTBOM=1)
    do i = 1, nptbom
!
        sens = 1.d0
!
!        RECUPERATION DE L'ANCIEN VECTEUR DE DIRECTION DE PROPAGATION
        if (.not. zl(jborl-1+iptbor(i))) then
            zr(jdirol-1+1+3*(iptbor(i)-1)) = zr(jbas-1+6*(iptbor(i)-1)+4)
            zr(jdirol-1+2+3*(iptbor(i)-1)) = zr(jbas-1+6*(iptbor(i)-1)+5)
            zr(jdirol-1+3+3*(iptbor(i)-1)) = zr(jbas-1+6*(iptbor(i)-1)+6)
        end if
!
        vdirol(1) = zr(jdirol-1+1+3*(iptbor(i)-1))
        vdirol(2) = zr(jdirol-1+2+3*(iptbor(i)-1))
        vdirol(3) = zr(jdirol-1+3+3*(iptbor(i)-1))
!
!C--     CAS 1: ON A DEUX POINTS DE BORD DANS UNE MEME MAILLE
        if (nptbom .gt. 1) then
            normal(1) = vectn(1+3*(i-1))
            normal(2) = vectn(2+3*(i-1))
            normal(3) = vectn(3+3*(i-1))
        end if
!
!C--     CAS 2: LA MAILLE N'A QU'UNE FACE DE BORD
        if (nbfacb .eq. 1) then
!          ON VERIFIE QUE LA FACE EST A PRENDRE EN COMPTE
!
!          N
            normal(1) = vectn(1)
            normal(2) = vectn(2)
            normal(3) = vectn(3)
!
!          N.VDIROLD
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            proj = ddot(b_n, normal, b_incx, vdirol, b_incy)
!
            if (proj .lt. 0) then
                signe = 0
            else
                signe = 1
            end if
!
!          NVDIR
            if (.not. zl(jborl-1+iptbor(i))) then
                zi(jnvdir-1+iptbor(i)) = signe
            else
                if (zi(jnvdir-1+iptbor(i)) .lt. signe) then
                    vecmax = .true.
                    zi(jnvdir-1+iptbor(i)) = signe
                else if (zi(jnvdir-1+iptbor(i)) .gt. signe) then
                    change = .false.
                end if
            end if
!
!C--     CAS 3: LA MAILLE A PLUSIEURS FACES DE BORD
        else if ((nbfacb .gt. 1) .and. (nptbom .eq. 1)) then
!          ON CHOISIT LA BONNE NORMALE
            do h = 1, nbfacb
!            N.VDIROLD
                proj = vectn( &
                       1+3*(h-1))*vdirol(1)+vectn(2+3*(h-1))*vdirol(2)+vectn(3+3*(h-1))*vdirol(3)
!
                if (proj .ge. maxi) then
                    maxi = proj
                    ind = h
                end if
!
            end do
!
            normal(1) = vectn(1+3*(ind-1))
            normal(2) = vectn(2+3*(ind-1))
            normal(3) = vectn(3+3*(ind-1))
        end if
!
!        SI ON A TROUVE UNE 'BONNE' FACE DE BORD
        if (change) then
!          CALCUL DE VDIR, LE NOUVEAU VECTEUR DE DIRECTION DE
!          PROPAGATION
            vnor(1) = zr(jbas-1+6*(iptbor(i)-1)+1)
            vnor(2) = zr(jbas-1+6*(iptbor(i)-1)+2)
            vnor(3) = zr(jbas-1+6*(iptbor(i)-1)+3)
!
            call provec(vnor, normal, vdir)
!
!          VERIFICATION QUE VDIR EST DANS LE BON SENS
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            proj = ddot(b_n, vdir, b_incx, vdirol, b_incy)
!
            if (proj .lt. 0) sens = -1.d0
!
!          NORMALISATION DE VDIR
            call normev(vdir, norm)
!
!          SI LE VECTEUR DE DIRECTION DE PROPAGATION
!          N'A PAS ENCORE ETE RECALCULE, ON LE REMPLACE DANS LA BASE
            if ((.not. zl(jborl-1+iptbor(i))) .or. (vecmax)) then
                do k = 1, 3
                    zr(jbas-1+6*(iptbor(i)-1)+k+3) = sens*vdir(k)
                    zl(jborl-1+iptbor(i)) = .true.
                end do
!          SINON ON L'AJOUTE (ON NORMALISE LE VECTEUR PAR LA SUITE, CE
!          QUI REVIENT A FAIRE UNE MOYENNE DES VECTEURS CALCULES)
            else
                do k = 1, 3
                    temp = zr(jbas-1+6*(iptbor(i)-1)+k+3)
                    zr(jbas-1+6*(iptbor(i)-1)+k+3) = temp+sens*vdir(k)
                end do
            end if
        end if
    end do
!
    call jedema()
end subroutine
