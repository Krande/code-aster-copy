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

subroutine gkmet3(nnoff, chfond, iadrgk, milieu, connex, &
                  iadgks, iadgki, abscur, num, typdis)

    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gmatc3.h"
#include "asterfort/gmatl3.h"
#include "asterfort/getvtx.h"
#include "asterfort/gsyste.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
    integer(kind=8)           :: nnoff, iadrgk, iadgks, iadgki, num
    character(len=16) :: typdis
    character(len=24) :: chfond, abscur
    aster_logical     :: milieu, connex

!      METHODE THETA-LAGRANGE ET G-LAGRANGE POUR LE CALCUL DE G(S)
!      K1(S) K2(S) ET K3(S)
!
! ENTREE
!
!     NNOFF    --> NOMBRE DE NOEUDS DU FOND DE FISSURE
!     CHFOND   --> COORDONNEES ET ABSCISSES CURVILIGNES DES NOEUDS
!                  DU FOND DE FISSURE
!     IADRGK   --> ADRESSE DE VALEURS DE GKTHI
!                 (G, K1, K2, K3 POUR LES CHAMPS THETAI)
!     MILIEU   --> .TRUE.  : ELEMENT QUADRATIQUE
!                  .FALSE. : ELEMENT LINEAIRE
!     CONNEX   --> .TRUE.  : SI FOND FERME
!                  .FALSE. : SI FOND OUVERT
!
!  SORTIE
!
!   IADGKS     --> ADRESSE DE VALEURS DE GKS
!                   (VALEUR DE G(S), K1(S), K2(S), K3(S), G_IRWIN(S))
!   IADGKI     --> ADRESSE DE VALEURS DE GKTHI
!                  (G, K1, K2, K3 POUR LES CHAMPS THETAI)
!   ABSCUR     --> VALEURS DES ABSCISSES CURVILIGNES S
!      NUM     --> 3 (LAGRANGE-LAGRANGE)
!              --> 4 (NOEUD-NOEUD)
! ......................................................................

    integer(kind=8)                            :: ifon, iadabs, ivect
    integer(kind=8)                            :: i, ibid
    real(kind=8)                       :: s1, s2, s3, sn2, sn1, sn
    real(kind=8), dimension(nnoff)     :: gthi, k1th, k2th, k3th
    real(kind=8), dimension(nnoff)     :: gs, k1s, k2s, k3s
    real(kind=8), dimension(nnoff)     :: gith, gis
    real(kind=8), dimension(nnoff)     :: g1th, g2th, g3th
    real(kind=8), dimension(nnoff)     :: g1s, g2s, g3s
    character(len=24)                  :: lissg, vect, matr

! ......................................................................

    call jemarq()

    call jeveuo(chfond, 'L', ifon)
    call jeveuo(abscur, 'E', iadabs)

    do i = 1, nnoff

        zr(iadabs-1+(i-1)+1) = zr(ifon-1+4*(i-1)+4)
        gthi(i) = zr(iadrgk-1+(i-1)*8+1)
        g1th(i) = zr(iadrgk-1+(i-1)*8+2)
        g2th(i) = zr(iadrgk-1+(i-1)*8+3)
        g3th(i) = zr(iadrgk-1+(i-1)*8+4)
        k1th(i) = zr(iadrgk-1+(i-1)*8+5)
        k2th(i) = zr(iadrgk-1+(i-1)*8+6)
        k3th(i) = zr(iadrgk-1+(i-1)*8+7)
        if (typdis .ne. 'COHESIF') then
            gith(i) = g1th(i)*g1th(i)+g2th(i)*g2th(i)+g3th(i)*g3th(i)
        else
            gith(i) = gthi(i)
        end if
    end do

!   CHOIX DU LISSAGE
    call getvtx('LISSAGE', 'LISSAGE_G', iocc=1, scal=lissg, nbret=ibid)

    if (lissg .eq. 'LAGRANGE_NO_NO') then
        num = 4

!       CALCUL DE LA MATRICE DU SYTEME LINAIRE : MATRICE LUMPEE
        vect = '&&METHO3.VECT'
        call gmatl3(nnoff, milieu, connex, &
                    abscur, vect)

        call jeveuo(vect, 'L', ivect)

!       RESOLUTION DU SYSTEME LINEAIRE : MATRICE DIAGONALE
        do i = 1, nnoff
            gs(i) = gthi(i)/zr(ivect+i-1)
            k1s(i) = k1th(i)/zr(ivect+i-1)
            k2s(i) = k2th(i)/zr(ivect+i-1)
            k3s(i) = k3th(i)/zr(ivect+i-1)
            gis(i) = gith(i)/(zr(ivect+i-1)*zr(ivect+i-1))
        end do

        if (.not. connex) then
!       CORRECTION DES VALEURS ASSOCIEES AU 1ER ET DERNIER CHAMPS THETA
            if (nnoff .gt. 2) then
                s1 = zr(iadabs-1+1)
                s2 = zr(iadabs-1+2)
                s3 = zr(iadabs-1+3)
                sn2 = zr(iadabs-1+nnoff-2)
                sn1 = zr(iadabs-1+nnoff-1)
                sn = zr(iadabs-1+nnoff)

                gs(1) = gs(2)+(s1-s2)*(gs(3)-gs(2))/(s3-s2)
                k1s(1) = k1s(2)+(s1-s2)*(k1s(3)-k1s(2))/(s3-s2)
                k2s(1) = k2s(2)+(s1-s2)*(k2s(3)-k2s(2))/(s3-s2)
                k3s(1) = k3s(2)+(s1-s2)*(k3s(3)-k3s(2))/(s3-s2)
                gis(1) = gis(2)+(s1-s2)*(gis(3)-gis(2))/(s3-s2)
                gs(nnoff) = gs(nnoff-1)+(sn-sn1)*(gs(nnoff-2)-gs(nnoff-1))/(sn2-sn1)
                k1s(nnoff) = k1s(nnoff-1)+(sn-sn1)*(k1s(nnoff-2)-k1s(nnoff-1))/(sn2-sn1)
                k2s(nnoff) = k2s(nnoff-1)+(sn-sn1)*(k2s(nnoff-2)-k2s(nnoff-1))/(sn2-sn1)
                k3s(nnoff) = k3s(nnoff-1)+(sn-sn1)*(k3s(nnoff-2)-k3s(nnoff-1))/(sn2-sn1)
                gis(nnoff) = gis(nnoff-1)+(sn-sn1)*(gis(nnoff-2)-gis(nnoff-1))/(sn2-sn1)

            end if
        end if

    else if (lissg .eq. 'LAGRANGE') then
        num = 3

!       CALCUL DE LA MATRICE DU SYTEME LINAIRE
        matr = '&&METHO3.MATRI'
        call gmatc3(nnoff, milieu, connex, &
                    abscur, matr)

!       X-FEM : CORRECTION VALEURS EXTREMITES (RESULTAT + PRECIS)
        if (.not. connex) then
            if (nnoff .ne. 2) then

                s1 = zr(iadabs-1+1)
                s2 = zr(iadabs-1+2)
                s3 = zr(iadabs-1+3)
                sn2 = zr(iadabs-1+nnoff-2)
                sn1 = zr(iadabs-1+nnoff-1)
                sn = zr(iadabs-1+nnoff)

                !            CORRECTION DANS LE CAS LINEAIRE
                if (.not. milieu) then
                    gthi(1) = gthi(2)*(s2-s1)/(s3-s1)
                    k1th(1) = k1th(2)*(s2-s1)/(s3-s1)
                    k2th(1) = k2th(2)*(s2-s1)/(s3-s1)
                    k3th(1) = k3th(2)*(s2-s1)/(s3-s1)
                    gith(1) = gith(2)*(s2-s1)/(s3-s1)
                    gthi(nnoff) = gthi(nnoff-1)*(sn-sn1)/(sn-sn2)
                    k1th(nnoff) = k1th(nnoff-1)*(sn-sn1)/(sn-sn2)
                    k2th(nnoff) = k2th(nnoff-1)*(sn-sn1)/(sn-sn2)
                    k3th(nnoff) = k3th(nnoff-1)*(sn-sn1)/(sn-sn2)
                    gith(nnoff) = gith(nnoff-1)*(sn-sn1)/(sn-sn2)

                    !            CORRECTION DANS LE CAS QUADRATIQUE
                else if (milieu) then
                    gthi(1) = gthi(2)/4.d0
                    k1th(1) = k1th(2)/4.d0
                    k2th(1) = k2th(2)/4.d0
                    k3th(1) = k3th(2)/4.d0
                    gith(1) = gith(2)/4.d0
                    gthi(nnoff) = gthi(nnoff-1)/4.d0
                    k1th(nnoff) = k1th(nnoff-1)/4.d0
                    k2th(nnoff) = k2th(nnoff-1)/4.d0
                    k3th(nnoff) = k3th(nnoff-1)/4.d0
                    gith(nnoff) = gith(nnoff-1)/4.d0
                end if

            end if
        end if
!       SYSTEME LINEAIRE:  MATR*GS = GTHI
        call gsyste(matr, nnoff, nnoff, gthi, gs)

!       SYSTEME LINEAIRE:  MATR*K1S = K1TH
        call gsyste(matr, nnoff, nnoff, k1th, k1s)

!       SYSTEME LINEAIRE:  MATR*K2S = K2TH
        call gsyste(matr, nnoff, nnoff, k2th, k2s)

!       SYSTEME LINEAIRE:  MATR*K3S = K3TH
        call gsyste(matr, nnoff, nnoff, k3th, k3s)

!       SYSTEMES LINEAIRES POUR GIRWIN
        call gsyste(matr, nnoff, nnoff, g1th, g1s)
        call gsyste(matr, nnoff, nnoff, g2th, g2s)
        call gsyste(matr, nnoff, nnoff, g3th, g3s)

        do i = 1, nnoff
            gis(i) = g1s(i)*g1s(i)+g2s(i)*g2s(i)+g3s(i)*g3s(i)
        end do

    end if

!   ----------------------------------------------------------------
!                              RECOPIES
!   ----------------------------------------------------------------
    if (typdis .ne. 'COHESIF') then
        do i = 1, nnoff
            zr(iadgks-1+(i-1)*5+1) = gs(i)
            zr(iadgks-1+(i-1)*5+2) = k1s(i)
            zr(iadgks-1+(i-1)*5+3) = k2s(i)
            zr(iadgks-1+(i-1)*5+4) = k3s(i)
            zr(iadgks-1+(i-1)*5+5) = gis(i)
        end do
    else if (typdis .eq. 'COHESIF') then
        do i = 1, nnoff
            zr(iadgks-1+(i-1)*5+1) = gs(i)
            k1s(i) = sqrt(k1s(i))
            zr(iadgks-1+(i-1)*5+2) = k1s(i)
            if (g2th(i) .ge. 0.d0) k2s(i) = sqrt(abs(k2s(i)))
            if (g2th(i) .lt. 0.d0) k2s(i) = -sqrt(abs(k2s(i)))
            zr(iadgks-1+(i-1)*5+3) = k2s(i)
            if (g3th(i) .ge. 0.d0) k3s(i) = sqrt(abs(k3s(i)))
            if (g3th(i) .lt. 0.d0) k3s(i) = -sqrt(abs(k3s(i)))
            zr(iadgks-1+(i-1)*5+4) = k3s(i)
            zr(iadgks-1+(i-1)*5+5) = gs(i)
        end do
    end if

    do i = 1, nnoff
        zr(iadgki-1+(i-1)*5+1) = zr(iadrgk-1+(i-1)*8+1)
        zr(iadgki-1+(i-1)*5+2) = zr(iadrgk-1+(i-1)*8+5)
        zr(iadgki-1+(i-1)*5+3) = zr(iadrgk-1+(i-1)*8+6)
        zr(iadgki-1+(i-1)*5+4) = zr(iadrgk-1+(i-1)*8+7)
    end do

!!   CALCUL DES ANGLES DE PROPAGATION DE FISSURE LOCAUX BETA
!    do i = 1, nnoff
!        betas(i) = 0.0d0
!        if ( abs(k2s(i)) .ge. 1.e-12 ) betas(i) = 2.0d0*atan2(0.25d0*(k1s(i)/k2s(i) - &
!                                                  sign(1.0d0, k2s(i))*sqrt((k1s(i)/k2s(i)) &
!                                                  **2.0d0+8.0d0)), 1.0d0)
!        zr(iadgks-1+(i-1)*6+6) = betas(i)
!    enddo

!!    LISSAGE PAR MOYENNE GLISSANTE SI METHODE COHESIVE
!     if(typdis.eq.'COHESIF') then
!         if(nnoff.gt.2) then
!             zr(iadgks-1+6) = (betas(1)+betas(2)+betas(3))/3.d0
!             zr(iadgks-1+(nnoff-1)*6+6)=(betas(nnoff-2)+betas(nnoff-1)+betas(nnoff))/3.d0
!             do i= 2,nnoff-1
!                 beta = (betas(i-1)+betas(i)+betas(i+1))/3.d0
!                 zr(iadgks-1+(i-1)*6+6)=beta
!             end do
!         else if(nnoff.eq.2) then
!             zr(iadgks-1+6) = (betas(1)+betas(2))/2.d0
!             zr(iadgks-1+6+6) = (betas(1)+betas(2))/2.d0
!        endif
!     endif

    call jedetr('&&METHO3.MATRI')
    call jedetr('&&METHO3.VECT')

    call jedema()

end subroutine
