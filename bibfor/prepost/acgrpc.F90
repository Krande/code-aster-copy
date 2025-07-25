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
subroutine acgrpc(nbordr, kwork, sompgw, jrwork, tspaq, &
                  ipg, nommet, forcri, nompar, vanocr, &
                  respc, vnmax)
    implicit none
#include "jeveux.h"
#include "asterc/loisem.h"
#include "asterc/lor8em.h"
#include "asterc/r8pi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/acgrcr.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedisp.h"
#include "asterfort/utmess.h"
#include "asterfort/vecnuv.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbordr, kwork
    integer(kind=8) :: sompgw, jrwork, tspaq, ipg
    character(len=16) :: nommet, forcri
    character(len=8) :: nompar(35)
    real(kind=8) :: respc(24), vnmax(6), vanocr(23)
!
!
! BUT: POUR LA FATIGUE A AMPLITUDE CONSTANTE
!      DETERMINER LE PLAN DES MAX DES TAU_MAX ET CALCULER DES GRANDEURS
!
!
! REMARQUE: CETTE SUBROUTINE EST APPLICABLE POUR UN NOEUD OU IPG EGALE
!           A 1 ET SOMPGW = SOMNOW,JVECPG = JVECNO
! ----------------------------------------------------------------------
! ARGUMENTS :
!     NBORDR  : IN  : NOMBRE DE NUMEROS D'ORDRE.
!     KWORK   : IN  : KWORK = 0 ON TRAITE LA 1ERE MAILLE DU PAQUET DE
!                               MAILLES ;
!                     KWORK = 1 ON TRAITE LA IEME (I>1) MAILLE DU PAQUET
!                               MAILLES.
!     SOMPGW  : IN  : SOMME DES POINTS DE GAUSS DES N MAILLES PRECEDANT
!                     LA MAILLE COURANTE.
!     JRWORK  : IN  : ADRESSE DU VECTEUR DE TRAVAIL CONTENANT
!                     L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                     ATTACHES A CHAQUE POINT DE GAUSS DES MAILLES
!                     DU <<PAQUET>> DE MAILLES.
!     TSPAQ   : IN  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE MAILLES
!                     COURANT.
!     IPG     : IN  : IEME POINT DE GAUSS.
!    NOMMET     IN    NOM DE METHOD D'APPROCHEMENT DE CERCLE ("CERCLE
!                     EXACT" ET "CERCLE APPROCHE")
!    VALA       IN    VALEUR DU PARAMETRE a ASSOCIE AU CRITERE.
!    COEFPA     IN    COEFFICIENT DE PASSAGE CISAILLEMENT - UNIAXIAL.
!   VRSESU      OUT   TABLEAU DES RESULTATS (GRANDEURS ET DOMMAGE).
!                     POUR L'INSTANT, LA DIMENSION DE VRESU EST 24
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, n
    integer(kind=8) :: nbvec, dim, jvpg2
    integer(kind=8) :: jvecn2, jvecu2, jvecv2
    integer(kind=8) :: ideb, ngam
    integer(kind=8) :: tneces, tdisp(1), jvecno, jnorm2
    integer(kind=8) :: tab2(18), jvectn, jvectu, jvectv
    integer(kind=8) :: vali(2), tnecno, jnorma
!
    real(kind=8) :: dgam, dphi, tab1(18)
    real(kind=8) :: epsilo, gamma, pi
    real(kind=8) :: gammam, phim, dgam2, dphi2, phi0
    real(kind=8) :: nxm(2), nym(2), nzm(2)
!
!
!-----------------------------------------------------------------------
    data tab1/180.0d0, 60.0d0, 30.0d0, 20.0d0, 15.0d0, 12.857d0,&
     &             11.25d0, 10.588d0, 10.0d0, 10.0d0, 10.0d0, 10.588d0,&
     &             11.25d0, 12.857d0, 15.0d0, 20.0d0, 30.0d0, 60.0d0/
!
    data tab2/1, 3, 6, 9, 12, 14, 16, 17, 18, 18, 18, 17, 16, 14,&
     &           12, 9, 6, 3/
!
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------
!
!234567
!
!
    epsilo = 1.0d-7
    pi = r8pi()
! PROJECTION DE L'HISTORIQUE DU STRESS ET STRAIN DANS UN PLAN.
!
! CONSTRUCTION DU VECTEUR CONTENANT DELTA_TAU_MAX
! CONSTRUCTION DU VECTEUR CONTENANT LA VALEUR DU POINTEUR PERMETTANT
!              DE RETROUVER LE VECTEUR NORMAL ASSOCIE A DELTA_TAU_MAX
!
!      call wkvect('&&ACGRPC.DTAU_MAX', 'V V R', 209, jdtaum)
!      call wkvect('&&ACGRPC.RESU_N', 'V V I', 209, jresun)
! !
! CONSTRUCTION DU VECTEUR NORMAL SUR UNE DEMI SPHERE
! CONSTRUCTION DU VECTEUR U DANS LE PLAN TANGENT, SUR UNE DEMI SPHERE
! CONSTRUCTION DU VECTEUR V DANS LE PLAN TANGENT, SUR UNE DEMI SPHERE
!
    call wkvect('&&ACGRPC.VECT_NORMA', 'V V R', 630, jvectn)
    call wkvect('&&ACGRPC.VECT_TANGU', 'V V R', 630, jvectu)
    call wkvect('&&ACGRPC.VECT_TANGV', 'V V R', 630, jvectv)
!
    tneces = 209*nbordr*2
    tnecno = 209*nbordr
    call jedisp(1, tdisp)
    tdisp(1) = (tdisp(1)*loisem())/lor8em()
    if (tdisp(1) .lt. tneces) then
        vali(1) = tdisp(1)
        vali(2) = tneces
        call utmess('F', 'PREPOST5_8', ni=2, vali=vali)
    else
        call wkvect('&&ACGRPC.VECTNO', 'V V R', tneces, jvecno)
        call wkvect('&&ACGRPC.VECT_NOR', 'V V R', tnecno, jnorma)
    end if
!
    dgam = 10.0d0
!
    n = 0
    k = 1
    ideb = 1
    dim = 627
    do j = 1, 18
        gamma = (j-1)*dgam*r8dgrd()
        dphi = tab1(j)*r8dgrd()
        phi0 = dphi/2.0d0
        ngam = tab2(j)
!
        call vecnuv(ideb, ngam, gamma, phi0, dphi, &
                    n, k, dim, zr(jvectn), zr(jvectu), &
                    zr(jvectv))
!
    end do
!
!
    do i = 1, 24
        respc(i) = 0.0d0
    end do
!
!!!!IDENTIFIER LE PLAN DE MAXIMUM DE GRANDEUR CRITIQUE
    nbvec = 209
!
    call acgrcr(nbvec, jvectn, jvectu, jvectv, nbordr, &
                kwork, sompgw, jrwork, tspaq, ipg, &
                nommet, jvecno, jnorma, forcri, nompar, &
                vanocr, respc, vnmax)
!
!
!!!REFINEMENT DE PLAN CRITIQUE
    call wkvect('&&ACGRPC.VECT_NORMA2', 'V V R', 27, jvecn2)
    call wkvect('&&ACGRPC.VECT_TANGU2', 'V V R', 27, jvecu2)
    call wkvect('&&ACGRPC.VECT_TANGV2', 'V V R', 27, jvecv2)
!
!
    call wkvect('&&ACGRPC.VECTPG2', 'V V R', 18*nbordr, jvpg2)
    call wkvect('&&ACGRPC.VECTNO2', 'V V R', 9*nbordr, jnorm2)
!
    dim = 27
!
    do k = 1, 2
!
        nxm(k) = vnmax(1+(k-1)*3)
        nym(k) = vnmax(2+(k-1)*3)
        nzm(k) = vnmax(3+(k-1)*3)
!
        gammam = atan2(sqrt(abs(1.0d0-nzm(k)**2)), nzm(k))
        if (gammam .lt. 0.0d0) then
            gammam = gammam+pi
        end if
!
        if ((abs(nym(k)) .lt. epsilo) .and. (abs(nxm(k)) .lt. epsilo)) then
            phim = 0.0d0
        else
            phim = atan2(abs(nym(k)), nxm(k))
        end if
        if (phim .lt. 0.0d0) then
            phim = phim+pi
        end if
!
        if (abs(gammam) .lt. epsilo) then
            gamma = 5.0d0*r8dgrd()
            dphi2 = 60.0d0*r8dgrd()
            phi0 = 0.0d0
            n = 0
!
            call vecnuv(1, 6, gamma, phi0, dphi2, &
                        n, 1, dim, zr(jvecn2), zr(jvecu2), &
                        zr(jvecv2))
!
            gamma = 0.0d0
            phi0 = pi
!
            call vecnuv(1, 1, gamma, phi0, dphi2, &
                        n, 1, dim, zr(jvecn2), zr(jvecu2), &
                        zr(jvecv2))
!
            nbvec = 7
!
            call acgrcr(nbvec, jvecn2, jvecu2, jvecv2, nbordr, &
                        kwork, sompgw, jrwork, tspaq, ipg, &
                        nommet, jvpg2, jnorm2, forcri, nompar, &
                        vanocr, respc, vnmax)
!
        else
            dgam2 = 2.0d0*r8dgrd()
            dphi2 = dgam2/sin(gammam)
            n = 0
            do j = 1, 3
                gamma = gammam+(j-2)*dgam2
!
                call vecnuv(1, 3, gamma, phim, dphi2, &
                            n, 2, dim, zr(jvecn2), zr(jvecu2), &
                            zr(jvecv2))
!
            end do
!
            nbvec = 9
!
            call acgrcr(nbvec, jvecn2, jvecu2, jvecv2, nbordr, &
                        kwork, sompgw, jrwork, tspaq, ipg, &
                        nommet, jvpg2, jnorm2, forcri, nompar, &
                        vanocr, respc, vnmax)
!
        end if
!
    end do
!
!
    call jedetr('&&ACGRPC.DTAU_MAX')
    call jedetr('&&ACGRPC.RESU_N')
    call jedetr('&&ACGRPC.VECT_NORMA')
    call jedetr('&&ACGRPC.VECT_TANGU')
    call jedetr('&&ACGRPC.VECT_TANGV')
    call jedetr('&&ACGRPC.VECTNO')
    call jedetr('&&ACGRPC.VECT_NOR')
!
!     call jedetr('&&ACGRPC.VECT_NORMA1')
!     call jedetr('&&ACGRPC.VECT_TANGU1')
!     call jedetr('&&ACGRPC.VECT_TANGV1')
!     call jedetr('&&ACGRPC.VECTPG1')
!
    call jedetr('&&ACGRPC.VECT_NORMA2')
    call jedetr('&&ACGRPC.VECT_TANGU2')
    call jedetr('&&ACGRPC.VECT_TANGV2')
    call jedetr('&&ACGRPC.VECTPG2')
    call jedetr('&&ACGRPC.VECTNO2')
!
!
end subroutine
