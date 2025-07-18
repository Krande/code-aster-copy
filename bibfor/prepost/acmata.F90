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

subroutine acmata(nbordr, kwork, sompgw, jrwork, tspaq, &
                  ipg, nommet, vrespc)
    implicit none
#include "jeveux.h"
#include "asterc/loisem.h"
#include "asterc/lor8em.h"
#include "asterc/r8pi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedisp.h"
#include "asterfort/jemarq.h"
#include "asterfort/raycir.h"
#include "asterfort/taurlo.h"
#include "asterfort/utmess.h"
#include "asterfort/vecnuv.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbordr, kwork
    integer(kind=8) :: sompgw, jrwork, tspaq, ipg
    character(len=16) :: nommet
    real(kind=8) :: vrespc(24)
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
!   vrespc      OUT   TABLEAU DES RESULTATS (GRANDEURS ET DOMMAGE).
!                     POUR L'INSTANT, LA DIMENSION DE VRESU EST 24
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, n
    integer(kind=8) :: nbvec, dim, mnmax(2), jvpg1, jvpg2
    integer(kind=8) :: jvecn2, jvecu2, jvecv2, jvecn1, jvecu1, jvecv1
    integer(kind=8) :: adrs, decal, tab2(18), vali(2), ideb, ngam, jresun
    integer(kind=8) :: tneces, tdisp(1), jvecno, tnecno, jnorma, dectau
    integer(kind=8) :: jdtaum, jvectn, jvectu, jvectv
!
    real(kind=8) :: dgam, pi, dphi, tab1(18)
    real(kind=8) :: epsilo, gamma
    real(kind=8) :: gammam, phim, dgam2, dphi2, phi0, dtaum(2)
    real(kind=8) :: nxm(2), nym(2), nzm(2)
    real(kind=8) :: sixx, siyy, sizz, sixy, sixz, siyz, fxm(2), fym(2)
    real(kind=8) :: fzm(2), epsxx, epsyy, epszz, epsxy, epsxz, epsyz
    real(kind=8) :: norm(2), normax(2), snorm(2), epsxm(2), epsym(2)
    real(kind=8) :: epszm(2), epnorm(2), epnmax(2), sepnmx(2), normoy(2)
    real(kind=8) :: epnmoy(2)
    real(kind=8) :: phydro, phydrm
!
!     ------------------------------------------------------------------
!
!234567
!
!-----------------------------------------------------------------------
    data tab1/180.0d0, 60.0d0, 30.0d0, 20.0d0, 15.0d0, 12.857d0,&
     &             11.25d0, 10.588d0, 10.0d0, 10.0d0, 10.0d0, 10.588d0,&
     &             11.25d0, 12.857d0, 15.0d0, 20.0d0, 30.0d0, 60.0d0/
!
    data tab2/1, 3, 6, 9, 12, 14, 16, 17, 18, 18, 18, 17, 16, 14,&
     &           12, 9, 6, 3/
!
    pi = r8pi()
!-----------------------------------------------------------------------
!
    call jemarq()
!
! PROJECTION DE L'HISTORIQUE DU STRESS ET STRAIN DANS UN PLAN.
!
! CONSTRUCTION DU VECTEUR CONTENANT DELTA_TAU_MAX
! CONSTRUCTION DU VECTEUR CONTENANT LA VALEUR DU POINTEUR PERMETTANT
!              DE RETROUVER LE VECTEUR NORMAL ASSOCIE A DELTA_TAU_MAX
!
    call wkvect('&&ACMATA.DTAU_MAX', 'V V R', 209, jdtaum)
    call wkvect('&&ACMATA.RESU_N', 'V V I', 209, jresun)
!
! CONSTRUCTION DU VECTEUR NORMAL SUR UNE DEMI SPHERE
! CONSTRUCTION DU VECTEUR U DANS LE PLAN TANGENT, SUR UNE DEMI SPHERE
! CONSTRUCTION DU VECTEUR V DANS LE PLAN TANGENT, SUR UNE DEMI SPHERE
!
    call wkvect('&&ACMATA.VECT_NORMA', 'V V R', 630, jvectn)
    call wkvect('&&ACMATA.VECT_TANGU', 'V V R', 630, jvectu)
    call wkvect('&&ACMATA.VECT_TANGV', 'V V R', 630, jvectv)
!
    tneces = 209*nbordr*2
    tnecno = 209*nbordr
!
    call jedisp(1, tdisp)
    tdisp(1) = (tdisp(1)*loisem())/lor8em()
    if (tdisp(1) .lt. tneces) then
        vali(1) = tdisp(1)
        vali(2) = tneces
        call utmess('F', 'PREPOST5_8', ni=2, vali=vali)
!
    else
        call wkvect('&&ACMATA.VECT_NORMA3', 'V V R', tneces, jvecno)
        call wkvect('&&ACMATA.VECT_NORMA4', 'V V R', tnecno, jnorma)
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
    call wkvect('&&ACMATA.VECT_NORMA1', 'V V R', 27, jvecn1)
    call wkvect('&&ACMATA.VECT_TANGU1', 'V V R', 27, jvecu1)
    call wkvect('&&ACMATA.VECT_TANGV1', 'V V R', 27, jvecv1)
    call wkvect('&&ACMATA.VECT_NORMA2', 'V V R', 27, jvecn2)
    call wkvect('&&ACMATA.VECT_TANGU2', 'V V R', 27, jvecu2)
    call wkvect('&&ACMATA.VECT_TANGV2', 'V V R', 27, jvecv2)
!
    call wkvect('&&ACMATA.VECTPG1', 'V V R', 18*nbordr, jvpg1)
    call wkvect('&&ACMATA.VECTPG2', 'V V R', 18*nbordr, jvpg2)
!
    epsilo = 1.0d-7
    pi = r8pi()
!
!
! PROJECTION DE L'HISTORIQUE DU CISAILLEMENT DANS UN PLAN.
!
    nbvec = 209
    dectau = 0
!
    call taurlo(nbvec, jvectn, jvectu, jvectv, nbordr, &
                kwork, sompgw, jrwork, tspaq, ipg, &
                dectau, jvecno, jnorma)
!
!
! CALCMAX DES DELTA_TAU MAX ET DU VECTEUR NORMAL ASSOCIE POUR
! LE PE GAUSS COURANT DE LA MAILLE COURANTE.
!
! 1/ RA ZERO DU VECTEUR DE TRAVAIL CONTENANT LES VALEURS DE
!    DAU POUR UN POINT DE GAUSS ET DU VECTEUR DE TRAVAIL
!    PANT DE POINTER SUR LE VECTEUR NORMAL ASSOCIE.
!
! 2/ CDU RAYON CIRCONSCRIT
!
    call raycir(jvecno, jdtaum, jresun, nbordr, nbvec, &
                nommet)
!
! 3/ CDU 1ER MAX DES DELTA_TAU ET DU VECTEUR NORMAL ASSOCIE
!
    dtaum(1) = 0.0d0
    dtaum(2) = 0.0d0
    mnmax(1) = 1
    mnmax(2) = 1
!
    do i = 1, nbvec
        if (zr(jdtaum+(i-1)) .gt. epsilo) then
            if ((zr(jdtaum+(i-1))-dtaum(1))/zr(jdtaum+(i-1)) .gt. epsilo) then
                dtaum(2) = dtaum(1)
                mnmax(2) = mnmax(1)
                dtaum(1) = zr(jdtaum+(i-1))
                mnmax(1) = i
            end if
            if (((zr(jdtaum+(i-1))-dtaum(2))/zr(jdtaum+(i-1)) .gt. epsilo) .and. &
                (i .ne. mnmax(1))) then
                dtaum(2) = zr(jdtaum+(i-1))
                mnmax(2) = i
            end if
        end if
    end do
!
! 4/ 1-ER RAFFINEMENT CONCERNANT LA DETERMINATION DU VECTEUR NORMAL
!    EAX DES DELTA_TAU (DETERMINATION DU VECTEUR NORMAL A 2
!    DPRES).
!
    phydro = 0.0d0
    phydrm = 0.0d0
    dim = 27
!
    do k = 1, 2
        norm(k) = 0.0d0
        normax(k) = 0.0d0
        snorm(k) = 0.0d0
        epnorm(k) = 0.0d0
        epnmax(k) = 0.0d0
        sepnmx(k) = 0.0d0
!
        nxm(k) = zr(jvectn+(mnmax(k)-1)*3)
        nym(k) = zr(jvectn+(mnmax(k)-1)*3+1)
        nzm(k) = zr(jvectn+(mnmax(k)-1)*3+2)
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
            call taurlo(nbvec, jvecn2, jvecu2, jvecv2, nbordr, &
                        kwork, sompgw, jrwork, tspaq, ipg, &
                        dectau, jvpg2, jnorma)
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
            call taurlo(nbvec, jvecn2, jvecu2, jvecv2, nbordr, &
                        kwork, sompgw, jrwork, tspaq, ipg, &
                        dectau, jvpg2, jnorma)
        end if
!
! 4-1/E A ZERO DU VECTEUR DE TRAVAIL CONTENANT LES VALEURS DE
!     TAU POUR UN POINT DE GAUSS ET DU VECTEUR DE TRAVAIL
!     TANT DE POINTER SUR LE VECTEUR NORMAL ASSOCIE.
!
! 4-2/L DU RAYON CIRCONSCRIT
!
        call raycir(jvpg2, jdtaum, jresun, nbordr, nbvec, &
                    nommet)
!
! 4-3/L DU 2EME MAX DES DELTA_TAU ET DU VECTEUR NORMAL ASSOCIE
!
        dtaum(k) = 0.0d0
        mnmax(k) = 1
!
        do i = 1, nbvec
            if (zr(jdtaum+(i-1)) .gt. dtaum(k)) then
                dtaum(k) = zr(jdtaum+(i-1))
                mnmax(k) = i
            end if
        end do
!
! 5/ 2-EXIME RAFFINEMENT CONCERNANT LA DETERMINATION DU VECTEUR NORMAL
!    EAX DES DELTA_TAU (DETERMINATION DU VECTEUR NORMAL A 1
!    DRES).
!
        nxm(k) = zr(jvecn2+(mnmax(k)-1)*3)
        nym(k) = zr(jvecn2+(mnmax(k)-1)*3+1)
        nzm(k) = zr(jvecn2+(mnmax(k)-1)*3+2)
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
            gamma = 1.0d0*r8dgrd()
            dphi2 = 60.0d0*r8dgrd()
            phi0 = 0.0d0
            n = 0
!
            call vecnuv(1, 6, gamma, phi0, dphi2, &
                        n, 1, dim, zr(jvecn1), zr(jvecu1), &
                        zr(jvecv1))
!
            gamma = 0.0d0
            phi0 = pi
!
            call vecnuv(1, 1, gamma, phi0, dphi2, &
                        n, 1, dim, zr(jvecn1), zr(jvecu1), &
                        zr(jvecv1))
!
            nbvec = 7
            call taurlo(nbvec, jvecn1, jvecu1, jvecv1, nbordr, &
                        kwork, sompgw, jrwork, tspaq, ipg, &
                        dectau, jvpg1, jnorma)
        else
            dgam2 = 1.0d0*r8dgrd()
            dphi2 = dgam2/sin(gammam)
            n = 0
            do j = 1, 3
                gamma = gammam+(j-2)*dgam2
!
                call vecnuv(1, 3, gamma, phim, dphi2, &
                            n, 2, dim, zr(jvecn1), zr(jvecu1), &
                            zr(jvecv1))
!
            end do
!
            nbvec = 9
            call taurlo(nbvec, jvecn1, jvecu1, jvecv1, nbordr, &
                        kwork, sompgw, jrwork, tspaq, ipg, &
                        dectau, jvpg1, jnorma)
        end if
!
! 5-1/E A ZERO DU VECTEUR DE TRAVAIL CONTENANT LES VALEURS DE
!     TAU POUR UN POINT DE GAUSS ET DU VECTEUR DE TRAVAIL
!     TANT DE POINTER SUR LE VECTEUR NORMAL ASSOCIE.
!
!
! 5-2/L DU RAYON CIRCONSCRIT
!
        call raycir(jvpg1, jdtaum, jresun, nbordr, nbvec, &
                    nommet)
!
! 5-3/L DU 2EME MAX DES DELTA_TAU ET DU VECTEUR NORMAL ASSOCIE
!
        dtaum(k) = 0.0d0
        mnmax(k) = 1
!
        do i = 1, nbvec
            if (zr(jdtaum+(i-1)) .gt. dtaum(k)) then
                dtaum(k) = zr(jdtaum+(i-1))
                mnmax(k) = i
            end if
        end do
!
        nxm(k) = zr(jvecn1+(mnmax(k)-1)*3)
        nym(k) = zr(jvecn1+(mnmax(k)-1)*3+1)
        nzm(k) = zr(jvecn1+(mnmax(k)-1)*3+2)
        gammam = atan2(sqrt(abs(1.0d0-nzm(k)**2)), nzm(k))
!
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
            gamma = 0.5d0*r8dgrd()
            dphi2 = 60.0d0*r8dgrd()
            phi0 = 0.0d0
            n = 0
!
            call vecnuv(1, 6, gamma, phi0, dphi2, &
                        n, 1, dim, zr(jvecn1), zr(jvecu1), &
                        zr(jvecv1))
!
            gamma = 0.0d0
            phi0 = pi
!
            call vecnuv(1, 1, gamma, phi0, dphi2, &
                        n, 1, dim, zr(jvecn1), zr(jvecu1), &
                        zr(jvecv1))
!
            nbvec = 7
            call taurlo(nbvec, jvecn1, jvecu1, jvecv1, nbordr, &
                        kwork, sompgw, jrwork, tspaq, ipg, &
                        dectau, jvpg1, jnorma)
        else
            dgam2 = 0.5d0*r8dgrd()
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
            call taurlo(nbvec, jvecn2, jvecu2, jvecv2, nbordr, &
                        kwork, sompgw, jrwork, tspaq, ipg, &
                        dectau, jvpg2, jnorma)
        end if
!
!
! 5/ 3-IEME RAFFINEMENT CONCERNANT LA DETERMINATION DU VECTEUR NORMAL
!    EAX DES DELTA_TAU (DETERMINATION DU VECTEUR NORMAL A 1
!    DRES)
!
! 5-2/L DU RAYON CIRCONSCRIT
!
        call raycir(jvpg2, jdtaum, jresun, nbordr, nbvec, &
                    nommet)
!
! 5-3/L DU 2EME MAX DES DELTA_TAU ET DU VECTEUR NORMAL ASSOCIE
!
        dtaum(k) = 0.0d0
        mnmax(k) = 1
!
        do i = 1, nbvec
            if (zr(jdtaum+(i-1)) .gt. dtaum(k)) then
                dtaum(k) = zr(jdtaum+(i-1))
                mnmax(k) = i
            end if
        end do
!
        nxm(k) = zr(jvecn2+(mnmax(k)-1)*3)
        nym(k) = zr(jvecn2+(mnmax(k)-1)*3+1)
        nzm(k) = zr(jvecn2+(mnmax(k)-1)*3+2)
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
            gamma = 0.25d0*r8dgrd()
            dphi2 = 60.0d0*r8dgrd()
            phi0 = 0.0d0
            n = 0
!
            call vecnuv(1, 6, gamma, phi0, dphi2, &
                        n, 1, dim, zr(jvecn1), zr(jvecu1), &
                        zr(jvecv1))
!
            gamma = 0.0d0
            phi0 = pi
!
            call vecnuv(1, 1, gamma, phi0, dphi2, &
                        n, 1, dim, zr(jvecn1), zr(jvecu1), &
                        zr(jvecv1))
!
            nbvec = 7
            call taurlo(nbvec, jvecn1, jvecu1, jvecv1, nbordr, &
                        kwork, sompgw, jrwork, tspaq, ipg, &
                        dectau, jvpg1, jnorma)
        else
            dgam2 = 1.0d0*r8dgrd()
            dphi2 = dgam2/sin(gammam)
            n = 0
            do j = 1, 3
                gamma = gammam+(j-2)*dgam2
!
                call vecnuv(1, 3, gamma, phim, dphi2, &
                            n, 2, dim, zr(jvecn1), zr(jvecu1), &
                            zr(jvecv1))
!
            end do
!
            nbvec = 9
            call taurlo(nbvec, jvecn1, jvecu1, jvecv1, nbordr, &
                        kwork, sompgw, jrwork, tspaq, ipg, &
                        dectau, jvpg1, jnorma)
        end if
!
! CALCLA CONTRAINTE NORMALE MAXIMALE SUR LE PLAN CRITIQUE,
! DE LRAINTE NORMALE MOYENNE SUR LE PLAN CRITIQUE,
! DE LRMATION NORMALE MAXIMALE SUR LE PLAN CRITIQUE,
! DE LRMATION NORMALE MOYENNE SUR LE PLAN CRITIQUE.
!
!          CALL RCVALE(NOMMAT,'ELAS',0,'        ',R8B,1,'E       ',
!      &               VALE,ICODRE,0)
!          IF (ICODRE .EQ. 1) THEN
!             CALL UTMESS('F','PREPOST_11')
!          ENDIF
!          CALL RCVALE(NOMMAT,'ELAS',0,'        ',R8B,1,'NU      ',
!      &               VALNU,ICODRE,0)
!          IF (ICODRE .EQ. 1) THEN
!             CALL UTMESS('F','PREPOST_12')
!          ENDIF
!          C1 = (1+VALNU)/VALE
!          C2 = VALNU/VALE
!
!
! 5/ 3-IEME RAFFINEMENT CONCERNANT LA DETERMINATION DU VECTEUR NORMAL
!    EAX DES DELTA_TAU (DETERMINATION DU VECTEUR NORMAL A 1
!    DRES)
!
! 5-2/L DU RAYON CIRCONSCRIT
!
        call raycir(jvpg1, jdtaum, jresun, nbordr, nbvec, &
                    nommet)
!
! 5-3/L DU 2EME MAX DES DELTA_TAU ET DU VECTEUR NORMAL ASSOCIE
!
        dtaum(k) = 0.0d0
        mnmax(k) = 1
!
        do i = 1, nbvec
            if (zr(jdtaum+(i-1)) .gt. dtaum(k)) then
                dtaum(k) = zr(jdtaum+(i-1))
                mnmax(k) = i
            end if
        end do
!
        nxm(k) = zr(jvecn1+(mnmax(k)-1)*3)
        nym(k) = zr(jvecn1+(mnmax(k)-1)*3+1)
        nzm(k) = zr(jvecn1+(mnmax(k)-1)*3+2)
!
        do i = 1, nbordr
            decal = 18
!           ADR = (J-1)*TSPAQ+KWORK*SOMPGW*6+(IPG-1)*6
            adrs = (i-1)*tspaq+kwork*sompgw*decal+(ipg-1)*decal
            sixx = zr(jrwork+adrs+0)
            siyy = zr(jrwork+adrs+1)
            sizz = zr(jrwork+adrs+2)
            sixy = zr(jrwork+adrs+3)
            sixz = zr(jrwork+adrs+4)
            siyz = zr(jrwork+adrs+5)
!
            epsxx = zr(jrwork+adrs+6)
            epsyy = zr(jrwork+adrs+7)
            epszz = zr(jrwork+adrs+8)
            epsxy = zr(jrwork+adrs+9)
            epsxz = zr(jrwork+adrs+10)
            epsyz = zr(jrwork+adrs+11)
!
!
! CALCLA PRESSION HYDROSTATIQUE MAXIMALE = Max_t(1/3 Tr[SIG])
!
            if (k .lt. 2) then
!
! ON C PHYDRM QU'UNE FOIS, PARCE QUE LA PRESSION HYDROSTATIQUE
! EST ANTE PAR RAPPORT AU vect_n.
!
                phydro = (sixx+siyy+sizz)/3.0d0
!
                if (phydro .gt. phydrm) then
                    phydrm = phydro
                end if
            end if
!
!             EPSXX = C1*SIXX - C2*(SIXX + SIYY + SIZZ)
!             EPSYY = C1*SIYY - C2*(SIXX + SIYY + SIZZ)
!             EPSZZ = C1*SIZZ - C2*(SIXX + SIYY + SIZZ)
!             EPSXY = C1*SIXY
!             EPSXZ = C1*SIXZ
!             EPSYZ = C1*SIYZ
!
! CALCvect_F = [SIG].vect_n
!
            fxm(k) = sixx*nxm(k)+sixy*nym(k)+sixz*nzm(k)
            fym(k) = sixy*nxm(k)+siyy*nym(k)+siyz*nzm(k)
            fzm(k) = sixz*nxm(k)+siyz*nym(k)+sizz*nzm(k)
!
! CALCNORM = vect_F.vect_n
!
            norm(k) = fxm(k)*nxm(k)+fym(k)*nym(k)+fzm(k)*nzm(k)
!
            if (abs(norm(k)) .gt. normax(k)) then
                normax(k) = norm(k)
            end if
!
            snorm(k) = snorm(k)+norm(k)
!
! CALCvect_EPS = [EPS].vect_n
!
            epsxm(k) = epsxx*nxm(k)+epsxy*nym(k)+epsxz*nzm(k)
            epsym(k) = epsxy*nxm(k)+epsyy*nym(k)+epsyz*nzm(k)
            epszm(k) = epsxz*nxm(k)+epsyz*nym(k)+epszz*nzm(k)
!
! CALCEPSILON NORMALE = vect_EPS.vect_n
!
            epnorm(k) = epsxm(k)*nxm(k)+epsym(k)*nym(k)+epszm(k)*nzm(k)
!
            if (abs(epnorm(k)) .gt. epnmax(k)) then
                epnmax(k) = epnorm(k)
            end if
!
            sepnmx(k) = sepnmx(k)+epnorm(k)
        end do
!
        normoy(k) = snorm(k)/nbordr
        epnmoy(k) = sepnmx(k)/nbordr
!
!
    end do
!
! CONSON D'UN CHAM_ELEM SIMPLE PUIS D'UN CHAM_ELEM CONTENANT
! POURE POINT DE GAUSS DE CHAQUE MAILLE MAX DE DTAU_MAX ET LE
! VECTRMAL ASSOCIE.
!
!
    vrespc(1) = dtaum(1)
    vrespc(2) = nxm(1)
    vrespc(3) = nym(1)
    vrespc(4) = nzm(1)
    vrespc(5) = normax(1)
    vrespc(6) = normoy(1)
    vrespc(7) = epnmax(1)
    vrespc(8) = epnmoy(1)
    vrespc(9) = 0.0d0
    vrespc(10) = 0.0d0
    vrespc(11) = 0.0d0
    vrespc(12) = dtaum(2)
    vrespc(13) = nxm(2)
    vrespc(14) = nym(2)
    vrespc(15) = nzm(2)
    vrespc(16) = normax(2)
    vrespc(17) = normoy(2)
    vrespc(18) = epnmax(2)
    vrespc(19) = epnmoy(2)
    vrespc(20) = 0.0d0
    vrespc(21) = 0.0d0
    vrespc(22) = 0.0d0
    vrespc(23) = 0.0d0
    vrespc(24) = 0.0d0
!
!
    call jedetr('&&ACMATA.DTAU_MAX')
    call jedetr('&&ACMATA.RESU_N')
!
    call jedetr('&&ACMATA.VECT_NORMA')
    call jedetr('&&ACMATA.VECT_TANGU')
    call jedetr('&&ACMATA.VECT_TANGV')
!
    call jedetr('&&ACMATA.VECT_NORMA1')
    call jedetr('&&ACMATA.VECT_TANGU1')
    call jedetr('&&ACMATA.VECT_TANGV1')
    call jedetr('&&ACMATA.VECT_NORMA2')
    call jedetr('&&ACMATA.VECT_TANGU2')
    call jedetr('&&ACMATA.VECT_TANGV2')
    call jedetr('&&ACMATA.VECTPG1')
    call jedetr('&&ACMATA.VECTPG2')
!
    call jedetr('&&ACMATA.VECT_NORMA3')
    call jedetr('&&ACMATA.VECT_NORMA4')
    call jedema()
!
!
!     call jedetr('&&ACMATA.VECTNOD')
!     call jedetr('&&ACMATA.NORMALD')
!
end subroutine
