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
subroutine avplcr(nbvec, vectn, vectu, vectv, nbordr, &
                  kwork, somnow, vwork, tdisp, tspaq, &
                  i, nomcri, nomfor, grdvie, forvie, &
                  fordef, fatsoc, proaxe, nommat, vala, &
                  coefpa, post, cudomx, nxm, nym, &
                  nzm)
! aslint: disable=W1306,W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterc/r8dgrd.h"
#include "asterfort/avcipr.h"
#include "asterfort/avgrdo.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedetr.h"
#include "asterfort/vecnuv.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbordr, kwork, i, nbvec
    integer(kind=8) :: somnow, tdisp, tspaq
    aster_logical :: fordef, post
    real(kind=8) :: vectn(3*nbvec), vectu(3*nbvec), vectv(3*nbvec)
    real(kind=8) :: vwork(tdisp), fatsoc
    character(len=16) :: nomcri, proaxe, nomfor, forvie, grdvie
    character(len=8) :: nommat
    real(kind=8) :: vala, coefpa
    real(kind=8) :: cudomx, nxm(2), nym(2), nzm(2)
!
! BUT:  POUR LA FATIGUE A AMPLITUDE VARIABLE
!       DETERMINER LE PLAN CRITIQUE OU DOMMAGE EST MAX
! ---------------------------------------------------------------------
! ARGUMENTS :
!  NBVEC   IN  I  : NOMBRE MAX DE VECTEUR(209 POUR LA VERSION ACTUELLE)
!  VECTN    IN  R  : VECTEUR CONTENANT LES COMPOSANTES DES
!                    VECTEURS NORMAUX.
!  VECTU    IN  R  : VECTEUR CONTENANT LES COMPOSANTES DES
!                    VECTEURS u DU PLAN DE CISAILLEMENT.
!  VECTV    IN  R  : VECTEUR CONTENANT LES COMPOSANTES DES
!                    VECTEURS v DU PLAN DE CISAILLEMENT.
!  NBORDR   IN  I  : NOMBRE DE NUMEROS D'ORDRE.
!  KWORK    IN  I  : KWORK = 0 ON TRAITE LA 1ERE MAILLE DU PAQUET
!                              MAILLES OU LE 1ER NOEUD DU PAQUET DE
!                              NOEUDS;
!                    KWORK = 1 ON TRAITE LA IEME (I>1) MAILLE DU PAQUET
!                              MAILLES OU LE IEME NOEUD DU PAQUET
!                              DE NOEUDS.
!  SOMMW    IN  I  : SOMME DES POINTS DE GAUSS OU DES NOEUDS DES N
!                    MAILLES PRECEDANT LA MAILLE COURANTE.
!  VWORK    IN  R  : VECTEUR DE TRAVAIL CONTENANT
!                    L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                    ATTACHES A CHAQUE POINT DE GAUSS OU NOEUD DES
!                    MAILLE OU NOEUD DU <<PAQUET>> DE MAILLES OU
!                    DE NOEUDS.
!  TDISP    IN  I  : DIMENSION DU VECTEUR VWORK
!  TSPAQ    IN  I  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE MAILLES
!                    OU DE NOEUDS COURANT.
!  I        IN  I  : IEME POINT DE GAUSS OU IEME NOEUD.
!  NOMCRI   IN  K16: NOM DU CRITERE D'ENDOMMAGEMENT PAR FATIGUE.
!  FATSOC   IN  R  : COEFFICIENT PERMETTANT D'UTILISER LES MEMES
!                    ROUTINES POUR LE TRAITEMENT DES CONTRAINTES ET
!                    DES DEFORMATIONS.
!  PROAXE    IN   K16: TYPE DE PROJECTION (UN OU DEUX AXES).
!  NOMMAT   IN   K  : NOM DU MATERIAU.
!  VALA     IN   R  : VALEUR DU PARAMETRE a ASSOCIE AU CRITERE.
!  COEFPA   IN   R  : COEFFICIENT DE PASSAGE CISAILLEMENT - UNIAXIAL.
!  VNORMX   OUT  I  : NUMERO DU VECTEUR NORMAL ASSOCIE AU MAX DES CUMULS
!                     DE DOMMAGE.
!  CUDOMX   OUT  R  : VALEUR DU MAX DES CUMULS DE DOMMAGE.
! REMARQUE : CETTE ROUTINE SERT POUR LE TRAITEMENT DES POINTS DE GAUSS
!            ET DES NOEUDS.
! ----------------------------------------------------------------------
    integer(kind=8) :: ncycl(nbvec), nbvec1, nval, ibid
!    integer :: omin(nbvec*(nbordr+2)), omax(nbvec*(nbordr+2))
    integer(kind=8) :: jomin, jomax, jvmin, jvmax
    integer(kind=8) :: vnormx(2), ideb, ifin, n, k, dim, j, kp, nbp
    integer(kind=8) :: nbplan, vnorm(2)
!    real(kind=8) :: vmin(nbvec*(nbordr+2)), vmax(nbvec*(nbordr+2))
    real(kind=8) :: pseuil, gammam, phim, dphi2, epsilo, gamma
    real(kind=8) :: vecn2(3*nbvec), vecu2(3*nbvec), vecv2(3*nbvec)
    real(kind=8) :: vecn1(3*nbvec), vecu1(3*nbvec), vecv1(3*nbvec)
    real(kind=8) :: dgam2, pi, phi0, cudom1, cudom2
    real(kind=8) :: prec
    character(len=8) :: method
!     --------------------------
    epsilo = 1.0d-7
    pi = r8pi()
!
    prec = 100.d0*r8prem()
!
    nbvec1 = 209
!
    method = 'RAINFLOW'
!
    call getvr8(' ', 'DELTA_OSCI', scal=pseuil, nbret=nval)
!
! CONSTRUCTION DU VECTEUR NORMAL SUR UNE DEMI SPHERE
!
!
    call wkvect('&&AVPLCR.VECT_OMIN', 'V V I', nbvec*(nbordr+2), jomin)
    call wkvect('&&AVPLCR.VECT_OMAX', 'V V I', nbvec*(nbordr+2), jomax)
    call wkvect('&&AVPLCR.VECT_VMIN', 'V V R', nbvec*(nbordr+2), jvmin)
    call wkvect('&&AVPLCR.VECT_VMAX', 'V V R', nbvec*(nbordr+2), jvmax)
!
    call avcipr(nbvec1, vectn, vectu, vectv, nbordr, &
                kwork, somnow, vwork, tdisp, tspaq, &
                i, nomcri, nomfor, fordef, fatsoc, &
                proaxe, pseuil, method, ncycl, jvmin, &
                jvmax, jomin, jomax)
!
! REMPACER PAR SUBROUTINE AVGRDO
!
    call avgrdo(nbvec1, nbordr, vectn, vwork, tdisp, &
                kwork, somnow, tspaq, i, nommat, &
                nomcri, nomfor, grdvie, forvie, vala, &
                coefpa, ncycl, jvmin, jvmax, jomin, &
                jomax, post, cudomx, vnorm, nbplan)
!
!
! 9. PREMIER RAFFINEMENT CONCERNANT LA DETERMINATION DU VECTEUR NORMAL
!    CORRESPONDANT AU MAX DES CUMULS DE DOMMAGE.
!
    if ((post) .and. (nbplan .gt. 2)) then
        write (6, *) 'IL EXISTE  PLUS DE 2 PLANS DU MAX DOMMAGE'
    end if
!
!      IF (NBPLAN .EQ. 2) THEN
!
    do kp = 1, 2
        nxm(kp) = vectn((vnorm(kp)-1)*3+1)
        nym(kp) = vectn((vnorm(kp)-1)*3+2)
        nzm(kp) = vectn((vnorm(kp)-1)*3+3)
    end do
!
    do kp = 1, 2
!
        nxm(kp) = vectn((vnorm(kp)-1)*3+1)
        nym(kp) = vectn((vnorm(kp)-1)*3+2)
        nzm(kp) = vectn((vnorm(kp)-1)*3+3)
!
        gammam = atan2(sqrt(abs(1.0d0-nzm(kp)**2)), nzm(kp))
        if (gammam .lt. 0.0d0) then
            gammam = gammam+pi
        end if
!
        if ((abs(nym(kp)) .lt. epsilo) .and. (abs(nxm(kp)) .lt. epsilo)) then
            phim = 0.0d0
        else
            phim = atan2(abs(nym(kp)), nxm(kp))
        end if
        if (phim .lt. 0.0d0) then
            phim = phim+pi
        end if
!
        if (abs(gammam) .lt. epsilo) then
            gamma = 5.0d0*r8dgrd()
            dphi2 = 60.0d0*r8dgrd()
            ideb = 1
            ifin = 6
            n = 0
            k = 1
            dim = 27
            phi0 = 0.0d0
!
            call vecnuv(ideb, ifin, gamma, phi0, dphi2, &
                        n, k, dim, vecn2, vecu2, &
                        vecv2)
            gamma = 0.0d0
            phi0 = pi
            ideb = 1
            ifin = 1
            k = 1
!
            call vecnuv(ideb, ifin, gamma, phi0, dphi2, &
                        n, k, dim, vecn2, vecu2, &
                        vecv2)
!
! 9.1 PROJECTION DE L'HISTORIQUE DU CISAILLEMENT SUR UN PLAN
!
            nbvec1 = 7
!
            call avcipr(nbvec1, vecn2, vecu2, vecv2, nbordr, &
                        kwork, somnow, vwork, tdisp, tspaq, &
                        i, nomcri, nomfor, fordef, fatsoc, &
                        proaxe, pseuil, method, ncycl, jvmin, &
                        jvmax, jomin, jomax)
!
!
        else
            dgam2 = 2.0d0*r8dgrd()
            dphi2 = dgam2/sin(gammam)
            n = 0
            k = 2
            dim = 27
            ideb = 1
            ifin = 3
            do j = 1, 3
                gamma = gammam+(j-k)*dgam2
                call vecnuv(ideb, ifin, gamma, phim, dphi2, &
                            n, k, dim, vecn2, vecu2, &
                            vecv2)
            end do
!
            nbvec1 = 9
!
            call avcipr(nbvec1, vecn2, vecu2, vecv2, nbordr, &
                        kwork, somnow, vwork, tdisp, tspaq, &
                        i, nomcri, nomfor, fordef, fatsoc, &
                        proaxe, pseuil, method, ncycl, jvmin, &
                        jvmax, jomin, jomax)
!
        end if
!
! REMPACER PAR SUBROUTINE AVGRDO
!
        call avgrdo(nbvec1, nbordr, vecn2, vwork, tdisp, &
                    kwork, somnow, tspaq, i, nommat, &
                    nomcri, nomfor, grdvie, forvie, vala, &
                    coefpa, ncycl, jvmin, jvmax, jomin, &
                    jomax, post, cudomx, vnormx, ibid)
!
!
!
! 10. SECOND RAFFINEMENT CONCERNANT LA DETERMINATION DU VECTEUR NORMAL
!     CORRESPONDANT AU MAX DES CUMULS DE DOMMAGE.
!        C
        nxm(kp) = vecn2((vnormx(kp)-1)*3+1)
        nym(kp) = vecn2((vnormx(kp)-1)*3+2)
        nzm(kp) = vecn2((vnormx(kp)-1)*3+3)
!
        gammam = atan2(sqrt(abs(1.0d0-nzm(kp)**2)), nzm(kp))
        if (gammam .lt. 0.0d0) then
            gammam = gammam+pi
        end if
!
        if ((abs(nym(kp)) .lt. epsilo) .and. (abs(nxm(kp)) .lt. epsilo)) then
            phim = 0.0d0
        else
            phim = atan2(abs(nym(kp)), nxm(kp))
        end if
        if (phim .lt. 0.0d0) then
            phim = phim+pi
        end if
!
        if (abs(gammam) .lt. epsilo) then
            gamma = 5.0d0*r8dgrd()
            dphi2 = 60.0d0*r8dgrd()
            ideb = 1
            ifin = 6
            n = 0
            k = 1
            dim = 27
            phi0 = 0.0d0
            call vecnuv(ideb, ifin, gamma, phi0, dphi2, &
                        n, k, dim, vecn1, vecu1, &
                        vecv1)
!
            gamma = 0.0d0
            phi0 = pi
            ideb = 1
            ifin = 1
            k = 1
            call vecnuv(ideb, ifin, gamma, phi0, dphi2, &
                        n, k, dim, vecn1, vecu1, &
                        vecv1)
!
! 10.1 PROJECTION DE L'HISTORIQUE DU CISAILLEMENT SUR UN PLAN
!
            nbvec1 = 7
!
!
            call avcipr(nbvec1, vecn1, vecu1, vecv1, nbordr, &
                        kwork, somnow, vwork, tdisp, tspaq, &
                        i, nomcri, nomfor, fordef, fatsoc, &
                        proaxe, pseuil, method, ncycl, jvmin, &
                        jvmax, jomin, jomax)
!
        else
            dgam2 = 1.0d0*r8dgrd()
            dphi2 = dgam2/sin(gammam)
            n = 0
            k = 2
            dim = 27
            ideb = 1
            ifin = 3
            do j = 1, 3
                gamma = gammam+(j-k)*dgam2
                call vecnuv(ideb, ifin, gamma, phim, dphi2, &
                            n, k, dim, vecn1, vecu1, &
                            vecv1)
            end do
!
            nbvec1 = 9
!
            call avcipr(nbvec1, vecn1, vecu1, vecv1, nbordr, &
                        kwork, somnow, vwork, tdisp, tspaq, &
                        i, nomcri, nomfor, fordef, fatsoc, &
                        proaxe, pseuil, method, ncycl, jvmin, &
                        jvmax, jomin, jomax)
        end if
!
! REMPACER PAR SUBROUTINE AVGRDO
!
        call avgrdo(nbvec1, nbordr, vecn1, vwork, tdisp, &
                    kwork, somnow, tspaq, i, nommat, &
                    nomcri, nomfor, grdvie, forvie, vala, &
                    coefpa, ncycl, jvmin, jvmax, jomin, &
                    jomax, post, cudomx, vnormx, ibid)
!
!
! 11. 3E RAFFINEMENT CONCERNANT LA DETERMINATION DU VECTEUR NORMAL
!     CORRESPONDANT AU MAX DES CUMULS DE DOMMAGE.
!        C
        nxm(kp) = vecn1((vnormx(kp)-1)*3+1)
        nym(kp) = vecn1((vnormx(kp)-1)*3+2)
        nzm(kp) = vecn1((vnormx(kp)-1)*3+3)
!
        gammam = atan2(sqrt(abs(1.0d0-nzm(kp)**2)), nzm(kp))
        if (gammam .lt. 0.0d0) then
            gammam = gammam+pi
        end if
!
        if ((abs(nym(kp)) .lt. epsilo) .and. (abs(nxm(kp)) .lt. epsilo)) then
            phim = 0.0d0
        else
            phim = atan2(abs(nym(kp)), nxm(kp))
        end if
        if (phim .lt. 0.0d0) then
            phim = phim+pi
        end if
!
        if (abs(gammam) .lt. epsilo) then
            gamma = 5.0d0*r8dgrd()
            dphi2 = 60.0d0*r8dgrd()
            ideb = 1
            ifin = 6
            n = 0
            k = 1
            dim = 27
            phi0 = 0.0d0
            call vecnuv(ideb, ifin, gamma, phi0, dphi2, &
                        n, k, dim, vecn2, vecu2, &
                        vecv2)
!
            gamma = 0.0d0
            phi0 = pi
            ideb = 1
            ifin = 1
            k = 1
            call vecnuv(ideb, ifin, gamma, phi0, dphi2, &
                        n, k, dim, vecn2, vecu2, &
                        vecv2)
!
! 11.1 PROJECTION DE L'HISTORIQUE DU CISAILLEMENT SUR UN PLAN
!
            nbvec1 = 7
!
!
            call avcipr(nbvec1, vecn2, vecu2, vecv2, nbordr, &
                        kwork, somnow, vwork, tdisp, tspaq, &
                        i, nomcri, nomfor, fordef, fatsoc, &
                        proaxe, pseuil, method, ncycl, jvmin, &
                        jvmax, jomin, jomax)
!
        else
            dgam2 = 0.5d0*r8dgrd()
            dphi2 = dgam2/sin(gammam)
            n = 0
            k = 2
            dim = 27
            ideb = 1
            ifin = 3
            do j = 1, 3
                gamma = gammam+(j-k)*dgam2
                call vecnuv(ideb, ifin, gamma, phim, dphi2, &
                            n, k, dim, vecn2, vecu2, &
                            vecv2)
            end do
!
            nbvec1 = 9
!
            call avcipr(nbvec1, vecn2, vecu2, vecv2, nbordr, &
                        kwork, somnow, vwork, tdisp, tspaq, &
                        i, nomcri, nomfor, fordef, fatsoc, &
                        proaxe, pseuil, method, ncycl, jvmin, &
                        jvmax, jomin, jomax)
        end if
!
! REMPACER PAR SUBROUTINE AVGRDO
!
        call avgrdo(nbvec1, nbordr, vecn2, vwork, tdisp, &
                    kwork, somnow, tspaq, i, nommat, &
                    nomcri, nomfor, grdvie, forvie, vala, &
                    coefpa, ncycl, jvmin, jvmax, jomin, &
                    jomax, post, cudomx, vnormx, ibid)
!
! 12. 4E RAFFINEMENT CONCERNANT LA DETERMINATION DU VECTEUR NORMAL
!     CORRESPONDANT AU MAX DES CUMULS DE DOMMAGE.
!        C
        nxm(kp) = vecn2((vnormx(kp)-1)*3+1)
        nym(kp) = vecn2((vnormx(kp)-1)*3+2)
        nzm(kp) = vecn2((vnormx(kp)-1)*3+3)
!
        gammam = atan2(sqrt(abs(1.0d0-nzm(kp)**2)), nzm(kp))
        if (gammam .lt. 0.0d0) then
            gammam = gammam+pi
        end if
!
        if ((abs(nym(kp)) .lt. epsilo) .and. (abs(nxm(kp)) .lt. epsilo)) then
            phim = 0.0d0
        else
            phim = atan2(abs(nym(kp)), nxm(kp))
        end if
        if (phim .lt. 0.0d0) then
            phim = phim+pi
        end if
!
        if (abs(gammam) .lt. epsilo) then
            gamma = 5.0d0*r8dgrd()
            dphi2 = 60.0d0*r8dgrd()
            ideb = 1
            ifin = 6
            n = 0
            k = 1
            dim = 27
            phi0 = 0.0d0
            call vecnuv(ideb, ifin, gamma, phi0, dphi2, &
                        n, k, dim, vecn1, vecu1, &
                        vecv1)
!
            gamma = 0.0d0
            phi0 = pi
            ideb = 1
            ifin = 1
            k = 1
            call vecnuv(ideb, ifin, gamma, phi0, dphi2, &
                        n, k, dim, vecn1, vecu1, &
                        vecv1)
!
! 12.1 PROJECTION DE L'HISTORIQUE DU CISAILLEMENT SUR UN PLAN
!
            nbvec1 = 7
!
!
            call avcipr(nbvec1, vecn1, vecu1, vecv1, nbordr, &
                        kwork, somnow, vwork, tdisp, tspaq, &
                        i, nomcri, nomfor, fordef, fatsoc, &
                        proaxe, pseuil, method, ncycl, jvmin, &
                        jvmax, jomin, jomax)
!
        else
            dgam2 = 0.25d0*r8dgrd()
            dphi2 = dgam2/sin(gammam)
            n = 0
            k = 2
            dim = 27
            ideb = 1
            ifin = 3
            do j = 1, 3
                gamma = gammam+(j-k)*dgam2
                call vecnuv(ideb, ifin, gamma, phim, dphi2, &
                            n, k, dim, vecn1, vecu1, &
                            vecv1)
            end do
!
            nbvec1 = 9
!
            call avcipr(nbvec1, vecn1, vecu1, vecv1, nbordr, &
                        kwork, somnow, vwork, tdisp, tspaq, &
                        i, nomcri, nomfor, fordef, fatsoc, &
                        proaxe, pseuil, method, ncycl, jvmin, &
                        jvmax, jomin, jomax)
        end if
!
! REMPACER PAR SUBROUTINE AVGRDO
!
        call avgrdo(nbvec1, nbordr, vecn1, vwork, tdisp, &
                    kwork, somnow, tspaq, i, nommat, &
                    nomcri, nomfor, grdvie, forvie, vala, &
                    coefpa, ncycl, jvmin, jvmax, jomin, &
                    jomax, post, cudomx, vnormx, nbp)
!  VECTEUR NORMAL ASSOCIE AUX PLAN CRITIQUE  TROUVE
!
        nxm(kp) = vecn1((vnormx(kp)-1)*3+1)
        nym(kp) = vecn1((vnormx(kp)-1)*3+2)
        nzm(kp) = vecn1((vnormx(kp)-1)*3+3)
!
        if (kp .eq. 1) cudom1 = cudomx
        if (kp .eq. 2) cudom2 = cudomx
!
    end do
!
!      ENDIF
    if (abs(cudom1-cudom2) .lt. prec) then
        if ((post) .and. (nbplan .eq. 2)) then
            write (6, *) 'IL EXISTE  2 PLANS DU DOMMAGE MAXIMUM'
        end if
!
    end if
!
    if ((cudom1-cudom2) .gt. prec) then
        if ((post) .and. (nbplan .eq. 2)) then
            write (6, *) 'IL EXISTE  1 PLAN DU DOMMAGE MAXIMUM'
        end if
!
        nxm(2) = nxm(1)
        nym(2) = nym(1)
        nzm(2) = nzm(1)
        cudomx = cudom1
    end if
!
    if ((cudom2-cudom1) .gt. prec) then
        if ((post) .and. (nbplan .eq. 2)) then
            write (6, *) 'IL EXISTE  1 PLAN DU DOMMAGE MAXIMUM'
        end if
!
        nxm(1) = nxm(2)
        nym(1) = nym(2)
        nzm(1) = nzm(2)
        cudomx = cudom2
    end if
!
    call jedetr('&&AVPLCR.VECT_OMIN')
    call jedetr('&&AVPLCR.VECT_OMAX')
    call jedetr('&&AVPLCR.VECT_VMIN')
    call jedetr('&&AVPLCR.VECT_VMAX')
!
!
!
end subroutine
