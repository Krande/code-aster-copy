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
subroutine lcdpeq(vind, vinf, rela_comp, nbcomm, cpmono, &
                  nmat, nvi, sig, detot, epsd, &
                  materf, pgl)
!
! aslint: disable=W1306
    implicit none
!     DEFORMATION PLASTIQUE EQUIVALENTE CUMULEE MACROSCOPIQUE
!     POUR LE MONOCRISTAL
!     IN  VIND   :  VARIABLES INTERNES A T
!     IN  VINF   :  VARIABLES INTERNES A T+DT
!          COMP   :  NOM MODELE DE COMPORTEMENT
!          NBCOMM :  INDICES DES COEF MATERIAU
!          CPMONO :  NOMS DES LOIS MATERIAU PAR FAMILLE
!          NMAT   :  DIMENSION MATER
!          VIND   :  VARIABLES INTERNES A T
!          SIG    :  CONTRAINTES A T
!          DETOT  :  INCREMENT DE  DEFORMATION TOTALE OU DF
!          EPSD   :  DEFORMATION TOTALE A T OU F A T
!     VAR  NVI    :  NOMBRE DE VARIABLES INTERNES
!          VINF   :  VARIABLES INTERNES A T+DT
!          MATERF :  COEF MATERIAU
!     ----------------------------------------------------------------
#include "asterfort/lcgrla.h"
#include "asterfort/lcloca.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/lcnrte.h"
#include "asterfort/matinv.h"
#include "asterfort/pk2sig.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    integer(kind=8) :: nvi, nmat, nbcomm(nmat, 3), nbphas, i, iphas, indfv, nuvi, ifa
    integer(kind=8) :: ifl, is, nbfsys, nbsys, nsfv, indpha, indcp, numirr, ns, indtau
    integer(kind=8) :: iei, is3, iv3, iv, irr2
    real(kind=8) :: vind(nvi), vinf(nvi), dvin(nvi), sig(6), granb(6)
    real(kind=8) :: epseq, fv, sigg(6), mus(6), ng(3), lg(3), pgl(3, 3)
    real(kind=8) :: id(3, 3), f(3, 3), fpm(3, 3), fp(3, 3), fe(3, 3), detp
    real(kind=8) :: detot(*), epsd(*), pk2(6), devi(6), endoc, dp, xi, qm(3, 3)
    real(kind=8) :: materf(nmat, 2), rhoirr(12), tau(60)
    real(kind=8) :: rhosat, phisat, dz, roloop(12), fivoid(12), sdp
    character(len=16) :: rela_comp, loca, necoul, nomfam
    character(len=24) :: cpmono(5*nmat+1)
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    blas_int :: b_incx, b_incy, b_n
    common/polycr/irr, decirr, nbsyst, decal, gdef
    data id/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
!
    if (rela_comp .eq. 'MONOCRISTAL') then
        nvi = nvi+3
        if (gdef .eq. 1) then
            nvi = nvi+9
        end if
    end if
!
    epseq = 0.d0
    nbsys = 0
!
    if (rela_comp .eq. 'MONOCRISTAL') then
!
        nbfsys = nbcomm(nmat, 2)
!        NSFV : debut de la famille IFA dans les variables internes
        nsfv = 6
        do ifa = 1, nbfsys
            ifl = nbcomm(ifa, 1)
!            NUECOU=NINT(MATERF(IFL,2))
            nomfam = cpmono(5*(ifa-1)+1) (1:16)
            necoul = cpmono(5*(ifa-1)+3) (1:16)
            call lcmmsg(nomfam, nbsys, 0, pgl, mus, &
                        ng, lg, 0, qm)
            if (necoul .eq. 'MONO_DD_CC_IRRA') then
                b_n = to_blas_int(12)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vind(nsfv+3*nbsys+1), b_incx, rhoirr, b_incy)
                xi = materf(ifl+20, 2)
                irr = 1
                irr2 = 1
            else if (necoul .eq. 'MONO_DD_CFC_IRRA') then
                b_n = to_blas_int(12)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vind(nsfv+3*nbsys+1), b_incx, roloop, b_incy)
                b_n = to_blas_int(12)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vind(nsfv+3*nbsys+13), b_incx, fivoid, b_incy)
                irr = 1
                irr2 = 2
                iei = nbcomm(ifa, 3)
                rhosat = materf(iei+8, 2)
                phisat = materf(iei+9, 2)
                xi = materf(iei+10, 2)
                dz = materf(iei+11, 2)
            else
                irr = 0
                irr2 = 0
            end if
!
            if (irr2 .eq. 1) then
                do is = 1, 12
!                 VARIABLES INTERNES PAR SYSTEME DE GLISSEMENT
                    nuvi = nsfv+3*(is-1)+3
                    dp = vinf(nuvi)
                    rhoirr(is) = rhoirr(is)*exp(-xi*dp)
                end do
                b_n = to_blas_int(12)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, rhoirr, b_incx, vinf(nsfv+3*nbsys+1), b_incy)
!
            end if
!
            if (irr2 .eq. 2) then
                do is = 1, 12
!                 SOMME SUR COPLA(S)
                    sdp = 0.d0
                    do iv = 1, 12
                        is3 = (is-1)/3
                        iv3 = (iv-1)/3
!                    VARIABLES INTERNES PAR SYSTEME DE GLISSEMENT
                        nuvi = nsfv+3*(iv-1)+3
                        dp = vinf(nuvi)
                        if (is3 .eq. iv3) then
                            sdp = sdp+dp
                        end if
                    end do
                    roloop(is) = rhosat+(roloop(is)-rhosat)*exp(-xi*sdp)
                    fivoid(is) = phisat+(fivoid(is)-phisat)*exp(-dz*sdp)
                end do
                b_n = to_blas_int(12)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, roloop, b_incx, vinf(nsfv+3*nbsys+1), b_incy)
                b_n = to_blas_int(12)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, fivoid, b_incx, vinf(nsfv+3*nbsys+13), b_incy)
            end if
!
            nsfv = nsfv+nbsys*3
        end do
!
!
        indtau = nsfv
        if (irr2 .eq. 1) indtau = indtau+12
        if (irr2 .eq. 2) indtau = indtau+24
!        CISSIONS TAU_S
        ns = 0
        do ifa = 1, nbfsys
            ifl = nbcomm(ifa, 1)
            nomfam = cpmono(5*(ifa-1)+1) (1:16)
            call lcmmsg(nomfam, nbsys, 0, pgl, mus, &
                        ng, lg, 0, qm)
            do is = 1, nbsys
!              CALCUL DE LA SCISSION REDUITE =
!              PROJECTION DE SIG SUR LE SYSTEME DE GLISSEMENT
!              TAU      : SCISSION REDUITE TAU=SIG:MUS
                call lcmmsg(nomfam, nbsys, is, pgl, mus, &
                            ng, lg, 0, qm)
                tau(ns+is) = 0.d0
                do i = 1, 6
                    tau(ns+is) = tau(ns+is)+sig(i)*mus(i)
                end do
            end do
            ns = ns+nbsys
        end do
        b_n = to_blas_int(ns)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, tau, b_incx, vinf(indtau+1), b_incy)
!
!
        if (gdef .eq. 1) then
!           ICI CONTRAIREMENT A LCMMON, NVI EST LE NOMBRE TOTAL DE V.I
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vinf(nvi-3-18+1), b_incx, fp, b_incy)
            call matinv('S', 3, fp, fpm, detp)
            f = matmul(reshape(detot(1:9), (/3, 3/)), reshape(epsd(1:9), (/3, 3/)))
            fe = matmul(f, fpm)
!           CALCUL DES CONTRAINTES DE KIRCHOFF
            b_n = to_blas_int(6)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, sig, b_incx, pk2, b_incy)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            call dscal(b_n, sqrt(2.d0), pk2(4), b_incx)
            call pk2sig(3, fe, 1.d0, pk2, sig, &
                        1)
!           LES RACINE(2) ATTENDUES PAR NMCOMP :-)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            call dscal(b_n, sqrt(2.d0), sig(4), b_incx)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -1.d0, id, b_incx, fe, &
                       b_incy)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, fe, b_incx, vinf(nvi-3-18+10), b_incy)
            call lcgrla(fp, devi)
            b_n = to_blas_int(6)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, devi, b_incx, vinf, b_incy)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            call dscal(b_n, sqrt(2.d0), devi(4), b_incx)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -1.d0, id, b_incx, fp, &
                       b_incy)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, fp, b_incx, vinf(nvi-3-18+1), b_incy)
            epseq = lcnrte(devi)
        else
!           V.I. 1 A 6 REPRéSENTE LA DEFORMATION VISCOPLASTIQUE MACRO
            epseq = 0.d0
            do i = 1, 6
                dvin(i) = vinf(i)-vind(i)
                epseq = epseq+dvin(i)*dvin(i)
            end do
            epseq = sqrt(2.0d0/3.0d0*epseq)
        end if
        vinf(nvi-1) = vind(nvi-1)+epseq
!
    else if (rela_comp(1:8) .eq. 'POLYCRIS') then
!
!        V.I. 1 A 6 REPRéSENTE LA DEFORMATION VISCOPLASTIQUE MACRO
        epseq = 0.d0
        do i = 1, 6
            dvin(i) = vinf(i)-vind(i)
            epseq = epseq+dvin(i)*dvin(i)
        end do
        epseq = sqrt(2.0d0/3.0d0*epseq)
        vinf(7) = vind(7)+epseq
!        LOCALISATION
!        RECUPERATION DU NOMBRE DE PHASES
        nbphas = nbcomm(1, 1)
        loca = cpmono(1) (1:16)
!        CALCUL DE  B
        do i = 1, 6
            granb(i) = 0.d0
        end do
        do i = 1, 6
            do iphas = 1, nbphas
                indfv = nbcomm(1+iphas, 3)
                fv = materf(indfv, 2)
                granb(i) = granb(i)+fv*vinf(7+6*(iphas-1)+i)
            end do
        end do
        nuvi = nvi-6*nbphas-1
        do iphas = 1, nbphas
            indfv = nbcomm(1+iphas, 3)
!         RECUPERER L'ORIENTATION DE LA PHASE ET LA PROPORTION
            fv = materf(indfv, 2)
            call lcloca(materf(1, 2), nmat, nbcomm, nbphas, sig, &
                        vinf, iphas, granb, loca, sigg)
            do i = 1, 6
                vinf(nuvi+6*(iphas-1)+i) = sigg(i)
            end do
        end do
!
!        IRRADIATION
        nsfv = 7+6*nbphas
        numirr = 0
        do iphas = 1, nbphas
            indpha = nbcomm(1+iphas, 1)
            nbfsys = nbcomm(indpha, 1)
            indcp = nbcomm(1+iphas, 2)
            do ifa = 1, nbfsys
                necoul = cpmono(indcp+5*(ifa-1)+3) (1:16)
!
                if (necoul .eq. 'MONO_DD_CC_IRRA') then
                    nbsys = 12
                    b_n = to_blas_int(12)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, vind(decirr+numirr+1), b_incx, rhoirr, b_incy)
                    ifl = nbcomm(indpha+ifa, 1)
                    xi = materf(ifl+20, 2)
                    do is = 1, nbsys
!                    VARIABLES INTERNES PAR SYSTEME DE GLISSEMENT
                        nuvi = nsfv+3*(is-1)+3
                        dp = vinf(nuvi)
                        if (irr .eq. 1) then
                            rhoirr(is) = rhoirr(is)*exp(-xi*dp)
                        end if
                    end do
                    b_n = to_blas_int(12)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, rhoirr, b_incx, vinf(decirr+numirr+1), b_incy)
                    numirr = numirr+nbsys
                end if
!
                if (necoul .eq. 'MONO_DD_CFC_IRRA') then
                    nbsys = 12
                    b_n = to_blas_int(12)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, vind(decirr+numirr+1), b_incx, roloop, b_incy)
                    b_n = to_blas_int(12)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, vind(decirr+numirr+13), b_incx, fivoid, b_incy)
                    iei = nbcomm(indpha+ifa, 3)
                    rhosat = materf(iei+8, 2)
                    phisat = materf(iei+9, 2)
                    xi = materf(iei+10, 2)
                    dz = materf(iei+11, 2)
                    do is = 1, nbsys
!                    SOMME SUR COPLA(S)
                        sdp = 0.d0
                        do iv = 1, 12
!                       VARIABLES INTERNES PAR SYSTEME DE GLISSEMENT
                            nuvi = nsfv+3*(iv-1)+3
                            dp = vinf(nuvi)
                            is3 = (is-1)/3
                            iv3 = (iv-1)/3
!                       PARTIE POSITIVE DE ALPHA
                            if (is3 .eq. iv3) then
                                sdp = sdp+dp
                            end if
                        end do
                        roloop(is) = rhosat+(roloop(is)-rhosat)*exp(-xi*sdp)
                        fivoid(is) = phisat+(fivoid(is)-phisat)*exp(-dz*sdp)
                    end do
                    b_n = to_blas_int(12)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, roloop, b_incx, vinf(decirr+numirr+1), b_incy)
                    b_n = to_blas_int(12)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, fivoid, b_incx, vinf(decirr+numirr+13), b_incy)
                    numirr = numirr+nbsys+nbsys
                end if
!
                nsfv = nsfv+nbsys*3
            end do
        end do
    end if
!
    if (epseq .eq. 0.d0) then
        vinf(nvi) = 0.d0
    else
        vinf(nvi) = 1.d0
    end if
!
! --    DEBUT TRAITEMENT DE VENDOCHAB --
! --    CALCUL DES CONTRAINTES SUIVANT QUE LE MATERIAU EST
! --    ENDOMMAGE OU PAS
!
    if (rela_comp(1:9) .eq. 'VENDOCHAB') then
! --    DEBUT TRAITEMENT DE VENDOCHAB --
! --    CALCUL DE DSDE SUIVANT QUE LE MATERIAU EST ENDOMMAGE OU PAS
        endoc = (1.0d0-vinf(9))
        materf(1, 1) = materf(1, 1)*endoc
    end if
!
    if (rela_comp(1:8) .eq. 'HAYHURST') then
! --    DEBUT TRAITEMENT DE HAYHURST --
! --    CALCUL DE DSDE SUIVANT QUE LE MATERIAU EST
! --    ENDOMMAGE OU PAS
        endoc = (1.0d0-vinf(11))
        materf(1, 1) = materf(1, 1)*endoc
! --    FIN   TRAITEMENT DE HAYHURST --
    end if
!
end subroutine
