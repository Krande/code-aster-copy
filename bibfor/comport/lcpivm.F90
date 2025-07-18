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
subroutine lcpivm(fami, kpg, ksp, mate, compor, &
                  carcri, instam, instap, fm, df, &
                  vim, option, taup, vip, dtaudf, &
                  iret)
!
!
! aslint: disable=
    implicit none
#include "asterf_types.h"
#include "asterfort/calcdp.h"
#include "asterfort/ecpuis.h"
#include "asterfort/gdsmci.h"
#include "asterfort/gdsmhy.h"
#include "asterfort/gdsmin.h"
#include "asterfort/gdsmtg.h"
#include "asterfort/lcpima.h"
#include "asterfort/lcpitg.h"
#include "asterfort/nmcri6.h"
#include "asterfort/rcfonc.h"
#include "asterfort/zerofr.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    integer(kind=8) :: mate, iret, kpg, ksp
    character(len=16) :: compor, option
    character(len=*) :: fami
    real(kind=8) :: instam, instap
    real(kind=8) :: df(3, 3), fm(3, 3)
    real(kind=8) :: vim(8), vip(8)
    real(kind=8) :: taup(6), dtaudf(6, 3, 3)
    real(kind=8) :: carcri(*)
!
! ----------------------------------------------------------------------
!       INTEGRATION DE LA LOI DE COMPORTEMENT PLASTIQUE ISOTROPE
!              EN GRANDES DEFORMATIONS DE TYPE SIMO-MIEHE
!              AINSI QUE SA VERSION VISQUEUSE (LOI SINH)
! ----------------------------------------------------------------------
!
! IN  MATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR : COMPORTEMENT
! IN  CARCRI : PARAMETRES POUR L INTEGRATION DE LA LOI DE COMMPORTEMENT
!                 CARCRI(1) = NOMBRE D ITERATIONS
!                 CARCRI(3) = PRECISION SUR LA CONVERGENCE
! IN  INSTAM : INSTANT PRECEDENT
! IN  INSTAP : INSTANT COURANT
! IN  DF     : INCREMENT DU GRADIENT DE LA TRANSFORMATION
! IN  FM     : GRADIENT DE LA TRANSFORMATION A L INSTANT PRECEDENT
! IN  VIM    : VARIABLES INTERNES A L INSTANT DU CALCUL PRECEDENT
!                 VIM(1)=P (DEFORMATION PLASTIQUE CUMULEE)
!                 VIM(2)=INDICATEUR DE PLASTICITE
!                          0 : ELASTIQUE  1: PLASTIQUE
!                 TRE/3 AVEC E=(ID-BE)/2.D0 EST STOCKE DANS :
!                 VIP(3) POUR LA PLASTICITE OU
!                 VIP(1) POUR L ELASTICITE
! IN  OPTION : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT TAUP   : CONTRAINTES DE KIRCHHOFF A L'INSTANT ACTUEL
! OUT VIP    : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DTAUDF : DERIVEE DE TAU PAR RAPPORT A DF  * (F)T
! OUT IRET   : CODE RETOUR DE  L'INTEGRATION DE LA LDC
!               IRET=0 => PAS DE PROBLEME
!               IRET=1 => DJ<0 ET INTEGRATION IMPOSSIBLE
! ----------------------------------------------------------------------
!  COMMON MATERIAU POUR VON MISES
!
    integer(kind=8) :: jprol, jvale, nbval, itmx
    real(kind=8) :: pm, young, nu, mu, unk, troisk, cother, sigy
    real(kind=8) :: sigm0, epsi0, dt, coefm, rpm, pente, apui, npui
    character(len=1) :: poum
    real(kind=8) :: xap, precr, rprim
!
    common/lcpim/&
     &          pm, young, nu, mu, unk, troisk, cother,&
     &          sigm0, epsi0, dt, coefm, rpm, pente,&
     &          apui, npui, sigy, jprol, jvale, nbval
! ----------------------------------------------------------------------
! COMMON GRANDES DEFORMATIONS SIMO - MIEHE
!
    integer(kind=8) :: ind(3, 3), ind1(6), ind2(6)
    real(kind=8) :: kr(6), rac2, rc(6), id(6, 6)
    real(kind=8) :: bem(6), betr(6), dvbetr(6), eqbetr, trbetr
    real(kind=8) :: jp, dj, jm, dfb(3, 3), mutrbe, tauteq
    real(kind=8) :: djdf(3, 3), dbtrdf(6, 3, 3)
!
    common/rconm6/mutrbe, tauteq
!
    common/gdsmc/&
     &            bem, betr, dvbetr, eqbetr, trbetr,&
     &            jp, dj, jm, dfb,&
     &            djdf, dbtrdf,&
     &            kr, id, rac2, rc, ind, ind1, ind2
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    aster_logical :: resi, rigi, elas
    integer(kind=8) :: i, ij, line, n
    real(kind=8) :: dp, seuil
    real(kind=8) :: rp, pentep, airerp
    real(kind=8) :: em(6), ep(6), trtau, dvbe(6)
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
!
! 1 - INITIALISATION
!-----------------------------------------------------------------------
!
!    DONNEES DE CONTROLE DE L'ALGORITHME
    resi = option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL'
    rigi = option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL'
    elas = option(11:14) .eq. 'ELAS'
    call gdsmin()
!
!    LECTURE DES VARIABLES INTERNES (DEFORMATION PLASTIQUE CUMULEE ET
!                                   -DEFORMATION ELASTIQUE)
    pm = vim(1)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vim(3), b_incx, em, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    call dscal(b_n, rac2, em(4), b_incx)
!
!    CALCUL DES ELEMENTS CINEMATIQUES
    call gdsmci(fm, df, em)
!
!    CARACTERISTIQUES MATERIAU
    if (resi) then
        poum = '+'
    else
        poum = '-'
    end if
    call lcpima(fami, kpg, ksp, poum, mate, &
                compor, instam, instap, carcri, taup, &
                vim)
!
! 2 - RESOLUTION
!-----------------------------------------------------------------------
    if (resi) then
        seuil = mu*eqbetr-rpm
!
        if (seuil .le. 0.d0) then
            dp = 0.d0
            line = 0
!
        else
            line = 1
            if (compor .eq. 'VMIS_ISOT_LINE') then
                dp = seuil/(pente+mu*trbetr)
!
            else if (compor .eq. 'VMIS_ISOT_PUIS') then
                tauteq = mu*eqbetr
                mutrbe = mu*trbetr
                call ecpuis(young, sigy, apui, 1.d0/npui, pm, &
                            0.d0, rp, rprim)
                xap = (tauteq-rp)/mutrbe
                precr = carcri(3)*sigy
                itmx = nint(carcri(1))
!
                call zerofr(0, 'DEKKER', nmcri6, 0.d0, xap, &
                            precr, itmx, dp, iret, n)
                if (iret .ne. 0) goto 999
                call ecpuis(young, sigy, apui, 1.d0/npui, pm, &
                            dp, rp, rprim)
                pente = rprim
            else if (compor .eq. 'VMIS_ISOT_TRAC') then
                call rcfonc('E', 1, jprol, jvale, nbval, &
                            e=young*trbetr/3, nu=nu, p=pm, rp=rp, rprim=pente, &
                            airerp=airerp, sieleq=mu*eqbetr, dp=dp)
            else
! CAS VISQUEUX : CALCUL DE DP PAR RESOLUTION DE
!  FPLAS - (R'+MU TR BEL)DP - PHI(DP) = 0
                call calcdp(carcri, seuil, dt, pente, mu*trbetr, &
                            sigm0, epsi0, coefm, dp, iret)
! DANS LE CAS NON LINEAIRE ON VERFIE QUE L ON A LA BONNE PENTE
                if (compor(10:14) .eq. '_TRAC') then
                    call rcfonc('V', 1, jprol, jvale, nbval, &
                                p=pm+dp, rp=rp, rprim=pentep)
                    do i = 1, nbval
                        if (abs(pente-pentep) .le. 1.d-3) then
                            goto 20
                        else
                            pente = pentep
                            seuil = mu*eqbetr-(rp-pente*dp)
                            call calcdp(carcri, seuil, dt, pente, mu*trbetr, &
                                        sigm0, epsi0, coefm, dp, iret)
                            call rcfonc('V', 1, jprol, jvale, nbval, &
                                        p=vim(1)+dp, rp=rp, rprim=pentep)
                        end if
                    end do
20                  continue
                end if
            end if
        end if
!
! 4 - MISE A JOUR DES CHAMPS
! 4.1 - CONTRAINTE
!
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dvbetr, b_incx, dvbe, b_incy)
        if (line .eq. 1) then
            b_n = to_blas_int(6)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1-dp*trbetr/eqbetr, dvbe, b_incx)
        end if
!
        trtau = (troisk*(jp**2-1)-3.d0*cother*(jp+1.d0/jp))/2.d0
!
        do ij = 1, 6
            taup(ij) = mu*dvbe(ij)+trtau/3.d0*kr(ij)
        end do
!
! 4.2 - CORRECTION HYDROSTATIQUE A POSTERIORI
!
        do ij = 1, 6
            ep(ij) = (kr(ij)-jp**(2.d0/3.d0)*(dvbe(ij)+trbetr/3.d0*kr(ij)))/2.d0
        end do
        call gdsmhy(jp, ep)
!
! 4.3 - P, DEFORMATION ELASTIQUE ET INDICE DE PLASTICITE
!
        vip(1) = pm+dp
        vip(2) = line
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, ep, b_incx, vip(3), b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, 1.d0/rac2, vip(6), b_incx)
    end if
!
! 5 - CALCUL DE LA MATRICE TANGENTE
!
    if (rigi) then
        if (.not. resi) then
            dp = 0.d0
            line = nint(vim(2))
            b_n = to_blas_int(6)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, dvbetr, b_incx, dvbe, b_incy)
        end if
!
        if (elas) line = 0
!
        call gdsmtg()
        call lcpitg(compor, df, line, dp, dvbe, &
                    dtaudf)
    end if
999 continue
end subroutine
