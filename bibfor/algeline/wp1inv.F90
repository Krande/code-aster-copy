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

subroutine wp1inv(lmasse, lamor, lraide, tolf, nitf, &
                  mxresf, nbfreq, neq, resufi, resufr, &
                  resufk, vecpro, solveu)
    implicit none
#include "jeveux.h"
#include "asterfort/detrsd.h"
#include "asterfort/freqom.h"
#include "asterfort/ggubsc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mcmult.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/mtdscr.h"
#include "asterfort/preres.h"
#include "asterfort/resoud.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: lmasse, lamor, lraide, nitf, nbfreq, neq
    integer(kind=8) :: mxresf
    integer(kind=8) :: resufi(mxresf, *)
    complex(kind=8) :: vecpro(neq, *)
    real(kind=8) :: tolf, resufr(mxresf, *)
    character(len=*) :: resufk(mxresf, *)
    character(len=19) :: solveu
!     CALCUL DES VECTEURS PROPRES COMPLEXES DU SYSTEME QUADRATIQUE
!                         2
!                        L (M) Y + L (C) Y + (K) Y = 0
!     PAR LA METHODE D'ITERATION INVERSE
!     -----------------------------------------------------------------
!
!     -----------------------------------------------------------------
    integer(kind=8) :: lmat(3), ibid
    character(len=1) :: typcst(3), kbid
    character(len=8) :: nomddl
    character(len=19) :: matdyn, k19bid, matass, chcine, criter, matpre
    character(len=24) :: ndynam, nmat(3)
    complex(kind=8) :: rp1, rp, rnorm, rpp
    complex(kind=8) :: czero, cun
    complex(kind=8) :: rmasse, ramor, rraide
    real(kind=8) :: const(6)
!     -----------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: icomb, ieq, imode, iter, jter, lacc1, lacc2
    integer(kind=8) :: ldynam
    real(kind=8) :: dseed, err, err2
    integer(kind=8) :: iret
    complex(kind=8), pointer :: yn_associe_a_xn(:) => null()
!-----------------------------------------------------------------------
    data nomddl/'        '/
!     -----------------------------------------------------------------
    call jemarq()
    czero = dcmplx(0.d0, 0.d0)
    cun = dcmplx(1.d0, 0.d0)
!
!     INIT. OBJETS ASTER
    chcine = ' '
    criter = ' '
    matpre = ' '
    k19bid = ' '
!
!     --- CREATION DES VECTEURS DE TRAVAIL ---
    AS_ALLOCATE(vc=yn_associe_a_xn, size=neq)
    call wkvect('&&WP1INV.XN_MOINS_1     ', 'V V C', neq, lacc1)
    call wkvect('&&WP1INV.YN_MOINS_1     ', 'V V C', neq, lacc2)
!CC   CALL WKVECT('&&WP1INV.M_ORTHOGONALISE','V V C',NEQ*NBFREQ,LMORTH)
!
!     --- CREATION DE LA MATRICE DYNAMIQUE A VALEUR COMPLEXE ---
    neq = zi(lmasse+2)
    matdyn = '&&WP1INV.MATR_DYNA'
    call mtdefs(matdyn, zk24(zi(lmasse+1)), 'V', 'C')
    call mtdscr(matdyn)
    ndynam = matdyn(1:19)//'.&INT'
    call jeveuo(ndynam, 'E', ldynam)
    matass = zk24(zi(ldynam+1))
!
!      --- DEFINITION DES TYPES DE CONSTANTES ET DES MATRICES ---
    lmat(1) = lmasse
    lmat(2) = lamor
    lmat(3) = lraide
    nmat(1) = zk24(zi(lmat(1)+1))
    nmat(2) = zk24(zi(lmat(2)+1))
    nmat(3) = zk24(zi(lmat(3)+1))
    do icomb = 1, 3
        typcst(icomb) = 'C'
    end do
    const(5) = -1.d0
    const(6) = -0.d0
!
!
    do imode = 1, nbfreq
!
        rp = dcmplx(resufr(imode, 3), resufr(imode, 2))
!
!        --- FACTORISATION DE LA MATRICE DYNAMIQUE ---
        const(1) = dble(-rp*rp)
        const(2) = dimag(-rp*rp)
        const(3) = dble(-rp)
        const(4) = dimag(-rp)
        call mtcmbl(3, typcst, const, nmat, ndynam, &
                    nomddl, ' ', 'ELIM=')
        call preres(solveu, 'V', ibid, matpre, matass, &
                    ibid, 2)
!
!        --- CHOIX D'UN VECTEUR INITIAL POUR LA METHODE ---
        dseed = 123457.d0
        call ggubsc(dseed, neq, vecpro(1, imode))
!CC      CALL WP1ORT(NEQ,VECPRO,VECPRO(1,IMODE),ZC(LMORTH),IMODE)
        do ieq = 1, neq
            yn_associe_a_xn(ieq) = rp*vecpro(ieq, imode)
        end do
!
!        --- METHODE D'ITERATION INVERSE ---
        rpp = rp
        do jter = 1, nitf
            iter = jter
!
!           --- M-NORMALISATION DES VECTEURS XN-1 ET YN-1 ---
            call mcmult('ZERO', lmasse, vecpro(1, imode), zc(lacc2), 1, &
                        .false._1)
!
!           --- RENORMALISATION ---
            rnorm = czero
            do ieq = 1, neq
                rnorm = rnorm+dconjg(vecpro(ieq, imode))*zc(lacc2+ieq-1)
            end do
            rnorm = sign(1.d0, dble(rnorm))*cun/sqrt(abs(dble(rnorm)))
            do ieq = 1, neq
                vecpro(ieq, imode) = vecpro(ieq, imode)*rnorm
                yn_associe_a_xn(ieq) = yn_associe_a_xn(ieq)*rnorm
            end do
!
!           --- CONSTITUTION DU SECOND MEMBRE POUR CALCULER XN ---
            call mcmult('ZERO', lamor, vecpro(1, imode), zc(lacc2), 1, &
                        .false._1)
            do ieq = 1, neq
                zc(lacc1+ieq-1) = rp*vecpro(ieq, imode)
            end do
            call mcmult('CUMU', lmasse, zc(lacc1), zc(lacc2), 1, &
                        .false._1)
            call mcmult('CUMU', lmasse, yn_associe_a_xn, zc(lacc2), 1, &
                        .false._1)
!
!           --- RESOLUTION ---
            call resoud(matass, k19bid, solveu, chcine, 1, &
                        k19bid, k19bid, kbid, [0.d0], zc(lacc2), &
                        criter, .false._1, 0, iret)
!
!CC         --- ORTHOGONALISATION DU VECTEUR AVEC LES PRECEDENTS ---
!CC         CALL WP1ORT(NEQ,VECPRO,ZC(LACC2),ZC(LMORTH),IMODE)
!
!           --- CALCUL DE YN ---
            do ieq = 1, neq
                yn_associe_a_xn(ieq) = rp*zc(lacc2+ieq-1)+vecpro(ieq, imode)
            end do
!
!           --- CALCUL DE LA VALEUR PROPRE ---
            call mcmult('ZERO', lmasse, zc(lacc2), zc(lacc1), 1, &
                        .false._1)
            rmasse = czero
            do ieq = 0, neq-1
                rmasse = rmasse+dconjg(zc(lacc2+ieq))*zc(lacc1+ieq)
            end do
            call mcmult('ZERO', lraide, zc(lacc2), zc(lacc1), 1, &
                        .false._1)
            rraide = czero
            do ieq = 0, neq-1
                rraide = rraide+dconjg(zc(lacc2+ieq))*zc(lacc1+ieq)
            end do
            call mcmult('ZERO', lamor, zc(lacc2), zc(lacc1), 1, &
                        .false._1)
            ramor = czero
            do ieq = 0, neq-1
                ramor = ramor+dconjg(zc(lacc2+ieq))*zc(lacc1+ieq)
            end do
!
            rp1 = -ramor+sqrt(ramor*ramor-4.d0*rmasse*rraide)
            rp1 = rp1/(2.d0*rmasse)
            rp1 = dcmplx(dble(rp1), abs(dimag(rp1)))
!
!           --- CONVERGENCE GLOBALE EN VALEUR RELATIVE ---
            err = abs(rpp-rp1)/abs(rp1)
!
!           --- ON ARCHIVE XN  ---
            do ieq = 1, neq
                vecpro(ieq, imode) = zc(lacc2+ieq-1)
            end do
            rpp = rp1
!
            if (err .lt. tolf) then
!             --- ERREUR SUR CHAQUE PARTIE EN ABSOLU ---
                err2 = max(abs(dble(rpp-rp1)), abs(dimag(rpp-rp1)))
                if (err2 .lt. tolf) goto 300
            end if
        end do
!
        iter = -nitf
!
300     continue
!
!        --- ARCHIVAGE DES VALEURS ----
        resufr(imode, 2) = dimag(rp1)*dimag(rp1)
        resufr(imode, 1) = freqom(resufr(imode, 2))
        resufr(imode, 3) = -dble(rp1)/abs(rp1)
        resufr(imode, 15) = err
        resufi(imode, 4) = iter
        resufk(imode, 2) = 'INVERSE_C'
!
!        --- M-NORMALISATION DES VECTEURS TROUVEE ---
        call mcmult('ZERO', lmasse, vecpro(1, imode), zc(lacc2), 1, &
                    .false._1)
!
!        --- RENORMALISATION ---
        rnorm = czero
        do ieq = 1, neq
            rnorm = rnorm+dconjg(vecpro(ieq, imode))*zc(lacc2+ieq-1)
        end do
        rnorm = sign(1.d0, dble(rnorm))*cun/sqrt(abs(dble(rnorm)))
        do ieq = 1, neq
            vecpro(ieq, imode) = vecpro(ieq, imode)*rnorm
        end do
!CC      LMORT = LMORTH + (IMODE-1)*NEQ
!CC      CALL MCMULT('ZERO',LMASSE,VECPRO(1,IMODE),'C',ZC(LMORT),1)
!
    end do
!
!     --- MENAGE ---
    AS_DEALLOCATE(vc=yn_associe_a_xn)
    call jedetr('&&WP1INV.XN_MOINS_1     ')
    call jedetr('&&WP1INV.YN_MOINS_1     ')
    call detrsd('MATR_ASSE', '&&WP1INV.MATR_DYNA')
!
    call jedema()
end subroutine
