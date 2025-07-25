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
subroutine wp1mul(lmasse, lamor, lraide, ptorig, tolf, &
                  nitf, nbfreq, mxresf, nprec, resufi, &
                  resufr, solveu)
    implicit none
#include "jeveux.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/mtdscr.h"
#include "asterfort/preres.h"
#include "asterfort/wkvect.h"
#include "asterfort/wp1dft.h"
    integer(kind=8) :: mxresf, nprec
    integer(kind=8) :: lmasse, lamor, lraide, nitf, nbfreq, resufi(mxresf, *)
    complex(kind=8) :: ptorig(3, *)
    real(kind=8) :: tolf, resufr(mxresf, *)
    character(len=19) :: solveu
!     CALCUL DES VALEURS PROPRES COMPLEXES DU SYSTEME QUADRATIQUE
!                         2
!                        L (M) Y + L (C) Y + (K) Y = 0
!     PAR RECHERCHE DES ZEROS DU POLYNOME CARACTERISTIQUE PAR UNE
!     METHODE COURBE D'INTERPOLATION, METHODE A 3 POINTS DE MULLER AVEC
!     DEFLATION
!     ------------------------------------------------------------------
! OUT RESUFR : R : ZERO DU POLYNOME CARACTERISTIQUE
!            POUR IMODE = 1, NBFREQ
!          (IMODE,1) : NUMERO D'ORDRE DU ZERO DU POLYNOME
!          (IMODE,2) : PARTIE IMAGINAIRE DU ZERO DU POLYNOME
!          (IMODE,3) : PARTIE REELLE DU ZERO DU POLYNOME
!     ------------------------------------------------------------------
!     REMARQUE: LES MATRICES QUE L'ON TRAITE SONT SYMETRIQUES ET DONC
!     LES VALEURS PROPRES COMPLEXES SE PRESENTENT PAR PAIRES CONJUGUEES,
!     ON NE RETIENT QUE CELLE A PARTIE IMAGINAIRE POSITIVE ET L'ON
!     ELIMINE L'AUTRE PAR DEFLATION.
!     ------------------------------------------------------------------
!
!
    character(len=1) :: typcst(3), base
    character(len=8) :: nomddl
    character(len=19) :: matdyn, matpre
    character(len=24) :: nmat(3), ndynam
    complex(kind=8) :: res0, res1, res2, h0, h1, lambda, delta, zz, g0, gg, gg1
    complex(kind=8) :: gg2, z0, z1, z2
    integer(kind=8) :: idet0, idet1, idet2, ibid
    real(kind=8) :: det0, det1, det2, rn1, rn2, err, errz, const(6)
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icomb, imode, istu, iter, ldynam, lzero
!-----------------------------------------------------------------------
    data nomddl/'        '/
!     ------------------------------------------------------------------
!
    call jemarq()
    matdyn = '&&WP1MUL.MAT.DYNA'
!
!
!     --- CREATION DE LA MATRICE DYNAMIQUE A VALEUR COMPLEXE ---
    call mtdefs(matdyn, zk24(zi(lmasse+1)), 'V', 'C')
    call mtdscr(matdyn)
    ndynam = matdyn(1:19)//'.&INT'
    call jeveuo(matdyn(1:19)//'.&INT', 'E', ldynam)
!
!      --- DEFINITION DES TYPES DE CONSTANTES ET DES MATRICES ---
    nmat(1) = zk24(zi(lmasse+1))
    nmat(2) = zk24(zi(lamor+1))
    nmat(3) = zk24(zi(lraide+1))
    do icomb = 1, 3
        typcst(icomb) = 'C'
    end do
    const(5) = 1.d0
    const(6) = 0.d0
!
!     --- CREATION D'UN OBJET DE TAVAIL POUR CONTENIR LES ZEROS ---
    call wkvect('&&WP1MUL.ZERO.POLYNOME', 'V V C', nbfreq, lzero)
!
!     --- BOUCLE SUR LE NOMBRE DE MODES DEMANDE ----
    base = 'V'
    matpre = ' '
    do imode = 1, nbfreq
!
        z2 = ptorig(3, imode)
        const(1) = dble(z2*z2)
        const(2) = dimag(z2*z2)
        const(3) = dble(z2)
        const(4) = dimag(z2)
        call mtcmbl(3, typcst, const, nmat, ndynam, &
                    nomddl, ' ', 'ELIM=')
        call preres(solveu, base, ibid, matpre, matdyn, &
                    ibid, 2)
        call wp1dft(ldynam, imode, zc(lzero), z2, res2, &
                    det2, idet2, istu)
!
        z1 = ptorig(2, imode)
        const(1) = dble(z1*z1)
        const(2) = dimag(z1*z1)
        const(3) = dble(z1)
        const(4) = dimag(z1)
        call mtcmbl(3, typcst, const, nmat, ndynam, &
                    nomddl, ' ', 'ELIM=')
        call preres(solveu, base, ibid, matpre, matdyn, &
                    ibid, 2)
        call wp1dft(ldynam, imode, zc(lzero), z1, res1, &
                    det1, idet1, istu)
!
!         --- BOUCLE JUSQU'A LA CONVERGENCE ---
        z0 = ptorig(1, imode)
        do i = 1, nitf
!
            const(1) = dble(z0*z0)
            const(2) = dimag(z0*z0)
            const(3) = dble(z0)
            const(4) = dimag(z0)
            call mtcmbl(3, typcst, const, nmat, ndynam, &
                        nomddl, ' ', 'ELIM=')
            call preres(solveu, base, ibid, matpre, matdyn, &
                        ibid, 2)
            call wp1dft(ldynam, imode, zc(lzero), z0, res0, &
                        det0, idet0, istu)
!
!           --- CALCUL DES COEFFICIENTS DE L'EQUATION DU 2ND DEGRE --
            h0 = z0-z1
            h1 = z1-z2
            lambda = h0/h1
            delta = 1.d0+lambda
            g0 = res2/res0*det2/det0*10.d0**(idet2-idet0)*lambda*lambda-res1/res0*det1/det0*10&
                 &.d0**(idet1-idet0)*delta*delta+lambda+delta
            gg = res2/res0*det2/det0*10.d0**(idet2-idet0)*lambda-res1/res0*det1/det0*10.d0**(id&
                 &et1-idet0)*delta+1.d0
            gg = gg*4.d0*delta*lambda
            gg1 = g0+sqrt(g0*g0-gg)
            gg2 = g0-sqrt(g0*g0-gg)
            rn1 = abs(gg1)
            rn2 = abs(gg2)
            if (rn1 .ge. rn2) then
                zz = -2.d0*delta/gg1
            else
                zz = -2.d0*delta/gg2
            end if
!
!           --- CORRECTION DE LA VALEUR PROPRE POUR CONVERGER VERS LA
!           --- FREQUENCE POSITIVE (SOLUTION CONJUGUEE)
            zz = z0+zz*h0
            zz = dcmplx(dble(zz), abs(dimag(zz)))
!
!           --- CALCUL DE L'ERREUR  ---
            errz = abs(z0)
            err = abs(zz-z0)
            err = sqrt(err/errz)
            if (err .ge. tolf) then
!
!              --- INCREMENTATION PAR PERMUTATION DES OBJETS ---
                res2 = res1
                res1 = res0
                det2 = det1
                det1 = det0
                idet2 = idet1
                idet1 = idet0
                z2 = z1
                z1 = z0
                z0 = zz
            else
                iter = i
                goto 120
            end if
        end do
!
!         --- FIN DES ITERATIONS ---
        iter = -nitf
120     continue
!
        zc(lzero+imode-1) = zz
        resufr(imode, 2) = dimag(zz)
        resufr(imode, 3) = dble(zz)
        resufr(imode, 14) = err
        resufi(imode, 2) = iter
!
    end do
!
! --- MENAGE
    call detrsd('MATR_ASSE', '&&WP1MUL.MAT.DYNA')
    call jedetr('&&WP1MUL.ZERO.POLYNOME')
!
    call jedema()
end subroutine
