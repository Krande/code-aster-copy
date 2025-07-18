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
subroutine cmphii(ck, cm, ndim, nbmod, niter, &
                  xcrit, ceigen, cmod, ndimax, cmat1, &
                  cmat2, cvec, ific)
! aslint: disable=W1306
    implicit none
!
!***********************************************************************
!    P. RICHARD                                 DATE 31/07/91
!-----------------------------------------------------------------------
!  BUT:  < COMPLEXE MODES PROBLEME HERMITIEN ITERATION INVERSE >
!
!   CALCULER LES PREMIERS MODES PROPRES D'UN PROBLEME
!       AUX VALEURS PROPRES A MATRICES RAIDEURS ET MASSE COMPLEXES
!       HERMITTIENNES STOCKEES TRIANGULAIRES SUPERIEURES
!
!                     CK*X= L CM*X
!
!    METHODE D'ITERATION INVERSE
!
!-----------------------------------------------------------------------
!
! CK       /I/: MATRICE RAIDEUR DU PROBLEME
! CM       /I/: MATRICE MASSE DU PROBLEME
! NDIM     /I/: DIMENSION DES MATRICES
! NBMOD    /I/: NOMBRE DE MODES PROPRES DESIRE
! NITER    /I/: NOMBRE MAX D'ITERATIONS PAR MODE
! XCRIT    /I/: TOLERANCE DE COLINEARITE RELATIVE (CRITERE CONVERGENCE)
! CEIGEN   /O/: VALEURS PROPRES COMPLEXES DU PROBLEME
! CMOD     /O/: MODES PROPRES COMPLEXES SOLUTIONS
! NDIMAX   /I/: NOMBRE DE DDL GENERALISES DES MODES >=NDIM
! CMAT1    /M/: MATRICE COMPLEXE DE TRAVAIL
! CMAT2    /M/: MATRICE COMPLEXE DE TRAVAIL
! CVEC     /M/: VECTEUR COMPLEXE DE TRAVAIL
! IFIC     /I/: NUMERO UNITE LOGIQUE POUR MESSAGE
!
!-----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "asterfort/cmatve.h"
#include "asterfort/cschmi.h"
#include "asterfort/ctescv.h"
#include "asterfort/cvalea.h"
#include "asterfort/rrldc.h"
#include "asterfort/sesqui.h"
#include "asterfort/trldc.h"
#include "asterfort/utmess.h"
#include "blas/zcopy.h"
    integer(kind=8) :: vali(2), nbmod, ndim, ndimax
    complex(kind=8) :: ck(*), cm(*), ceigen(nbmod)
    complex(kind=8) :: cmod(ndimax, nbmod), cprod, cmod0(ndim)
    complex(kind=8) :: cmat1(*), cmat2(ndim, ndim), cvec(ndim), cvec0(ndim)
    aster_logical :: convok
    integer(kind=8) :: i, idiag, ific, ipivo, iv, ivdiag, j
    integer(kind=8) :: k, niter
    real(kind=8) :: valr(3), xcrit, xer
    character(len=6) :: valk
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
!
    valk = 'CMPHII'
    call utmess('I', 'ALGELINE7_2', sk=valk)
    call utmess('I', 'ALGELINE7_3')
!
!      RECOPIE DE LA MATRICE DE RAIDEUR
    b_n = to_blas_int(ndim*(ndim+1)/2)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call zcopy(b_n, ck, b_incx, cmat1, b_incy)
!
!    FACTORISATION DE LA MATRICE DE RAIDEUR
    call trldc(cmat1, ndim, ipivo)
!    GESTION DU PIVOT NUL
    if (ipivo .ne. 0) then
        vali(1) = ipivo
        call utmess('F', 'ALGORITH12_53', si=vali(1))
    end if
!
!
!   CALCUL DE L'INVERSE DE LA MATRICE DE MASSE
    do iv = 1, ndim
        ivdiag = iv*(iv-1)/2+1
        do i = 1, ndim
            if (i .le. iv) then
                cmat2(i, iv) = cm(ivdiag+iv-i)
            else
                idiag = i*(i-1)/2+1
                cmat2(i, iv) = dconjg(cm(idiag+i-iv))
            end if
        end do
    end do
    call rrldc(cmat1, ndim, cmat2, ndim)
!
    do iv = 1, ndim
        cvec(iv) = dcmplx(0.d0, 0.d0)
        cvec0(iv) = dcmplx(0.d0, 0.d0)
        cmod0(iv) = dcmplx(0.d0, 0.d0)
    end do
!
!
!   INITIALISATION ALEATOIRE DES VECTEURS PROPRES DE DEPART
    call cvalea(ndim, cmod, ndimax, nbmod)
!
!
!     DEBUT DE LA BOUCLE D'ITERATION SUR LES MODES
!
    do j = 1, nbmod
!
!       INITIALISATION DES CRITERES D'ARRET
        k = 0
        convok = .true.
!
!   BOUCLE D'ITERATION SUR CHAQUE MODES
100     continue
!
        k = k+1
!
!    PRODUIT MATRICIEL INV(M)*VECTEUR
        call cmatve(cmat2, cmod(1, j), cvec, ndim)
!
!    CALCUL DE L'ERREUR COLINEARITE ET REECOPIE
!    DE CVEC DANS CMOD
        call ctescv(cvec, cmod(1, j), cvec0, cmod0, ndim, &
                    xer)
!
!      RECOPIE DU VECTEUR DE L'ITERATION PRECEDENTE
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, cmod(1, j), b_incx, cmod0, b_incy)
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, cvec, b_incx, cvec0, b_incy)
!
!   ORTHORMALISATION PAR RAPPORT MATRICE DE MASSE
        call cschmi(cm, ndim, cmod(1, j), cmod, ndimax, &
                    j-1)
!
!
!   CALCUL VALEURS PROPRES PAR COEF RAYLEIGH
!    En fait on calcul explicitement CMOD*inv(M)*K*CMOD
        call sesqui(ck, cmod(1, j), ndim, ceigen(j))
        call sesqui(cm, cmod(1, j), ndim, cprod)
        ceigen(j) = ceigen(j)/cprod
!
!         TEST SUR LA PRECISION
        if (xer .le. xcrit) convok = .false.
!
        if (k .lt. niter .and. convok) goto 100
!
!
!     IMPRESSION DES FREQUENCES PROPRES
        vali(1) = j
        vali(2) = k
        valr(1) = xer
        valr(2) = dble(ceigen(j))
        valr(3) = dimag(ceigen(j))
        call utmess('I', 'ALGELINE7_4', ni=2, vali=vali, nr=3, &
                    valr=valr)
!
    end do
!
    write (ific, *) '     '
    write (ific, *) '     '
!
end subroutine
