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
subroutine cmphdi(ck, cm, ndim, nbmod, niter, &
                  xcrit, ceigen, cmod, ndimax, cmat1, &
                  cmat2, cvect, cvect1, alpha, beta, &
                  lambd1, lambd2, interv)
! aslint: disable=W1306
    implicit none
!
!***********************************************************************
!    B. GUIGON     P. RICHARD                  DATE 06/04/92
!-----------------------------------------------------------------------
!  BUT:  < COMPLEXE MODE PROBLEME HERMITIEN DICHOTOMIE INVERSE >
!
!   CALCULER LES N PREMIERS MODES PROPRES DE LA BANDE LAMBD1, LAMBD2
!   D'UN PROBLEME AUX VALEURS PROPRES A MATRICES RAIDEUR ET MASSE
!   COMPLEXES HERMITIENNES STOCKEES TRIANGULAIRES SUPERIEURES
!                     CK*X= L CM*X
!
!    METHODE D'ITERATION INVERSE ASSOCIEE A UN SHIFT DE CHACUNE DES
!    VALEURS PROPRES ET A UNE SEPARATION PREALABLE DE CES VALEURS
!
!-----------------------------------------------------------------------
!
! CK       /I/: MATRICE RAIDEUR DU PROBLEME
! CM       /I/: MATRICE MASSE DU PROBLEME
! NDIM     /I/: DIMENSION DES MATRICES
! NBMOD    /M/: NOMBRE DE MODES PROPRES DESIRE/EXISTANT
! NITER    /I/: NOMBRE MAX D'ITERATIONS PAR MODE
! XCRIT    /I/: TOLERANCE DE COLINEARITE RELATIVE (CRITERE CONVERGENCE)
! CEIGEN   /O/: VALEURS PROPRES COMPLEXES DU PROBLEME
! CMOD     /O/: MODES PROPRES COMPLEXES SOLUTIONS
! NDIMAX   /I/: NOMBRE DE DDL GENERALISES DES MODES >=NDIM
! CMAT1    /M/: MATRICE COMPLEXE DE TRAVAIL
! CMAT2    /M/: MATRICE COMPLEXE DE TRAVAIL
! CVECT    /M/: VECTEUR COMPLEXE DE TRAVAIL
! CVECT1   /M/: VECTEUR COMPLEXE DE TRAVAIL
! ALPHA    /M/: VECTEUR REEL DE TRAVAIL
! BETA     /M/: VECTEUR REEL DE TRAVAIL
! LAMBD1   /I/: BORNE INFERIEURE DE L'INTERVALLE DE RECHERCHE
! LAMBD2   /I/: BORNE SUPERIEURE DE L'INTERVALLE DE RECHERCHE
! INTERV   /I/: LONGUEUR MAXIMALE D'UN INTERVALLE CONTENANT UNE VP
!
!-----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "asterfort/cmatve.h"
#include "asterfort/ctescv.h"
#include "asterfort/cvalea.h"
#include "asterfort/cvnorm.h"
#include "asterfort/rrldc.h"
#include "asterfort/sepavp.h"
#include "asterfort/sesqui.h"
#include "asterfort/trldc.h"
#include "asterfort/utmess.h"
#include "blas/zcopy.h"
    integer(kind=8) :: ndim, nbmod, niter, ndimax
    complex(kind=8) :: ck(*), cm(*)
    complex(kind=8) :: ceigen(nbmod), cmod(ndimax, nbmod), cmod0(ndim)
    complex(kind=8) :: cmat1(*), cmat2(ndim, ndim)
    complex(kind=8) :: cvect(ndim), cvect1(ndim), cvec0(ndim)
    real(kind=8) :: alpha(ndim+1), beta(ndim+1), interv
    real(kind=8) :: lambd1, lambd2, xcrit
    integer(kind=8) :: i, j, ct, ipivo
    integer(kind=8) :: vali(2)
    complex(kind=8) :: cshift
    real(kind=8) :: ecart, valr(3)
    aster_logical :: sortie
    integer(kind=8) :: idiag, iretou, iv, ivdiag
    character(len=6) :: valk
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
!
    valk = 'CMPHDI'
    call utmess('I', 'ALGELINE7_2', sk=valk)
!
!
!        SEPARATION DES VALEURS PROPRES
!
    call sepavp(ck, cm, cmat1, ndim, alpha, &
                beta, nbmod, lambd1, lambd2, interv)
!
    call utmess('I', 'ALGELINE7_3')
!
!        INITIALISATION DES VECTEURS POUR LES ITERATIONS
!
    do iv = 1, ndim
        cvect(iv) = dcmplx(0.d0, 0.d0)
        cvec0(iv) = dcmplx(0.d0, 0.d0)
        cmod0(iv) = dcmplx(0.d0, 0.d0)
    end do
!
    call cvalea(ndim, cmod, ndimax, nbmod)
!
!        BOUCLE SUR LES MODES
!
    do j = 1, nbmod
!
!        INITIALISATION DE LA MATRICE SHIFTEE, ET COPIE DANS LES
!        MATRICES DE TRAVAIL
!
        cshift = dcmplx((alpha(j)+beta(j))/2.d0, 0.d0)
        do i = 1, ndim*(ndim+1)/2
            cmat1(i) = ck(i)-cshift*cm(i)
        end do
!
!        CALCUL DE LA MATRICE CM*(CK-SHIFT*CM)**-1
!
        call trldc(cmat1, ndim, ipivo)
        if (ipivo .ne. 0) then
            vali(1) = ipivo
            call utmess('F', 'ALGORITH12_53', si=vali(1))
        end if
!
!
!   CALCUL DE LA MATRICE INVERSE DU PROBLEME
!
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
!        INITIALISATION DES VARIABLES DE LA BOUCLE
!
        sortie = .false.
        call cvnorm(cm, cmod(1, j), ndim, iretou)
        if (iretou .eq. 1) then
            call utmess('F', 'ALGORITH2_22')
        end if
!
        ct = 0
!
!         ITERATION INVERSE PROPREMENT DITE
!
30      continue
        if (sortie) goto 40
!      RECOPIE DU VECTEUR DE L'ITERATION PRECEDENTE
        ct = ct+1
        call cmatve(cmat2, cmod(1, j), cvect1, ndim)
        call cvnorm(cm, cvect1, ndim, iretou)
        if (iretou .eq. 1) then
            call utmess('F', 'ALGORITH2_22')
        end if
        call ctescv(cvect1, cmod(1, j), cvec0, cmod0, ndim, &
                    ecart)
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, cmod(1, j), b_incx, cmod0, b_incy)
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, cvect1, b_incx, cvec0, b_incy)
        if (ecart .le. xcrit) sortie = .true.
        if (ct .ge. niter) sortie = .true.
        goto 30
40      continue
!
!         CALCUL DE LA VALEUR PROPRE PAR LE COEFFICIENT DE RAYLEIGH
!
        call sesqui(ck, cmod(1, j), ndim, ceigen(j))
!
        vali(1) = j
        vali(2) = ct
        valr(1) = ecart
        valr(2) = dble(ceigen(j))
        valr(3) = dimag(ceigen(j))
        call utmess('I', 'ALGELINE7_4', ni=2, vali=vali, nr=3, &
                    valr=valr)
    end do
!
end subroutine
