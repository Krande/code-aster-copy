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
subroutine unista(h, ldh, v, ldv, ddlsta, &
                  n, vectp, csta, beta, etat, &
                  ldynfa, ddlexc, redem)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
! ROUTINE PRINCIPALE - ALGORITHME D'OPTIMISATION SOUS CONTRAINTES
! POUR LES ETUDES DE STABILITE
!
! ----------------------------------------------------------------------
!
! IN  H        : MATRICE REDUITE
! IN  LDH      : NOMBRE DE COEFFICIENTS DE H
! IN  V        : MATRICE DE CHANGEMENT DE BASE
! IN  LDV      : NOMBRE DE COEFFICIENTS DE V
! IN  DDLSTA   : POSITION DES DDL_STAB
! IN  N        : DIMENSION ESPACE GLOBAL
! IN/OUT VECTP : MODE DE STABILITE
! OUT CSTA     : VALEUR CRITERE DE STABILITE
! IN  BETA     : PLUS GRANDE VALEUR PROPRE NEGATIVE
! IN  ETAT     : =0 TTES LES VP SONT POSITIVES
!                =1 AU MOINS UNE VP EST NEGATIVE
! IN  LDYNFA   : DESCRIPTEUR OPERATEUR TANGENT
! IN  DDLEXC   : POSITION DDL IMPOSES
! IN  REDEM    : NOMBRE REDEMARRAGES METHODE SORENSEN
!
! ----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     %-----------------%
!     | ARRAY ARGUMENTS |
!     %-----------------%
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedetr.h"
#include "asterfort/mgauss.h"
#include "asterfort/mppsta.h"
#include "asterfort/mrmult.h"
#include "asterfort/r8inir.h"
#include "asterfort/wkvect.h"
#include "blas/dnrm2.h"
#include "blas/dscal.h"
    integer(kind=8) :: n, ldh, ldv
    integer(kind=8) :: ddlsta(n), ddlexc(n)
    real(kind=8) :: h(ldh, ldh), v(ldv, ldh)
    real(kind=8) :: vectp(ldv)
!
!     %------------------%
!     | SCALAR ARGUMENTS |
!     %------------------%
!
    integer(kind=8) :: etat, ldynfa
    integer(kind=8) :: redem
    real(kind=8) :: beta, csta
!
!
!     %------------------------%
!     | LOCAL SCALARS & ARRAYS |
!     %------------------------%
!
    integer(kind=8) :: i, j, iret, indico, proj
    integer(kind=8) :: vectt, xsol, vect2
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: q, b
    real(kind=8) :: gama, det
    real(kind=8) :: vtest, err
    real(kind=8) :: zero, one
    parameter(one=1.0d+0, zero=0.0d+0)
    character(len=4) :: cara
    blas_int :: b_incx, b_n
!
    call infniv(ifm, niv)
!
    call wkvect('&&UNISTA.VECT.TEM1', 'V V R', ldv, vectt)
    call wkvect('&&UNISTA.VECT.TEM2', 'V V R', ldv, xsol)
    call wkvect('&&UNISTA.VECT.TEM3', 'V V R', ldv, vect2)
    call wkvect('&&UNISTA.VECT.TEM4', 'V V R', ldh*ldh, q)
    call wkvect('&&UNISTA.VECT.TEM5', 'V V R', ldh*ldh, b)
!
    call r8inir(ldv, 0.d0, zr(vectt), 1)
    call r8inir(ldv, 0.d0, zr(xsol), 1)
    call r8inir(ldv, 0.d0, zr(vect2), 1)
    call r8inir(ldh*ldh, 0.d0, zr(q), 1)
    call r8inir(ldh*ldh, 0.d0, zr(b), 1)
!
    indico = 0
    proj = 0
!
    if (etat .eq. 0) then
!
!     1ER APPEL A LA METHODE DES PUISSANCES
!
        call mppsta(h, ldh, v, ldv, ddlsta, &
                    n, zr(vectt), ddlexc, indico, proj)
!
    end if
!
    if (etat .eq. 1) then
!
!     MISE EN PLACE SHIFT POUR VALEUR NEGATIVE
!
        gama = abs(beta)+1.d0
!
        do i = 1, ldh
            do j = 1, ldh
                zr(q+(j-1)*ldh+i-1) = -h(i, j)
            end do
        end do
        do i = 1, ldh
            zr(q+(i-1)*(ldh+1)) = gama+zr(q+(i-1)*(ldh+1))
        end do
!
!     2ND APPEL A LA METHODE DES PUISSANCES
!
        proj = 1
!
        call mppsta(zr(q), ldh, v, ldv, ddlsta, &
                    n, zr(vectt), ddlexc, indico, proj)
!
    end if
!
    b_n = to_blas_int(ldv)
    b_incx = to_blas_int(1)
    err = dnrm2(b_n, zr(vectt), b_incx)
    b_n = to_blas_int(ldv)
    b_incx = to_blas_int(1)
    call dscal(b_n, one/err, zr(vectt), b_incx)
    call mrmult('ZERO', ldynfa, zr(vectt), zr(xsol), 1, &
                .true._1)
!
    vtest = 0.d0
    do i = 1, ldv
        vtest = vtest+zr(vectt+i-1)*zr(xsol+i-1)
        if (ddlsta(i) .eq. 0 .and. proj .eq. 1) then
            if (zr(vectt+i-1) .lt. zero) then
                ASSERT(.false.)
            end if
        end if
    end do
!
    write (ifm, *) 'VAL1_STAB : ', vtest
!
    if (vtest .lt. 0.d0 .or. etat .eq. 0) then
        do i = 1, ldv
            vectp(i) = zr(vectt+i-1)
        end do
        csta = vtest
        goto 300
    end if
!
!     CALCUL DU CRITERE
!
    do i = 1, ldh
        do j = 1, ldh
            if (i .eq. j) then
                zr(b+(i-1)*ldh+j-1) = 1.d0
            else
                zr(b+(i-1)*ldh+j-1) = 0.d0
            end if
        end do
    end do
!
    cara = 'NFSP'
!
    call mgauss(cara, h, zr(b), ldh, ldh, &
                ldh, det, iret)
!
    proj = 0
!
    call mppsta(zr(b), ldh, v, ldv, ddlsta, &
                n, zr(vect2), ddlexc, indico, proj)
!
    b_n = to_blas_int(ldv)
    b_incx = to_blas_int(1)
    err = dnrm2(b_n, zr(vect2), b_incx)
    b_n = to_blas_int(ldv)
    b_incx = to_blas_int(1)
    call dscal(b_n, one/err, zr(vect2), b_incx)
    call mrmult('ZERO', ldynfa, zr(vect2), zr(xsol), 1, &
                .true._1)
    vtest = 0.d0
    do i = 1, ldv
        vtest = vtest+zr(vect2+i-1)*zr(xsol+i-1)
    end do
!
!      WRITE (IFM,*) 'VAL_PROPRE_MAX : ',VTEST
!
    do i = 1, ldh
        do j = 1, ldh
            zr(q+(i-1)*ldh+j-1) = -zr(b+(i-1)*ldh+j-1)-zr(b+(j-1)*ldh+i-1)
        end do
    end do
    do i = 1, ldh
!        Q(I,I) = 2*VTEST*(1.D0/4.D0)+Q(I,I)
        zr(q+(i-1)*ldh+i-1) = 2*vtest*(1.5d0/2.d0)+zr(q+(i-1)*ldh+i-1)
    end do
!
    if (redem .eq. 0) then
        do i = 1, ldv
            vectp(i) = zr(vectt+i-1)
        end do
    end if
!
    indico = 1
    proj = 1
!
    call mppsta(zr(q), ldh, v, ldv, ddlsta, &
                n, vectp, ddlexc, indico, proj)
!
    b_n = to_blas_int(ldv)
    b_incx = to_blas_int(1)
    err = dnrm2(b_n, vectp, b_incx)
    b_n = to_blas_int(ldv)
    b_incx = to_blas_int(1)
    call dscal(b_n, one/err, vectp, b_incx)
    call mrmult('ZERO', ldynfa, vectp, zr(xsol), 1, &
                .true._1)
    vtest = 0.d0
    do i = 1, ldv
        vtest = vtest+vectp(i)*zr(xsol+i-1)
        if (ddlsta(i) .eq. 0) then
            if (vectp(i) .lt. zero) then
                ASSERT(.false.)
            end if
        end if
    end do
!
    write (ifm, *) 'VAL2_STAB : ', vtest
    write (ifm, 9070)
    write (ifm, 9080)
!
    csta = vtest
!
300 continue
!
    redem = redem+1
!
9070 format(72(' '))
9080 format(72('-'))
!
! ----------------------------------------------
!
! --- MENAGE
    call jedetr('&&UNISTA.VECT.TEM1')
    call jedetr('&&UNISTA.VECT.TEM2')
    call jedetr('&&UNISTA.VECT.TEM3')
    call jedetr('&&UNISTA.VECT.TEM4')
    call jedetr('&&UNISTA.VECT.TEM5')
!
end subroutine
