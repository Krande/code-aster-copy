! --------------------------------------------------------------------
! Copyright (C) LAPACK
! Copyright (C) 2007 - 2025 - EDF R&D - www.code-aster.org
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
! ===============================================================
! THIS LAPACK 2.0 ROUTINE IS DEPRECATED
! DO NOT USE IT : YOU SHOULD PREFER UP-TO-DATE LAPACK ROUTINE
!
! BUT DO NOT REMOVE IT :
! THE PRESENT ROUTINE IS MANDATORY FOR ARPACK LIBRARY
! WHICH STICKS TO LAPACK 2.0 VERSION
! ==============================================================
subroutine ar_dlahqr(wantt, wantz, n, ilo, ihi, &
                     h, ldh, wr, wi, iloz, &
                     ihiz, z, ldz, info)
!
!     SUBROUTINE LAPACK DE MISE A JOUR DES VALEURS PROPRES ET DE LA
!     DECOMPOSITION DE SCHUR DEJA CALCULEES PAR DHSEQR.
!-----------------------------------------------------------------------
!  -- LAPACK AUXILIARY ROUTINE (VERSION 2.0) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     OCTOBER 31, 1992
!
!  PURPOSE
!  =======
!
!  DLAHQR IS AN AUXILIARY ROUTINE CALLED BY DHSEQR TO UPDATE THE
!  EIGENVALUES AND SCHUR DECOMPOSITION ALREADY COMPUTED BY DHSEQR, BY
!  DEALING WITH THE HESSENBERG SUBMATRIX IN ROWS AND COLUMNS ILO TO IHI.
!
!  ARGUMENTS
!  =========
!
!  WANTT   (INPUT) LOGICAL
!          = .TRUE. : THE FULL SCHUR FORM T IS REQUIRED,
!          = .FALSE.: ONLY EIGENVALUES ARE REQUIRED.
!
!  WANTZ   (INPUT) LOGICAL
!          = .TRUE. : THE MATRIX OF SCHUR VECTORS Z IS REQUIRED,
!          = .FALSE.: SCHUR VECTORS ARE NOT REQUIRED.
!
!  N       (INPUT) INTEGER
!          THE ORDER OF THE MATRIX H.  N >= 0.
!
!  ILO     (INPUT) INTEGER
!  IHI     (INPUT) INTEGER
!          IT IS ASSUMED THAT H IS ALREADY UPPER QUASI-TRIANGULAR IN
!          ROWS AND COLUMNS IHI+1:N, AND THAT H(ILO,ILO-1) = 0 (UNLESS
!          ILO = 1). DLAHQR WORKS PRIMARILY WITH THE HESSENBERG
!          SUBMATRIX IN ROWS AND COLUMNS ILO TO IHI, BUT APPLIES
!          TRANSFORMATIONS TO ALL OF H IF WANTT IS .TRUE..
!          1 <= ILO <= MAX(1,IHI), IHI <= N.
!
!  H       (INPUT/OUTPUT) REAL*8 ARRAY, DIMENSION (LDH,N)
!          ON ENTRY, THE UPPER HESSENBERG MATRIX H.
!          ON EXIT, IF WANTT IS .TRUE., H IS UPPER QUASI-TRIANGULAR IN
!          ROWS AND COLUMNS ILO:IHI, WITH ANY 2-BY-2 DIAGONAL BLOCKS IN
!          STANDARD FORM. IF WANTT IS .FALSE., THE CONTENTS OF H ARE
!          UNSPECIFIED ON EXIT.
!
!  LDH     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY H. LDH >= MAX(1,N).
!
!  WR      (OUTPUT) REAL*8 ARRAY, DIMENSION (N)
!  WI      (OUTPUT) REAL*8 ARRAY, DIMENSION (N)
!          THE REAL AND IMAGINARY PARTS, RESPECTIVELY, OF THE COMPUTED
!          EIGENVALUES ILO TO IHI ARE STORED IN THE CORRESPONDING
!          ELEMENTS OF WR AND WI. IF TWO EIGENVALUES ARE COMPUTED AS A
!          COMPLEX CONJUGATE PAIR, THEY ARE STORED IN CONSECUTIVE
!          ELEMENTS OF WR AND WI, SAY THE I-TH AND (I+1)TH, WITH
!          WI(I) > 0 AND WI(I+1) < 0. IF WANTT IS .TRUE., THE
!          EIGENVALUES ARE STORED IN THE SAME ORDER AS ON THE DIAGONAL
!          OF THE SCHUR FORM RETURNED IN H, WITH WR(I) = H(I,I), AND, IF
!          H(I:I+1,I:I+1) IS A 2-BY-2 DIAGONAL BLOCK,
!          WI(I) = SQRT(H(I+1,I)*H(I,I+1)) AND WI(I+1) = -WI(I).
!
!  ILOZ    (INPUT) INTEGER
!  IHIZ    (INPUT) INTEGER
!          SPECIFY THE ROWS OF Z TO WHICH TRANSFORMATIONS MUST BE
!          APPLIED IF WANTZ IS .TRUE..
!          1 <= ILOZ <= ILO, IHI <= IHIZ <= N.
!
!  Z       (INPUT/OUTPUT) REAL*8 ARRAY, DIMENSION (LDZ,N)
!          IF WANTZ IS .TRUE., ON ENTRY Z MUST CONTAIN THE CURRENT
!          MATRIX Z OF TRANSFORMATIONS ACCUMULATED BY DHSEQR, AND ON
!          EXIT Z HAS BEEN UPDATED, TRANSFORMATIONS ARE APPLIED ONLY TO
!          THE SUBMATRIX Z(ILOZ:IHIZ,ILO:IHI).
!          IF WANTZ IS .FALSE., Z IS NOT REFERENCED.
!
!  LDZ     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY Z. LDZ >= MAX(1,N).
!
!  INFO    (OUTPUT) INTEGER
!          = 0: SUCCESSFUL EXIT
!          > 0: DLAHQR FAILED TO COMPUTE ALL THE EIGENVALUES ILO TO IHI
!               IN A TOTAL OF 30*(IHI-ILO+1) ITERATIONS, IF INFO = I,
!               ELEMENTS I+1:IHI OF WR AND WI CONTAIN THOSE EIGENVALUES
!               WHICH HAVE BEEN SUCCESSFULLY COMPUTED.
!
!-----------------------------------------------------------------------
! ASTER INFORMATION
! 14/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            DISPARITION DE DLAMCH ET DE DLABAD,
!            UTILISATION DE R8PREM(), R8MIEM() ET ISBAEM(),
!            REMPLACEMENT DE 3 RETURN PAR GOTO 1000,
!            MODIFICATION DES APPELS BLAS (ROUTINE ASTER BL...),
!            IMPLICIT NONE.
! INTRINSIC FUNCTIONS
!            ABS, MAX, MIN.
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     .. SCALAR ARGUMENTS ..
#include "asterf_types.h"
#include "asterc/isbaem.h"
#include "asterc/matfpe.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/ar_dlanv2.h"
#include "asterfort/ar_dlarfg.h"
#include "blas/dcopy.h"
#include "blas/dlanhs.h"
#include "blas/drot.h"
    aster_logical :: wantt, wantz
    integer(kind=8) :: ihi, ihiz, ilo, iloz, info, ldh, ldz, n
!     ..
!     .. ARRAY ARGUMENTS ..
    real(kind=8) :: h(ldh, *), wi(*), wr(*), z(ldz, *)
!     ..
!     .. PARAMETERS ..
    real(kind=8) :: zero
    parameter(zero=0.0d+0)
    real(kind=8) :: dat1, dat2
    parameter(dat1=0.75d+0, dat2=-0.4375d+0)
!     ..
!     .. LOCAL SCALARS ..
    integer(kind=8) :: i, i1, i2, itn, its, j, k, l, m, nh, nr, nz
    real(kind=8) :: cs, h00, h10, h11, h12, h21, h22, h33, h33s, h43h34, h44
    real(kind=8) :: h44s, s, smlnum, sn, sum, t1, t2, t3, tst1, ulp, unfl, v1
    real(kind=8) :: v2, v3
! DUE TO CRS512       REAL*8 OVFL
!     ..
!     .. LOCAL ARRAYS ..
    real(kind=8) :: v(3), work(1)
    blas_int :: b_incx, b_incy, b_n
    blas_int :: b_lda
!     ..
!     .. EXTERNAL FUNCTIONS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
    call matfpe(-1)
!
! DUE TO CRS513
    work(1) = zero
    info = 0
!
!     QUICK RETURN IF POSSIBLE
!
    if (n .eq. 0) goto 1000
    if (ilo .eq. ihi) then
        wr(ilo) = h(ilo, ilo)
        wi(ilo) = zero
        goto 1000
    end if
!
    nh = ihi-ilo+1
    nz = ihiz-iloz+1
!
!     SET MACHINE-DEPENDENT CONSTANTS FOR THE STOPPING CRITERION.
!     IF NORM(H) <= SQRT(OVFL), OVERFLOW SHOULD NOT OCCUR.
!
    unfl = r8miem()
! DUE TO CRS512      OVFL = ONE / UNFL
    ulp = r8prem()*0.5d0*isbaem()
    smlnum = unfl*(nh/ulp)
!
!     I1 AND I2 ARE THE INDICES OF THE FIRST ROW AND LAST COLUMN OF H
!     TO WHICH TRANSFORMATIONS MUST BE APPLIED. IF EIGENVALUES ONLY ARE
!     BEING COMPUTED, I1 AND I2 ARE SET INSIDE THE MAIN LOOP.
!
    if (wantt) then
        i1 = 1
        i2 = n
    end if
!
!     ITN IS THE TOTAL NUMBER OF QR ITERATIONS ALLOWED.
!
    itn = 30*nh
!
!     THE MAIN LOOP BEGINS HERE. I IS THE LOOP INDEX AND DECREASES FROM
!     IHI TO ILO IN STEPS OF 1 OR 2. EACH ITERATION OF THE LOOP WORKS
!     WITH THE ACTIVE SUBMATRIX IN ROWS AND COLUMNS L TO I.
!     EIGENVALUES I+1 TO IHI HAVE ALREADY CONVERGED. EITHER L = ILO OR
!     H(L,L-1) IS NEGLIGIBLE SO THAT THE MATRIX SPLITS.
!
    i = ihi
10  continue
    l = ilo
    if (i .lt. ilo) goto 150
!
!     PERFORM QR ITERATIONS ON ROWS AND COLUMNS ILO TO I UNTIL A
!     SUBMATRIX OF ORDER 1 OR 2 SPLITS OFF AT THE BOTTOM BECAUSE A
!     SUBDIAGONAL ELEMENT HAS BECOME NEGLIGIBLE.
!
    do its = 0, itn
!
!        LOOK FOR A SINGLE SMALL SUBDIAGONAL ELEMENT.
!
        do k = i, l+1, -1
            tst1 = abs(h(k-1, k-1))+abs(h(k, k))
            b_lda = to_blas_int(ldh)
            b_n = to_blas_int(i-l+1)
            if (tst1 .eq. zero) tst1 = dlanhs('1', b_n, h(l, l), b_lda, work)
            if (abs(h(k, k-1)) .le. max(ulp*tst1, smlnum)) goto 30
        end do
30      continue
        l = k
        if (l .gt. ilo) then
!
!           H(L,L-1) IS NEGLIGIBLE
!
            h(l, l-1) = zero
        end if
!
!        EXIT FROM LOOP IF A SUBMATRIX OF ORDER 1 OR 2 HAS SPLIT OFF.
!
        if (l .ge. i-1) goto 140
!
!        NOW THE ACTIVE SUBMATRIX IS IN ROWS AND COLUMNS L TO I. IF
!        EIGENVALUES ONLY ARE BEING COMPUTED, ONLY THE ACTIVE SUBMATRIX
!        NEED BE TRANSFORMED.
!
        if (.not. wantt) then
            i1 = l
            i2 = i
        end if
!
        if (its .eq. 10 .or. its .eq. 20) then
!
!           EXCEPTIONAL SHIFT.
!
            s = abs(h(i, i-1))+abs(h(i-1, i-2))
            h44 = dat1*s
            h33 = h44
            h43h34 = dat2*s*s
        else
!
!           PREPARE TO USE WILKINSON'S DOUBLE SHIFT
!
            h44 = h(i, i)
            h33 = h(i-1, i-1)
            h43h34 = h(i, i-1)*h(i-1, i)
        end if
!
!        LOOK FOR TWO CONSECUTIVE SMALL SUBDIAGONAL ELEMENTS.
!
        do m = i-2, l, -1
!
!           DETERMINE THE EFFECT OF STARTING THE DOUBLE-SHIFT QR
!           ITERATION AT ROW M, AND SEE IF THIS WOULD MAKE H(M,M-1)
!           NEGLIGIBLE.
!
            h11 = h(m, m)
            h22 = h(m+1, m+1)
            h21 = h(m+1, m)
            h12 = h(m, m+1)
            h44s = h44-h11
            h33s = h33-h11
            v1 = (h33s*h44s-h43h34)/h21+h12
            v2 = h22-h11-h33s-h44s
            v3 = h(m+2, m+1)
            s = abs(v1)+abs(v2)+abs(v3)
            v1 = v1/s
            v2 = v2/s
            v3 = v3/s
            v(1) = v1
            v(2) = v2
            v(3) = v3
            if (m .eq. l) goto 50
            h00 = h(m-1, m-1)
            h10 = h(m, m-1)
            tst1 = abs(v1)*(abs(h00)+abs(h11)+abs(h22))
            if (abs(h10)*(abs(v2)+abs(v3)) .le. ulp*tst1) goto 50
        end do
50      continue
!
!        DOUBLE-SHIFT QR STEP
!
        do k = m, i-1
!
!           THE FIRST ITERATION OF THIS LOOP DETERMINES A REFLECTION G
!           FROM THE VECTOR V AND APPLIES IT FROM LEFT AND RIGHT TO H,
!           THUS CREATING A NONZERO BULGE BELOW THE SUBDIAGONAL.
!
!           EACH SUBSEQUENT ITERATION DETERMINES A REFLECTION G TO
!           RESTORE THE HESSENBERG FORM IN THE (K-1)TH COLUMN, AND THUS
!           CHASES THE BULGE ONE STEP TOWARD THE BOTTOM OF THE ACTIVE
!           SUBMATRIX. NR IS THE ORDER OF G.
!
            nr = min(3, i-k+1)
            if (k .gt. m) then
                b_n = to_blas_int(nr)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, h(k, k-1), b_incx, v, b_incy)
            end if
            call ar_dlarfg(nr, v(1), v(2), 1, t1)
            if (k .gt. m) then
                h(k, k-1) = v(1)
                h(k+1, k-1) = zero
                if (k .lt. i-1) h(k+2, k-1) = zero
            else if (m .gt. l) then
                h(k, k-1) = -h(k, k-1)
            end if
            v2 = v(2)
            t2 = t1*v2
            if (nr .eq. 3) then
                v3 = v(3)
                t3 = t1*v3
!
!              APPLY G FROM THE LEFT TO TRANSFORM THE ROWS OF THE MATRIX
!              IN COLUMNS K TO I2.
!
                do j = k, i2
                    sum = h(k, j)+v2*h(k+1, j)+v3*h(k+2, j)
                    h(k, j) = h(k, j)-sum*t1
                    h(k+1, j) = h(k+1, j)-sum*t2
                    h(k+2, j) = h(k+2, j)-sum*t3
                end do
!
!              APPLY G FROM THE RIGHT TO TRANSFORM THE COLUMNS OF THE
!              MATRIX IN ROWS I1 TO MIN(K+3,I).
!
                do j = i1, min(k+3, i)
                    sum = h(j, k)+v2*h(j, k+1)+v3*h(j, k+2)
                    h(j, k) = h(j, k)-sum*t1
                    h(j, k+1) = h(j, k+1)-sum*t2
                    h(j, k+2) = h(j, k+2)-sum*t3
                end do
!
                if (wantz) then
!
!                 ACCUMULATE TRANSFORMATIONS IN THE MATRIX Z
!
                    do j = iloz, ihiz
                        sum = z(j, k)+v2*z(j, k+1)+v3*z(j, k+2)
                        z(j, k) = z(j, k)-sum*t1
                        z(j, k+1) = z(j, k+1)-sum*t2
                        z(j, k+2) = z(j, k+2)-sum*t3
                    end do
                end if
            else if (nr .eq. 2) then
!
!              APPLY G FROM THE LEFT TO TRANSFORM THE ROWS OF THE MATRIX
!              IN COLUMNS K TO I2.
!
                do j = k, i2
                    sum = h(k, j)+v2*h(k+1, j)
                    h(k, j) = h(k, j)-sum*t1
                    h(k+1, j) = h(k+1, j)-sum*t2
                end do
!
!              APPLY G FROM THE RIGHT TO TRANSFORM THE COLUMNS OF THE
!              MATRIX IN ROWS I1 TO MIN(K+3,I).
!
                do j = i1, i
                    sum = h(j, k)+v2*h(j, k+1)
                    h(j, k) = h(j, k)-sum*t1
                    h(j, k+1) = h(j, k+1)-sum*t2
                end do
!
                if (wantz) then
!
!                 ACCUMULATE TRANSFORMATIONS IN THE MATRIX Z
!
                    do j = iloz, ihiz
                        sum = z(j, k)+v2*z(j, k+1)
                        z(j, k) = z(j, k)-sum*t1
                        z(j, k+1) = z(j, k+1)-sum*t2
                    end do
                end if
            end if
        end do
!
    end do
!
!     FAILURE TO CONVERGE IN REMAINING NUMBER OF ITERATIONS
!
    info = i
    goto 1000
!
140 continue
!
    if (l .eq. i) then
!
!        H(I,I-1) IS NEGLIGIBLE: ONE EIGENVALUE HAS CONVERGED.
!
        wr(i) = h(i, i)
        wi(i) = zero
    else if (l .eq. i-1) then
!
!        H(I-1,I-2) IS NEGLIGIBLE: A PAIR OF EIGENVALUES HAVE CONVERGED.
!
!        TRANSFORM THE 2-BY-2 SUBMATRIX TO STANDARD SCHUR FORM,
!        AND COMPUTE AND STORE THE EIGENVALUES.
!
        call ar_dlanv2(h(i-1, i-1), h(i-1, i), h(i, i-1), h(i, i), wr(i-1), &
                       wi(i-1), wr(i), wi(i), cs, sn)
!
        if (wantt) then
!
!           APPLY THE TRANSFORMATION TO THE REST OF H.
!
            if (i2 .gt. i) then
                b_n = to_blas_int(i2-i)
                b_incx = to_blas_int(ldh)
                b_incy = to_blas_int(ldh)
                call drot(b_n, h(i-1, i+1), b_incx, h(i, i+1), b_incy, &
                          cs, sn)
            end if
            b_n = to_blas_int(i-i1-1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call drot(b_n, h(i1, i-1), b_incx, h(i1, i), b_incy, &
                      cs, sn)
        end if
        if (wantz) then
!
!           APPLY THE TRANSFORMATION TO Z.
!
            b_n = to_blas_int(nz)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call drot(b_n, z(iloz, i-1), b_incx, z(iloz, i), b_incy, &
                      cs, sn)
        end if
    end if
!
!     DECREMENT NUMBER OF REMAINING ITERATIONS, AND RETURN TO START OF
!     THE MAIN LOOP WITH NEW VALUE OF I.
!
    itn = itn-its
    i = l-1
    goto 10
!
150 continue
1000 continue
!
    call matfpe(1)
!
!     END OF DLAHQR
!
end subroutine
