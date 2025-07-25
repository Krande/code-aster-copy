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
subroutine mgauss(cara, a, b, dim, nordre, &
                  nb, det, iret)
!
!
! aslint: disable=W1306
    implicit none
!
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/mgausw.h"
#include "asterfort/utmess.h"
#include "blas/dgesv.h"
#include "blas/dgesvx.h"
    character(len=*) :: cara
    integer(kind=8) :: dim, nb, nordre, iret, idb
    real(kind=8) :: a(dim, dim), b(dim, nb), det
!
! ----------------------------------------------------------------------
!  RESOLUTION PAR FACTORISATION LU D'UN SYSTEME LINEAIRE
! ----------------------------------------------------------------------
!  VARIABLES D'ENTREE
!  K4  CARA : CHAINE DE CARACTERES PRECISANT CERTAINES CARACTERISTIQUES
!             DE LA RESOLUTION :
!
!       CARA(1:1) = /'N' : SYSTEME "NORMAL"    : A*X=B
!                   /'T' : SYSTEME "TRANSPOSE" : A'*X=B
!
!       CARA(2:2) = /'F' : ON S'ARRETE EN ERREUR <F> EN CAS DE PROBLEME
!                   /'C' : SI PROBLEME, ON "CONTINUE" AVEC IRET > 0
!
!       CARA(3:3) = /'S' : ON VEUT ALLER SUREMENT : LA METHODE DE
!                          RESOLUTION EST EN PRINCIPE PLUS PRECISE
!                          ET SI INFO=2, ON IMPRIME DES INFORMATIONS
!                          SUR LA QUALITE DE LA SOLUTION:
!                          REMARQUE: LA MATRICE 'A' N'EST PAS MODIFIEE
!                         (POUR EVITER UNE MODIFICATION DE CETTE MATRICE
!                         ,ON FAIT PASSER DANS LA ROUTINE DE RESOLUTION
!                         UNE COPIE DE 'A', QUI ELLE PEUT ETRE MODIFIEE
!                         SI UN EQUILIBRAGE EST NECESSAIRE)
!                   /'V' : ON VEUT ALLER VITE, ON NE CONTROLE PAS
!                          LA QUALITE DE LA SOLUTION
!                   /'W' : ON VEUT ALLER ENCORE PLUS VITE ET UTILISER
!                          L'ANCIENNE ROUTINE MGAUSS DE G. RATEAU
!                    ATTENTION : SI 'W', LA MATRICE A EST MODIFIEE
!
!
!       CARA(4:4) = /'D' : ON VEUT AUSSI LE DETERMINANT DE A
!                   /'P' : ON NE VEUT PAS CALCULER LE DETERMINANT DE A
!                   ATTENTION : LE CALCUL DU DETERMINANT PEUT PROVOQUER
!                                ASSEZ FACILEMENT DES "OVERFLOW"
!
!  REAL*8      A(DIM, DIM) : MATRICE CARREE PLEINE
!                            A N'EST PAS MODIFIEE SAUF SI CARA(3:3)='W'
!  REAL*8      B(DIM, NB)  : SECONDS MEMBRES
!  INTEGER     DIM         : DIMENSION DE A
!  INTEGER     NORDRE      : RANG DE A (NORDRE <= DIM)
!  INTEGER     NB          : NOMBRE DE SECONDS MEMBRES
!
!  VARIABLES DE SORTIE
!  REAL*8      B(DIM, NB)  : A-1 * B
!  REAL*8      DET         : DETERMINANT DE A
!  INTEGER     IRET        : CODE RETOUR
!                              IRET = 0 : OK
!                              IRET > 0 : PB
!
! ----------------------------------------------------------------------
    integer(kind=8) :: n, nrhs, ldb, ldx, ifm, niv, i, j, lda, ldaf
    integer(kind=4) :: ipiv4(dim), inf4, iwork4(dim)
    integer(kind=8) :: vali(2)
    real(kind=8) :: af(dim, dim), r(dim), c(dim), x(dim, nb), rcond
    real(kind=8) :: work(4*dim), ferr(nb), berr(nb), aa(dim, dim)
    real(kind=8) :: detr, detc
    character(len=1) :: fact, equed, trans2
    character(len=4) :: cara2
    character(len=24) :: valk(2)
    aster_logical :: ltrans, lstop, ldet, lret
    blas_int :: b_lda, b_ldb, b_n, b_nrhs
    blas_int :: b_ldaf, b_ldx
!----------------------------------------------------------------------
    call matfpe(-1)
!
    cara2 = cara
    ASSERT((cara2(1:1) .eq. 'N') .or. (cara2(1:1) .eq. 'T'))
    ASSERT((cara2(2:2) .eq. 'F') .or. (cara2(2:2) .eq. 'C'))
    ASSERT((cara2(3:3) .eq. 'V') .or. (cara2(3:3) .eq. 'S') .or. (cara2(3:3) .eq. 'W'))
    ASSERT((cara2(4:4) .eq. 'D') .or. (cara2(4:4) .eq. 'P'))
!
    ltrans = (cara2(1:1) .eq. 'T')
    lstop = (cara2(2:2) .eq. 'F')
    ldet = (cara2(4:4) .eq. 'D')
!
    det = r8nnem()
!
!
!     -- 1. ON VEUT ALLER SUREMENT (QUITTE A PERDRE DU TEMPS):
!     ---------------------------------------------------------
    if (cara2(3:3) .eq. 'S') then
! ---   DEFINITION DES PARAMETRES D'ENTREE POUR L'APPEL A LA ROUTINE
!       LAPACK : DGESVX
!       FACT : PERMET D'EQUILIBRER LA MATRICE (SI BESOIN)
        fact = 'E'
!       TRANS2 : CARACTERE PERMETTANT DE TRANSPOSER A
        trans2 = 'N'
        if (ltrans) trans2 = 'T'
!       N : ORDRE DE LA MATRICE A
        n = nordre
!       NRHS : NOMBRE DE COLONNES DE X
        nrhs = nb
        lda = dim
        ldaf = dim
        ldb = dim
        ldx = dim
!
!       SAUVEGARDE DE LA MATRICE 'A' DANS UNE MATRICE DE TRAVAIL: 'AA'
        aa = a
!
! --- RESOLUTION
        b_ldx = to_blas_int(ldx)
        b_ldb = to_blas_int(ldb)
        b_ldaf = to_blas_int(ldaf)
        b_lda = to_blas_int(lda)
        b_n = to_blas_int(n)
        b_nrhs = to_blas_int(nrhs)
        call dgesvx(fact, trans2, b_n, b_nrhs, aa, &
                    b_lda, af, b_ldaf, ipiv4, equed, &
                    r, c, b, b_ldb, x, &
                    b_ldx, rcond, ferr, berr, work, &
                    iwork4, inf4)
        iret = inf4
!
!       -- RECOPIE DE X DANS B :
        do i = 1, n
            do j = 1, nb
                b(i, j) = x(i, j)
            end do
        end do
!
        if (ldet) then
            det = 1.d0
            detr = 1.d0
            detc = 1.d0
            do i = 1, n
                if (ipiv4(i) .ne. i) det = (-1.d0)*det
                det = det*af(i, i)
                detr = detr*r(i)
                detc = detc*c(i)
            end do
            if (equed .eq. 'R') then
                det = det/detr
            else if (equed .eq. 'C') then
                det = det/detc
            else if (equed .eq. 'B') then
                det = det/(detr*detc)
            end if
        end if
!
!
!     -- 2. ON VEUT ALLER VITE :
!     ---------------------------------------------------------
    else if (cara2(3:3) .eq. 'V') then
! ---   DEFINITION DES PARAMETRES D'ENTREE POUR L'APPEL A LA ROUTINE
!       LAPACK : DGESV
!       N : ORDRE DE LA MATRICE A
        n = nordre
!       NRHS : NOMBRE DE COLONNES DE X
        nrhs = nb
        lda = dim
        ldb = dim
!
        if (ltrans) then
            do i = 1, n
                do j = 1, n
                    af(j, i) = a(i, j)
                end do
            end do
        else
            do i = 1, n
                do j = 1, n
                    af(i, j) = a(i, j)
                end do
            end do
        end if
!
!       ---   RESOLUTION
        b_ldb = to_blas_int(ldb)
        b_lda = to_blas_int(lda)
        b_n = to_blas_int(n)
        b_nrhs = to_blas_int(nrhs)
        call dgesv(b_n, b_nrhs, af, b_lda, ipiv4, &
                   b, b_ldb, inf4)
        iret = inf4
        if (ldet) then
            det = 1.d0
            do i = 1, n
                if (ipiv4(i) .ne. i) det = (-1.d0)*det
                det = det*af(i, i)
            end do
        end if
!
!
!     -- 3. ON VEUT ALLER ENCORE PLUS VITE :
!     ---------------------------------------------------------
    else if (cara2(3:3) .eq. 'W') then
        n = nordre
        if (ltrans) then
            do i = 1, n
                do j = 1, n
                    af(j, i) = a(i, j)
                end do
            end do
        end if
        if (ldet) then
            det = 1.d0
        else
            det = 0.d0
        end if
        if (ltrans) then
            call mgausw(af, b, dim, nordre, nb, &
                        det, lret)
        else
            call mgausw(a, b, dim, nordre, nb, &
                        det, lret)
        end if
        iret = 0
        if (.not. lret) iret = 1
    end if
!
!
!     -- 4. EN CAS DE PROBLEME : IRET > 0
!     ---------------------------------------
    if (iret .gt. 0) then
        if (lstop) then
            if (cara2(3:3) .eq. 'S') then
                if (iret .eq. n+1) then
                    call utmess('F', 'CALCULEL3_79')
                else
                    vali(1) = iret
                    vali(2) = iret
                    valk(1) = ' '
                    valk(2) = ' '
                    call utmess('F', 'CALCULEL6_15', nk=2, valk=valk, ni=2, &
                                vali=vali)
                end if
            else
                call utmess('F', 'CALCULEL3_79')
            end if
        else
!         -- ON CONTINUE
        end if
    end if
!
!
!     -- 5. IMPRESSIONS (SI LE MOT CLE 'INFO' = 2): RCOND, BERR, FERR
!     ---------------------------------------------------------------
    if (cara2(3:3) .eq. 'S') then
        call infniv(ifm, niv)
!       JMP : TROP D'IMPRESSIONS EN INFO=2
        idb = 0
        if (iret .ne. 0 .or. niv .le. 1 .or. idb .eq. 0) goto 110
        write (ifm, 1001) 'DEBUT DE MGAUSS'
        if (equed .eq. 'N') then
            write (ifm, *) 'L''EQUILIBRAGE DE LA MATRICE ''A'' '// &
                ' N''A PAS ETE NECESSAIRE'
        else if (equed .eq. 'R') then
            write (ifm, *) 'LA MATRICE ''A'' A ETE EQUILIBREE SOUS LA'//&
     &      ' FORME : DIAG(R)*A'
        else if (equed .eq. 'C') then
            write (ifm, *) 'LA MATRICE ''A'' A ETE EQUILIBREE SOUS LA'//&
     &      ' FORME : A*DIAG(C)'
        else if (equed .eq. 'B') then
            write (ifm, *) 'LA MATRICE ''A'' A ETE EQUILIBREE SOUS LA'//&
     &      ' FORME : DIAG(R)*A*DIAG(C)'
        end if
        write (ifm, *) 'ESTIMATION DE LA VALEUR DU CONDITIONNEMENT '// &
            'DE A :', rcond
!       L'ERREUR ARRIERE (BACKWARD ERROR) EST L'ERREUR FAITE EN
!       ASSIMILANT LA MATRICE 'A' AU PRODUIT 'LU'
        write (ifm, *) 'ERREUR ARRIERE : ', berr
!       L'ERREUR AVANT (FORWARD ERROR) EST LA DIFFERENCE NORMALISEE
!       ENTRE LA VALEUR CALCULEE X ET SA VALEUR EXACTE
        write (ifm, *) 'ERREUR AVANT : ', ferr
        write (ifm, 1001) 'FIN DE MGAUSS'
    end if
!
!
110 continue
!
1001 format(10('='), a, 10('='))
!
    call matfpe(1)
!
end subroutine
