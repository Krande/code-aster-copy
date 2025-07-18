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
subroutine gcpc(m, in, ip, ac, inpc, &
                perm, ippc, acpc, bf, xp, &
                r, rr, p, irep, niter, &
                epsi, criter, solveu, matas, istop, &
                iret)
!     RESOLUTION D'UN SYSTEME LINEAIRE SYMETRIQUE PAR UNE METHODE DE
!     GRADIENT CONJUGUE PRECONDITIONNE
!               LA MATRICE EST STOCKEE SOUS FORME COMPACTE (IN,IP,AC)
!    -------------------------------------------------------------------
!    . M             -->   NOMBRE DE COLONNES DE LA MATRICE
!    . IN            -->   POINTEUR DE FIN DE COLONNE DE LA MATRICE
!    . IP            -->   TABLEAU DES NUMEROS DE LIGNE
!    . AC            -->   TABLEAU DES COEFFICIENTS DE LA MATRICE
!
!    . INPC          -->   IDEM IN POUR MATRICE DE PRECOND.
!    . IPPC          -->   IDEM IP POUR MATRICE DE PRECOND.
!    . ACPC          -->   IDEM AC POUR MATRICE DE PRECOND.
!    . BF            -->   VECTEUR SECOND MEMBRE
!    . XP           <-->   VECTEUR SOLUTION
!    . R            <--    VECTEUR RESIDU
!    . RR           <--    DIRECTION DE DESCENTE AVANT CONJUGAISON
!    . P            <--    DIRECTION DE DESCENTE APRES CONJUGAISON
!    -------------------------------------------------------------------
!    . IREP          -->    0  XP INITIAL MIS A ZERO
!                           1  XP INITIAL DONNEE DE GCPC
!    -------------------------------------------------------------------
!    . NITER         -->   NOMBRE MAXIMUM D'ITERATIONS
!    . EPSI          -->   CRITERE DE CONVERGENCE
!    . CRITER        -->   SD_CRITER (CRITERES DE CONVERGENCE)
!    -------------------------------------------------------------------
!    . SOLVEU        -->   SD_SOLVEUR (POUR LDLT_SP)
!    . MATASS        -->   MATRICE ASSEMBLEE DU SYSTEME (POUR LDLT_SP)
!     ------------------------------------------------------------------
!     - PRECAUTIONS D'EMPLOI:  XP PEUT ETRE EVENTUELLEMENT CONFONDU
!                              AVEC BF SI MEME ARGUMENT
!     ------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "jeveux.h"
#include "asterc/matfpe.h"
#include "asterc/r8prem.h"
#include "asterfort/amumph.h"
#include "asterfort/assert.h"
#include "asterfort/crsvfm.h"
#include "asterfort/detrsd.h"
#include "asterfort/gcax.h"
#include "asterfort/gcldm1.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dnrm2.h"
#include "blas/dscal.h"
#include "asterfort/isParallelMatrix.h"
!
    integer(kind=4) :: ip(*), ippc(*)
    integer(kind=8) :: m, in(m), inpc(m), irep, niter, perm(m)
    real(kind=8) :: ac(m), acpc(m), bf(m), xp(m), r(m), rr(m), p(m), epsi
    character(len=19) :: criter, matas, solveu
    integer(kind=8) :: istop, iret
!
! -----------------------------------------------------------------
    real(kind=8) :: zero, bnorm, anorm, epsix, anormx, rrri, gama, rrrim1
    real(kind=8) :: paraaf, anorxx, rau, valr(2), blreps
    integer(kind=8) :: ifm, niv, jcri, jcrr, jcrk, iter, ier, vali, pcpiv, redmpi
    character(len=24) :: precon, solvbd, usersm, renum
    character :: prec, rank
    complex(kind=8) :: cbid
    integer(kind=8), pointer :: slvi(:) => null()
    character(len=24), pointer :: slvk(:) => null()
    real(kind=8), pointer :: slvr(:) => null()
    real(kind=8), pointer :: xtrav(:) => null()
    real(kind=8), pointer :: ytrav(:) => null()
    aster_logical :: l_parallel_matrix
    blas_int :: b_incx, b_incy, b_n
! -----------------------------------------------------------------
!
    cbid = (0.d0, 0.d0)
    call jemarq()
!
    call matfpe(-1)
    iter = 0
!
!-----RECUPERATION DU NIVEAU D'IMPRESSION
    call infniv(ifm, niv)
!
! --- GCPC interdit avec un ParallelMesh
    l_parallel_matrix = isParallelMatrix(matas)
    if (l_parallel_matrix) then
        call utmess('F', 'ALGELINE4_44')
    end if
!
!-----PARAMETRE D'AFFICHAGE DE LA DECROISSANCE DU RESIDU
!     (SI ON GAGNE PARAAF * 100%)
    paraaf = 0.1d0
!
!-----INITS DIVERS
    iret = 0
    zero = 0.d0
    ASSERT(irep .eq. 0 .or. irep .eq. 1)
!
!-----RECUPERATION DU PRECONDITIONNEUR
!  -- CREATION DE LA SD SOLVEUR MUMPS SIMPLE PRECISION/LOW_RANK
!  -- (A DETRUIRE A LA SORTIE)
    call jeveuo(solveu//'.SLVK', 'L', vk24=slvk)
    call jeveuo(solveu//'.SLVI', 'L', vi=slvi)
    call jeveuo(solveu//'.SLVR', 'L', vr=slvr)
    redmpi = slvi(1)
    precon = slvk(2)
    usersm = slvk(9)
    pcpiv = slvi(7)
    blreps = slvr(4)
    solvbd = slvk(3)
    renum = slvk(4)
    if (precon == 'LDLT_SP') then
        prec = 'S'
    else if (precon == 'LDLT_DP') then
        prec = 'D'
    end if
    if (abs(blreps) < r8prem()) then
        rank = 'F'
    else
        rank = 'L'
    end if
    if ((precon == 'LDLT_SP') .or. (precon == 'LDLT_DP')) then
        call crsvfm(solvbd, matas, prec, rank, pcpiv, &
                    usersm, blreps, renum, redmpi)
    end if
!-----Pour tenir compte de la renumerotation de la matrice de preconditionnement (LDLT):
    if (precon .eq. 'LDLT_INC') then
        AS_ALLOCATE(vr=xtrav, size=m)
        AS_ALLOCATE(vr=ytrav, size=m)
    end if
!
!-----CALCULS PRELIMINAIRES
!
!      ---- CALCUL DE NORME DE BF
    b_n = to_blas_int(m)
    b_incx = to_blas_int(1)
    bnorm = dnrm2(b_n, bf, b_incx)
    if (bnorm .eq. zero) then
        call r8inir(m, zero, xp, 1)
!        WRITE (IFM,*)'>>>>>>> SECOND MEMBRE = 0 DONC SOLUTION = 0 '
        goto 80
    end if
!
    if (irep .eq. 0) then
!       ---- INITIALISATION X1 = 0    ===>   CALCUL DE R1 = A*X0 - B
        call r8inir(m, zero, xp, 1)
        b_n = to_blas_int(m)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, bf, b_incx, r, b_incy)
        b_n = to_blas_int(m)
        b_incx = to_blas_int(1)
        call dscal(b_n, -1.d0, r, b_incx)
        anorm = bnorm
        epsix = epsi*anorm
        if (niv .eq. 2) write (ifm, 101) anorm, epsix, epsi
    else
!       ---- INITIALISATION PAR X PRECEDENT: CALCUL DE R1 = A*X1 - B
        call gcax(m, in, ip, ac, xp, &
                  r)
        b_n = to_blas_int(m)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, bf, b_incx, r, &
                   b_incy)
        b_n = to_blas_int(m)
        b_incx = to_blas_int(1)
        anorm = dnrm2(b_n, r, b_incx)
        epsix = epsi*anorm
        if (niv .eq. 2) write (ifm, 102) anorm, epsix, epsi
    end if
!
    call jeexin(criter//'.CRTI', ier)
    if (ier .eq. 0) then
        if (criter .ne. ' ') then
            call wkvect(criter//'.CRTI', 'V V I', 1, jcri)
            call wkvect(criter//'.CRTR', 'V V R8', 1, jcrr)
            call wkvect(criter//'.CRDE', 'V V K16', 2, jcrk)
            zk16(jcrk) = 'ITER_GCPC'
            zk16(jcrk+1) = 'RESI_GCPC'
        else
            jcri = 0
        end if
    else
        call jeveuo(criter//'.CRTI', 'E', jcri)
        call jeveuo(criter//'.CRTR', 'E', jcrr)
    end if
!
! ---- ITERATIONS
    anormx = anorm
    anorxx = anorm
!
    do iter = 1, niter
!       ---- PRECONDITIONNEMENT DU RESIDU:
!                                             ZK = (LDLT)-1. RK
!                                                   RK <--- R()
!                                                  ZK <--- RR()
        if (precon .eq. 'LDLT_INC') then
            call gcldm1(m, inpc, ippc, acpc, r, &
                        rr, perm, xtrav, ytrav)
        else if ((precon .eq. 'LDLT_SP') .or. (precon .eq. 'LDLT_DP')) then
            b_n = to_blas_int(m)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, r, b_incx, rr, b_incy)
!         ON PASSE ' ' AU LIEU DE VCINE, DEJA PRIS EN COMPTE DANS RESGRA
            call amumph('RESOUD', solvbd, matas, rr, [cbid], &
                        ' ', 1, ier, .true._1)
        else
            ASSERT(.false.)
        end if
!
!                                             RRRI <--- (RK,ZK)
        b_n = to_blas_int(m)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        rrri = ddot(b_n, r, b_incx, rr, b_incy)
!       ---- NOUVELLE DIRECTION DE DESCENTE:
!                                    BETAK = (RK,ZK)/(RK-1,ZK-1)
!                                               BETAK <--- GAMA
!                                        PK = BETAK * PK-1 + ZK
!                                                   PK <--- P()
        if (iter .gt. 1) then
            gama = rrri/rrrim1
            b_n = to_blas_int(m)
            b_incx = to_blas_int(1)
            call dscal(b_n, gama, p, b_incx)
            b_n = to_blas_int(m)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, rr, b_incx, p, &
                       b_incy)
        else
            b_n = to_blas_int(m)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, rr, b_incx, p, b_incy)
        end if
        rrrim1 = rrri
!
!       ---- NOUVEAUX RESIDU ET DEPLACEMENT:
!                       ZZK = A.PK ET ALPHAK = -(RK,ZK)/(PK,ZZK)
!                                       XK+1 = XK + ALPHAK * PK
!                                      RK+1 = RK + ALPHAK * ZZK
!                                                 ZZK <--- RR()
!                                                 XK  <--- XP()
        call gcax(m, in, ip, ac, p, &
                  rr)
        b_n = to_blas_int(m)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        rau = -rrri/ddot(b_n, p, b_incx, rr, b_incy)
        b_n = to_blas_int(m)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, rau, p, b_incx, xp, &
                   b_incy)
        b_n = to_blas_int(m)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, rau, rr, b_incx, r, &
                   b_incy)
!
!       ---- CALCUL TEST D'ARRET ET AFFICHAGE
        b_n = to_blas_int(m)
        b_incx = to_blas_int(1)
        anorm = dnrm2(b_n, r, b_incx)
        if (anorm .le. anormx*paraaf) then
            if (niv .eq. 2) write (*, 104) iter, anorm, anorm/anorxx
            anormx = anorm
        end if
        if (niv .eq. 3) write (ifm, 104) iter, anorm, anorm/anorxx
!
!       --- TEST DE CONVERGENCE
        if (anorm .lt. epsix) then
            if (niv .eq. 2) write (ifm, 103) anorxx, anorm, anorm/anorxx
            if (niv .eq. 2) write (ifm, 105) iter
            if (jcri .ne. 0) then
                zi(jcri) = iter
                zr(jcrr) = anorm
            end if
            goto 80
        end if
    end do
!
!
!
!        ---  NON CONVERGENCE
    vali = iter
    valr(1) = anorm/anorxx
    valr(2) = epsi
    if (istop == 0) then
!            ERREUR <F>
        select case (precon)
        case ('LDLT_INC')
            call utmess('F', 'ALGELINE4_3', si=vali, nr=2, valr=valr)
        case ('LDLT_SP', 'LDLT_DP')
            call utmess('F', 'ALGELINE4_6', si=vali, nr=2, valr=valr)
        case default
            ASSERT(.false.)
        end select
    else if (istop == 2) then
!            ON CONTINUE EN RETOURNANT UN CODE D'ERREUR IRET=1
        iret = 1
        goto 80
    else
        ASSERT(.false.)
    end if
!    -----------
101 format(/'   * GCPC   NORME DU RESIDU =', d11.4,&
     &       '  (INITIALISATION PAR X = ZERO)', /,&
     &'   *        NORME DU RESIDU A ATTEINDRE EN ABS/RELA=',&
     &d11.4, d11.4,/)
102 format(/'   * GCPC   NORME DU RESIDU =', d11.4,&
     &       '  (INITIALISATION PAR X PRECEDENT)', /,&
     & '   *        NORME DU RESIDU A ATTEINDRE EN ABS/RELA=',&
     & d11.4, d11.4)
103 format('   * NORME DU RESIDU INITIAL/FINAL/RELATIF=',&
     &         d11.4, d11.4, d11.4)
104 format('   * ITERATION', i5, ' NORME DU RESIDU EN ABS/RELA =',&
     &         d11.4, d11.4)
105 format(1x, /, 2x, 32('*')/'  * CONVERGENCE EN ', i4,&
     &       ' ITERATIONS'/2x, 32('*'),/)
!    -----------
80  continue
!
! --  DESTRUCTION DE LA SD SOLVEUR MUMPS SIMPLE PRECISION
    if ((precon .eq. 'LDLT_SP') .or. (precon .eq. 'LDLT_DP')) then
        call detrsd('SOLVEUR', solvbd)
!       ON STOCKE LE NOMBRE D'ITERATIONS DU GCPC
        call jeveuo(solveu//'.SLVI', 'E', vi=slvi)
        slvi(5) = iter
    end if
!
    call matfpe(1)
!
    if (precon .eq. 'LDLT_INC') then
        AS_DEALLOCATE(vr=xtrav)
        AS_DEALLOCATE(vr=ytrav)
    end if
    call jedema()
!
end subroutine
