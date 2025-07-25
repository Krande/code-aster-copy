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
subroutine cfgcpr(resoco, matass, solveu, neq, nbliai, &
                  search, alpha)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/calatm.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/r8inir.h"
#include "asterfort/resoud.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    character(len=24) :: resoco
    character(len=16) :: search
    integer(kind=8) :: neq, nbliai
    character(len=19) :: matass, solveu
    real(kind=8) :: alpha
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (RESOLUTION - GCP)
!
! PROJECTION DU PAS D'AVANCEMENT
!
! ----------------------------------------------------------------------
!
!
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  SOLVEU : SD SOLVEUR
! IN  MATASS : NOM DE LA MATRICE DU PREMIER MEMBRE ASSEMBLEE
! IN  NBLIAI : NOMBRE DE LIAISONS DE CONTACT
! IN  NEQ    : NOMBRE D'EQUATIONS
! IN  SEARCH : RECHERCHE LINEAIRE
!              'ADMISSIBLE'
!              'NON_ADMISSIBLE'
! I/O ALPHA  : COEFFICIENT DE RECHERCHE LINEAIRE
!
!
!
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iliai, jdecal, nbddl
    complex(kind=8) :: c16bid
    character(len=19) :: k19bla
    character(len=24) :: apcoef, apddl, appoin
    integer(kind=8) :: japcoe, japddl, japptr
    character(len=19) :: direct
    integer(kind=8) :: jdirec
    character(len=19) :: ddeplc, ddepl0, ddelt
    character(len=24) :: secmbr, cncin0
    integer(kind=8) :: jsecmb
    character(len=19) :: mu
    integer(kind=8) :: jmu
    integer(kind=8) :: iret
    real(kind=8), pointer :: vddelt(:) => null()
    real(kind=8), pointer :: ddep0(:) => null()
    real(kind=8), pointer :: ddepc(:) => null()
    blas_int :: b_incx, b_incy, b_n
    c16bid = dcmplx(0.d0, 0.d0)
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    k19bla = ' '
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    appoin = resoco(1:14)//'.APPOIN'
    apcoef = resoco(1:14)//'.APCOEF'
    apddl = resoco(1:14)//'.APDDL'
    direct = resoco(1:14)//'.DIRE'
    mu = resoco(1:14)//'.MU'
!
    call jeveuo(appoin, 'L', japptr)
    call jeveuo(apcoef, 'L', japcoe)
    call jeveuo(apddl, 'L', japddl)
    call jeveuo(direct, 'L', jdirec)
    call jeveuo(mu, 'E', jmu)
!
! --- ACCES AUX CHAMPS DE TRAVAIL
! --- DDEPL0: INCREMENT DE SOLUTION SANS CORRECTION DU CONTACT
! --- DDEPLC: INCREMENT DE SOLUTION APRES CORRECTION DU CONTACT
! --- DDELT : INCREMENT DE SOLUTION ITERATION DE CONTACT
!
    ddepl0 = resoco(1:14)//'.DEL0'
    ddeplc = resoco(1:14)//'.DELC'
    ddelt = resoco(1:14)//'.DDEL'
    call jeveuo(ddepl0(1:19)//'.VALE', 'L', vr=ddep0)
    call jeveuo(ddeplc(1:19)//'.VALE', 'E', vr=ddepc)
    call jeveuo(ddelt(1:19)//'.VALE', 'E', vr=vddelt)
!
! --- ACCES AUX CHAMPS DE TRAVAIL
!
    secmbr = resoco(1:14)//'.SECM'
    cncin0 = resoco(1:14)//'.CIN0'
    call jeveuo(secmbr(1:19)//'.VALE', 'E', jsecmb)
!
! --- RECALCUL DE ALPHA POUR UNE SOLUTION ADMISSIBLE
!
    if (search .eq. 'ADMISSIBLE') then
        do iliai = 1, nbliai
            if (zr(jdirec-1+iliai) .lt. 0.d0) then
                alpha = min(alpha, -zr(jmu+iliai-1)/zr(jdirec-1+iliai))
            end if
        end do
        if (niv .eq. 2) then
            write (ifm, 9050) alpha
        end if
    end if
!
! --- MISE A JOUR DE MU
!
    b_n = to_blas_int(nbliai)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, alpha, zr(jdirec), b_incx, zr(jmu), &
               b_incy)
!
! --- DESACTIVATION DE MU POUR UNE SOLUTION NON-ADMISSIBLE
!
    if (search .eq. 'NON_ADMISSIBLE') then
        do iliai = 1, nbliai
            if (zr(jmu-1+iliai) .lt. 0.d0) then
                zr(jmu-1+iliai) = 0.d0
            end if
        end do
    end if
!
! --- RECALCUL D'UN SOLUTION
!
    if (search .eq. 'NON_ADMISSIBLE') then
!
        call r8inir(neq, 0.d0, zr(jsecmb), 1)
        call r8inir(neq, 0.d0, vddelt, 1)
!
! ----- SECOND MEMBRE: [A]T.{MU}
!
        do iliai = 1, nbliai
            jdecal = zi(japptr+iliai-1)
            nbddl = zi(japptr+iliai)-zi(japptr+iliai-1)
            call calatm(neq, nbddl, zr(jmu-1+iliai), zr(japcoe+jdecal), zi(japddl+jdecal), &
                        zr(jsecmb))
        end do
!
! ----- RESOLUTION [K].{DDELT} = [A]T.{MU} -> {DDELT}
!
        call resoud(matass, k19bla, solveu, cncin0, 0, &
                    secmbr, ddelt, 'V', [0.d0], [c16bid], &
                    k19bla, .true._1, 0, iret)
!
! ----- RECOPIE DE LA SPOLUTION SANS CONTACT
!
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, ddep0, b_incx, ddepc, b_incy)
        alpha = 1.d0
    end if
!
9050 format(' <CONTACT><CALC> PAS D''AVANCEMENT APRES PROJECTION : ',&
&       1pe12.5)
!
    call jedema()
!
end subroutine
