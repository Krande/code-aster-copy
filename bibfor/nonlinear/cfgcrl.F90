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
subroutine cfgcrl(resoco, neq, nbliai, matass, solveu, &
                  alpha)
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
#include "asterfort/utmess.h"
#include "blas/ddot.h"
    character(len=24) :: resoco
    integer(kind=8) :: neq, nbliai
    character(len=19) :: matass, solveu
    real(kind=8) :: alpha
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (RESOLUTION - GCP)
!
! RECHERCHE LINEAIRE
!
! ----------------------------------------------------------------------
!
!
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  SOLVEU : SD SOLVEUR
! IN  MATASS : NOM DE LA MATRICE DU PREMIER MEMBRE ASSEMBLEE
! IN  NBLIAI : NOMBRE DE LIAISONS DE CONTACT
! IN  NEQ    : NOMBRE D'EQUATIONS
! OUT ALPHA  : COEFFICIENT DE RECHERCHE LINEAIRE
!
!
!
!
    integer(kind=8) :: ifm, niv
    real(kind=8) :: numer, denom
    integer(kind=8) :: iliai, jdecal, nbddl
    complex(kind=8) :: c16bid
    character(len=19) :: k19bla
    character(len=24) :: apcoef, apddl, appoin
    integer(kind=8) :: japcoe, japddl, japptr
    character(len=19) :: sgradp, sgrprp, direct
    integer(kind=8) :: jsgrap, jsgprp, jdirec
    character(len=24) :: secmbr, ddelt, cncin0
    integer(kind=8) :: jsecmb, jddelt
    integer(kind=8) :: iret
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
    sgradp = resoco(1:14)//'.SGDP'
    sgrprp = resoco(1:14)//'.SGPP'
!
    call jeveuo(appoin, 'L', japptr)
    call jeveuo(apcoef, 'L', japcoe)
    call jeveuo(apddl, 'L', japddl)
    call jeveuo(direct, 'L', jdirec)
    call jeveuo(sgradp, 'L', jsgrap)
    call jeveuo(sgrprp, 'L', jsgprp)
!
! --- ACCES AUX CHAMPS DE TRAVAIL
!
    secmbr = resoco(1:14)//'.SECM'
    cncin0 = resoco(1:14)//'.CIN0'
    ddelt = resoco(1:14)//'.DDEL'
    call jeveuo(ddelt(1:19)//'.VALE', 'E', jddelt)
    call jeveuo(secmbr(1:19)//'.VALE', 'E', jsecmb)
!
! --- INITIALISATIONS DES VECTEURS DE TRAVAIL
!
    call r8inir(neq, 0.d0, zr(jsecmb), 1)
    call r8inir(neq, 0.d0, zr(jddelt), 1)
!
! --- SECOND MEMBRE: [A]T .{DIRECP}
!
    do iliai = 1, nbliai
        jdecal = zi(japptr+iliai-1)
        nbddl = zi(japptr+iliai)-zi(japptr+iliai-1)
        call calatm(neq, nbddl, zr(jdirec+iliai-1), zr(japcoe+jdecal), zi(japddl+jdecal), &
                    zr(jsecmb))
    end do
!
! --- RESOLUTION [K].{DDELT} = [A]T .{DIRECP} -> {DDELT}
!
    call resoud(matass, k19bla, solveu, cncin0, 0, &
                secmbr, ddelt, 'V', [0.d0], [c16bid], &
                k19bla, .true._1, 0, iret)
!
! --- PRODUIT SCALAIRE  NUMER = <DIRECP>.{DIRECP}
!
    b_n = to_blas_int(nbliai)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    numer = ddot(b_n, zr(jsgprp), b_incx, zr(jsgrap), b_incy)
!
! --- PRODUIT SCALAIRE  DENOM = <DIRECP>.[A].[K]-1.[A]T .{DIRECP}
!
    call jeveuo(ddelt(1:19)//'.VALE', 'L', jddelt)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    denom = ddot(b_n, zr(jddelt), b_incx, zr(jsecmb), b_incy)
!
    if (denom .lt. 0.d0) then
        call utmess('A', 'CONTACT_7')
    end if
!
! --- COEFFICIENT DE RECHERCHE LINEAIRE
!
    alpha = numer/denom
!
! --- AFFICHAGE
!
    if (niv .eq. 2) then
        write (ifm, 9040) alpha
    end if
!
9040 format(' <CONTACT><CALC> PAS D''AVANCEMENT INITIAL : ', 1pe12.5)
!
    call jedema()
!
end subroutine
