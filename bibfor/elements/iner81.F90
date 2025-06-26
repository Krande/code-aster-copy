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
subroutine iner81(nomres, classe, basmod, nommat)
    implicit none
!
!  BUT : CALCUL DES FORCES D'INERTIES SUR BASE MODALE
!
!-----------------------------------------------------------------------
!
! NOMRES /I/ : NOM K19 DE LA MATRICE CARREE RESULTAT
! CLASSE /I/ : CLASSE DE LA BASE JEVEUX DE L'OBJET RESULTAT
! BASMOD /I/ : NOM UT DE LA BASE MODALE DE PROJECTION
! NOMMAT /I/ : NOM K8 DE LA MATRICE A PROJETER
!
!
!
!
#include "jeveux.h"
#include "asterfort/copmod.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/mtexis.h"
#include "asterfort/pteddl.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ia, iad, idbase, ieq
    integer(kind=8) :: ier, if, ldref, ldres, lmat, ltvec1
    integer(kind=8) :: ltvec2, ltvec3, mxddl, nbdef, neq
!-----------------------------------------------------------------------
    parameter(mxddl=6)
    character(len=8), parameter :: nomddl(mxddl) = ['DX      ', 'DY      ', 'DZ      ', &
                                                    'DRX     ', 'DRY     ', 'DRZ     ']
    character(len=1) :: classe
    character(len=6), parameter :: pgc = 'INER81'
    character(len=19) :: nommat
    character(len=8) :: basmod
    character(len=14) :: num
    character(len=24) :: nomres
    character(len=24) :: valk
    complex(kind=8) :: cbid
    integer(kind=8), pointer :: deeq(:) => null()
    blas_int :: b_incx, b_incy, b_n
    cbid = dcmplx(0.d0, 0.d0)
!
! --- CREATION DU .REFE
!
    call jemarq()
    call wkvect(nomres(1:18)//'_REFE', 'G V K24', 2, ldref)
    zk24(ldref) = basmod
    zk24(ldref+1) = nommat
!
! --- NOMBRE TOTAL DE MODES ET DEFORMEES
!
    call dismoi('NB_MODES_TOT', basmod, 'RESULTAT', repi=nbdef)
!
!
! --- ALLOCATION DE LA MATRICE RESULTAT
!
    call wkvect(nomres(1:18)//'_VALE', classe//' V R', 3*nbdef, ldres)
!
! --- CONTROLE D'EXISTENCE DE LA MATRICE
!
    call mtexis(nommat(1:8), ier)
    if (ier .eq. 0) then
        valk = nommat(1:8)
        call utmess('F', 'ALGORITH12_39', sk=valk)
    end if
!
! --- ALLOCATION DESCRIPTEUR DE LA MATRICE
!
    call mtdscr(nommat(1:8))
    call jeveuo(nommat(1:19)//'.&INT', 'E', lmat)
!
! --- RECUPERATION NUMEROTATION ET NB EQUATIONS
!
    call dismoi('NB_EQUA', nommat(1:8), 'MATR_ASSE', repi=neq)
    call dismoi('NOM_NUME_DDL', nommat(1:8), 'MATR_ASSE', repk=num)
!
! --- ALLOCATION VECTEURS DE TRAVAIL
!
    call wkvect('&&'//pgc//'.VECT1', 'V V R', neq, ltvec1)
    call wkvect('&&'//pgc//'.VECT2', 'V V R', neq, ltvec2)
    call wkvect('&&'//pgc//'.VECT3', 'V V I', mxddl*neq, ltvec3)
    call pteddl('NUME_DDL', num, mxddl, nomddl, neq, &
                tabl_equa=zi(ltvec3))
!
    call jeveuo(num//'.NUME.DEEQ', 'L', vi=deeq)
    call wkvect('&&'//pgc//'.BASEMO', 'V V R', nbdef*neq, idbase)
    call copmod(basmod, numer=num, bmodr=zr(idbase))
!
! --- CALCUL DES FORCES D'INERTIES
!
    do if = 1, 3
!
!     --- MODE RIGIDE EN DX , DY , DZ
!
        ia = (if-1)*neq
        do ieq = 0, neq-1
            zr(ltvec1+ieq) = zi(ltvec3+ia+ieq)
        end do
!
!     --- MULTIPLICATION DU MODE RIGIDE PAR LA MATRICE MASSE
!
        call mrmult('ZERO', lmat, zr(ltvec1), zr(ltvec2), 1, &
                    .true._1)
!
!     --- PROJECTION SUR LES MODES PROPRES ET LES DEFORMEES NON MODALES
!
        iad = (if-1)*nbdef
        do i = 1, nbdef
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(idbase+(i-1)*neq), b_incx, zr(ltvec1), b_incy)
            call zerlag(neq, deeq, vectr=zr(ltvec1))
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            zr(ldres+iad+i-1) = ddot(b_n, zr(ltvec1), b_incx, zr(ltvec2), b_incy)
        end do
!
    end do
!
! --- DESTRUCTION VECTEURS DE TRAVAIL
!
    call jedetr('&&'//pgc//'.BASEMO')
    call jedetr('&&'//pgc//'.VECT1')
    call jedetr('&&'//pgc//'.VECT2')
    call jedetr('&&'//pgc//'.VECT3')
!
    call jedema()
end subroutine
