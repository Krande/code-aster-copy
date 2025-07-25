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
subroutine ajlagr(rigid, masse, masinv)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mtcmbi.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/mtdscr.h"
#include "asterfort/pteddl.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: rigid, masse, masinv
!     AJOUTE LES "LAGRANGE" DANS LA MATRICE DE MASSE A PARTIR DES
!     DONNEES STOCKEES DANS LA MATRICE DE RAIDEUR.
!
! IN  : RIGID  : NOM DE LA MATRICE DE RAIDEUR
! IN  : MASSE  : NOM DE LA MATRICE DE MASSE
! OUT : MASINV : NOM DE LA MATRICE DE MASSE AVEC LES LAGRANGES
!-----------------------------------------------------------------------
    integer(kind=8) :: neq, hbloc, nbbloc, mxddl
    real(kind=8) :: zero, un, mmax, kmax, coef, lcoef(2)
    character(len=1) :: typmat, typma2, typcst(2)
    character(len=8) :: raid, mass, masi, nomddl, matrer
    character(len=14) :: numddl, nu2ddl
    character(len=19) :: rigi2, mass2, matre2, masin2
    complex(kind=8) :: cun, cmmax, ckmax, ccoef
    character(len=24) :: nmat(4), nmati
    character(len=24) :: valk(2)
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, imati, imatm, imatr, imtrer, j, jconl
    integer(kind=8) :: jmass, jraid, nbmat
    integer(kind=8), pointer :: lagr(:) => null()
    integer(kind=8), pointer :: smde(:) => null()
    character(len=24), pointer :: refa1(:) => null()
    character(len=24), pointer :: refa2(:) => null()
!
!-----------------------------------------------------------------------
    call jemarq()
    zero = 0.d0
    un = 1.d0
    cun = dcmplx(un, zero)
!
    raid = rigid
    call mtdscr(rigid)
    rigi2 = rigid
    call jeveuo(rigi2//'.&INT', 'E', imatr)
    if (zi(imatr+3) .eq. 1) then
        typmat = 'R'
    else if (zi(imatr+3) .eq. 2) then
        typmat = 'C'
    else
        call utmess('F', 'ALGORITH_3')
    end if
    call jeveuo(raid//'           .REFA', 'L', vk24=refa1)
    numddl = refa1(2) (1:14)
!
    mass = masse
    call mtdscr(masse)
    mass2 = masse
    nmat(2) = mass2//'.&INT'
    call jeveuo(nmat(2), 'E', imatm)
    if (zi(imatm+3) .eq. 1) then
        typma2 = 'R'
    else if (zi(imatm+3) .eq. 2) then
        typma2 = 'C'
    else
        call utmess('F', 'ALGORITH_3')
    end if
    call jeveuo(mass//'           .REFA', 'L', vk24=refa2)
    nu2ddl = refa2(2) (1:14)
!
    if (typma2 .ne. typmat) then
        valk(1) = typmat
        valk(2) = typma2
        call utmess('F', 'ALGORITH14_77', nk=2, valk=valk)
    end if
    if (nu2ddl .ne. numddl) then
        valk(1) = numddl
        valk(2) = nu2ddl
        call utmess('F', 'ALGORITH14_78', nk=2, valk=valk)
    end if
!
!
    call jeveuo(numddl//'.SMOS.SMDE', 'L', vi=smde)
    neq = smde(1)
    hbloc = smde(2)
    nbbloc = smde(3)
    ASSERT(nbbloc .eq. 1)
!
!     --- DETERMINATION DU COEFFICIENT DE CONDITIONNEMENT ---
    if (typmat .eq. 'R') then
        mmax = zero
        kmax = zero
        do i = 1, nbbloc
            call jeveuo(jexnum(raid//'           .VALM', i), 'L', jraid)
            call jeveuo(jexnum(mass//'           .VALM', i), 'L', jmass)
            do j = 0, hbloc-1
                mmax = max(zr(jmass+j), mmax)
                kmax = max(zr(jraid+j), kmax)
            end do
            call jelibe(jexnum(mass//'           .VALM', i))
            call jelibe(jexnum(raid//'           .VALM', i))
        end do
        coef = mmax/kmax
    else
        cmmax = zero
        ckmax = zero
        do i = 1, nbbloc
            call jeveuo(jexnum(raid//'           .VALM', i), 'L', jraid)
            call jeveuo(jexnum(mass//'           .VALM', i), 'L', jmass)
            do j = 0, hbloc-1
                cmmax = max(abs(zc(jmass+j)), abs(cmmax))
                ckmax = max(abs(zc(jraid+j)), abs(ckmax))
            end do
            call jelibe(jexnum(mass//'           .VALM', i))
            call jelibe(jexnum(raid//'           .VALM', i))
        end do
        ccoef = cmmax/ckmax
    end if
!
    matrer = '&&RIGIL'
    call mtdefs(matrer, rigid, 'V', typmat)
    call mtdscr(matrer)
    matre2 = matrer
    nmat(1) = matre2//'.&INT'
    call jeveuo(nmat(1), 'E', imtrer)
!
    call mtcmbi(typmat, imatr, coef, ccoef, imtrer)
!
    masi = masinv
    call mtdefs(masinv, rigid, 'V', typmat)
    call mtdscr(masinv)
    masin2 = masinv
    nmati = masin2//'.&INT'
    call jeveuo(nmati, 'E', imati)
!
    nbmat = 2
    nomddl = ' '
    lcoef(1) = 1.d0
    lcoef(2) = 1.d0
    typcst(1) = typmat
    typcst(2) = typmat
    call mtcmbl(nbmat, typcst, lcoef, nmat, nmati, &
                nomddl, ' ', 'ELIM=')
!
    nomddl = 'LAGR    '
    mxddl = 1
    AS_ALLOCATE(vi=lagr, size=neq*mxddl)
    call pteddl('NUME_DDL', numddl, mxddl, nomddl, neq, &
                tabl_equa=lagr)
    call jeveuo(masi//'           .CONL', 'E', jconl)
    if (typmat .eq. 'R') then
        do i = 0, neq-1
            if (lagr(1+i) .ne. 0) then
                zr(jconl+i) = mmax
            else
                zr(jconl+i) = un
            end if
        end do
    else
        do i = 0, neq-1
            if (lagr(1+i) .ne. 0) then
                zc(jconl+i) = cmmax
            else
                zc(jconl+i) = cun
            end if
        end do
    end if
!
!
! --- MENAGE
!
    AS_DEALLOCATE(vi=lagr)
    call detrsd('MATR_ASSE', '&&RIGIL')
!
    call jedema()
end subroutine
