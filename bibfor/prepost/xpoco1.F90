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
subroutine xpoco1(dirma, nbma, dirno, nbno, ma1, &
                  ma2, jnivgr)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
!
    character(len=8) :: ma1, ma2
    integer(kind=8) :: nbma, dirma(nbma), nbno, dirno(nbno), jnivgr
!
!   COPIE DANS LE MAILLAGE MA2 DES MAILLES ET DES NOEUDS DU MAILLAGE MA1
!   CONTENUS DANS LES TABLEAUX D'INDIRECTION DIRMA ET DIRNO
!
!   IN
!       DIRMA : TABLEAU DE CORRESPONDANCE DES NUMEROS DE MAILLES
!       DIRNO : TABLEAU DE CORRESPONDANCE DES NUMEROS DE NOEUDS
!       NBMA  : LONGUEUR DE DIRMA
!       NBNO  : LONGUEUR DE DIRNO
!       MA1   : NOM DU MAILLAGE SAIN
!      JNIVGR : ADRESSE DU VECTEUR DE REMPLISSAGE DES GROUP_MA DE MAXFEM
!
!   OUT
!       MA2   : NOM DU MAILLAGE FISSURE
!
!
!
    integer(kind=8) :: i, j, ino1, iret, nbgn, iagno
    integer(kind=8) :: iacon1, n, iacon2
    integer(kind=8) :: ino2, nbgm2, i1, i2, iagma1, iagma2, n1, n2, ima
    character(len=7) :: ch7
    character(len=8) :: noma2, nono2
    character(len=24) :: nogma
    integer(kind=8), pointer :: typm1(:) => null()
    integer(kind=8), pointer :: typm2(:) => null()
    real(kind=8), pointer :: coo1(:) => null()
    real(kind=8), pointer :: coo2(:) => null()
!
    call jemarq()
!
!     RECUP DES .TYPMAIL, .COORDO DU MAILLAGE 1 ET 2
    call jeveuo(ma1//'.TYPMAIL', 'L', vi=typm1)
    call jeveuo(ma2//'.TYPMAIL', 'E', vi=typm2)
    call jeveuo(ma1//'.COORDO    .VALE', 'L', vr=coo1)
    call jeveuo(ma2//'.COORDO    .VALE', 'E', vr=coo2)
!
    call jeexin(ma2//'.GROUPENO', iret)
    nbgn = 0
    if (iret .gt. 0) call jelira(ma2//'.GROUPENO', 'NUTIOC', nbgn)
!
    call jeexin(ma2//'.GROUPEMA', iret)
    nbgm2 = 0
    if (iret .gt. 0) call jelira(ma2//'.GROUPEMA', 'NUTIOC', nbgm2)
!
!     ---------------------------------------------------------------
!     COPIE DES VECTEURS
!     ---------------------------------------------------------------
!
!     .NOMMAI ET .TYPMAIL
    do i = 1, nbma
        if (dirma(i) .ne. 0) then
            call codent(i, 'G', ch7)
            noma2 = 'M'//ch7
            call jecroc(jexnom(ma2//'.NOMMAI', noma2))
            typm2(dirma(i)) = typm1(i)
        end if
    end do
!
!     .NOMNOE
    do i = 1, nbno
        if (dirno(i) .ne. 0) then
            call codent(i, 'G', ch7)
            nono2 = 'N'//ch7
            call jecroc(jexnom(ma2//'.NOMNOE', nono2))
        end if
    end do
!
!     .COORDO
    do i = 1, nbno
        if (dirno(i) .ne. 0) then
            do j = 1, 3
                coo2(3*(dirno(i)-1)+j) = coo1(3*(i-1)+j)
            end do
        end if
    end do
!
!     .CONNEX
    do i = 1, nbma
        if (dirma(i) .ne. 0) then
            call jeveuo(jexnum(ma1//'.CONNEX', i), 'L', iacon1)
            call jelira(jexnum(ma1//'.CONNEX', i), 'LONMAX', n)
            call jeecra(jexnum(ma2//'.CONNEX', dirma(i)), 'LONMAX', n)
            call jeveuo(jexnum(ma2//'.CONNEX', dirma(i)), 'E', iacon2)
            do j = 1, n
                ino1 = zi(iacon1-1+j)
                ino2 = dirno(ino1)
                zi(iacon2-1+j) = ino2
            end do
        end if
    end do
!
!     .GROUPENO
    do i = 1, nbgn
        call jeveuo(jexnum(ma2//'.GROUPENO', i), 'E', iagno)
        call jelira(jexnum(ma2//'.GROUPENO', i), 'LONUTI', n)
        do j = 1, n
            if (dirno(zi(iagno-1+j)) .ne. 0) then
                zi(iagno-1+j) = dirno(zi(iagno-1+j))
            end if
        end do
    end do
!
!     .GROUPEMA
    do i2 = 1, nbgm2
        call jenuno(jexnum(ma2//'.GROUPEMA', i2), nogma)
        call jenonu(jexnom(ma1//'.GROUPEMA', nogma), i1)
        call jeveuo(jexnum(ma1//'.GROUPEMA', i1), 'L', iagma1)
        call jelira(jexnum(ma1//'.GROUPEMA', i1), 'LONUTI', n1)
        call jeveuo(jexnum(ma2//'.GROUPEMA', i2), 'E', iagma2)
        call jelira(jexnum(ma2//'.GROUPEMA', i2), 'LONUTI', n2)
        do i = 1, n1
            ima = zi(iagma1-1+i)
            if (dirma(ima) .ne. 0) then
!           NIVEAU DE REMPLISSAGE DU GROUP_MA
                zi(jnivgr-1+i2) = zi(jnivgr-1+i2)+1
                zi(iagma2-1+zi(jnivgr-1+i2)) = dirma(ima)
            end if
        end do
        ASSERT(zi(jnivgr-1+i2) .le. n2)
    end do
!
    call jedema()
end subroutine
