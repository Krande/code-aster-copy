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
subroutine asmasu(ma1, ma2, mag)
    implicit none
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/infniv.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
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
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: ma1, ma2, mag
!
!     OPERATEUR: ASSE_MAILLAGE / CAS DE L ASSEMBLAGE DE MAILLAGES
!     AVEC SUPERPOSITION
!
!-----------------------------------------------------------------------
!
    character(len=1) :: kkk
    character(len=19) :: coordo
    character(len=24) :: nogma, nogmab, nogno, nognob
    integer(kind=8) :: nbma, nbm1, nbm2, nbno, nbn1, nbn2, nbgma, nbgm1, nbgm2
    integer(kind=8) :: i1, icompt, ino, l1, l2, l3, i, n, ncoor, k, ifm, niv
    integer(kind=8) :: iadime
    integer(kind=8) :: iagma1, iagma2, iagmax
    integer(kind=8) :: iacon1, iacon2, iaconx
    integer(kind=8) :: iagno1, iagno2, iagnox
    integer(kind=8) :: iatyp1, iatyp2, iatypx
    integer(kind=8) :: nbgno, nbgn1, nbgn2, ii, igeomr, iadesc, ibid
    integer(kind=8) :: iatyma, iavale, iret, iret1, iret2
    real(kind=8), pointer :: coo1(:) => null()
    real(kind=8), pointer :: coo2(:) => null()
    integer(kind=8), pointer :: dim1(:) => null()
    integer(kind=8), pointer :: dim2(:) => null()
!
!     ------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
!     --OBJET .DIME :
!     ---------------
    call jeveuo(ma1//'.DIME', 'L', vi=dim1)
    call jeveuo(ma2//'.DIME', 'L', vi=dim2)
    call wkvect(mag//'.DIME', 'G V I', 6, iadime)
!CC SOMME POUR : 1 LE NB DE NOEUDS,
!CC              2       DE NOEUDS LAGRANGES,
!CC              3       DE MAILLES,
!CC              4       DE SUPER MAILLES
!CC              5       DU MAJORANT DE SUPER MAILLES
    do i = 1, 5
        zi(iadime-1+i) = dim1(i)+dim2(i)
    end do
!
    ncoor = max(dim1(6), dim2(6))
    zi(iadime-1+6) = ncoor
!
    nbma = zi(iadime-1+3)
    nbm1 = dim1(3)
    nbm2 = dim2(3)
!
    nbno = zi(iadime-1+1)
    nbn1 = dim1(1)
    nbn2 = dim2(1)
!
!
!     --OBJET .CONNEX:
!     -----------------
    if (nbma .gt. 0) then
        call jecrec(mag//'.CONNEX', 'G V I', 'NU', 'CONTIG', 'VARIABLE', &
                    nbma)
        call jeecra(mag//'.CONNEX', 'NUTIOC', nbma)

        l1 = 0
        l2 = 0
        if (nbm1 .gt. 0) call jelira(ma1//'.CONNEX', 'LONT', l1)
        if (nbm2 .gt. 0) call jelira(ma2//'.CONNEX', 'LONT', l2)
        l3 = l1+l2
        call jeecra(mag//'.CONNEX', 'LONT', l3)
        do i = 1, nbm1
            call jeveuo(jexnum(ma1//'.CONNEX', i), 'L', iacon1)
            call jelira(jexnum(ma1//'.CONNEX', i), 'LONMAX', n)
            call jeecra(jexnum(mag//'.CONNEX', i), 'LONMAX', n)
            call jeveuo(jexnum(mag//'.CONNEX', i), 'E', iaconx)
            do ii = 1, n
                zi(iaconx-1+ii) = zi(iacon1-1+ii)
            end do
        end do
        do i = 1, nbm2
            i1 = i+nbm1
            call jeveuo(jexnum(ma2//'.CONNEX', i), 'L', iacon2)
            call jelira(jexnum(ma2//'.CONNEX', i), 'LONMAX', n)
            call jeecra(jexnum(mag//'.CONNEX', i1), 'LONMAX', n)
            call jeveuo(jexnum(mag//'.CONNEX', i1), 'E', iaconx)
            do ii = 1, n
                zi(iaconx-1+ii) = zi(iacon2-1+ii)+nbn1
            end do
        end do
    end if
!
!     -- CREATION DU CHAMP .COORDO :
!     ------------------------------
    coordo = mag//'.COORDO'
!
    call jenonu(jexnom('&CATA.GD.NOMGD', 'GEOM_R'), igeomr)
    call wkvect(coordo//'.DESC', 'G V I', 3, iadesc)
    call jeecra(coordo//'.DESC', 'DOCU', ibid, 'CHGO')
    zi(iadesc-1+1) = igeomr
!     -- TOUJOURS 3 COMPOSANTES X, Y ET Z
    zi(iadesc-1+2) = -3
!     -- 14 = 2**1 + 2**2 + 2**3
    zi(iadesc-1+3) = 14
!
    call jeveuo(ma1//'.COORDO    .VALE', 'L', vr=coo1)
    call jeveuo(ma2//'.COORDO    .VALE', 'L', vr=coo2)
    call wkvect(coordo//'.VALE', 'G V R', 3*nbno, iavale)
!     -- COORDONNEES DES NOEUDS :
    do ino = 1, nbn1
        do k = 1, 3
            zr(iavale-1+3*(ino-1)+k) = coo1(3*(ino-1)+k)
        end do
    end do
    do ino = 1, nbn2
        do k = 1, 3
            zr(iavale-1+3*(nbn1+ino-1)+k) = coo2(3*(ino-1)+k)
        end do
    end do
!
!
!     --OBJET .TYPMAIL:
!     -----------------
    if (nbma .gt. 0) then
        call wkvect(mag//'.TYPMAIL', 'G V I', nbma, ibid)
        do i = 1, nbm1
            call jeveuo(ma1//'.TYPMAIL', 'L', iatyma)
            iatyp1 = iatyma-1+i
            call jeveuo(mag//'.TYPMAIL', 'E', iatyma)
            iatypx = iatyma-1+i
            zi(iatypx) = zi(iatyp1)
        end do
        do i = 1, nbm2
            i1 = i+nbm1
            call jeveuo(ma2//'.TYPMAIL', 'L', iatyma)
            iatyp2 = iatyma-1+i
            call jeveuo(mag//'.TYPMAIL', 'E', iatyma)
            iatypx = iatyma-1+i1
            zi(iatypx) = zi(iatyp2)
        end do
    end if
!
!
!     --OBJET .GROUPEMA:
!     -----------------
    call jeexin(ma1//'.GROUPEMA', iret1)
    call jeexin(ma2//'.GROUPEMA', iret2)
    nbgm1 = 0
    nbgm2 = 0
    if (iret1 .gt. 0) call jelira(ma1//'.GROUPEMA', 'NUTIOC', nbgm1)
    if (iret2 .gt. 0) call jelira(ma2//'.GROUPEMA', 'NUTIOC', nbgm2)
    nbgma = nbgm1+nbgm2
    if (nbgma .gt. 0) then
        call jecreo(mag//'.PTRNOMMAI', 'G N K24')
        call jeecra(mag//'.PTRNOMMAI', 'NOMMAX', nbgma)
        call jecrec(mag//'.GROUPEMA', 'G V I', 'NO '//mag//'.PTRNOMMAI', 'DISPERSE', 'VARIABLE', &
                    nbgma)
        do i = 1, nbgm1
            call jeveuo(jexnum(ma1//'.GROUPEMA', i), 'L', iagma1)
            call jelira(jexnum(ma1//'.GROUPEMA', i), 'LONUTI', n)
            call jenuno(jexnum(ma1//'.GROUPEMA', i), nogma)
            call jecroc(jexnom(mag//'.GROUPEMA', nogma))
            call jeecra(jexnum(mag//'.GROUPEMA', i), 'LONMAX', max(1, n))
            call jeecra(jexnum(mag//'.GROUPEMA', i), 'LONUTI', n)
            call jeveuo(jexnum(mag//'.GROUPEMA', i), 'E', iagmax)
            do ii = 1, n
                zi(iagmax-1+ii) = zi(iagma1-1+ii)
            end do
        end do
        icompt = 0
        do i = 1, nbgm2
            call jeveuo(jexnum(ma2//'.GROUPEMA', i), 'L', iagma2)
            call jelira(jexnum(ma2//'.GROUPEMA', i), 'LONUTI', n)
            call jenuno(jexnum(ma2//'.GROUPEMA', i), nogma)
            call jeexin(jexnom(mag//'.GROUPEMA', nogma), iret)
            if (iret .gt. 0) then
                call utmess('A', 'MODELISA2_21', sk=nogma)
                nogmab = nogma
                ii = lxlgut(nogmab(1:7))
                do k = ii+1, 7
                    nogmab(k:k) = '_'
                end do
                do k = 0, 9
                    call codent(k, 'G', kkk)
                    nogmab(8:8) = kkk
                    call jeexin(jexnom(mag//'.GROUPEMA', nogmab), iret)
                    if (iret .eq. 0) goto 723
                end do
723             continue
                write (ifm, *) ' LE GROUP_MA '//nogma//' DU MAILLAGE ' &
                    //ma2//' EST RENOMME '//nogmab//' DANS '//mag
                nogma = nogmab
            end if
            icompt = icompt+1
            i1 = nbgm1+icompt
            call jecroc(jexnom(mag//'.GROUPEMA', nogma))
            call jeecra(jexnum(mag//'.GROUPEMA', i1), 'LONMAX', max(1, n))
            call jeecra(jexnum(mag//'.GROUPEMA', i1), 'LONUTI', n)
            call jeveuo(jexnum(mag//'.GROUPEMA', i1), 'E', iagmax)
            do ii = 1, n
                zi(iagmax-1+ii) = zi(iagma2-1+ii)+nbm1
            end do
        end do
    end if
!
!
!     --OBJET .GROUPENO:
!     -----------------
    call jeexin(ma1//'.GROUPENO', iret1)
    call jeexin(ma2//'.GROUPENO', iret2)
    nbgn1 = 0
    nbgn2 = 0
    if (iret1 .gt. 0) call jelira(ma1//'.GROUPENO', 'NUTIOC', nbgn1)
    if (iret2 .gt. 0) call jelira(ma2//'.GROUPENO', 'NUTIOC', nbgn2)
    nbgno = nbgn1+nbgn2
    if (nbgno .gt. 0) then
        call jecreo(mag//'.PTRNOMNOE', 'G N K24')
        call jeecra(mag//'.PTRNOMNOE', 'NOMMAX', nbgno)
        call jecrec(mag//'.GROUPENO', 'G V I', 'NO '//mag//'.PTRNOMNOE', 'DISPERSE', 'VARIABLE', &
                    nbgno)
        do i = 1, nbgn1
            call jeveuo(jexnum(ma1//'.GROUPENO', i), 'L', iagno1)
            call jelira(jexnum(ma1//'.GROUPENO', i), 'LONUTI', n)
            call jenuno(jexnum(ma1//'.GROUPENO', i), nogma)
            call jecroc(jexnom(mag//'.GROUPENO', nogma))
            call jeecra(jexnum(mag//'.GROUPENO', i), 'LONMAX', max(1, n))
            call jeecra(jexnum(mag//'.GROUPENO', i), 'LONUTI', n)
            call jeveuo(jexnum(mag//'.GROUPENO', i), 'E', iagnox)
            do ii = 1, n
                zi(iagnox-1+ii) = zi(iagno1-1+ii)
            end do
        end do
        icompt = 0
        do i = 1, nbgn2
            call jeveuo(jexnum(ma2//'.GROUPENO', i), 'L', iagno2)
            call jelira(jexnum(ma2//'.GROUPENO', i), 'LONUTI', n)
            call jenuno(jexnum(ma2//'.GROUPENO', i), nogno)
            call jeexin(jexnom(mag//'.GROUPENO', nogno), iret)
            if (iret .gt. 0) then
                call utmess('A', 'MODELISA2_22', sk=nogno)
                nognob = nogno
                ii = lxlgut(nognob(1:7))
                do k = ii+1, 7
                    nognob(k:k) = '_'
                end do
                do k = 0, 9
                    call codent(k, 'G', kkk)
                    nognob(8:8) = kkk
                    call jeexin(jexnom(mag//'.GROUPENO', nognob), iret)
                    if (iret .eq. 0) goto 823
                end do
823             continue
                write (ifm, *) ' LE GROUP_NO '//nogno//' DU MAILLAGE ' &
                    //ma2//' EST RENOMME '//nognob//' DANS '//mag
                nogno = nognob
            end if
            icompt = icompt+1
            i1 = nbgn1+icompt
            call jecroc(jexnom(mag//'.GROUPENO', nogno))
            call jeecra(jexnum(mag//'.GROUPENO', i1), 'LONMAX', max(1, n))
            call jeecra(jexnum(mag//'.GROUPENO', i1), 'LONUTI', n)
            call jeveuo(jexnum(mag//'.GROUPENO', i1), 'E', iagnox)
            do ii = 1, n
                zi(iagnox-1+ii) = zi(iagno2-1+ii)+nbn1
            end do
        end do
    end if
!
    call jedema()
end subroutine
