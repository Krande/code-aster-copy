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
subroutine asmael(ma1, ma2, mag)
    implicit none
#include "jeveux.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/ssdmte.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: ma1, ma2, mag
!
!     OPERATEUR: ASSE_MAILLAGE / CAS DE L ASSEMBLAGE DE MAILLAGES
!     POUR LA SOUS-STRUCTURATION
!
!-----------------------------------------------------------------------
!
    character(len=8) ::  nono, nosma, nomacr
    character(len=24) :: valk(3), nogma, nogno
    real(kind=8) :: x, y, z, drefe, dij
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, i1noe, i1nol, iacon1, iacon2, iaconx
    integer(kind=8) :: iacoo2, iacoor, iadime, iadimp
    integer(kind=8) :: iagma1, iagma2, iagmax, iagno1, iagno2, iagnox, iancnf
    integer(kind=8) :: ianmcr, ianon2, iaparr
    integer(kind=8) :: iasup1, iasup2, iasupm, iatyma, iatyp1, iatyp2, iatypx
    integer(kind=8) :: ibid, icompt, ii, iret, iret1, iret2
    integer(kind=8) :: itrou, j, l1, l2, l3, n, nbgm1
    integer(kind=8) :: nbgm2, nbgma, nbgn1, nbgn2, nbgno, nbl1, nbm1
    integer(kind=8) :: nbm2, nbma, nbn1, nbn2, nbno, nbnoe, nbnol
    integer(kind=8) :: nbsm1, nbsm2, nbsma, ncoor, ier
    integer(kind=8), pointer :: dim1(:) => null()
    integer(kind=8), pointer :: dim2(:) => null()
    real(kind=8), pointer :: par1(:) => null()
    real(kind=8), pointer :: par2(:) => null()
    integer(kind=8), pointer :: desm(:) => null()
    character(len=8), pointer :: nmc1(:) => null()
    character(len=8), pointer :: nmc2(:) => null()
    aster_logical :: lcolle, lcolle2
!-----------------------------------------------------------------------
    call jemarq()
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
    if (dim1(6) .ne. dim2(6)) then
        call utmess('A', 'SOUSTRUC_1')
    end if
!
    ncoor = dim1(6)
    zi(iadime-1+6) = ncoor
!
    nbsma = zi(iadime-1+4)
    nbsm1 = dim1(4)
    nbsm2 = dim2(4)
!
    nbma = zi(iadime-1+3)
    nbm1 = dim1(3)
    nbm2 = dim2(3)
!
    nbno = zi(iadime-1+1)
    nbn1 = dim1(1)
    nbn2 = dim2(1)
!
    nbl1 = dim1(2)
!
!
!     --OBJET .NOMACR :
!     -----------------
    if (nbsma .gt. 0) then
        call wkvect(mag//'.NOMACR', 'G V K8', nbsma, ianmcr)
        if (nbsm1 .gt. 0) call jeveuo(ma1//'.NOMACR', 'L', vk8=nmc1)
        if (nbsm2 .gt. 0) call jeveuo(ma2//'.NOMACR', 'L', vk8=nmc2)
        do i = 1, nbsm1
            zk8(ianmcr-1+i) = nmc1(i)
        end do
        do i = 1, nbsm2
            zk8(ianmcr-1+nbsm1+i) = nmc2(i)
        end do
    end if
!
!
!     --OBJET .DIME_2 (V):
!     -----------------
    if (nbsma .gt. 0) then
        call wkvect(mag//'.DIME_2', 'V V I', 4*nbsma, iadimp)
        i1noe = 0
        i1nol = 0
        do i = 1, nbsma
            nomacr = zk8(ianmcr-1+i)
            call jeveuo(nomacr//'.DESM', 'L', vi=desm)
            nbnoe = desm(2)
            nbnol = desm(8)+desm(9)
            zi(iadimp-1+4*(i-1)+1) = nbnoe
            zi(iadimp-1+4*(i-1)+2) = nbnol
            zi(iadimp-1+4*(i-1)+3) = i1noe
            zi(iadimp-1+4*(i-1)+4) = i1nol
            i1noe = i1noe+nbnoe
            i1nol = i1nol+nbnol
        end do
    end if
!
!
!     --OBJET .PARA_R :
!     -----------------
    if (nbsma .gt. 0) then
        call wkvect(mag//'.PARA_R', 'G V R', 14*nbsma, iaparr)
        if (nbsm1 .gt. 0) call jeveuo(ma1//'.PARA_R', 'L', vr=par1)
        if (nbsm2 .gt. 0) call jeveuo(ma2//'.PARA_R', 'L', vr=par2)
        do i = 1, 14*nbsm1
            zr(iaparr-1+i) = par1(i)
        end do
        do i = 1, 14*nbsm2
            zr(iaparr-1+nbsm1+i) = par2(i)
        end do
    end if
!
!
!     --OBJET .SUPMAIL:
!     -----------------
    if (nbsma .gt. 0) then
        call jecrec(mag//'.SUPMAIL', 'G V I', 'NO', 'DISPERSE', 'VARIABLE', &
                    nbsma)
        do i = 1, nbsm1
            call jeveuo(jexnum(ma1//'.SUPMAIL', i), 'L', iasup1)
            call jelira(jexnum(ma1//'.SUPMAIL', i), 'LONMAX', n)
            call jenuno(jexnum(ma1//'.SUPMAIL', i), nosma)
            call jecroc(jexnom(mag//'.SUPMAIL', nosma))
            call jeecra(jexnum(mag//'.SUPMAIL', i), 'LONMAX', n)
            call jeveuo(jexnum(mag//'.SUPMAIL', i), 'E', iasupm)
            do ii = 1, n
                if (zi(iasup1-1+ii) .le. nbn1) then
                    zi(iasupm-1+ii) = zi(iasup1-1+ii)
                else
                    zi(iasupm-1+ii) = zi(iasup1-1+ii)+nbn2
                end if
            end do
        end do
        do i = 1, nbsm2
            i1 = i+nbsm1
            call jeveuo(jexnum(ma2//'.SUPMAIL', i), 'L', iasup2)
            call jelira(jexnum(ma2//'.SUPMAIL', i), 'LONMAX', n)
            call jenuno(jexnum(ma2//'.SUPMAIL', i), nosma)
            call jeexin(jexnom(mag//'.SUPMAIL', nosma), iret)
            if (iret .gt. 0) then
                call utmess('F', 'SOUSTRUC_2', sk=nosma)
            end if
            call jecroc(jexnom(mag//'.SUPMAIL', nosma))
            call jeecra(jexnum(mag//'.SUPMAIL', i1), 'LONMAX', n)
            call jeveuo(jexnum(mag//'.SUPMAIL', i1), 'E', iasupm)
            do ii = 1, n
                if (zi(iasup2-1+ii) .le. nbn2) then
                    zi(iasupm-1+ii) = zi(iasup2-1+ii)+nbn1
                else
                    zi(iasupm-1+ii) = zi(iasup2-1+ii)+nbn1+nbl1
                end if
            end do
        end do
    end if
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
                call utmess('A', 'SOUSTRUC_4', sk=nogma)
                goto 32
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
32          continue
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
                call utmess('A', 'SOUSTRUC_5', sk=nogno)
                goto 34
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
34          continue
        end do
    end if
!
!
!     --OBJET .COORDO_2 (V):
!     ----------------------
    call wkvect(mag//'.COORDO_2', 'V V R', 3*nbno, iacoo2)
    call jeveuo(ma1//'.COORDO    .VALE', 'L', iacoor)
    do i = 1, 3*nbn1
        zr(iacoo2-1+i) = zr(iacoor-1+i)
    end do
    call jeveuo(ma2//'.COORDO    .VALE', 'L', iacoor)
    do i = 1, 3*nbn2
        zr(iacoo2-1+3*nbn1+i) = zr(iacoor-1+i)
    end do
!
!
!     --OBJET .NOEUD_CONF ET .NOMNOE_2 (V):
!     -------------------------------------
    call wkvect(mag//'.NOEUD_CONF', 'V V I', nbno, iancnf)
    call wkvect(mag//'.NOMNOE_2', 'V V K8', nbno, ianon2)
    do i = 1, nbno
        zi(iancnf-1+i) = i
    end do
    call jeexin(ma1//'.NOMNOE', iret)
    if (iret .eq. 0) then
        call jecreo(ma1//'.NOMNOE', 'V N K8')
        call jeecra(ma1//'.NOMNOE', 'NOMMAX', nbn1)
        do i = 1, nbn1
            nono = int_to_char8(i)
            nono = 'N'//nono(1:7)
            call jecroc(jexnom(ma1//'.NOMNOE', nono))
        end do
    end if
    lcolle = .true.
    do i = 1, nbn1
        nono = int_to_char8(i, lcolle, ma1, 'NOEUD')
        zk8(ianon2-1+i) = nono
    end do
    call jeexin(ma2//'.NOMNOE', iret)
    if (iret .eq. 0) then
        call jecreo(ma2//'.NOMNOE', 'V N K8')
        call jeecra(ma2//'.NOMNOE', 'NOMMAX', nbn2)
        do i = 1, nbn2
            nono = int_to_char8(i)
            nono = 'N'//nono(1:7)
            call jecroc(jexnom(ma2//'.NOMNOE', nono))
        end do
    end if
    lcolle = .false.
    call jeexin(ma1//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
    lcolle2 = .false.
    call jeexin(ma2//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle2 = .true.
    end if
    do i = 1, nbn2
        nono = int_to_char8(i, lcolle2, ma2, 'NOEUD')
        zk8(ianon2-1+nbn1+i) = nono
        itrou = char8_to_int(nono, lcolle, ma1, 'NOEUD')
        if (itrou .gt. 0) then
            zi(iancnf-1+nbn1+i) = itrou
        end if
    end do
!
!
!     --ON VERIFIE QUE LES NOEUDS CONFONDUS NE SONT PAS TROP DISTANTS:
!     ----------------------------------------------------------------
    drefe = 0.0d0
    do i = 1, nbno
        x = zr(iacoo2-1+3*(i-1)+1)-zr(iacoo2-1+1)
        y = zr(iacoo2-1+3*(i-1)+2)-zr(iacoo2-1+2)
        z = zr(iacoo2-1+3*(i-1)+3)-zr(iacoo2-1+3)
        drefe = max(drefe, sqrt(x**2+y**2+z**2))
    end do
    do i = 1, nbno
        j = zi(iancnf-1+i)
        if (j .ne. i) then
            x = zr(iacoo2-1+3*(i-1)+1)-zr(iacoo2-1+3*(j-1)+1)
            y = zr(iacoo2-1+3*(i-1)+2)-zr(iacoo2-1+3*(j-1)+2)
            z = zr(iacoo2-1+3*(i-1)+3)-zr(iacoo2-1+3*(j-1)+3)
            dij = sqrt(x**2+y**2+z**2)
            if (dij .gt. 1.0d-6*drefe) then
                valk(1) = zk8(ianon2-1+i)
                valk(2) = ma1
                valk(3) = ma2
                call utmess('A', 'SOUSTRUC_6', nk=3, valk=valk)
            end if
        end if
    end do
!
!
!     --ON "TERMINE" LE MAILLAGE:
!     ---------------------------
    call ssdmte(mag)
!
    call jedema()
end subroutine
