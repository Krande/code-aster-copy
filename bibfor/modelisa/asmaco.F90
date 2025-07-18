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
subroutine asmaco(ma1, ma2, mag)
! aslint: disable=W1501
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/getvtx.h"
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
#include "asterfort/ltnotb.h"
#include "asterfort/lxlgut.h"
#include "asterfort/tbliva.h"
#include "asterfort/tri.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: ma1, ma2, mag
!     OPERATEUR: ASSE_MAILLAGE / CAS DE L ASSEMBLAGE DE MAILLAGES
!     AVEC COLLAGE DE DEUX GROUPES DE MAILLES
!
!-----------------------------------------------------------------------
!
    character(len=1) :: kkk
    character(len=8) :: kind, k8b
    character(len=8) :: noma, nono
    character(len=19) :: coordo, nomt19
    character(len=24) :: para, valk(2), cgpm1, cgpm2, nogma, nogmab, nogno
    character(len=24) :: nognob
    integer(kind=8) :: nbma, nbm1, nbm2, nbno, nbn1, nbn2, nbgma, nbgm1, nbgm2
    integer(kind=8) :: nbngm1, nbngm2, nbngm, nno1, nno2, ianode, nnodif
    integer(kind=8) :: i1, icompt, ino, l1, l2, l3, i, n, ncoor, k, ifm, niv, j
    integer(kind=8) :: iadime
    integer(kind=8) :: iagma1, iagma2, iagmax
    integer(kind=8) :: iacon1, iacon2, iaconx
    integer(kind=8) :: iagno1, iagno2, iagnox
    integer(kind=8) :: iatyp1, iatyp2, iatypx
    integer(kind=8) :: nbgno, nbgn1, nbgn2, ii, jj, igeomr, iadesc, ibid
    integer(kind=8) :: iatyma, iacoo1, iacoo2, iavale, iret, iret1, iret2
    integer(kind=8) :: iamam1, iamam2, nbpar
    integer(kind=8) :: ilgma, ilgm2, decal
    aster_logical :: match, elim
    real(kind=8) :: prec1, prec2, prec, dist, x1, y1, z1, x2, y2, z2, r8b, armin
    real(kind=8) :: xi1, yi1, zi1
    complex(kind=8) :: c16b
    integer(kind=8), pointer :: dim1(:) => null()
    integer(kind=8), pointer :: dim2(:) => null()
!
!     ------------------------------------------------------------------
!
    call jemarq()
    r8b = 0.d0
    call infniv(ifm, niv)
!CC   ------------------------------------------------------------------
!CC RECUPERATION DE L'ARETE MINIMUM DES MAILLAGE
!CC   ------------------------------------------------------------------
    call jeexin(ma1//'           .LTNT', iret)
    if (iret .ne. 0) then
        call ltnotb(ma1, 'CARA_GEOM', nomt19)
        nbpar = 0
        para = 'AR_MIN                  '
        call tbliva(nomt19, nbpar, ' ', [ibid], [r8b], &
                    [c16b], k8b, k8b, [r8b], para, &
                    k8b, ibid, armin, c16b, k8b, &
                    iret)
        if (iret .ne. 0) then
            call utmess('F', 'MODELISA2_13')
        end if
        prec1 = armin*1.d-06
    else
        prec1 = 1.d-10
    end if
    if (prec1 .le. 0.d0) then
        call utmess('F', 'MODELISA2_14')
    end if
!CC   ------------------------------------------------------------------
    call jeexin(ma2//'           .LTNT', iret)
    if (iret .ne. 0) then
        call ltnotb(ma2, 'CARA_GEOM', nomt19)
        nbpar = 0
        para = 'AR_MIN                  '
        call tbliva(nomt19, nbpar, ' ', [ibid], [r8b], &
                    [c16b], k8b, k8b, [r8b], para, &
                    k8b, ibid, armin, c16b, k8b, &
                    iret)
        if (iret .ne. 0) then
            call utmess('F', 'MODELISA2_13')
        end if
        prec2 = armin*1.d-06
    else
        prec2 = 1.d-10
    end if
    if (prec2 .le. 0.d0) then
        call utmess('F', 'MODELISA2_14')
    end if
    prec = min(prec1, prec2)
!CC   ------------------------------------------------------------------
!CC RECUPERATION DES 2 GROUP_MA A COLLER
!CC   ------------------------------------------------------------------
    call getvtx('COLLAGE', 'GROUP_MA_1', iocc=1, scal=cgpm1, nbret=ibid)
    call getvtx('COLLAGE', 'GROUP_MA_2', iocc=1, scal=cgpm2, nbret=ibid)
    elim = .false.
    call jeexin(jexnom(ma1//'.GROUPEMA', cgpm1), iret1)
    if (iret1 .eq. 0) then
        valk(1) = cgpm1
        valk(2) = ma1
        call utmess('F', 'MODELISA2_15', nk=2, valk=valk)
    end if
    call jeexin(jexnom(ma2//'.GROUPEMA', cgpm2), iret2)
    if (iret2 .eq. 0) then
        valk(1) = cgpm2
        valk(2) = ma2
        call utmess('F', 'MODELISA2_16', nk=2, valk=valk)
    end if
!CC   ------------------------------------------------------------------
!CC VERIFICATION QUE LES 2 GROUP_MA A COLLER ONT LE MM NOMBRE DE MAILLES
!CC   ------------------------------------------------------------------
    call jelira(jexnom(ma1//'.GROUPEMA', cgpm1), 'LONUTI', nbngm1)
    call jelira(jexnom(ma2//'.GROUPEMA', cgpm2), 'LONUTI', nbngm2)
    nbngm = nbngm1
    if (nbngm1 .ne. nbngm2) then
        valk(1) = cgpm1
        valk(2) = cgpm2
        call utmess('F', 'MODELISA2_17', nk=2, valk=valk)
    end if
    call jeveuo(jexnom(ma1//'.GROUPEMA', cgpm1), 'L', iagma1)
    call jeveuo(jexnom(ma2//'.GROUPEMA', cgpm2), 'L', iagma2)
!CC   ------------------------------------------------------------------
!CC CREATION DU VECTEUR DES NOEUDS A APPARIER DANS CHACUN DES MAILLAGES
!CC VECTEUR '.NODE' :  - DE    1 A   NNO1 : NOEUDS DU MAILLAGE 2
!CC                    - DE NNO1 A 2*NNO1 : NOEUDS DU MAILLAGE 1
!CC   ------------------------------------------------------------------
    nno1 = 0
    nno2 = 0
    do i = 1, nbngm
        call jelira(jexnum(ma1//'.CONNEX', zi(iagma1+i-1)), 'LONMAX', ii)
        call jelira(jexnum(ma2//'.CONNEX', zi(iagma2+i-1)), 'LONMAX', jj)
        nno1 = nno1+ii
        nno2 = nno2+jj
    end do
    if (nno1 .ne. nno2) then
        call utmess('F', 'MODELISA2_18')
    end if
    call wkvect('&&ASMACO'//'.NODE', 'V V I', nno1*2, ianode)
    nno1 = 0
    do i = 1, nbngm
        call jelira(jexnum(ma1//'.CONNEX', zi(iagma1+i-1)), 'LONMAX', ii)
        call jeveuo(jexnum(ma1//'.CONNEX', zi(iagma1+i-1)), 'L', iagno1)
        do j = 1, ii
            nno1 = nno1+1
            zi(ianode+nno1-1) = zi(iagno1+j-1)
        end do
    end do
    call tri(zi(ianode), zi, 0, nno1)
    nnodif = 1
    zi(ianode+nno1) = zi(ianode)
    do i = 2, nno1
        if (zi(ianode+i-1) .ne. zi(ianode+nno1+nnodif-1)) then
            nnodif = nnodif+1
            zi(ianode+nno1+nnodif-1) = zi(ianode+i-1)
        end if
    end do
    nno2 = 0
    call jeveuo(ma1//'.COORDO    .VALE', 'L', iacoo1)
    call jeveuo(ma2//'.COORDO    .VALE', 'L', iacoo2)
!
    do k = 1, nnodif
        x1 = zr(iacoo1+3*(zi(ianode+nno1+k-1)-1)-1+1)
        y1 = zr(iacoo1+3*(zi(ianode+nno1+k-1)-1)-1+2)
        z1 = zr(iacoo1+3*(zi(ianode+nno1+k-1)-1)-1+3)
!
!     TEST DE LA PRESENCE DE NOEUDS GEOMETRIQUEMENT CONFONDUS
!     DANS LA ZONE A COLLER
!     SI C EST LE CAS, L APPARIEMENT SE PASSERA MAL ENSUITE
!     PUISQU IL EST FAIT SUR UN CRITERE DE DISTANCE NULLE
!
!     NORMALEMENT FAIRE LE TEST SUR LA ZONE A COLLER DU MAILLAGE 1
!     SUFFIT, PUISQU ON VERIFIE ENSUITE QUE LA ZONE DU MAILLAGE 2
!     LUI EST TOPOLOGIQUEMENT ET GEOMETRIQUMENT IDENTIQUE
!
        do i = 1, k-1
            xi1 = zr(iacoo1+3*(zi(ianode+nno1+i-1)-1)-1+1)
            yi1 = zr(iacoo1+3*(zi(ianode+nno1+i-1)-1)-1+2)
            zi1 = zr(iacoo1+3*(zi(ianode+nno1+i-1)-1)-1+3)
            dist = (x1-xi1)**2+(y1-yi1)**2+(z1-zi1)**2
            if (dist .le. prec) then
                call utmess('F', 'MODELISA2_97')
            end if
        end do
!
        match = .false.
        do i = 1, nbngm
            call jelira(jexnum(ma2//'.CONNEX', zi(iagma2+i-1)), 'LONMAX', ii)
            call jeveuo(jexnum(ma2//'.CONNEX', zi(iagma2+i-1)), 'L', iagno2)
            do j = 1, ii
                nno2 = nno2+1
                x2 = zr(iacoo2+3*(zi(iagno2+j-1)-1)-1+1)
                y2 = zr(iacoo2+3*(zi(iagno2+j-1)-1)-1+2)
                z2 = zr(iacoo2+3*(zi(iagno2+j-1)-1)-1+3)
                dist = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
                if (dist .le. prec) then
                    match = .true.
                    goto 1060
                end if
            end do
        end do
        if (.not. match) then
            call utmess('F', 'MODELISA2_19')
        end if
1060    continue
        zi(ianode+k-1) = zi(iagno2+j-1)
    end do
!CC   ------------------------------------------------------------------
!CC   --OBJET .DIME :
!CC   ------------------------------------------------------------------
    call jeveuo(ma1//'.DIME', 'L', vi=dim1)
    call jeveuo(ma2//'.DIME', 'L', vi=dim2)
    call wkvect(mag//'.DIME', 'G V I', 6, iadime)
!CC SOMME POUR : 1 LE NB DE NOEUDS - ON LAISSE LES NOEUDS EN DOUBLE
!CC                                  INUTILISES DANS LA CONNECTIVITE
!CC              2       DE NOEUDS LAGRANGES,
!CC              3       DE MAILLES - ON SOUSTRAIT LES MAILLES APPARIES
!CC                                   POUR LES 2 MAILLAGES
!CC              4       DE SUPER MAILLES
!CC              5       DU MAJORANT DE SUPER MAILLES
    zi(iadime-1+1) = dim1(1)+dim2(1)
    zi(iadime-1+2) = dim1(2)+dim2(2)
    if (elim) then
!CC   SI ELIM : ON SUPPRIME LES MAILLES DES 2 GROUPES CGPM1 ET CGPM2
        zi(iadime-1+3) = dim1(3)+dim2(3)-2*nbngm
    else
!CC   SINON, SEULES LES MAILLES DE CGPM1 SONT SUPPRIMEES
        zi(iadime-1+3) = dim1(3)+dim2(3)-nbngm
    end if
!
    zi(iadime-1+4) = dim1(4)+dim2(4)
    zi(iadime-1+5) = dim1(5)+dim2(5)
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
!CC   ------------------------------------------------------------------
!CC   --OBJET .CONNEX :
!CC   ON NE RETIENT QUE LES MAILLES HORS DES 2 GROUP_MA CGPM1 ET CGPM2
!CC   POUR LES MAILLES DU MAILLAGE 2 CONTENANT DES NOEUDS DE CGPM2,
!CC   ON SUBSTITUE DANS LEUR CONNECTIVITE LES NOMS DE NOEUDS DU
!CC   MAILLAGE 1 QUI LEUR ONT ETE APPARIES
!CC   ------------------------------------------------------------------
    if (nbma .gt. 0) then
        call jecrec(mag//'.CONNEX', 'G V I', 'NU', 'CONTIG', 'VARIABLE', &
                    nbma)
        call jeecra(mag//'.CONNEX', 'NUTIOC', nbma)

        call wkvect('&&ASMACO'//'.MAM1', 'V V I', nbm1*2, iamam1)
        call wkvect('&&ASMACO'//'.MAM2', 'V V I', nbm2*2, iamam2)
        do i = 1, nbm1
            zi(iamam1+i-1) = i
        end do
        do i = 1, nbm2
            zi(iamam2+i-1) = i
        end do
        do i = 1, nbngm
            zi(iamam1+zi(iagma1+i-1)-1) = 0
            if (elim) zi(iamam2+zi(iagma2+i-1)-1) = 0
        end do
        ii = 0
        do i = 1, nbm1
            if (zi(iamam1+i-1) .eq. 0) then
                zi(iamam1+nbm1+i-1) = 0
            else
                ii = ii+1
                zi(iamam1+nbm1+i-1) = ii
            end if
        end do
        ii = 0
        do i = 1, nbm2
            if (zi(iamam2+i-1) .eq. 0) then
                zi(iamam2+nbm2+i-1) = 0
            else
                ii = ii+1
                zi(iamam2+nbm2+i-1) = ii
            end if
        end do
        call tri(zi(iamam1), zi, 0, nbm1)
        call tri(zi(iamam2), zi, 0, nbm2)
        l1 = 0
        l2 = 0
        if (nbm1 .gt. 0) call jelira(ma1//'.CONNEX', 'LONT', l1)
        if (nbm2 .gt. 0) call jelira(ma2//'.CONNEX', 'LONT', l2)
        l3 = l1+l2
        call jeecra(mag//'.CONNEX', 'LONT', l3)
        do i = 1, nbm1-nbngm
            call jeveuo(jexnum(ma1//'.CONNEX', zi(iamam1+nbngm+i-1)), 'L', iacon1)
            call jelira(jexnum(ma1//'.CONNEX', zi(iamam1+nbngm+i-1)), 'LONMAX', n)
            call jeecra(jexnum(mag//'.CONNEX', i), 'LONMAX', n)
            call jeveuo(jexnum(mag//'.CONNEX', i), 'E', iaconx)
            do ii = 1, n
                zi(iaconx-1+ii) = zi(iacon1-1+ii)
            end do
        end do
        if (elim) then
            decal = nbngm
        else
            decal = 0
        end if
        do i = 1, nbm2-decal
            i1 = i+nbm1-nbngm
            call jeveuo(jexnum(ma2//'.CONNEX', zi(iamam2+decal+i-1)), 'L', iacon2)
            call jelira(jexnum(ma2//'.CONNEX', zi(iamam2+decal+i-1)), 'LONMAX', n)
            call jeecra(jexnum(mag//'.CONNEX', i1), 'LONMAX', n)
            call jeveuo(jexnum(mag//'.CONNEX', i1), 'E', iaconx)
            do ii = 1, n
                match = .false.
                do jj = 1, nnodif
                    if (zi(iacon2+ii-1) .eq. zi(ianode+jj-1)) then
                        match = .true.
                        goto 423
                    end if
                end do
423             continue
                if (match) then
                    zi(iaconx+ii-1) = zi(ianode+nno1+jj-1)
                else
                    zi(iaconx+ii-1) = zi(iacon2+ii-1)+nbn1
                end if
            end do
        end do
    end if
!CC   ------------------------------------------------------------------
!CC   --OBJET .COORDO :
!CC   ------------------------------------------------------------------
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
    call jeveuo(ma1//'.COORDO    .VALE', 'L', iacoo1)
    call jeveuo(ma2//'.COORDO    .VALE', 'L', iacoo2)
    call wkvect(coordo//'.VALE', 'G V R', 3*nbno, iavale)
!     -- COORDONNEES DES NOEUDS :
    do 51, ino = 1, nbn1
    do k = 1, 3
        zr(iavale-1+3*(ino-1)+k) = zr(iacoo1-1+3*(ino-1)+k)
    end do
51  end do
    do 52, ino = 1, nbn2
    do k = 1, 3
        zr(iavale-1+3*(nbn1+ino-1)+k) = zr(iacoo2-1+3*(ino-1)+k)
    end do
52  end do
!CC   ------------------------------------------------------------------
!CC   --OBJET .TYPMAIL :
!CC   ------------------------------------------------------------------
    if (nbma .gt. 0) then
        call wkvect(mag//'.TYPMAIL', 'G V I', nbma, ibid)
        do i = 1, nbm1-nbngm
            call jeveuo(ma1//'.TYPMAIL', 'L', iatyma)
            iatyp1 = iatyma-1+zi(iamam1+nbngm+i-1)
            call jeveuo(mag//'.TYPMAIL', 'E', iatyma)
            iatypx = iatyma-1+i
            zi(iatypx) = zi(iatyp1)
        end do
        if (elim) then
            decal = nbngm
        else
            decal = 0
        end if
        do i = 1, nbm2-decal
            i1 = i+nbm1-nbngm
            call jeveuo(ma2//'.TYPMAIL', 'L', iatyma)
            iatyp2 = iatyma-1+zi(iamam2+decal+i-1)
            call jeveuo(mag//'.TYPMAIL', 'E', iatyma)
            iatypx = iatyma-1+i1
            zi(iatypx) = zi(iatyp2)
        end do
    end if
!CC   ------------------------------------------------------------------
!CC   --OBJET .GROUPEMA:
!CC   ON RECREE TOUS LES GROUP_MA DANS LE NOUVEAU MAILLAGE
!CC - SAUF : CELUI QUI SERT A REALISER LE COLLAGE DANS LE MAILLAGE 1
!CC - SAUF : SON EQUIVALENT DANS LE MAILLAGE 2 SI ELIM=.TRUE.
!CC - SAUF : LES GROUPES QUI SE RETROUVENT VIDES (TOUTES MAILLES SUPPR.)
!CC   DANS CE DERNIER CAS, LA COLLECTION EST INUTILEMENT SURDIMENSIONNEE
!CC   A NBGMA (PAS GRAVE) EN TENANT COMPTE DES GMA DE COLLAGE SUPPRIMES
!CC   MAIS PAS DES GMA VIDES NON CREES
!CC - ICOMPT COMPTE LES GROUP_MA EFFECTIVEMENT CREES (PAR JECROC)
!CC   ------------------------------------------------------------------
    call jeexin(ma1//'.GROUPEMA', iret1)
    call jeexin(ma2//'.GROUPEMA', iret2)
    nbgm1 = 0
    nbgm2 = 0
    if (iret1 .gt. 0) call jelira(ma1//'.GROUPEMA', 'NUTIOC', nbgm1)
    if (iret2 .gt. 0) call jelira(ma2//'.GROUPEMA', 'NUTIOC', nbgm2)
    if (elim) then
        nbgma = nbgm1-1+nbgm2-1
    else
        nbgma = nbgm1-1+nbgm2
    end if
!
    if (nbgma .gt. 0) then
        call jecreo(mag//'.PTRNOMMAI', 'G N K24')
        call jeecra(mag//'.PTRNOMMAI', 'NOMMAX', nbgma)
        call jecrec(mag//'.GROUPEMA', 'G V I', 'NO '//mag//'.PTRNOMMAI', 'DISPERSE', 'VARIABLE', &
                    nbgma)
        icompt = 0
        do i = 1, nbgm1
            call jeveuo(jexnum(ma1//'.GROUPEMA', i), 'L', iagma1)
            call jelira(jexnum(ma1//'.GROUPEMA', i), 'LONUTI', n)
            call jenuno(jexnum(ma1//'.GROUPEMA', i), nogma)
            if (nogma .ne. cgpm1) then
!
                ilgma = 0
                do ii = 1, n
                    if (zi(iamam1+nbm1+zi(iagma1-1+ii)-1) .ne. 0) then
                        ilgma = ilgma+1
                    end if
                end do
!
                if (ilgma .ne. 0) then
                    icompt = icompt+1
                    call jecroc(jexnom(mag//'.GROUPEMA', nogma))
                    call jeecra(jexnum(mag//'.GROUPEMA', icompt), 'LONMAX', max(1, ilgma))
                    call jeecra(jexnum(mag//'.GROUPEMA', icompt), 'LONUTI', ilgma)
                    call jeveuo(jexnum(mag//'.GROUPEMA', icompt), 'E', iagmax)
                    ilgm2 = 0
                    do ii = 1, n
                        if (zi(iamam1+nbm1+zi(iagma1-1+ii)-1) .ne. 0) then
                            ilgm2 = ilgm2+1
                            zi(iagmax-1+ilgm2) = zi(iamam1+nbm1+zi( &
                                                    iagma1-1+ii)-1)
                        end if
                    end do
                else
                    valk(1) = nogma
                    valk(2) = ma1
                    call utmess('A', 'MODELISA7_97', nk=2, valk=valk)
                end if
            end if
        end do
        do i = 1, nbgm2
            call jeveuo(jexnum(ma2//'.GROUPEMA', i), 'L', iagma2)
            call jelira(jexnum(ma2//'.GROUPEMA', i), 'LONUTI', n)
            call jenuno(jexnum(ma2//'.GROUPEMA', i), nogma)
            if ((nogma .ne. cgpm2) .or. (.not. elim)) then
                call jeexin(jexnom(mag//'.GROUPEMA', nogma), iret)
!
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
723                 continue
                    write (ifm, *) ' LE GROUP_MA '//nogma//' DU MAILLAGE '&
     &             //ma2//' EST RENOMME '//nogmab//' DANS '//mag
                    nogma = nogmab
                end if
!
                ilgma = 0
                do ii = 1, n
                    if (zi(iamam2+nbm2+zi(iagma2-1+ii)-1) .ne. 0) then
                        ilgma = ilgma+1
                    end if
                end do
!
                if (ilgma .ne. 0) then
                    icompt = icompt+1
                    call jecroc(jexnom(mag//'.GROUPEMA', nogma))
                    call jeecra(jexnum(mag//'.GROUPEMA', icompt), 'LONMAX', max(1, ilgma))
                    call jeecra(jexnum(mag//'.GROUPEMA', icompt), 'LONUTI', ilgma)
                    call jeveuo(jexnum(mag//'.GROUPEMA', icompt), 'E', iagmax)
                    ilgm2 = 0
                    do ii = 1, n
                        if (zi(iamam2+nbm2+zi(iagma2-1+ii)-1) .ne. 0) then
                            ilgm2 = ilgm2+1
                            zi(iagmax-1+ilgm2) = zi(iamam2+nbm2+zi( &
                                                    iagma2-1+ii)-1)+nbm1-nbngm
                        end if
                    end do
                else
                    valk(1) = nogma
                    valk(2) = ma2
                    call utmess('A', 'MODELISA7_97', nk=2, valk=valk)
                end if
            end if
        end do
    end if
!CC   ------------------------------------------------------------------
!CC   --OBJET .GROUPENO:
!CC   LES GROUP_NO SONT CONSERVES TELS QUELS, DANS LA MESURE OU ON NE
!CC   SUPPRIME PAS DE NOEUDS. POUR LES GROUP_NO DU MAILLAGE 2, SI DES
!CC   NOEUDS FONT PARTI DE CEUX APPARIES, ON LES SUBSTITUE PAR LEUR
!CC   HOMOLOGUE DU MAILLAGE 1, PLUTOT QUE DE LAISSER LE(S) NOEUD
!CC   DESORMAIS ORPHELIN DANS LE GROUPE.
!CC   ------------------------------------------------------------------
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
            call jenuno(jexnum(ma1//'.GROUPENO', i), nogno)
            call jecroc(jexnom(mag//'.GROUPENO', nogno))
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
                match = .false.
                do jj = 1, nnodif
                    if (zi(iagno2+ii-1) .eq. zi(ianode+jj-1)) then
                        match = .true.
                        goto 826
                    end if
                end do
826             continue
                if (match) then
                    zi(iagnox+ii-1) = zi(ianode+nno1+jj-1)
                else
                    zi(iagnox+ii-1) = zi(iagno2+ii-1)+nbn1
                end if
            end do
        end do
    end if
!
    call jedema()
end subroutine
