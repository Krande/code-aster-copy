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
subroutine nuadrf(nuag1, nuag2, ic1, ic2, dref)
    implicit none
#include "jeveux.h"
#include "asterc/r8gaem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    character(len=19) :: nuag1, nuag2
    integer(kind=8) :: ic1, ic2
    real(kind=8) :: dref(*)
!
!  BUT : CALCULER POUR TOUS LES POINTS DE NUAG2 UNE DISTANCE
!        DE REFERENCE POUR QU L'INTERPOLATION NE "CAPTE"
!        QU'UN NOMBRE LIMITE DE POINTS NE NUAG1 :
!        ON CHERCHE UNE INTERPOLATION LA PLUS LOCALE POSSIBLE
!
!        EN UN POINT DONNE IP2 DE NUAG2, ON ASSOCIE LA DISTANCE DREF
!        DREF EST TELLE QUE :
!                LA BOULE (IP2,SQRT(DREF)) CONTIENNE :
!                 . 2 POINTS EN 1D (NON CONFONDUS)
!                 . 3 POINTS EN 2D (NON ALIGNES)
!                 . 4 POINTS EN 3D (NON COPLANAIRES)
!
! IN/JXIN  NUAG1   : NUAGE A PROJETER
! IN/JXIN  NUAG2   : NUAGE A EVALUER
! IN       IC1     : NUMERO DE LA CMP DANS NUAG1
! IN       IC2     : NUMERO DE LA CMP DANS NUAG2
! OU       DREF    : VECTEUR QUI CONTIENDRA LES DISTANCE**2 CHERCHEES
!                    DIMENSION : NP2 = NOMBRE DE POINTS DE NUAG2
! VARIABLES LOCALES :
    integer(kind=8) :: inual1, inual2
    integer(kind=8) :: np1, np2, nx1, nx2, nc1, nc2, ip1, ip2, im1, im2, im3, im4
    real(kind=8) :: x2, y2, z2, x1, y1, z1, xm1, ym1, zm1
    real(kind=8) :: d, dm0, dm, l2, s, s2, v, l
    real(kind=8) :: m1m2(3), m1m3(3), m1p1(3), n2(3), n(3), epsabs
    real(kind=8), pointer :: vdm0(:) => null()
    real(kind=8), pointer :: nuax1(:) => null()
    real(kind=8), pointer :: nuax2(:) => null()
    integer(kind=8), pointer :: nuai1(:) => null()
    integer(kind=8), pointer :: nuai2(:) => null()
!
! DEB-------------------------------------------------------------------
    call jemarq()
!
    epsabs = sqrt(1.d0/r8gaem())
!
    call jeveuo(nuag1//'.NUAI', 'L', vi=nuai1)
    call jeveuo(nuag2//'.NUAI', 'L', vi=nuai2)
    call jeveuo(nuag1//'.NUAX', 'L', vr=nuax1)
    call jeveuo(nuag2//'.NUAX', 'L', vr=nuax2)
    call jeveuo(nuag1//'.NUAL', 'L', inual1)
    call jeveuo(nuag2//'.NUAL', 'L', inual2)
!
    np1 = nuai1(1)
    np2 = nuai2(1)
    nx1 = nuai1(2)
    nx2 = nuai2(2)
    nc1 = nuai1(3)
    nc2 = nuai2(3)
!
!
!     RECHERCHE DE LA PLUS GRANDE DISTANCE**2 ENTRE CHAQUE IP2
!     ET L'ENSEMBLE DES IP1 :
!     ------------------------------------------------------
    AS_ALLOCATE(vr=vdm0, size=np2)
    do 1, ip2 = 1, np2
        if (.not. zl(inual2-1+(ip2-1)*nc2+ic2)) goto 1
!
!       -- DM0 EST LA PLUS GRANDE DISTANCE**2 ENTRE IP2 ET
!          L'ENSEMBLE DES IP1
        dm0 = 0.d0
!
        if (nx1 .eq. 1) then
            x2 = nuax2((ip2-1)*nx2+1)
            do ip1 = 1, np1
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 2
                x1 = nuax1((ip1-1)*nx1+1)
                dm0 = max(dm0, (x2-x1)**2)
2               continue
            end do
        else if (nx1 .eq. 2) then
            x2 = nuax2((ip2-1)*nx2+1)
            y2 = nuax2((ip2-1)*nx2+2)
            do ip1 = 1, np1
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 3
                x1 = nuax1((ip1-1)*nx1+1)
                y1 = nuax1((ip1-1)*nx1+2)
                dm0 = max(dm0, (x2-x1)**2+(y2-y1)**2)
3               continue
            end do
        else if (nx1 .eq. 3) then
            x2 = nuax2((ip2-1)*nx2+1)
            y2 = nuax2((ip2-1)*nx2+2)
            z2 = nuax2((ip2-1)*nx2+3)
            do ip1 = 1, np1
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 4
                x1 = nuax1((ip1-1)*nx1+1)
                y1 = nuax1((ip1-1)*nx1+2)
                z1 = nuax1((ip1-1)*nx1+3)
                dm0 = max(dm0, (x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
4               continue
            end do
        end if
!
        if (dm0 .eq. 0.d0) goto 9994
        vdm0(ip2) = dm0
1   end do
!
!
    if (nx1 .eq. 1) then
!     ------------
        do ip2 = 1, np2
            if (.not. zl(inual2-1+(ip2-1)*nc2+ic2)) goto 10
            x2 = nuax2((ip2-1)*nx2+1)
!
!         -- IM1 EST L'INDICE DU POINT LE + PROCHE DE IP2
            im1 = 0
            dm = vdm0(ip2)
            do ip1 = 1, np1
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 12
                x1 = nuax1((ip1-1)*nx1+1)
                d = (x1-x2)**2
                if (d .le. dm) then
                    dm = d
                    im1 = ip1
                end if
12              continue
            end do
            if (im1 .eq. 0) goto 9995
            xm1 = nuax1((im1-1)*nx1+1)
!
!         -- IM2 EST L'INDICE DU POINT LE + PROCHE DE IP2
!            ET DIFFERENT DE IM1
            im2 = 0
            dm = vdm0(ip2)
            do ip1 = 1, np1
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 13
                x1 = nuax1((ip1-1)*nx1+1)
                if ((x1-xm1)**2 .lt. epsabs) goto 13
                d = (x1-x2)**2
                if (d .le. dm) then
                    dm = d
                    im2 = ip1
                end if
13              continue
            end do
            if (im2 .eq. 0) goto 9996
            dref(ip2) = dm
!
10          continue
        end do
        goto 999
!
    else if (nx1 .eq. 2) then
!     ------------------
        do ip2 = 1, np2
            if (.not. zl(inual2-1+(ip2-1)*nc2+ic2)) goto 20
            x2 = nuax2((ip2-1)*nx2+1)
            y2 = nuax2((ip2-1)*nx2+2)
!
!
!         -- IM1 EST L'INDICE DU POINT LE + PROCHE DE IP2
            im1 = 0
            dm = vdm0(ip2)
            do ip1 = 1, np1
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 22
                x1 = nuax1((ip1-1)*nx1+1)
                y1 = nuax1((ip1-1)*nx1+2)
                d = (x2-x1)**2+(y2-y1)**2
                if (d .le. dm) then
                    dm = d
                    im1 = ip1
                end if
22              continue
            end do
            if (im1 .eq. 0) goto 9995
            xm1 = nuax1((im1-1)*nx1+1)
            ym1 = nuax1((im1-1)*nx1+2)
!
!         -- IM2 EST L'INDICE DU POINT LE + PROCHE DE IP2
!            ET DIFFERENT DE IM1
            im2 = 0
            dm = vdm0(ip2)
            do ip1 = 1, np1
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 23
                x1 = nuax1((ip1-1)*nx1+1)
                y1 = nuax1((ip1-1)*nx1+2)
                if ((x1-xm1)**2+(y1-ym1)**2 .le. epsabs) goto 23
                d = (x2-x1)**2+(y2-y1)**2
                if (d .le. dm) then
                    dm = d
                    im2 = ip1
                end if
23              continue
            end do
            if (im2 .eq. 0) goto 9996
!
!         VECTEUR M1M2 :
            m1m2(1) = nuax1((im1-1)*nx1+1)-nuax1((im2-1)* &
                                                 nx1+1)
            m1m2(2) = nuax1((im1-1)*nx1+2)-nuax1((im2-1)* &
                                                 nx1+2)
            l2 = m1m2(1)**2+m1m2(2)**2
!
!         -- IM3 EST L'INDICE DU POINT M3 LE + PROCHE DE P2
!            DIFFERENT DE M1 ET M2 ET TEL QUE M1 M2 M3 FORMENT UN PLAN
            im3 = 0
            dm = vdm0(ip2)
            do ip1 = 1, np1
!           IF ((IP1.EQ.IM1).OR.(IP1.EQ.IM2)) GOTO 24
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 24
                x1 = nuax1((ip1-1)*nx1+1)
                y1 = nuax1((ip1-1)*nx1+2)
!
!           SI LES POINTS M1 M2 ET P1 NE FORMENT PAS UN PLAN GOTO 24
                m1p1(1) = nuax1((im1-1)*nx1+1)-x1
                m1p1(2) = nuax1((im1-1)*nx1+2)-y1
                s = abs(m1m2(1)*m1p1(2)-m1m2(2)*m1p1(1))
                if (s .le. (1.d-3*l2)) goto 24
                d = (x2-x1)**2+(y2-y1)**2
                if (d .le. dm) then
                    dm = d
                    im3 = ip1
                end if
24              continue
            end do
            if (im3 .eq. 0) goto 9997
            dref(ip2) = dm
20          continue
        end do
        goto 999
!
!
    else if (nx1 .eq. 3) then
!     ------------------
        do ip2 = 1, np2
            if (.not. zl(inual2-1+(ip2-1)*nc2+ic2)) goto 30
            x2 = nuax2((ip2-1)*nx2+1)
            y2 = nuax2((ip2-1)*nx2+2)
            z2 = nuax2((ip2-1)*nx2+3)
!
!
!         -- IM1 EST L'INDICE DU POINT LE + PROCHE DE IP2
            im1 = 0
            dm = vdm0(ip2)
            do ip1 = 1, np1
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 32
                x1 = nuax1((ip1-1)*nx1+1)
                y1 = nuax1((ip1-1)*nx1+2)
                z1 = nuax1((ip1-1)*nx1+3)
                d = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
                if (d .le. dm) then
                    dm = d
                    im1 = ip1
                end if
32              continue
            end do
            if (im1 .eq. 0) goto 9995
            xm1 = nuax1((im1-1)*nx1+1)
            ym1 = nuax1((im1-1)*nx1+2)
            zm1 = nuax1((im1-1)*nx1+3)
!
!         -- IM2 EST L'INDICE DU POINT LE + PROCHE DE IP2
!            ET DIFFERENT DE IM1
            im2 = 0
            dm = vdm0(ip2)
            do ip1 = 1, np1
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 33
                x1 = nuax1((ip1-1)*nx1+1)
                y1 = nuax1((ip1-1)*nx1+2)
                z1 = nuax1((ip1-1)*nx1+3)
                if ((x1-xm1)**2+(y1-ym1)**2+(z1-zm1)**2 .le. epsabs) goto 33
                d = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
                if (d .le. dm) then
                    dm = d
                    im2 = ip1
                end if
33              continue
            end do
            if (im2 .eq. 0) goto 9996
!
!         -- VECTEUR M1M2 :
            m1m2(1) = nuax1((im2-1)*nx1+1)-nuax1((im1-1)* &
                                                 nx1+1)
            m1m2(2) = nuax1((im2-1)*nx1+2)-nuax1((im1-1)* &
                                                 nx1+2)
            m1m2(3) = nuax1((im2-1)*nx1+3)-nuax1((im1-1)* &
                                                 nx1+3)
            l2 = m1m2(1)**2+m1m2(2)**2+m1m2(3)**2
            l = sqrt(l2)
!
!         -- IM3 EST L'INDICE DU POINT M3 LE + PROCHE DE P2
!            DIFFERENT DE M1 ET M2 ET TEL QUE M1 M2 M3 FORMENT UN PLAN
            im3 = 0
            dm = vdm0(ip2)
            do ip1 = 1, np1
!           IF ((IP1.EQ.IM1).OR.(IP1.EQ.IM2)) GOTO 34
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 34
                x1 = nuax1((ip1-1)*nx1+1)
                y1 = nuax1((ip1-1)*nx1+2)
                z1 = nuax1((ip1-1)*nx1+3)
!
!           SI LES POINTS M1 M2 ET P1 NE FORMENT PAS UN PLAN GOTO 34
                m1p1(1) = nuax1((im1-1)*nx1+1)-x1
                m1p1(2) = nuax1((im1-1)*nx1+2)-y1
                m1p1(3) = nuax1((im1-1)*nx1+3)-z1
                n2(1) = m1m2(2)*m1p1(3)-m1m2(3)*m1p1(2)
                n2(2) = m1m2(3)*m1p1(1)-m1m2(1)*m1p1(3)
                n2(3) = m1m2(1)*m1p1(2)-m1m2(2)*m1p1(1)
                s2 = n2(1)**2+n2(2)**2+n2(3)**2
!
                if (s2 .le. (1.d-3*l2)**2) goto 34
                d = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
                if (d .le. dm) then
                    dm = d
                    im3 = ip1
                end if
34              continue
            end do
            if (im3 .eq. 0) goto 9997
!
!         -- VECTEUR M1M3 :
            m1m3(1) = nuax1((im3-1)*nx1+1)-nuax1((im1-1)* &
                                                 nx1+1)
            m1m3(2) = nuax1((im3-1)*nx1+2)-nuax1((im1-1)* &
                                                 nx1+2)
            m1m3(3) = nuax1((im3-1)*nx1+3)-nuax1((im1-1)* &
                                                 nx1+3)
!
!         -- N = M1M2 X M1M3 :
            n(1) = m1m2(2)*m1m3(3)-m1m2(3)*m1m3(2)
            n(2) = m1m2(3)*m1m3(1)-m1m2(1)*m1m3(3)
            n(3) = m1m2(1)*m1m3(2)-m1m2(2)*m1m3(1)
            s2 = n(1)**2+n(2)**2+n(3)**2
            s = sqrt(s2)
!
!         -- IM4 EST L'INDICE DU POINT LE + PROCHE DE P2
!            DIFFERENT DE M1 M2 M3 ET TEL QUE M1 M2 M3 M4
!            FORMENT UN VOLUME
            im4 = 0
            dm = vdm0(ip2)
            do ip1 = 1, np1
!           IF ((IP1.EQ.IM1).OR.(IP1.EQ.IM2).OR.(IP1.EQ.IM3)) GOTO 35
                if (.not. zl(inual1-1+(ip1-1)*nc1+ic1)) goto 35
                x1 = nuax1((ip1-1)*nx1+1)
                y1 = nuax1((ip1-1)*nx1+2)
                z1 = nuax1((ip1-1)*nx1+3)
!
!           SI LES POINTS M1 M2 M3 ET P1 NE FORMENT PAS
!           UN VOLUME GOTO 35
                m1p1(1) = nuax1((im1-1)*nx1+1)-x1
                m1p1(2) = nuax1((im1-1)*nx1+2)-y1
                m1p1(3) = nuax1((im1-1)*nx1+3)-z1
!
                v = abs(m1p1(1)*n(1)+m1p1(2)*n(2)+m1p1(3)*n(3))
                if (v .le. (1.d-3*s*l)) goto 35
!
                d = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
                if (d .le. dm) then
                    dm = d
                    im4 = ip1
                end if
35              continue
            end do
            if (im4 .eq. 0) goto 9998
!
            dref(ip2) = dm
30          continue
        end do
        goto 999
!
!
    else
        call utmess('F', 'UTILITAI2_54')
    end if
!
9994 continue
    call utmess('F', 'UTILITAI2_55')
!
9995 continue
    call utmess('F', 'UTILITAI2_56')
!
9996 continue
    call utmess('F', 'UTILITAI2_57')
!
9997 continue
    call utmess('F', 'UTILITAI2_58')
!
!
9998 continue
    call utmess('F', 'UTILITAI2_59')
!
!
999 continue
    AS_DEALLOCATE(vr=vdm0)
    call jedema()
end subroutine
