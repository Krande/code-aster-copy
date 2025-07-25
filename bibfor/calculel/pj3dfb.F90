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

subroutine pj3dfb(boite, maillz, geom1, geom2)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utimsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    real(kind=8) :: geom1(*), geom2(*)
    character(len=14) :: boite
    character(len=*) :: maillz
!     BUT :
!       CONSTRUIRE LA STRUCTURE DE DONNEES BOITE_3D QUI PERMET DE SAVOIR
!       QUELLES SONT LES MAILLES QUI SE TROUVENT DANS UNE BOITE(P,Q,R)
!
!  IN/JXOUT   BOITE      K14 : NOM DE LA SD BOITE_3D A CREER
!  IN         GEOM2(*)   R  : COORDONNEES DES NOEUDS DU MAILLAGE M2
!  IN         GEOM1(*)   R  : COORDONNEES DES NOEUDS DU MAILLAGE M1
!  IN         MAILLZ   K*  : OBJET '&&PJXXCO.TETR4' OU '&&PJXXCO.TRIA3'
!                               OU '&&PJXXCO.SEG2'
! ----------------------------------------------------------------------
!
!
    real(kind=8) :: stotal, dx, dy, ddx, ddy, rbig, xxmax, xxmin, xmax, xmin
    real(kind=8) :: yymax, yymin, ymax, ymin
    real(kind=8) :: dz, ddz, zmin, zmax, zzmin, zzmax
    character(len=24) :: maille
    integer(kind=8) :: p1, q1, r1, p2, q2, r2, p, q, r, nx, ny, nz, ndec, nno
    integer(kind=8) :: iatr3, ntr3, nno1, nno2, i, iposi, ifm, niv, vali(2)
    integer(kind=8) :: iabtdi, iabtvr, iabtnb, iabtlc, k, ino, ib, lont, iabtco
    integer(kind=8) :: nbtot, nbmax, nbmin, nbtet, nbmoy
    aster_logical :: dbg
    integer(kind=8), pointer :: lino1(:) => null()
    integer(kind=8), pointer :: lino2(:) => null()
!
! DEB ------------------------------------------------------------------
    call jemarq()
    call infniv(ifm, niv)
    dbg = ASTER_FALSE
    maille = maillz
    if (maille(10:14) .eq. 'TRIA3') then
        ndec = 4
        nno = 3
    else if (maille(10:14) .eq. 'TETR4') then
        ndec = 6
        nno = 4
    else if (maille(10:13) .eq. 'SEG2') then
        ndec = 3
        nno = 2
    else
        ASSERT(ASTER_FALSE)
    end if
    call jeveuo(maille, 'L', iatr3)
    ntr3 = zi(iatr3-1+1)
    rbig = r8maem()
    if (ntr3 .eq. 0) then
        call utmess('F', 'CALCULEL4_57')
    end if
!
    call jeveuo('&&PJXXCO.LINO1', 'L', vi=lino1)
    call jeveuo('&&PJXXCO.LINO2', 'L', vi=lino2)
    call jelira('&&PJXXCO.LINO1', 'LONMAX', nno1)
    call jelira('&&PJXXCO.LINO2', 'LONMAX', nno2)
!
!
!     1. : ON CALCULE XMIN,XMAX,YMIN,YMAX,NX,NY,NZ,DX,DY,DZ...
!     -------------------------------------------------------
    xmin = rbig
    ymin = rbig
    zmin = rbig
    xmax = -rbig
    ymax = -rbig
    zmax = -rbig
    do i = 1, nno1
        if (lino1(i) .eq. 0) cycle
        xmin = min(xmin, geom1(3*(i-1)+1))
        xmax = max(xmax, geom1(3*(i-1)+1))
        ymin = min(ymin, geom1(3*(i-1)+2))
        ymax = max(ymax, geom1(3*(i-1)+2))
        zmin = min(zmin, geom1(3*(i-1)+3))
        zmax = max(zmax, geom1(3*(i-1)+3))
    end do
    do i = 1, nno2
        if (lino2(i) .eq. 0) cycle
        xmin = min(xmin, geom2(3*(i-1)+1))
        xmax = max(xmax, geom2(3*(i-1)+1))
        ymin = min(ymin, geom2(3*(i-1)+2))
        ymax = max(ymax, geom2(3*(i-1)+2))
        zmin = min(zmin, geom2(3*(i-1)+3))
        zmax = max(zmax, geom2(3*(i-1)+3))
    end do
!
!
    stotal = max((xmax-xmin), (ymax-ymin), (zmax-zmin))
    if (stotal .eq. 0.d0) then
        call utmess('F', 'CALCULEL4_58')
    end if
    dx = 2.d0*stotal/(dble(ntr3)**(1.d0/3.d0))
!
    dy = dx
    dz = dx
!
    nx = int((xmax-xmin)*1.05d0/dx)+1
    ny = int((ymax-ymin)*1.05d0/dy)+1
    nz = int((zmax-zmin)*1.05d0/dz)+1
    ASSERT(nx*ny*nz .ne. 0)
    ddx = (nx*dx-(xmax-xmin))/2.d0
    ddy = (ny*dy-(ymax-ymin))/2.d0
    ddz = (nz*dz-(zmax-zmin))/2.d0
    xmin = xmin-ddx
    xmax = xmax+ddx
    ymin = ymin-ddy
    ymax = ymax+ddy
    zmin = zmin-ddz
    zmax = zmax+ddz
!
    if (dbg) then
        write (ifm, *)
        write (ifm, *) '-----------------------------------------'
        write (ifm, *) ' MISE EN BOITES DES ELEMENTS DU MODELE_1'
        write (ifm, *) '-----------------------------------------'
        write (ifm, *)
        write (ifm, *) 'ZONE DE TRAVAIL : XMIN,XMAX=', xmin, xmax
        write (ifm, *) '                  YMIN,YMAX=', ymin, ymax
        write (ifm, *) '                  ZMIN,ZMAX=', zmin, ymax
        write (ifm, *)
        write (ifm, *) 'NOMBRE DE BOITES :'
        write (ifm, *) '  DANS LES DIRECTIONS X,Y,Z :', nx, ' ', ny, ' ', nz
        write (ifm, *) '  TOTAL                     :', nx*ny*nz
        write (ifm, *)
        write (ifm, *) 'DIMENSIONS DES BOITES LX,LY,LZ=', dx, dy, dz
    end if
!
!
!     2. : ALLOCATION DE LA SD BOITE_3D :
!     ---------------------------------------
    call wkvect(boite//'.BT3DDI', 'V V I', 3, iabtdi)
    call wkvect(boite//'.BT3DVR', 'V V R', 9, iabtvr)
    call wkvect(boite//'.BT3DNB', 'V V I', nx*ny*nz, iabtnb)
    call wkvect(boite//'.BT3DLC', 'V V I', 1+nx*ny*nz, iabtlc)
!
    zi(iabtdi-1+1) = nx
    zi(iabtdi-1+2) = ny
    zi(iabtdi-1+3) = nz
!
    zr(iabtvr-1+1) = xmin
    zr(iabtvr-1+2) = xmax
    zr(iabtvr-1+3) = ymin
    zr(iabtvr-1+4) = ymax
    zr(iabtvr-1+5) = zmin
    zr(iabtvr-1+6) = zmax
    zr(iabtvr-1+7) = dx
    zr(iabtvr-1+8) = dy
    zr(iabtvr-1+9) = dz
!
!
!
!     3. : ON COMPTE COMBIEN DE MAILLES SERONT CONTENUES
!             DANS CHAQUE BOITE(P,Q,R)
!     -------------------------------------------------------
    do i = 1, ntr3
        xxmin = rbig
        yymin = rbig
        zzmin = rbig
        xxmax = -rbig
        yymax = -rbig
        zzmax = -rbig
        do k = 1, nno
            ino = zi(iatr3+ndec*(i-1)+k)
            xxmin = min(xxmin, geom1(3*(ino-1)+1))
            xxmax = max(xxmax, geom1(3*(ino-1)+1))
            yymin = min(yymin, geom1(3*(ino-1)+2))
            yymax = max(yymax, geom1(3*(ino-1)+2))
            zzmin = min(zzmin, geom1(3*(ino-1)+3))
            zzmax = max(zzmax, geom1(3*(ino-1)+3))
        end do
        p1 = int((xxmin-xmin)/dx)+1
        p2 = int((xxmax-xmin)/dx)+1
        q1 = int((yymin-ymin)/dy)+1
        q2 = int((yymax-ymin)/dy)+1
        r1 = int((zzmin-zmin)/dz)+1
        r2 = int((zzmax-zmin)/dz)+1
        do p = p1, p2
            do q = q1, q2
                do r = r1, r2
                    zi(iabtnb-1+(r-1)*nx*ny+(q-1)*nx+p) = zi(iabtnb-1+(r-1)*nx*ny+(q-1)*nx+p &
                                                             )+1
                end do
            end do
        end do
!
    end do
!
!   3.2: calcul du nombre de tetraedres par boite :
!   -----------------------------------------------
    nbtot = 0
    nbmax = 0
    nbmin = ismaem()
    do p = 1, nx
        do q = 1, ny
            do r = 1, nz
                nbtet = zi(iabtnb-1+(r-1)*nx*ny+(q-1)*nx+p)
                if (dbg) write (ifm, *) 'P,Q,R,NBTET=', p, q, r, nbtet
                nbtot = nbtot+nbtet
                nbmin = min(nbmin, nbtet)
                nbmax = max(nbmax, nbtet)
            end do
        end do
    end do
    nbmoy = nbtot/(nx*ny*nz)
    if (dbg) then
        write (ifm, *)
        write (ifm, *) 'NOMBRE DE TETRAEDRES PAR BOITE:'
        write (ifm, *) '   EN MOYENNE :', nbmoy
        write (ifm, *) '   MIN        :', nbmin
        write (ifm, *) '   MAX        :', nbmax
    end if

    if (nbmax .gt. 20*nbmoy) then
        vali(1) = nbmoy
        vali(2) = nbmax
        call utmess('I', 'CALCULEL5_79', ni=2, vali=vali)
    end if
!
!
!
!   4. : ON REMPLIT .BT3DCO  ET .BT3DLC :
!   -------------------------------------------------------
    zi(iabtlc-1+1) = 0
    do ib = 1, nx*ny*nz
        zi(iabtlc-1+ib+1) = zi(iabtlc-1+ib)+zi(iabtnb-1+ib)
        zi(iabtnb-1+ib) = 0
    end do
!
!
    lont = zi(iabtlc-1+1+nx*ny*nz)
    call wkvect(boite//'.BT3DCO', 'V V I', lont, iabtco)
!
    do i = 1, ntr3
        xxmin = rbig
        yymin = rbig
        zzmin = rbig
        xxmax = -rbig
        yymax = -rbig
        zzmax = -rbig
        do k = 1, nno
            ino = zi(iatr3+ndec*(i-1)+k)
            xxmin = min(xxmin, geom1(3*(ino-1)+1))
            xxmax = max(xxmax, geom1(3*(ino-1)+1))
            yymin = min(yymin, geom1(3*(ino-1)+2))
            yymax = max(yymax, geom1(3*(ino-1)+2))
            zzmin = min(zzmin, geom1(3*(ino-1)+3))
            zzmax = max(zzmax, geom1(3*(ino-1)+3))
        end do
        p1 = int((xxmin-xmin)/dx)+1
        p2 = int((xxmax-xmin)/dx)+1
        q1 = int((yymin-ymin)/dy)+1
        q2 = int((yymax-ymin)/dy)+1
        r1 = int((zzmin-zmin)/dz)+1
        r2 = int((zzmax-zmin)/dz)+1
        do p = p1, p2
            do q = q1, q2
                do r = r1, r2
                    zi(iabtnb-1+(r-1)*nx*ny+(q-1)*nx+p) = zi(iabtnb-1+(r-1)*nx*ny+(q-1)*nx+p &
                                                             )+1
                    iposi = zi( &
                            iabtlc-1+(r-1)*nx*ny+(q-1)*nx+p)+zi(iabtnb-1+(r-1)*nx*ny+(q-1)*&
                            &nx+p &
                            )
                    ASSERT((iposi .ge. 1) .and. (iposi .le. lont))
                    zi(iabtco-1+iposi) = i
                end do
            end do
        end do
!
    end do
!
    if (dbg) call utimsd(ifm, 2, ASTER_FALSE, ASTER_TRUE, boite, &
                         1, ' ')
    call jedema()
end subroutine
