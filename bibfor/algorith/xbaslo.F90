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
subroutine xbaslo(noma, fiss, grlt, grln, ndim)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterc/r8dgrd.h"
#include "blas/ddot.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/assert.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/imprsd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/trigom.h"
    integer(kind=8) :: ndim
    character(len=8) :: noma, fiss
    character(len=19) :: grlt, grln
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (PREPARATION)
!
! CREATION D'UN CHAM_EL QUI CONTIENT LA BASE
! LOCALE AU POINT DU FOND DE FISSURE ASSOCIE - (2D/3D)
!
! ----------------------------------------------------------------------
!
!
! IN  FISS   : NOM DE LA FISSURE
! IN  NOMA   : NOM DU MAILLAGE
! IN  MODELE : NOM DE L'OBJET MODELE
! IN  GRLT   : CHAM_NO_S DES GRADIENTS DE LA LEVEL-SET TANGENTE
! IN  GRLN   : CHAM_NO_S DES GRADIENTS DE LA LEVEL-SET NORMALE
! IN  NDIM   : DIMENSION DU MAILLAGE
! OUT FISS   : FISSURE AVEC LE .BASLOC EN PLUS
!    .BASLOC CONTIENT ::
!      VR(1:NDIM)        : LE PROJETE EN FOND DE FISSURE
!      VR(NDIM+1:3*NDIM) : LA BASE LOCALE ASSOCIEE A PROJETE EN FOND DE FISSURE
!
!
    integer(kind=8) :: nbcmp
    character(len=8) :: licmp(9)
    character(len=24) :: coorn
    character(len=24) :: xfonfi, xbasfo
    integer(kind=8) :: ifon, npoint, ifm, niv, ier, ibid, ibas
    character(len=19) :: cnsbas, basloc
    integer(kind=8) :: iadrco, jgsl, jgtl
    integer(kind=8) :: long, nfon, nbno, ino, j
    integer(kind=8) :: nbfron, jnfon, ni, nf
    real(kind=8) :: xi1, yi1, zi1, xj1, yj1, zj1, xij, yij, zij, eps, d, norm2
    real(kind=8) :: xm, ym, zm, xim, yim, zim, s, dmin, xn, yn, zn, a(3)
    real(kind=8) :: u(3), v(3), un(3), vn(3)
    real(kind=8) :: ui(3), uf(3), angle_max, theta, cosi
    real(kind=8), pointer :: gn(:) => null()
    real(kind=8), pointer :: gsv(:) => null()
    real(kind=8), pointer :: gt(:) => null()
    aster_logical, pointer :: is_continu(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    data licmp/'X1', 'X2', 'X3',&
     &             'X4', 'X5', 'X6',&
     &             'X7', 'X8', 'X9'/
!
!  DEFINITION DU PARAMETRE POUR DETERMINER SI UN FRONT MULTIPLE EST DISCONTINU ::
!     ON PREND UN ANGLE MAXIMUM ENTRE LES VECTEURS DE PROPAGATION DES POINTS AUX EXTREMITES
!     EGAL A 10° AU MAXIMUM. CE SEUIL POURRA ETRE MODIFIE ULTERIEUREMENT...
    parameter(angle_max=10.)
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('XFEM', ifm, niv)
!
! --- CREATION DU CHAM_NO
!
    cnsbas = '&&OP0041.CNSBAS'
    nbcmp = ndim*3
    call cnscre(noma, 'NEUT_R', nbcmp, licmp, 'V', &
                cnsbas)
    call jeveuo(cnsbas//'.CNSV', 'E', vr=gsv)
    call jeveuo(cnsbas//'.CNSL', 'E', jgsl)
!
! --- ACCES AU MAILLAGE
!
    coorn = noma//'.COORDO    .VALE'
    call jeveuo(coorn, 'L', iadrco)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbno)
!
! --- ACCES AUX OBJETS: NOM DES POINTS DU FOND DE FISSURE
!
    xfonfi = fiss(1:8)//'.FONDFISS'
    xbasfo = fiss(1:8)//'.BASEFOND'
    call jeexin(xfonfi, ier)
!
    if (ier .eq. 0) then
!       LE FOND DE FISSURE N'EXISTE PAS (CAS D'UNE INTERFACE)
!       ON MET TOUT A ZERO ET ON SORT
        do ino = 1, nbno
            do j = 1, ndim
                gsv(3*ndim*(ino-1)+j) = 0.d0
                zl(jgsl-1+3*ndim*(ino-1)+j) = .true.
                gsv(3*ndim*(ino-1)+j+ndim) = 0.d0
                zl(jgsl-1+3*ndim*(ino-1)+j+ndim) = .true.
                gsv(3*ndim*(ino-1)+j+2*ndim) = 0.d0
                zl(jgsl-1+3*ndim*(ino-1)+j+2*ndim) = .true.
            end do
        end do
        goto 999
    end if
!
    call jeveuo(xfonfi, 'L', ifon)
    call jelira(xfonfi, 'LONMAX', long)
    nfon = long/4
    call jeveuo(xbasfo, 'L', ibas)
!
! EN CAS DE FOND MULTIPLES IL FAUT FAIRE ATTENTION A LA PROJECTION SUR LE FRONT
!   ON VERIFIE A MINIMA SI LES FRONT SONT CONTINUS
    if (ndim .eq. 3) then
        call jelira(fiss(1:8)//'.FONDMULT', 'LONMAX', long)
        nbfron = long/2
        AS_ALLOCATE(vl=is_continu, size=nfon-1)
        is_continu(:) = .true.
        call jeveuo(fiss(1:8)//'.FONDMULT', 'L', jnfon)
        do j = 1, nbfron-1
            ni = zi(jnfon-1+2*(j-1)+2)
            nf = zi(jnfon-1+2*(j+1-1)+1)
            ASSERT(nf .eq. ni+1)
            ui(1) = zr(ibas-1+6*(ni-1)+4)
            ui(2) = zr(ibas-1+6*(ni-1)+5)
            ui(3) = zr(ibas-1+6*(ni-1)+6)
            uf(1) = zr(ibas-1+6*(nf-1)+4)
            uf(2) = zr(ibas-1+6*(nf-1)+5)
            uf(3) = zr(ibas-1+6*(nf-1)+6)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            cosi = ddot(b_n, ui, b_incx, uf, b_incy)/(sqrt(ddot(b_n, ui, b_incx, ui, b_incy))*sq&
                   &rt(ddot(b_n, uf, b_incx, uf, b_incy)))
            theta = trigom('ACOS', cosi)
            is_continu(ni) = abs(theta) .le. (angle_max*r8dgrd())
        end do
    end if
!
! --- RÉCUPÉRATION DES GRADIENTS DE LST ET LSN
!
    call jeveuo(grlt//'.CNSV', 'L', vr=gt)
    call jeveuo(grlt//'.CNSL', 'L', jgtl)
    call jeveuo(grln//'.CNSV', 'L', vr=gn)
!
!     CALCUL DES PROJETÉS DES NOEUDS SUR LE FOND DE FISSURE
    eps = 1.d-12
    do ino = 1, nbno
        if (.not. zl(jgtl-1+ndim*(ino-1)+1)) goto 100
!       COORD DU NOEUD M DU MAILLAGE
        xm = zr(iadrco+(ino-1)*3+1-1)
        ym = zr(iadrco+(ino-1)*3+2-1)
        zm = zr(iadrco+(ino-1)*3+3-1)
!       INITIALISATION
        dmin = r8maem()
        u(:) = 0.d0
        v(:) = 0.d0
!       BOUCLE SUR PT DE FONFIS
        if (ndim .eq. 2) npoint = nfon
        if (ndim .eq. 3) npoint = nfon-1
        do j = 1, npoint
            if (ndim .eq. 2) then
!           COORD PT N
                xn = zr(ifon-1+4*(j-1)+1)
                yn = zr(ifon-1+4*(j-1)+2)
                zn = 0.d0
!           BASE AU PT N
                un(1) = zr(ibas-1+4*(j-1)+1)
                un(2) = zr(ibas-1+4*(j-1)+2)
                un(3) = 0.d0
                vn(1) = zr(ibas-1+4*(j-1)+3)
                vn(2) = zr(ibas-1+4*(j-1)+4)
                vn(3) = 0.d0
!           DISTANCE MN
                d = sqrt((xn-xm)*(xn-xm)+(yn-ym)*(yn-ym))
            else if (ndim .eq. 3) then
                if (.not. is_continu(j)) goto 200
!           COORD PT I, ET J
                xi1 = zr(ifon-1+4*(j-1)+1)
                yi1 = zr(ifon-1+4*(j-1)+2)
                zi1 = zr(ifon-1+4*(j-1)+3)
                xj1 = zr(ifon-1+4*(j-1+1)+1)
                yj1 = zr(ifon-1+4*(j-1+1)+2)
                zj1 = zr(ifon-1+4*(j-1+1)+3)
!           VECTEUR IJ ET IM
                xij = xj1-xi1
                yij = yj1-yi1
                zij = zj1-zi1
                xim = xm-xi1
                yim = ym-yi1
                zim = zm-zi1
!           PARAM S (PRODUIT SCALAIRE...)
                s = xij*xim+yij*yim+zij*zim
                norm2 = xij*xij+yij*yij+zij*zij
                s = s/norm2
!           SI N=P(M) SORT DU SEGMENT
                if ((s-1) .ge. eps) s = 1.d0
                if (s .le. eps) s = 0.d0
!           COORD DE N
                xn = s*xij+xi1
                yn = s*yij+yi1
                zn = s*zij+zi1
!           DISTANCE MN
                d = sqrt((xn-xm)*(xn-xm)+(yn-ym)*(yn-ym)+(zn-zm)*(zn-zm))
!           BASE AU PT N
                un(1) = s*zr(ibas-1+6*(j-1+1)+1)+(1-s)*zr(ibas-1+6*(j-1)+1)
                un(2) = s*zr(ibas-1+6*(j-1+1)+2)+(1-s)*zr(ibas-1+6*(j-1)+2)
                un(3) = s*zr(ibas-1+6*(j-1+1)+3)+(1-s)*zr(ibas-1+6*(j-1)+3)
                vn(1) = s*zr(ibas-1+6*(j-1+1)+4)+(1-s)*zr(ibas-1+6*(j-1)+4)
                vn(2) = s*zr(ibas-1+6*(j-1+1)+5)+(1-s)*zr(ibas-1+6*(j-1)+5)
                vn(3) = s*zr(ibas-1+6*(j-1+1)+6)+(1-s)*zr(ibas-1+6*(j-1)+6)
            end if
            if (d .lt. (dmin*(1-abs(r8prem())*100))) then
                dmin = d
                a(1) = xn
                a(2) = yn
                a(3) = zn
                u(1:3) = un(1:3)
                v(1:3) = vn(1:3)
            end if
200         continue
        end do
!       STOCKAGE DU PROJETÉ ET DES GRADIENTS
        do j = 1, ndim
            gsv(3*ndim*(ino-1)+j) = a(j)
            zl(jgsl-1+3*ndim*(ino-1)+j) = .true.
!            gsv(3*ndim*(ino-1)+j+ndim)=gt(ndim*(ino-1)+j)
            gsv(3*ndim*(ino-1)+j+ndim) = v(j)
            zl(jgsl-1+3*ndim*(ino-1)+j+ndim) = .true.
!            gsv(3*ndim*(ino-1)+j+2*ndim)=gn(ndim*(ino-1)+j)
            gsv(3*ndim*(ino-1)+j+2*ndim) = u(j)
            zl(jgsl-1+3*ndim*(ino-1)+j+2*ndim) = .true.
        end do
100     continue
    end do
!
999 continue
!
!     ENREGISTREMENT DU .BASLOC DANS LA SD FISS_XFEM
    basloc = fiss(1:8)//'.BASLOC'
    call cnscno(cnsbas, basloc(1:13)//'.NUMEQ', 'NON', 'G', basloc, &
                'F', ibid)
    call detrsd('CHAM_NO_S', cnsbas)
    if (ndim .eq. 3) AS_DEALLOCATE(vl=is_continu)
!
    if (niv .gt. 2) then
        call imprsd('CHAMP', basloc, ifm, 'FISSURE.BASLOC=')
    end if
!
    call jedema()
end subroutine
