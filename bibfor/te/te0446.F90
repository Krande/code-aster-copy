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
subroutine te0446(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/coqrep.h"
#include "asterfort/dxbsig.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/terefe.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "blas/dcopy.h"
    character(len=16) :: option, nomte
!
!   CALCUL DES OPTIONS DES ELEMENTS DE PLAQUE POUR LA MODELISATION DKTG
!   ET LA MODELISATION Q4GG
!     -----------------------------------------------------------------
!                            TRIANGLE  QUADRANGLE
!        KIRCHOFF  (MINCE)      DKT       DKQ
!
!                  (EPAIS)      Q4G       T3G
!
!        OPTIONS     FORC_NODA
!
    integer(kind=8) :: nnos, ipoids, ivf, idfdx, jgano
    integer(kind=8) :: jtab(7), jvDisp
    integer(kind=8) :: icompo, i, i1, i2, j, k, ivectu, ipg, npg
    integer(kind=8) :: jvSief, iretc
    integer(kind=8) :: nno, igeom
    integer(kind=8) :: ndim, iret, ind
    integer(kind=8) :: jcara
    real(kind=8) :: pgl(3, 3), xyzl(3, 4), bsigma(24)
    real(kind=8) :: effgt(32), effort(32)
    real(kind=8) :: effref, momref
    real(kind=8) :: alpha, beta, t2ev(4), t2ve(4), c, s
    real(kind=8) :: foref, moref
    aster_logical :: reactu
    blas_int :: b_incx, b_incy, b_n
!
    if (option .eq. 'FORC_NODA') then
!
! ---   RECUPERATION DES ADRESSES DANS ZR DES POIDS DES PG
!       DES FONCTIONS DE FORME DES VALEURS DES DERIVEES DES FONCTIONS
!       DE FORME ET DE LA MATRICE DE PASSAGE GAUSS -> NOEUDS
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
!
        call jevech('PGEOMER', 'L', igeom)
!
        if (nno .eq. 3) then
            call dxtpgl(zr(igeom), pgl)
        else if (nno .eq. 4) then
            call dxqpgl(zr(igeom), pgl)
        end if
!
        call utpvgl(nno, 3, pgl, zr(igeom), xyzl)
!
        call tecach('ONO', 'PCOMPOR', 'L', iretc, iad=icompo)
!
! --- CALCUL DES MATRICES DE CHANGEMENT DE REPERES
!
!     T2EV : LA MATRICE DE PASSAGE (2X2) : UTILISATEUR -> INTRINSEQUE
!     T2VE : LA MATRICE DE PASSAGE (2X2) : INTRINSEQUE -> UTILISATEUR
!
        call jevech('PCACOQU', 'L', jcara)
        alpha = zr(jcara+1)*r8dgrd()
        beta = zr(jcara+2)*r8dgrd()
        call coqrep(pgl, alpha, beta, t2ev, t2ve, &
                    c, s)
!
! --- VECTEUR DES EFFORTS GENERALISES AUX POINTS
! --- D'INTEGRATION DU REPERE LOCAL
        call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, &
                    itab=jtab)
!
! --- PASSAGE DU VECTEUR DES EFFORTS GENERALISES AUX POINTS
! --- D'INTEGRATION DU REPERE LOCAL AU REPERE INTRINSEQUE
        do ipg = 1, npg
            jvSief = jtab(1)+8*(ipg-1)
            b_n = to_blas_int(8)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvSief), b_incx, effort(8*(ipg-1)+1), b_incy)
        end do
        call dxefro(npg, t2ve, effort, effgt)
!
        reactu = .false.
        if (iretc .eq. 0) then
            if (zk16(icompo+2) (6:10) .eq. '_REAC') call utmess('A', 'ELEMENTS2_72')
            reactu = (zk16(icompo+2) .eq. 'PETIT_REAC' .or. zk16(icompo+2) .eq. 'GROT_GDEP')
        end if
!
        if (reactu) then
            call jevech('PDEPLAR', 'L', jvDisp)
            do i = 1, nno
                i1 = 3*(i-1)
                i2 = 6*(i-1)
                zr(igeom+i1) = zr(igeom+i1)+zr(jvDisp+i2)
                zr(igeom+i1+1) = zr(igeom+i1+1)+zr(jvDisp+i2+1)
                zr(igeom+i1+2) = zr(igeom+i1+2)+zr(jvDisp+i2+2)
            end do
            if (nno .eq. 3) then
                call dxtpgl(zr(igeom), pgl)
            else if (nno .eq. 4) then
                call dxqpgl(zr(igeom), pgl)
            end if
!
            call utpvgl(nno, 3, pgl, zr(igeom), xyzl)
        end if
!
! --- CALCUL DES EFFORTS INTERNES (I.E. SOMME_VOL(BT_SIG))
        call dxbsig(nomte, xyzl, pgl, effgt, bsigma, &
                    option)
!
! --- AFFECTATION DES VALEURS DE BSIGMA AU VECTEUR EN SORTIE
        call jevech('PVECTUR', 'E', ivectu)
!
        k = 0
        do i = 1, nno
            do j = 1, 6
                k = k+1
                zr(ivectu+k-1) = bsigma(k)
            end do
        end do
    else if (option .eq. 'REFE_FORC_NODA') then
!     -------------------------------------
!
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
        call jevech('PGEOMER', 'L', igeom)
!
        if (nno .eq. 3) then
            call dxtpgl(zr(igeom), pgl)
        else if (nno .eq. 4) then
            call dxqpgl(zr(igeom), pgl)
        end if
!
        call utpvgl(nno, 3, pgl, zr(igeom), xyzl)
!
        call terefe('EFFORT_REFE', 'MECA_COQUE', foref)
        call terefe('MOMENT_REFE', 'MECA_COQUE', moref)
!
        ind = 8
        do i = 1, nno
            do j = 1, 3
                effgt((i-1)*ind+j) = foref
                effgt((i-1)*ind+3+j) = moref
                effgt((i-1)*ind+7) = foref
                effgt((i-1)*ind+8) = foref
            end do
        end do
!
! ------ CALCUL DES EFFORTS INTERNES (I.E. SOMME_VOL(BT_SIG))
        call dxbsig(nomte, xyzl, pgl, effgt, bsigma, &
                    option)
!
! ------ AFFECTATION DES VALEURS DE BSIGMA AU VECTEUR EN SORTIE
        call jevech('PVECTUR', 'E', ivectu)
        k = 0
        do i = 1, nno
            effref = (abs(bsigma(k+1))+abs(bsigma(k+2))+abs(bsigma(k+3)))/3.d0
            momref = (abs(bsigma(k+4))+abs(bsigma(k+5))+abs(bsigma(k+6)))/3.d0
            ASSERT(abs(effref) .gt. r8prem())
            ASSERT(abs(momref) .gt. r8prem())
            do j = 1, 6
                k = k+1
                if (j .lt. 4) then
                    zr(ivectu+k-1) = effref
                else
                    zr(ivectu+k-1) = momref
                end if
            end do
        end do
    else
        ASSERT(.false.)
    end if
end subroutine
