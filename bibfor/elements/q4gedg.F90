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
subroutine q4gedg(xyzl, option, pgl, depl, edgl)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dsqbfb.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxqbm.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gquad4.h"
#include "asterfort/jquad4.h"
#include "asterfort/q4gbc.h"
    real(kind=8) :: xyzl(3, *), pgl(3, *), depl(*), edgl(*)
    character(len=16) :: option
!     EFFORTS ET DEFORMATIONS GENERALISES DE L'ELEMENT DE PLAQUE Q4GAMMA
!     ------------------------------------------------------------------
!     IN  XYZL   : COORDONNEES LOCALES DES QUATRE NOEUDS
!     IN  OPTION : NOM DE L'OPTION DE CALCUL
!     IN  PGL    : MATRICE DE PASSAGE GLOBAL - LOCAL
!     IN  DEPL   : DEPLACEMENTS
!     OUT EDGL   : EFFORTS OU DEFORMATIONS GENERALISES AUX NOEUDS DANS
!                  LE REPERE INTRINSEQUE A L'ELEMENT
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: multic, ne, k, j, i, ie
    real(kind=8) :: depf(12), depm(8)
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dmc(3, 2), dfc(3, 2)
    real(kind=8) :: dci(2, 2), dc(2, 2)
    real(kind=8) :: bf(3, 12), bm(3, 8), bc(2, 12)
    real(kind=8) :: bdm(3), bdf(3), bcdf(2), dcis(2)
    real(kind=8) :: vf(3), vm(3), vt(2)
    real(kind=8) :: vfm(3), vmf(3), vmc(3), vfc(3), vcm(2), vcf(2), caraq4(25)
    real(kind=8) :: t2iu(4), t2ui(4), t1ve(9), jacob(5), qsi, eta
    aster_logical :: coupmf
    character(len=4) :: fami
!     ------------------------------------------------------------------
!
    if (option(6:9) .eq. 'ELGA') then
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                         jgano=jgano)
        ne = npg
        fami = 'RIGI'
    else if (option(6:9) .eq. 'ELNO') then
        call elrefe_info(fami='NOEU', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                         jgano=jgano)
        ne = nno
        fami = 'NOEU'
    end if
!
!     ----- CALCUL DES MATRICES DE RIGIDITE DU MATERIAU EN FLEXION,
!           MEMBRANE ET CISAILLEMENT INVERSEES -------------------------
!
!     ----- CALCUL DES GRANDEURS GEOMETRIQUES SUR LE QUADRANGLE --------
    call gquad4(xyzl, caraq4)
!     ----- CARACTERISTIQUES DES MATERIAUX --------
    call dxmate(fami, df, dm, dmf, dc, &
                dci, dmc, dfc, nno, pgl, &
                multic, coupmf, t2iu, t2ui, t1ve)
!     ----- COMPOSANTES DEPLACEMENT MEMBRANE ET FLEXION ----------------
    do j = 1, nno
        do i = 1, 2
            depm(i+2*(j-1)) = depl(i+6*(j-1))
        end do
        depf(1+3*(j-1)) = depl(1+2+6*(j-1))
        depf(2+3*(j-1)) = depl(3+2+6*(j-1))
        depf(3+3*(j-1)) = -depl(2+2+6*(j-1))
    end do
    if (option(1:4) .eq. 'DEGE') then
        do ie = 1, ne
!
            qsi = zr(icoopg-1+ndim*(ie-1)+1)
            eta = zr(icoopg-1+ndim*(ie-1)+2)
!
!           ----- CALCUL DU JACOBIEN SUR LE QUADRANGLE -----------------
            call jquad4(xyzl, qsi, eta, jacob)
!           ----- CALCUL DE LA MATRICE BM ------------------------------
            call dxqbm(qsi, eta, jacob(2), bm)
!           ----- CALCUL DE LA MATRICE BF AU POINT QSI ETA -------------
            call dsqbfb(qsi, eta, jacob(2), bf)
!           ---- CALCUL DE LA MATRICE BC AU POINT QSI ETA --------------
            call q4gbc(qsi, eta, jacob(2), caraq4, bc)
!           ------ BCDF = BC.DEPF -------------------------------------
            bcdf(1) = 0.d0
            bcdf(2) = 0.d0
            do j = 1, 12
                bcdf(1) = bcdf(1)+bc(1, j)*depf(j)
                bcdf(2) = bcdf(2)+bc(2, j)*depf(j)
            end do
            do k = 1, 3
                bdf(k) = 0.d0
                bdm(k) = 0.d0
            end do
            do i = 1, 3
                do j = 1, 12
                    bdf(i) = bdf(i)+bf(i, j)*depf(j)
                end do
                do j = 1, 8
                    bdm(i) = bdm(i)+bm(i, j)*depm(j)
                end do
            end do
            do i = 1, 3
                edgl(i+8*(ie-1)) = bdm(i)
                edgl(i+3+8*(ie-1)) = bdf(i)
            end do
            edgl(7+8*(ie-1)) = bcdf(1)
            edgl(8+8*(ie-1)) = bcdf(2)
!           --- PASSAGE DE LA DISTORSION A LA DEFORMATION DE CIS. ------
            edgl(3+8*(ie-1)) = edgl(3+8*(ie-1))/2.d0
            edgl(6+8*(ie-1)) = edgl(6+8*(ie-1))/2.d0
            edgl(7+8*(ie-1)) = bcdf(1)/2.d0
            edgl(8+8*(ie-1)) = bcdf(2)/2.d0
        end do
    else
        do ie = 1, ne
!
            qsi = zr(icoopg-1+ndim*(ie-1)+1)
            eta = zr(icoopg-1+ndim*(ie-1)+2)
!
!           ----- CALCUL DU JACOBIEN SUR LE QUADRANGLE -----------------
            call jquad4(xyzl, qsi, eta, jacob)
!           ----- CALCUL DE LA MATRICE BM ------------------------------
            call dxqbm(qsi, eta, jacob(2), bm)
!           ----- CALCUL DE LA MATRICE BF AU POINT QSI ETA -------------
            call dsqbfb(qsi, eta, jacob(2), bf)
!           ---- CALCUL DE LA MATRICE BC AU POINT QSI ETA --------------
            call q4gbc(qsi, eta, jacob(2), caraq4, bc)
!           ------ VT = DC.BC.DEPF -------------------------------------
            vt(1) = 0.d0
            vt(2) = 0.d0
            bcdf(1) = 0.d0
            bcdf(2) = 0.d0
            do j = 1, 12
                bcdf(1) = bcdf(1)+bc(1, j)*depf(j)
                bcdf(2) = bcdf(2)+bc(2, j)*depf(j)
            end do
            vt(1) = dc(1, 1)*bcdf(1)+dc(1, 2)*bcdf(2)
            vt(2) = dc(2, 1)*bcdf(1)+dc(2, 2)*bcdf(2)
            do k = 1, 3
                bdf(k) = 0.d0
                bdm(k) = 0.d0
                vf(k) = 0.d0
                vm(k) = 0.d0
                vfm(k) = 0.d0
                vmf(k) = 0.d0
                vmc(k) = 0.0d0
                vfc(k) = 0.0d0
            end do
            vcm(1) = 0.0d0
            vcm(2) = 0.0d0
            vcf(1) = 0.0d0
            vcf(2) = 0.0d0
!           ------ VF = DF.BF.DEPF , VFM = DMF.BM.DEPM ----------------
!           ------ VM = DM.BM.DEPM , VMF = DMF.BF.DEPF ----------------
            do i = 1, 3
                do j = 1, 12
                    bdf(i) = bdf(i)+bf(i, j)*depf(j)
                end do
                do j = 1, 8
                    bdm(i) = bdm(i)+bm(i, j)*depm(j)
                end do
            end do
            do i = 1, 3
                do j = 1, 3
                    vf(i) = vf(i)+df(i, j)*bdf(j)
                    vfm(i) = vfm(i)+dmf(i, j)*bdm(j)
                    vm(i) = vm(i)+dm(i, j)*bdm(j)
                    vmf(i) = vmf(i)+dmf(i, j)*bdf(j)
                end do
            end do
!
            dcis(1) = dci(1, 1)*vt(1)+dci(1, 2)*vt(2)
            dcis(2) = dci(2, 1)*vt(1)+dci(2, 2)*vt(2)
!
            vmc(1) = dmc(1, 1)*dcis(1)+dmc(1, 2)*dcis(2)
            vmc(2) = dmc(2, 1)*dcis(1)+dmc(2, 2)*dcis(2)
            vmc(3) = dmc(3, 1)*dcis(1)+dmc(3, 2)*dcis(2)
!
            vcm(1) = dmc(1, 1)*vm(1)+dmc(2, 1)*vm(2)+dmc(3, 1)*vm(3)
            vcm(2) = dmc(1, 2)*vm(1)+dmc(2, 2)*vm(2)+dmc(3, 2)*vm(3)
!
            vfc(1) = dfc(1, 1)*dcis(1)+dfc(1, 2)*dcis(2)
            vfc(2) = dfc(2, 1)*dcis(1)+dfc(2, 2)*dcis(2)
            vfc(3) = dfc(3, 1)*dcis(1)+dfc(3, 2)*dcis(2)
!
            vcf(1) = dfc(1, 1)*vf(1)+dfc(2, 1)*vf(2)+dfc(3, 1)*vf(3)
            vcf(2) = dfc(1, 2)*vf(1)+dfc(2, 2)*vf(2)+dfc(3, 2)*vf(3)
!
            do i = 1, 3
                edgl(i+8*(ie-1)) = vm(i)+vmf(i)+vmc(i)
                edgl(i+3+8*(ie-1)) = vf(i)+vfm(i)+vfc(i)
            end do
            edgl(7+8*(ie-1)) = vt(1)+vcm(1)+vcf(1)
            edgl(8+8*(ie-1)) = vt(2)+vcm(2)+vcf(2)
        end do
    end if
!
end subroutine
