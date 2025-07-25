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

subroutine t3grig(nomte, xyzl, option, pgl, rig, &
                  ener)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/bsthpl.h"
#include "asterfort/dstbfb.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxtbm.h"
#include "asterfort/dxtloc.h"
#include "asterfort/dxtloe.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/r8inir.h"
#include "asterfort/t3gbc.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
    real(kind=8) :: xyzl(3, *), pgl(*), rig(*), ener(*)
    character(len=16) :: option, nomte
!
!     MATRICE DE RIGIDITE DE L'ELEMENT T3GAMMA (AVEC CISAILLEMENT)
!     ------------------------------------------------------------------
!     IN  XYZL   : COORDONNEES LOCALES DES QUATRE NOEUDS
!     IN  OPTION : OPTION RIGI_MECA OU EPOT_ELEM
!     IN  PGL    : MATRICE DE PASSAGE GLOBAL/LOCAL
!     OUT RIG    : MATRICE DE RIGIDITE
!     OUT ENER   : TERMES POUR ENER_POT (EPOT_ELEM)
!     ------------------------------------------------------------------
    integer(kind=8) :: multic
    real(kind=8) :: depl(18)
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: bfb(3, 9)
    real(kind=8) :: bc(2, 9)
    real(kind=8) :: bm(3, 6)
    real(kind=8) :: xab1(3, 6), xab2(3, 9), xab3(2, 9)
!         ---(9,9)---
    real(kind=8) :: kc(81)
!         -----(9,9) ----(9,9)
    real(kind=8) :: flexi(81), flex(81)
!          -----(6,6)
    real(kind=8) :: memb(36)
!                   -----(6,9)  -----(6,9)
    real(kind=8) :: mefl(54)
    real(kind=8) :: bsigth(24), enerth, carat3(25)
    real(kind=8) :: t2iu(4), t2ui(4), t1ve(9), qsi, eta
    aster_logical :: coupmf, indith
    integer(kind=8) :: i, jcoqu, jdepg, k
    real(kind=8) :: ctor, excent, zero
    real(kind=8) :: aire
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
!     ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                     jgano=jgano)
!
    zero = 0.0d0
    enerth = zero
!
    call jevech('PCACOQU', 'L', jcoqu)
    ctor = zr(jcoqu+3)
    excent = zr(jcoqu+4)
!
! --- ON NE CALCULE PAS ENCORE LA MATRICE DE RIGIDITE D'UN ELEMENT
! --- Q4G EXCENTRE, ON S'ARRETE EN ERREUR FATALE :
!     ------------------------------------------
    if (excent .ne. zero) then
        call utmess('F', 'ELEMENTS2_57')
    end if
!
    call r8inir(81, zero, kc, 1)
    call r8inir(81, zero, flex, 1)
    call r8inir(36, zero, memb, 1)
    call r8inir(54, zero, mefl, 1)
!
!     ----- CALCUL DES MATRICES DE RIGIDITE DU MATERIAU EN FLEXION,
!           MEMBRANE ET CISAILLEMENT INVERSEE --------------------------
    call dxmate('RIGI', df, dm, dmf, dc, &
                dci, dmc, dfc, nno, pgl, &
                multic, coupmf, t2iu, t2ui, t1ve)
!
!     ----- CALCUL DES GRANDEURS GEOMETRIQUES SUR LE TRIANGLE --------
    call gtria3(xyzl, carat3)
!
!     ------------------------------------------------------------------
!     CALCUL DE LA MATRICE DE RIGIDITE DE L'ELEMENT EN MEMBRANE
!     ------------------------------------------------------------------
!     ------ CALCUL DE LA MATRICE BM -----------------------------------
    call dxtbm(carat3(9), bm)
!     ------ CALCUL DU PRODUIT BMT.DM.BM -------------------------------
    call utbtab('ZERO', 3, 6, dm, bm, &
                xab1, memb)
    aire = carat3(8)
    do k = 1, 36
        memb(k) = memb(k)*aire
    end do
!
!     ------------------------------------------------------------------
!     CALCUL DES MATRICES DE RIGIDITE DE L'ELEMENT EN FLEXION ET
!     COUPLAGE MEMBRANE/FLEXION
!     ------------------------------------------------------------------
!
!     ------- CALCUL DE LA MATRICE BFB -------------------------------
    call dstbfb(carat3(9), bfb)
!
!     ------- CALCUL DU PRODUIT BFBT.DF.BFB --------------------------
    call utbtab('ZERO', 3, 9, df, bfb, &
                xab2, flex)
!
!        ---- CALCUL DE LA MATRICE BC ----------------------------------
    qsi = 1.d0/3.d0
    eta = qsi
    call t3gbc(xyzl, qsi, eta, bc)
!
!        ---- CALCUL DU PRODUIT BCT.DC.BC -----------------------------
    call utbtab('ZERO', 2, 9, dc, bc, &
                xab3, kc)
!
    do k = 1, 81
        flexi(k) = (flex(k)+kc(k))*aire
    end do
!
    if (option .eq. 'RIGI_MECA') then
        call dxtloc(flexi, memb, mefl, ctor, rig)
    else if (option .eq. 'EPOT_ELEM') then
        call jevech('PDEPLAR', 'L', jdepg)
        call utpvgl(3, 6, pgl, zr(jdepg), depl)
        call dxtloe(flex, memb, mefl, ctor, coupmf, &
                    depl, ener)
        call bsthpl(nomte, bsigth, indith)
        if (indith) then
            do i = 1, 18
                enerth = enerth+depl(i)*bsigth(i)
            end do
            ener(1) = ener(1)-enerth
        end if
    end if
!
end subroutine
