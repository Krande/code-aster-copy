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

subroutine q4gmas(xyzl, option, pgl, mas, ener)
    implicit none
#include "jeveux.h"
#include "asterfort/dialum.h"
#include "asterfort/dxqloc.h"
#include "asterfort/dxqloe.h"
#include "asterfort/dxroep.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gquad4.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
#include "asterfort/q4gbc.h"
#include "asterfort/q4gniw.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
    real(kind=8) :: xyzl(3, *), pgl(*), mas(*), ener(*)
    character(len=16) :: option
!     MATRICE MASSE DE L'ELEMENT DE PLAQUE Q4GAMMA (W LINEAIRE)
!     ------------------------------------------------------------------
!     IN  XYZL   : COORDONNEES LOCALES DES QUATRE NOEUDS
!     IN  OPTION : OPTION RIGI_MECA OU EPOT_ELEM
!     IN  PGL    : MATRICE DE PASSAGE GLOBAL/LOCAL
!     OUT MAS    : MATRICE DE RIGIDITE
!     OUT ENER   : TERMES POUR ENER_CIN (ECIN_ELEM)
!     ------------------------------------------------------------------
    integer(kind=8) :: i, j, k, int, jcoqu, jdepg, ii(8), jj(8), ll(16)
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: jvitg, iret
    real(kind=8) :: flex(12, 12), bc(2, 12)
    real(kind=8) :: memb(8, 8), amemb(64), mefl(8, 12)
    real(kind=8) :: wq4(12), depl(24), masloc(300), masglo(300), vite(24)
    real(kind=8) :: rho, epais, roe, ctor, excent, detj, wgt, zero, coefm
    real(kind=8) :: caraq4(25), jacob(5), qsi, eta
    character(len=3) :: stopz
!     ------------------------------------------------------------------
    data(ii(k), k=1, 8)/1, 10, 19, 28, 37, 46, 55, 64/
    data(jj(k), k=1, 8)/5, 14, 23, 32, 33, 42, 51, 60/
    data(ll(k), k=1, 16)/3, 7, 12, 16, 17, 21, 26, 30, 35, 39, 44, 48, 49, 53, 58, 62/
!     ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, &
                     jdfd2=idfd2, jgano=jgano)
!
    zero = 0.0d0
!
    call dxroep(rho, epais)
    roe = rho*epais
!
    call jevech('PCACOQU', 'L', jcoqu)
    ctor = zr(jcoqu+3)
    excent = zr(jcoqu+4)
!
! --- ON NE CALCULE PAS ENCORE LA MATRICE DE MASSE D'UN ELEMENT
! --- DE PLAQUE EXCENTRE, ON S'ARRETE EN ERREUR FATALE :
!     ------------------------------------------------
!
!     ----- CALCUL DES GRANDEURS GEOMETRIQUES SUR LE QUADRANGLE --------
    call gquad4(xyzl, caraq4)
!
    call r8inir(144, zero, flex, 1)
    call r8inir(96, zero, mefl, 1)
!
    do int = 1, npg
!
        qsi = zr(icoopg-1+ndim*(int-1)+1)
        eta = zr(icoopg-1+ndim*(int-1)+2)
!
!        ----- CALCUL DU JACOBIEN SUR LE QUADRANGLE -------------------
        call jquad4(xyzl, qsi, eta, jacob)
!
!        ---- CALCUL DE LA MATRICE BC ----------------------------------
        call q4gbc(qsi, eta, jacob(2), caraq4, bc)
!
!        ----- CALCUL DES FONCTIONS D'INTERPOLATION EN FLEXION ---------
        call q4gniw(qsi, eta, caraq4, bc, wq4)
!
        detj = jacob(1)
        wgt = zr(ipoids+int-1)*detj*roe
!
        do i = 1, 12
            do j = 1, 12
                flex(i, j) = flex(i, j)+wq4(i)*wq4(j)*wgt
            end do
        end do
    end do
!
!     ----- CALCUL DE LA MATRICE MASSE EN MEMBRANE ---------------------
    coefm = caraq4(21)*roe/9.d0
    amemb(:) = 0.d0
    do k = 1, 8
        amemb(ii(k)) = 1.d0
        amemb(jj(k)) = 0.25d0
    end do
    do k = 1, 16
        amemb(ll(k)) = 0.5d0
    end do
    do j = 1, 8
        do i = 1, 8
            memb(i, j) = coefm*amemb((j-1)*8+i)
        end do
    end do
!
    if ((option .eq. 'MASS_MECA') .or. (option .eq. 'M_GAMMA')) then
        call dxqloc(flex, memb, mefl, ctor, mas)
!
    else if (option .eq. 'MASS_MECA_DIAG' .or.&
 &         option .eq. 'MASS_MECA_EXPLI') then
        call dxqloc(flex, memb, mefl, ctor, masloc)
        wgt = caraq4(21)*roe
        call utpslg(4, 6, pgl, masloc, masglo)
        call dialum(4, 6, 24, wgt, masglo, &
                    mas)
!
    else if (option .eq. 'ECIN_ELEM') then
        stopz = 'ONO'
! IRET NE PEUT VALOIR QUE 0 (TOUT VA BIEN) OU 2 (CHAMP NON FOURNI)
        call tecach(stopz, 'PDEPLAR', 'L', iret, iad=jdepg)
        if (iret .eq. 0) then
            call utpvgl(4, 6, pgl, zr(jdepg), depl)
            call dxqloe(flex, memb, mefl, ctor, .false._1, &
                        depl, ener)
        else
            call tecach(stopz, 'PVITESR', 'L', iret, iad=jvitg)
            if (iret .eq. 0) then
                call utpvgl(4, 6, pgl, zr(jvitg), vite)
                call dxqloe(flex, memb, mefl, ctor, .false._1, &
                            vite, ener)
            else
                call utmess('F', 'ELEMENTS2_1', sk=option)
            end if
        end if
    end if
!
end subroutine
