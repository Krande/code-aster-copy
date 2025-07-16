! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! Copyright (C) 2022 - Anthony McDonald - anthony.mcdonald .at. wave-venture.com
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
subroutine te0435(option, nomte)
!
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/codere.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/jevecd.h"
#include "asterfort/mbxnlr.h"
#include "asterfort/mbgnlr.h"
#include "asterfort/nmpr3d_vect.h"
#include "asterfort/nmpr3d_matr.h"
#include "asterfort/nmprmb_matr.h"
#include "asterc/r8prem.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: GRILLE_MEMBRANE / GRILLE_EXCENTRE
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          RIGI_MECA
!          RIGI_MECA_PRSU_R
!          RIGI_MECA_PRSU_F
!          CHAR_MECA_PRSU_R
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: mxnpg, mxvect, mxmatr
    parameter(mxnpg=27, mxvect=3*9, mxmatr=3*9*3*9)
    integer(kind=8) :: mxpara
    parameter(mxpara=7)
!
    character(len=8) :: nompar(mxpara)
    real(kind=8) :: valpar(mxpara)
    integer(kind=8) :: ier
    real(kind=8) :: x, y, z, xf, yf, zf
!
    character(len=4) :: fami
    integer(kind=8) :: nddl, nno, nnos, npg, ndim, ncomp, nvari
    integer(kind=8) :: n, kpg, iret, cod(9)
    integer(kind=8) :: ipoids, ivf, idfde, jgano, jtab(7)
    integer(kind=8) :: igeom, icacoq, imate, icompo, icarcr
    integer(kind=8) :: iinstm, iinstp, icontm, ideplm, ideplp, ivarim, ivarix
    integer(kind=8) :: ivectu, icontp, ivarip, jcret, imatuu, imatun, icontx, i_pres
    integer(kind=8) :: i_temp, kdec, i, j, k, iddl, ino, ndofbynode
    real(kind=8) :: pres, pres_point(mxnpg)
    real(kind=8) :: matr(mxmatr), geom_reac(mxvect)
    real(kind=8) :: dff(2, 9), alpha, beta, h, preten
    aster_logical :: lNonLine, lLine
    aster_logical :: pttdef, grddef
    aster_logical :: lVect, lMatr, lVari, lSigm
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: defo_comp, rela_comp
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    lNonLine = (option(1:9) .eq. 'FULL_MECA') .or. (option(1:9) .eq. 'RAPH_MECA') .or. &
               (option(1:10) .eq. 'RIGI_MECA_') .and. (option(1:15) .ne. 'RIGI_MECA_PRSU_')
    lLine = option .eq. 'RIGI_MECA'
    cod = 0
!
! - NOMBRE DE COMPOSANTES DES TENSEURS
!
    ncomp = 3
    nddl = 3
!
! - FONCTIONS DE FORME ET POINTS DE GAUSS
!
    fami = 'RIGI'
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    ndofbynode = ndim+1
!
! - Get input fields
!
    lVect = ASTER_FALSE
    lVari = ASTER_FALSE
    lSigm = ASTER_FALSE
    lMatr = ASTER_FALSE
    call jevech('PGEOMER', 'L', igeom)
    if ((option(1:15) .ne. 'RIGI_MECA_PRSU_') .and. (option .ne. 'CHAR_MECA_PRSU_F')) then
        call jevech('PCACOQU', 'L', icacoq)
        call jevech('PMATERC', 'L', imate)
    end if
    if (lNonLine) then
        call jevech('PCOMPOR', 'L', icompo, vk16=compor)
        call jevech('PCARCRI', 'L', icarcr)
        call jevech('PINSTMR', 'L', iinstm)
        call jevech('PINSTPR', 'L', iinstp)
        call jevech('PCONTMR', 'L', icontm)
        call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                    itab=jtab)
        nvari = max(jtab(6), 1)*jtab(7)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PVARIMP', 'L', ivarix)
    end if
    if (option .ne. 'RIGI_MECA') then
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
    end if
!
    if (lNonLine) then
! ----- Select objects to construct from option name
        call behaviourOption(option, compor, lMatr, lVect, lVari, &
                             lSigm)
! ----- Properties of behaviour
        rela_comp = compor(RELA_NAME)
        defo_comp = compor(DEFO)
        pttdef = (defo_comp .eq. 'PETIT')
        grddef = (defo_comp .eq. 'GROT_GDEP')
        if (.not. pttdef .and. .not. grddef) then
            call utmess('F', 'MEMBRANE_2', sk=defo_comp)
        end if
    end if
    if (lLine) then
        pttdef = ASTER_TRUE
        grddef = ASTER_FALSE
        lMatr = ASTER_TRUE
    end if
!
! - PARAMETRES NECESSAIRE AU CALCUL DE LA MATRICE DE RIGITE POUR PRESSION SUIVEUSE
    if (option .eq. 'RIGI_MECA_PRSU_R') then
        call jevecd('PPRESSR', i_pres, 0.d0)
    end if
!
    if (option(10:16) .eq. '_PRSU_F') then
        call jevech('PPRESSF', 'L', i_pres)
        call jevech('PINSTR', 'L', i_temp)
        valpar(4) = zr(i_temp)
        nompar(4) = 'INST'
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(5) = 'XF'
        nompar(6) = 'YF'
        nompar(7) = 'ZF'
        do iddl = 1, nddl*nno
            geom_reac(iddl) = zr(igeom+iddl-1)+zr(ideplm+iddl-1)+zr(ideplp+iddl-1)
        end do
    end if
!
! - Get output fields
!
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        b_n = to_blas_int(npg*nvari)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
    end if
    if (lMatr) then
        call jevech('PMATUUR', 'E', imatuu)
    end if
    if (option(1:15) .eq. 'RIGI_MECA_PRSU_') then
        call jevech('PMATUNS', 'E', imatun)
    end if
    if (option .eq. 'RIGI_MECA_IMPLEX') then
        call jevech('PCONTXR', 'E', icontx)
        b_n = to_blas_int(npg*ncomp)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(icontm), b_incx, zr(icontx), b_incy)
    end if
!
!
! -----------------------------------------------------------------
! ---          IMPORTATION DES PARAMETRES MATERIAU              ---
! -----------------------------------------------------------------
!
! - DIRECTION DE REFERENCE POUR UN COMPORTEMENT ANISOTROPE
! - EPAISSEUR
! - PRECONTRAINTES
!
    if ((option(1:15) .ne. 'RIGI_MECA_PRSU_') .and. (option .ne. 'CHAR_MECA_PRSU_F')) then
        alpha = zr(icacoq+1)*r8dgrd()
        beta = zr(icacoq+2)*r8dgrd()
        h = zr(icacoq)
! ---   On empeche une epaisseur nulle ou negative
        if (h .lt. r8prem()) then
            call utmess('F', 'MEMBRANE_1')
        end if
        preten = zr(icacoq+3)/h
    end if
!
    do kpg = 1, npg
        do n = 1, nno
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do
        if (option .eq. 'RIGI_MECA_PRSU_R') then
            call nmprmb_matr(nno, npg, kpg, zr(ipoids+kpg), zr(ivf), &
                             dff, igeom, ideplm, ideplp, i_pres, &
                             imatun)
        elseif ((option .eq. 'RIGI_MECA_PRSU_F') .or. &
                (option .eq. 'CHAR_MECA_PRSU_F')) then
            kdec = (kpg-1)*nno
            x = 0.d0
            y = 0.d0
            z = 0.d0
            xf = 0.d0
            yf = 0.d0
            zf = 0.d0
            do ino = 1, nno
                x = x+zr(igeom+3*(ino-1)+1-1)*zr(ivf+kdec+ino-1)
                y = y+zr(igeom+3*(ino-1)+2-1)*zr(ivf+kdec+ino-1)
                z = z+zr(igeom+3*(ino-1)+3-1)*zr(ivf+kdec+ino-1)
                xf = xf+geom_reac(3*(ino-1)+1)*zr(ivf+kdec+ino-1)
                yf = yf+geom_reac(3*(ino-1)+2)*zr(ivf+kdec+ino-1)
                zf = zf+geom_reac(3*(ino-1)+3)*zr(ivf+kdec+ino-1)
            end do
            valpar(1) = x
            valpar(2) = y
            valpar(3) = z
            valpar(5) = xf
            valpar(6) = yf
            valpar(7) = zf
            call fointe('FM', zk8(i_pres), mxpara, nompar, valpar, &
                        pres, ier)
            pres_point(kpg) = pres
        else
            if (pttdef) then
                call mbxnlr(option, fami, nddl, nno, ncomp, &
                            kpg, ipoids, igeom, imate, ideplm, &
                            ideplp, ivectu, icontp, imatuu, dff, &
                            alpha, beta, lVect, lMatr)
            else if (grddef) then
                if (rela_comp(1:14) .eq. 'ELAS_MEMBRANE_') then
                    if ((abs(alpha) .gt. r8prem()) .or. (abs(beta) .gt. r8prem())) then
                        call utmess('A', 'MEMBRANE_6')
                    end if
                    call mbgnlr(lVect, lMatr, nno, ncomp, imate, &
                                icompo, dff, alpha, beta, h, &
                                preten, igeom, ideplm, ideplp, kpg, &
                                fami, ipoids, icontp, ivectu, imatuu)
                else
                    call utmess('F', 'MEMBRANE_3')
                end if
            end if
        end if
    end do
!
!
! - Second member
!
    if (option .eq. 'CHAR_MECA_PRSU_F') then
        call jevech('PVECTUR', 'E', ivectu)
        call nmpr3d_vect(nno, npg, ndofbynode, zr(ipoids), zr(ivf), &
                         zr(idfde), geom_reac, pres_point, zr(ivectu))
!
! - Tangent matrix
!
    else if (option .eq. 'RIGI_MECA_PRSU_F') then
        call nmpr3d_matr(nno, npg, zr(ipoids), zr(ivf), zr(idfde), &
                         geom_reac, pres_point, matr)
        call jevech('PMATUNS', 'E', imatun)
        k = 0
        do i = 1, nddl*nno
            do j = 1, nddl*nno
                k = k+1
                zr(imatun-1+k) = matr((j-1)*nddl*nno+i)
            end do
        end do
        ASSERT(k .eq. nddl*nno*nddl*nno)
    end if
!
    if (lSigm) then
        call codere(cod, npg, zi(jcret))
    end if
!
end subroutine
