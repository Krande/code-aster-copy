! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine te0249(option, nomte)
!
!     BUT: CALCUL DES MATRICES TANGENTES ELEMENTAIRES EN THERMIQUE
!          CORRESPONDANT AU TERME D'ECHANGE
!          SUR DES FACES D'ELEMENTS ISOPARAMETRIQUES 2D
!
!          OPTION : 'MTAN_THER_COEF_R'
!          OPTION : 'MTAN_THER_RAYO_R'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!   -------------------------------------------------------------------
!     ASTER INFORMATIONS:
!       04/04/02 (OB): CORRECTION BUG CALCUL TPG EN LUMPE
!       + MODIFS FORMELLES: IMPLICIT NONE, LAXI, LCOEF, ...
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
! PARAMETRES D'APPEL
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8t0.h"
#include "asterfort/assert.h"
#include "asterfort/connec.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/vff2dn.h"
    character(len=16) :: option, nomte
!
!
    real(kind=8) :: poids, r, nx, ny, theta, mrigt(9, 9), coorse(18), hech
    real(kind=8) :: sigma, epsil, tpg, tz0
    integer :: nno, nnos, jgano, ndim, kp, npg, ipoids, ivf, idfde, igeom
    integer :: c(6, 9), imattt, i, j, ij, l, li, lj, iray, itemp, ise, nse
    integer :: nnop2, iech, itemps, ibid
    aster_logical :: laxi, lcoef
    character(len=8) :: elrefe, alias8
!
!====
! 1.1 PREALABLES: RECUPERATION ADRESSES FONCTIONS DE FORMES...
!====
    tz0 = r8t0()
    call elref1(elrefe)
!
    if (lteatt('LUMPE','OUI')) then
        call teattr('S', 'ALIAS8', alias8, ibid)
        if (alias8(6:8) .eq. 'SE3') elrefe='SE2'
    endif
!
    call elrefe_info(elrefe=elrefe, fami='RIGI', ndim=ndim, nno=nno, nnos=nnos,&
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! INITS.
    if (option(11:14) .eq. 'COEF') then
        lcoef = .true.
        call jevech('PCOEFHR', 'L', iech)
        hech = zr(iech)
    else if (option(11:14).eq.'RAYO') then
        lcoef = .false.
        call jevech('PRAYONR', 'L', iray)
        call jevech('PTEMPEI', 'L', itemp)
        sigma = zr(iray)
        epsil = zr(iray+1)
    else
!C OPTION DE CALCUL INVALIDE
        ASSERT(.false.)
    endif
    if (lteatt('AXIS','OUI')) then
        laxi = .true.
    else
        laxi = .false.
    endif
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PMATTTR', 'E', imattt)
!
    theta = zr(itemps+2)
!
    call connec(nomte, nse, nnop2, c)
!
    do i = 1, nnop2
        do j = 1, nnop2
            mrigt(i,j) = 0.d0
        end do
    end do
!
! --- CALCUL ISO-P2 : BOUCLE SUR LES SOUS-ELEMENTS -------
!
    do ise = 1, nse
!
        do i = 1, nno
            do j = 1, 2
                coorse(2* (i-1)+j) = zr(igeom-1+2* (c(ise,i)-1)+j)
            end do
        end do
!
        do kp = 1, npg
            call vff2dn(ndim, nno, kp, ipoids, idfde,&
                        coorse, nx, ny, poids)
            if (laxi) then
                r = 0.d0
                do i = 1, nno
                    l = (kp-1)*nno + i
                    r = r + coorse(2* (i-1)+1)*zr(ivf+l-1)
                end do
                poids = poids*r
            endif
            ij = imattt - 1
            if (lcoef) then
                do i = 1, nno
                    li = ivf + (kp-1)*nno + i - 1
                    do j = 1, i
                        lj = ivf + (kp-1)*nno + j - 1
                        ij = ij + 1
                        mrigt(c(ise,i),c(ise,j)) = mrigt(&
                                                   c(ise, i),&
                                                   c( ise, j)) + poids*theta*zr(li)*zr(lj&
                                                   )* hech
                    end do
                end do
            else
                tpg = 0.d0
                do i = 1, nno
                    l = (kp-1)*nno + i
                    tpg = tpg + zr(itemp-1+c(ise,i))*zr(ivf+l-1)
                end do
                do i = 1, nno
                    li = ivf + (kp-1)*nno + i - 1
                    do j = 1, i
                        lj = ivf + (kp-1)*nno + j - 1
                        ij = ij + 1
                        mrigt(c(ise,i),c(ise,j)) = mrigt(&
                                                   c(ise, i),&
                                                   c( ise, j)) + poids*theta*zr(li)*zr(lj)* 4.d0*&
                                                   & sigma*epsil* (tpg+tz0&
                                                   )**3
                    end do
                end do
            endif
        end do
    end do
!
! MISE SOUS FORME DE VECTEUR
!
    ij = imattt - 1
    do i = 1, nnop2
        do j = 1, i
            ij = ij + 1
            zr(ij) = mrigt(i,j)
        end do
    end do
end subroutine
