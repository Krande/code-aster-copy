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

subroutine te0430(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/cargri.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/nmgrib.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecael.h"
#include "asterfort/terefe.h"
#include "asterfort/utmess.h"
#include "asterfort/lteatt.h"
#include "asterfort/verift.h"
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE CHARGEMENT :
!                                  - CHAR_MECA_EPSI_R
!                                  - CHAR_MECA_PESA_R
!                                  - CHAR_MECA_TEMP_R
!                                  - FORC_NODA
!                                  - REFE_FORC_NODA
!                          POUR LES GRILLES MEMBRANES EXCENTREES OU NON
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: codres(2)
    character(len=4) :: fami
    character(len=16) :: nomres(2)
    integer(kind=8) :: nddl, nno, npg, i, kpg, n, ndim, nnos, jgano
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate, icontm, ivectu, iret, ier
    integer(kind=8) :: ipesa, iepsin, iadzi, iazk24, itemps
    real(kind=8) :: dff(2, 8), vff(8), b(6, 8), p(3, 6), jac, epsthe, epsref
    real(kind=8) :: dir11(3), densit, pgl(3, 3), distn, vecn(3)
    real(kind=8) :: sig, rho(1), valres(2), b_max_rot
    aster_logical :: lexc
    character(len=8) :: nompar(4)
    real(kind=8) :: valpar(4)
    real(kind=8) :: xgau, ygau, zgau, exx
!
! - BOOLEEN POUR LES GRILLES EXCENTREES
!
    lexc = (lteatt('MODELI', 'GRC'))
!
! - FONCTIONS DE FORMES ET POINTS DE GAUSS
!
    fami = 'RIGI'
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! - PARAMETRES EN ENTREE
!
    call jevech('PGEOMER', 'L', igeom)
!
    if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', icontm)
!
    else if (option .eq. 'REFE_FORC_NODA') then
        call jevech('PMATERC', 'L', imate)
!
    else if (option .eq. 'CHAR_MECA_EPSI_R') then
        call jevech('PMATERC', 'L', imate)
        call jevech('PEPSINR', 'L', iepsin)
    else if (option .eq. 'CHAR_MECA_EPSI_F') then
        call jevech('PMATERC', 'L', imate)
        call jevech('PEPSINF', 'L', iepsin)
        call jevech('PINSTR', 'L', itemps)
!
    else if (option .eq. 'CHAR_MECA_PESA_R') then
        call jevech('PMATERC', 'L', imate)
        call jevech('PPESANR', 'L', ipesa)
!
    else if (option .eq. 'CHAR_MECA_TEMP_R') then
        call jevech('PMATERC', 'L', imate)
!
    end if
!
! - PARAMETRES EN SORTIE
!
    call jevech('PVECTUR', 'E', ivectu)
!
! - LECTURE DES CARACTERISTIQUES DE GRILLE ET
!   CALCUL DE LA DIRECTION D'ARMATURE
!
    call cargri(lexc, densit, distn, dir11)
!
! - SI EXCENTREE : RECUPERATION DE LA NORMALE ET DE L'EXCENTREMENT
!
    if (lexc) then
!
        if (nomte .eq. 'MEGCTR3') then
            call dxtpgl(zr(igeom), pgl)
        else if (nomte .eq. 'MEGCQU4') then
            call dxqpgl(zr(igeom), pgl)
        end if
!
        do i = 1, 3
            vecn(i) = distn*pgl(3, i)
        end do
        nddl = 6
!
    else
        nddl = 3
    end if
!
! - DEBUT DE LA BOUCLE SUR LES POINTS DE GAUSS
!
    do kpg = 1, npg
!
! - MISE SOUS FORME DE TABLEAU DES VALEURS DES FONCTIONS DE FORME
!   ET DES DERIVEES DE FONCTION DE FORME
!
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do
!
! - CALCUL DE LA MATRICE "B" : DEPL NODAL --> EPS11 ET DU JACOBIEN
!
        call nmgrib(nno, zr(igeom), dff, dir11, lexc, &
                    vecn, b, jac, p)
!
! - BRANCHEMENT DES DIFFERENTES OPTIONS
!
        if ((option .eq. 'FORC_NODA') .or. (option .eq. 'CHAR_MECA_TEMP_R') .or. &
            (option(1:15) .eq. 'CHAR_MECA_EPSI_')) then
!
! - FORC_NODA : IL SUFFIT DE RECOPIER SIGMA
!
            if (option .eq. 'FORC_NODA') then
                sig = zr(icontm+kpg-1)
!
! - CHAR_MECA_EPSI_* : SIG = E*EPSIN
!
            else if (option .eq. 'CHAR_MECA_EPSI_R') then
                nomres(1) = 'E'
                call rcvalb(fami, kpg, 1, '+', zi(imate), &
                            ' ', 'ELAS', 0, ' ', [0.d0], &
                            1, nomres, valres, codres, 1)
                sig = valres(1)*zr(iepsin+kpg-1)

            else if (option .eq. 'CHAR_MECA_EPSI_F') then
                nomres(1) = 'E'
                call rcvalb(fami, kpg, 1, '+', zi(imate), &
                            ' ', 'ELAS', 0, ' ', [0.d0], &
                            1, nomres, valres, codres, 1)

                nompar(1) = 'X'
                nompar(2) = 'Y'
                nompar(3) = 'Z'
                nompar(4) = 'INST'
                valpar(4) = zr(itemps)
                xgau = 0.d0
                ygau = 0.d0
                zgau = 0.d0
!
                do i = 1, nno
                    xgau = xgau+zr(ivf-1+i+nno*(kpg-1))*zr(igeom-1+1+3*(i-1))
                    ygau = ygau+zr(ivf-1+i+nno*(kpg-1))*zr(igeom-1+2+3*(i-1))
                    zgau = zgau+zr(ivf-1+i+nno*(kpg-1))*zr(igeom-1+3+3*(i-1))
                end do
!
                valpar(1) = xgau
                valpar(2) = ygau
                valpar(3) = zgau

                call fointe('FM', zk8(iepsin), 4, nompar, valpar, exx, ier)

                sig = valres(1)*exx
!
! - CHAR_MECA_TEMP_R : SIG = E*EPSTHE
!
            else if (option .eq. 'CHAR_MECA_TEMP_R') then
                call verift(fami, kpg, 1, '+', zi(imate), &
                            iret_=iret, epsth_=epsthe)
                if (iret .ne. 0) then
                    call tecael(iadzi, iazk24)
                    call utmess('S', 'CALCULEL2_81', sk=zk24(iazk24-1+3))
                end if
                nomres(1) = 'E'
                call rcvalb(fami, kpg, 1, '+', zi(imate), &
                            ' ', 'ELAS', 0, ' ', [0.d0], &
                            1, nomres, valres, codres, 1)
                sig = valres(1)*epsthe
            end if
!
            do n = 1, nno
                do i = 1, nddl
                    zr(ivectu+(n-1)*nddl+i-1) = zr(ivectu+(n-1)*nddl+i- &
                                                   1)+b(i, n)*sig*zr(ipoids+kpg-1)*jac*densit
                end do
            end do
!
! - REFE_FORC_NODA : ON CALCULE DES FORCES DE REFERENCE
!      (N'EST VALABLE QUE POUR LES GRILLES MEMBRANES)
!
        else if (option .eq. 'REFE_FORC_NODA') then
!
            call terefe('EPSI_REFE', 'GRILLE', epsref)
            if (epsref .eq. r8vide()) then
                ASSERT(.false.)
            end if
!
            nomres(1) = 'E'
            call rcvalb(fami, kpg, 1, '+', zi(imate), &
                        ' ', 'ELAS', 0, ' ', [0.d0], &
                        1, nomres, valres, codres, 1)
            sig = valres(1)*epsref
!
            do n = 1, nno
                do i = 1, 3
                    zr(ivectu+(n-1)*nddl+i-1) = zr(ivectu+(n-1)*nddl+i-1)+ &
                                                sig*sqrt(abs(jac))*densit/npg
                end do
                b_max_rot = 0.d0
                do i = 4, nddl
                    if (abs(b(i, n)) .gt. b_max_rot) b_max_rot = abs(b(i, n))
                end do
                do i = 4, nddl
                    zr(ivectu+(n-1)*nddl+i-1) = zr(ivectu+(n-1)*nddl+i-1)+ &
                                                b_max_rot*sig*sqrt(abs(jac))*densit/npg
                end do
            end do
!
! - CHAR_MECA_PESA_R
!
        else if (option .eq. 'CHAR_MECA_PESA_R') then
            call rcvalb(fami, kpg, 1, '+', zi(imate), &
                        ' ', 'ELAS', 0, ' ', [0.d0], &
                        1, 'RHO', rho, codres, 1)
            do n = 1, nno
                do i = 1, 3
                    zr(ivectu+(n-1)*nddl+i-1) = zr(ivectu+(n-1)*nddl+i-1)+ &
                                                rho(1)*zr(ipoids+kpg-1)*zr(ipesa)*zr(ipesa+i)* &
                                                vff(n)*densit*jac
                end do
            end do
        end if
!
! - FIN DE LA BOUCLE SUR LES POINTS DE GAUSS
    end do
!
end subroutine
