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

subroutine mbgchg(option, fami, nddl, nno, ncomp, kpg, imate, jvSief, &
                  ipoids, ipesa, igeom, ivectu, vff, dff, h, alpha, beta, preten)
!
    implicit none
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/mbpk2c.h"
#include "asterfort/mbvfie.h"
#include "asterfort/rcvalb.h"
#include "asterfort/subaco.h"
#include "asterfort/sumetr.h"
#include "asterfort/subacv.h"
!
    character(len=16) :: option
    character(len=4) :: fami
    integer(kind=8) :: nddl, nno, ncomp
    integer(kind=8) :: kpg
    integer(kind=8) :: ipoids, igeom, jvSief, imate, ipesa
    integer(kind=8) :: ivectu
    real(kind=8) :: vff(nno), dff(2, nno), h, preten, alpha, beta
! ----------------------------------------------------------------------
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE DE CHARGEMENT :
!                                  - FORC_NODA
!                                  - CHAR_MECA_PESA_R
!                          POUR LES MEMBRANES EN GRANDES DEFORMATIONS
! ----------------------------------------------------------------------
! IN  OPTION       OPTION DE CALCUL
! IN  FAMI         NOM DE LA FAMILLE DE POINTS DE GAUSS :
!                  'RIGI','MASS',..
! IN  NNO          NOMBRE DE NOEUDS
! IN  KPG          INCREMENT SUR LA BOUCLE DES PTS DE GAUSS
! IN  NPG          NOMBRE DE POINT DE GAUSS
! IN  IPOIDS       ADRESSE DANS ZR DU TABLEAU POIDS
! IN  IGEOM        ADRESSE DANS ZR DU TABLEAU PGEOMER
! IN  IMATE        ADRESSE DANS ZI DU TABLEAU PMATERC
! IN  IVECTU       ADRESSE DANS ZR DU TABLEAU PVECTUR
! IN  DFF          DERIVEE DES F. DE FORME
! IN  H            EPAISSEUR DE LA MEMBRANE
! IN ALPHA, BETA   ANGLES DEF. LA BASE DE L'ÉCRITURE DES CONTRAINTES
! IN  PRETEN       PRECONTRAINTES
!
! OUT ***          ***
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, n, c
    integer(kind=8) :: jvDisp
    integer(kind=8) :: codres(2)
    real(kind=8) :: posdef(3*nno)
    real(kind=8) :: covaini(3, 3), metrini(2, 2), jacini, cnvaini(3, 2), aini(2, 2)
    real(kind=8) :: covadef(3, 3), metrdef(2, 2), jacdef, cnvadef(3, 2), adef(2, 2)
    real(kind=8) :: sigpk2(2, 2)
    real(kind=8) :: vecfie(3*nno), sighca(3), sigout(3)
    real(kind=8) :: rho(1)
!
! - FORC_NODA
    if (option .eq. 'FORC_NODA') then
        call jevech('PDEPLAR', 'L', jvDisp)

!
! ---   CALCUL DES COORDONNEES COVARIANTES ET CONTRAVARIANTES DE LA SURFACE INITIALE
!
        call subaco(nno, dff, zr(igeom), covaini)
        call sumetr(covaini, metrini, jacini)
        call subacv(covaini, metrini, jacini, cnvaini, aini)

! ---   CALCUL DES COORDONNEES COVARIANTES ET CONTRAVARIANTES DE LA SURFACE DEFORMEE
!
        do n = 1, 3*nno
            posdef(n) = zr(igeom+n-1)+zr(jvDisp+n-1)
        end do

        call subaco(nno, dff, posdef, covadef)
        call sumetr(covadef, metrdef, jacdef)
        call subacv(covadef, metrdef, jacdef, cnvadef, adef)

! ---   ON EXTRAIT LES CONTRAINTES DE CAUCHY INTEGREES QUE L'ON TRANSFORME EN PKII (NON INTEGREES)
!
        do c = 1, ncomp
            sighca(c) = zr(jvSief+(kpg-1)*ncomp+c-1)
        end do

        call mbpk2c(1, alpha, beta, h, covaini, jacini, jacdef, sighca, sigout)

        sigpk2(1, 1) = sigout(1)
        sigpk2(2, 2) = sigout(2)
        sigpk2(1, 2) = sigout(3)
        sigpk2(2, 1) = sigout(3)

! ---   SI LA NORME EUCLIDIENNE DE SIGPK2 EST NULLE, ON APPLIQUE DES PRECONTRAINTES
        if (sqrt(sigpk2(1, 1)**2+2*sigpk2(1, 2)**2+sigpk2(2, 2)**2) .lt. 1.d-6) then
            sigpk2(1, 1) = preten
            sigpk2(2, 2) = preten
        end if

! ---   CALCUL DU VECTEUR FORCE INTERNE ELEMENTAIRE
        call mbvfie(nno, kpg, dff, sigpk2, ipoids, h, covadef, vecfie)

! ---   RANGEMENT DES RESULTATS
!
        do n = 1, 3*nno
            zr(ivectu+n-1) = zr(ivectu+n-1)+vecfie(n)*jacini
        end do

! - CHAR_MECA_PESA_R
    else if (option .eq. 'CHAR_MECA_PESA_R') then

        call subaco(nno, dff, zr(igeom), covaini)
        call sumetr(covaini, metrini, jacini)

        call rcvalb(fami, kpg, 1, '+', zi(imate), &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    1, 'RHO', rho, codres, 1)
        do n = 1, nno
            do i = 1, nddl
                zr(ivectu+(n-1)*nddl+i-1) = zr(ivectu+(n-1)*nddl+i-1)+ &
                                            rho(1)*zr(ipesa)*zr(ipesa+i)* &
                                            vff(n)*zr(ipoids+kpg-1)*h*jacini
            end do
        end do

    end if

end subroutine
