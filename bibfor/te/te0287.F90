! --------------------------------------------------------------------
! Copyright (C) 2019 Christophe Durand - www.code-aster.org
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine te0287(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/matrot.h"
#include "asterfort/ntfcma.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcdiff.h"
#include "asterfort/rcfode.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/utrcyl.h"
#include "asterfort/uttgel.h"
!
    character(len=16) :: nomte, option
! ----------------------------------------------------------------------
!
!    - FONCTION REALISEE:  CALCUL DES MATRICES TANGENTES ELEMENTAIRES
!                          OPTION : 'MASS_THER_TANG'
!                          ELEMENTS 3D ISO PARAMETRIQUES LUMPES
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! THERMIQUE NON LINEAIRE
!
    integer :: nbres
    parameter(nbres=3)
    integer :: icodre(nbres)
    character(len=2) :: typgeo
    character(len=32) :: phenom
    real(kind=8) :: rhocp, tpgi
    real(kind=8) :: dfdx(27), dfdy(27), dfdz(27), poids, r8bid
    integer :: igeom, imate
    integer :: jgano, nno, kp, i, j, ij, l, imattt, ifon(6)
    integer :: icomp, itempi, nnos, ndim
    integer :: npg2, ipoid2, ivf2, idfde2
    aster_logical :: aniso
!
!====
! 1.1 PREALABLES: RECUPERATION ADRESSES FONCTIONS DE FORMES...
!====
    call uttgel(nomte, typgeo)
    if ((lteatt('LUMPE', 'OUI')) .and. (typgeo .ne. 'PY')) then
        call elrefe_info(fami='NOEU', ndim=ndim, nno=nno, nnos=nnos, npg=npg2, &
                         jpoids=ipoid2, jvf=ivf2, jdfde=idfde2, jgano=jgano)
    else
        call elrefe_info(fami='MASS', ndim=ndim, nno=nno, nnos=nnos, npg=npg2, &
                         jpoids=ipoid2, jvf=ivf2, jdfde=idfde2, jgano=jgano)
    end if
!
!====
! 1.2 PREALABLES LIES AUX RECHERCHES DE DONNEES GENERALES
!====
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PTEMPEI', 'L', itempi)
    call jevech('PCOMPOR', 'L', icomp)
    call jevech('PMATTTR', 'E', imattt)
!
    if (zk16(icomp) (1:5) .eq. 'THER_') then
!====
! --- THERMIQUE
!====
        call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
        aniso = .false.
        if (phenom(1:12) .eq. 'THER_NL_ORTH') then
            aniso = .true.
        end if
        call ntfcma(zk16(icomp), zi(imate), aniso, ifon)
!
! ---   CALCUL DU DEUXIEME TERME
!
        do kp = 1, npg2
            l = (kp-1)*nno
            call dfdm3d(nno, kp, ipoid2, idfde2, zr(igeom), &
                        poids, dfdx, dfdy, dfdz)
!
! ---       EVALUATION DE LA CAPACITE CALORIFIQUE
! ---       PAS DE TRAITEMENT POUR L ORTHOTROPIE (RHOCP NON CONCERNE)
!
            tpgi = 0.d0
            do i = 1, nno
                tpgi = tpgi+zr(itempi+i-1)*zr(ivf2+l+i-1)
            end do
            call rcfode(ifon(1), tpgi, r8bid, rhocp)
!
! ---       CALCUL DE LA DEUXIEME COMPOSANTE DU TERME ELEMENTAIRE
!
            do i = 1, nno
                do j = 1, i
                    ij = (i-1)*i/2+j
                    zr(imattt+ij-1) = zr(imattt+ij-1)+poids*rhocp*zr(ivf2+l+i-1)*zr(ivf2+l+j&
                                      &-1)
                end do
            end do
        end do
!
    else if (zk16(icomp) (1:5) .eq. 'SECH_') then
!====
! --- SECHAGE
!====
        do kp = 1, npg2
            l = (kp-1)*nno
            call dfdm3d(nno, kp, ipoid2, idfde2, zr(igeom), &
                        poids, dfdx, dfdy, dfdz)
            do i = 1, nno
!
                do j = 1, i
                    ij = (i-1)*i/2+j
                    zr(imattt+ij-1) = zr(imattt+ij-1)+poids*(zr(ivf2+l+i-1)*zr(ivf2+l+j-1))
                end do
            end do
        end do
!
    end if
!
! FIN ------------------------------------------------------------------
end subroutine
