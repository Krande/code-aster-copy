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
subroutine te0296(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/matrot.h"
#include "asterfort/ntfcma.h"
#include "asterfort/ppgan2.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcdiff.h"
#include "asterfort/rcfode.h"
#include "asterfort/rcvalb.h"
#include "asterfort/runge6.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/utrcyl.h"
#include "asterfort/uttgel.h"
!
    character(len=16) :: nomte, option
! ----------------------------------------------------------------------
!
!    - FONCTION REALISEE:  CALCUL DES VECTEURS RESIDUS
!                          OPTION : 'MASS_THER_RESI'
!                          ELEMENTS 3D ISO PARAMETRIQUES
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! THERMIQUE NON LINEAIRE
!
!
!
    integer :: nbres
    parameter(nbres=3)
    integer :: icodre(nbres)
    character(len=2) :: typgeo
    character(len=32) :: phenom
    real(kind=8) :: beta, deltat, tpg, tpgm
    real(kind=8) :: dfdx(27), dfdy(27), dfdz(27), poids, hydrgm(27)
    real(kind=8) :: rbid, chal(1), hydrgp(27), err
    integer :: igeom, imate
    integer :: nno, kp, i, itemps, ifon(6), l, ndim
    integer :: ihydr, ihydrp, itempr
    integer :: isechf, jgano2
    integer :: icomp, itempi, iveres, nnos
    integer :: npg2, ipoid2, ivf2, idfde2
    aster_logical :: aniso
! ----------------------------------------------------------------------
! PARAMETER ASSOCIE AU MATERIAU CODE
! --- INDMAT : INDICE SAUVEGARDE POUR LE MATERIAU
!CC      PARAMETER        ( INDMAT = 8 )
! ----------------------------------------------------------------------
!====
! 1.1 PREALABLES: RECUPERATION ADRESSES FONCTIONS DE FORMES...
!====
    call uttgel(nomte, typgeo)
    if ((lteatt('LUMPE', 'OUI')) .and. (typgeo .ne. 'PY')) then
        call elrefe_info(fami='NOEU', ndim=ndim, nno=nno, nnos=nnos, npg=npg2, &
                         jpoids=ipoid2, jvf=ivf2, jdfde=idfde2, jgano=jgano2)
    else
        call elrefe_info(fami='MASS', ndim=ndim, nno=nno, nnos=nnos, npg=npg2, &
                         jpoids=ipoid2, jvf=ivf2, jdfde=idfde2, jgano=jgano2)
    end if
!
!====
! 1.2 PREALABLES LIES AUX RECHERCHES DE DONNEES GENERALES
!====
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PTEMPEI', 'L', itempi)
    call jevech('PCOMPOR', 'L', icomp)
    call jevech('PRESIDU', 'E', iveres)
!
    deltat = zr(itemps+1)
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
!----
!   INITIALISATION THER_HYDR
!----
        if (zk16(icomp) (1:9) .eq. 'THER_HYDR') then
            call jevech('PHYDRPM', 'L', ihydr)
            call jevech('PHYDRPP', 'E', ihydrp)
            call jevech('PTEMPER', 'L', itempr)
            call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                        ' ', 'THER_HYDR', 0, ' ', [0.d0], &
                        1, 'CHALHYDR', chal, icodre(1), 1)
            do kp = 1, npg2
                l = nno*(kp-1)
                hydrgm(kp) = 0.d0
                do i = 1, nno
                    hydrgm(kp) = hydrgm(kp)+zr(ihydr)*zr(ivf2+l+i-1)
                end do
            end do
        end if
!
! ---   CALCUL DU DEUXIEME TERME
!
        do kp = 1, npg2
            l = (kp-1)*nno
            call dfdm3d(nno, kp, ipoid2, idfde2, zr(igeom), &
                        poids, dfdx, dfdy, dfdz)
            tpg = 0.d0
            do i = 1, nno
                tpg = tpg+zr(itempi+i-1)*zr(ivf2+l+i-1)
            end do
!
! ---       RESOLUTION DE L EQUATION D HYDRATATION
!
            if (zk16(icomp) (1:9) .eq. 'THER_HYDR') then
                tpgm = 0.d0
                hydrgp(kp) = 0.d0
                do i = 1, nno
                    tpgm = tpgm+zr(itempr+i-1)*zr(ivf2+l+i-1)
                end do
                call runge6(ifon(3), deltat, tpg, tpgm, hydrgm(kp), &
                            hydrgp(kp), err)
            end if
!
            call rcfode(ifon(1), tpg, beta, rbid)
            if (zk16(icomp) (1:9) .eq. 'THER_HYDR') then
!
! ---           THERMIQUE NON LINEAIRE AVEC HYDRATATION
!
                do i = 1, nno
                    zr(iveres+i-1) = zr(iveres+i-1)+poids*((beta-chal(1)*hydrgp(kp))&
                                                          &*zr(ivf2+l+i-1))
                end do
            else
!
! ---           THERMIQUE NON LINEAIRE SEULE
!
                do i = 1, nno
                    zr(iveres+i-1) = zr(iveres+i-1)+poids*beta*zr(ivf2+l+i-1)
                end do
            end if
        end do
!
    else if (zk16(icomp) (1:5) .eq. 'SECH_') then
!====
! --- SECHAGE
!====
        if (zk16(icomp) (1:12) .eq. 'SECH_GRANGER' .or. zk16(icomp) (1:10) .eq. &
            'SECH_NAPPE') then
            call jevech('PTMPCHF', 'L', isechf)
        else
!          POUR LES AUTRES LOIS, PAS DE CHAMP DE TEMPERATURE
!          ISECHI ET ISECHF SONT FICTIFS
            isechf = itempi
        end if
!
! ---   CALCUL DU DEUXIEME TERME
!
        do kp = 1, npg2
            l = (kp-1)*nno
            call dfdm3d(nno, kp, ipoid2, idfde2, zr(igeom), &
                        poids, dfdx, dfdy, dfdz)
            tpg = 0.d0
            do i = 1, nno
                tpg = tpg+zr(itempi+i-1)*zr(ivf2+l+i-1)
            end do
            do i = 1, nno
                zr(iveres+i-1) = zr(iveres+i-1)+poids*(1.d0*zr(ivf2+l+i-1)*tpg)
            end do
        end do
    end if
!
    if (zk16(icomp) (1:9) .eq. 'THER_HYDR') call ppgan2(jgano2, 1, 1, hydrgp, zr(ihydrp))
!
! FIN ------------------------------------------------------------------
end subroutine
