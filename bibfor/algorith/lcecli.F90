! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine lcecli(fami, kpg, ksp, ndim, mate, &
                  option, lamb, saut, sigma, dsidep, &
                  vim, vip, r)
!
! person_in_charge: patrick.massin at edf.fr
!
! aslint: disable=W1306
    implicit none
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    integer :: mate, ndim, kpg, ksp
    real(kind=8) :: lamb(ndim), saut(ndim), sigma(6), dsidep(6, 6)
    real(kind=8) :: vim(*), vip(*), r, laug(ndim)
    character(len=16) :: option
    character(len=*) :: fami
!
!-----------------------------------------------------------------------
! LOI DE COMPORTEMENT COHESIVE : CZM_LIN_MIX
! Elements cohésifs X-FEM
!
! In : fami     => Schéma d'intégration
! In : kpg, ksp => numéro point et sous-point d'intégration
! In : ndim     => Dimension de l'espace
! In : mate     => Matériau
! In : option   => Option de calcul
! In : lamb     => Champ lambda (en base locale)
! In : saut     => Saut de deplacement (base locale)
! Out: sigma    => Contrainte cohésive (base locale)
! Out: dsidep   => Dérivée de la contrainte par rapport
!                  au multiplicateur augmenté (base locale)
! In : vim      => Variables internes
! Out: vip      => Variables internes actualisées
! Out: r        => Paramètre d'augmentation
!-----------------------------------------------------------------------
!
    aster_logical :: resi, rigi, elas
    integer :: i, j, diss, cass
    real(kind=8) :: sc, gc, lc, val(3), rtan, kapp, zero
    real(kind=8) :: na, ka, rk, ra, coef, coef2
    integer :: cod(3)
    character(len=16) :: nom(3)
    character(len=1) :: poum
    parameter(zero=0.d0)
!
! OPTION CALCUL DU RESIDU OU CALCUL DE LA MATRICE TANGENTE
!
    resi = option(1:9) .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
    rigi = option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RIGI_MECA'
    elas = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'RIGI_MECA_ELAS'
!
!
! les sauts et lambda sont deja ceux a l instant que l on veut calculer
!
!
! RECUPERATION DES PARAMETRES PHYSIQUES
!
    nom(1) = 'GC'
    nom(2) = 'SIGM_C'
    nom(3) = 'PENA_LAGR'
!
    if (option .eq. 'RIGI_MECA_TANG') then
        poum = '-'
    else
        poum = '+'
    end if
!
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'RUPT_FRAG', 0, ' ', [0.d0], &
                3, nom, val, cod, 2)
!
    gc = val(1)
    sc = val(2)
    r = val(3)*sc*sc/gc
    lc = 2*gc/sc
!
! INITIALISATION
!
    ka = max(sc, vim(1))
    rtan = 0.d0
!
! --- CALCUL DE LA FORCE COHESIVE AUGMENTEE
!
    do i = 1, ndim
        laug(i) = lamb(i)+r*saut(i)
    end do
!
! --- FORCE COHESIVE EQUIVALENTE
!
    do i = 2, ndim
        rtan = rtan+laug(i)**2
    end do
    na = sqrt(max(zero, laug(1))**2+rtan)
!
! --- MODULE TANGENT POUR LES REGIMES ELASTIQUES
!
    rk = sc/ka+(1.d0/(r*lc/sc-1.d0))*(sc/ka-1.d0)
!
!
! INITIALISATION COMPLEMENTAIRE POUR RIGI_MECA_TANG (SECANTE PENALISEE)
!
    if (.not. resi) then
!
        if (elas) then
            diss = 0
        else
            diss = nint(vim(2))
        end if
!
        cass = nint(vim(3))
!
        goto 500
    end if
!
! CALCUL DE LA CONTRAINTE
!
    call r8inir(6, 0.d0, sigma, 1)
!
!     CONTRAINTE DE CONTACT : PARTIE NEGATIVE ET NORMALE DE LA FORCE AUGMENTEE
!
    sigma(1) = min(zero, laug(1))
!
!     CONTRAINTE DE FISSURATION
!
    if ((na .ge. (lc*r)) .or. (ka .ge. (lc*r))) then
!
        diss = 0
        cass = 2
!
    else
!
        if (na .le. ka) then
!
            diss = 0
            if (ka .gt. sc) then
                cass = 1
            else
                cass = 0
            end if
            sigma(1) = sigma(1)+rk*max(zero, laug(1))
            do i = 2, ndim
                sigma(i) = sigma(i)+rk*laug(i)
            end do
!
        else
!
            diss = 1
            cass = 1
            ra = sc/na+(1.d0/(r*lc/sc-1.d0))*(sc/na-1.d0)
            sigma(1) = sigma(1)+ra*max(zero, laug(1))
            do i = 2, ndim
                sigma(i) = sigma(i)+ra*laug(i)
            end do
!
        end if
!
    end if
!
! ACTUALISATION DES VARIABLES INTERNES
!
! on conserve la distinction à suivre
    kapp = max(vim(1), na)
    vip(1) = kapp
    vip(2) = diss
    vip(3) = cass
!
! --- VARIABLE INTERNES DE POST-TRAITEMENT PAS REMPLIES POUR L INSTANT
!
! -- MATRICE TANGENTE
!
500 continue
    if (.not. rigi) goto 999
!
    call r8inir(36, 0.d0, dsidep, 1)
!
!    MATRICE TANGENTE DE CONTACT
!
    if (laug(1) .le. 0.d0) dsidep(1, 1) = dsidep(1, 1)+1.d0
!
! DANS LE CAS OU L'ELEMENT EST TOTALEMENT CASSE ON INTRODUIT UNE
! RIGIDITE ARTIFICIELLE DANS LA MATRICE TANGENTE POUR ASSURER
! LA CONVERGENCE
!
    if (cass .eq. 2) then
        if (laug(1) .gt. 0.d0) dsidep(1, 1) = dsidep(1, 1)-1.d-8
        do i = 2, ndim
            dsidep(i, i) = dsidep(i, i)-1.d-8
        end do
        goto 999
    end if
!
!    MATRICE TANGENTE DE FISSURATION
    if ((diss .eq. 0) .or. elas) then
!
        if (laug(1) .gt. 0.d0) dsidep(1, 1) = dsidep(1, 1)+rk
!
        do i = 2, ndim
            dsidep(i, i) = dsidep(i, i)+rk
        end do
!
    else
!
        coef = sc/na+(1.d0/(r*lc/sc-1.d0))*(sc/na-1.d0)
        coef2 = -sc*(1.d0+1.d0/(r*lc/sc-1.d0))/(na**3)
!
        if (laug(1) .le. 0.d0) then
!
            do i = 2, ndim
                dsidep(i, i) = dsidep(i, i)+coef+coef2*laug(i)*laug(i)
            end do
!
            if (ndim .eq. 3) then
                dsidep(2, 3) = dsidep(2, 3)+coef2*laug(2)*laug(3)
                dsidep(3, 2) = dsidep(3, 2)+coef2*laug(3)*laug(2)
            end if
!
        else
!
            do i = 1, ndim
                dsidep(i, i) = dsidep(i, i)+coef+coef2*laug(i)*laug(i)
            end do
!
            do j = 1, ndim-1
                do i = j+1, ndim
                    dsidep(j, i) = dsidep(j, i)+coef2*laug(j)*laug(i)
                    dsidep(i, j) = dsidep(i, j)+coef2*laug(i)*laug(j)
                end do
            end do
!
        end if
!
    end if
999 continue
!
!
end subroutine
