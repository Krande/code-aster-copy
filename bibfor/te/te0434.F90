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

subroutine te0434(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/mbgchg.h"
#include "asterfort/mbxchg.h"
#include "asterfort/rccoma.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE CHARGEMENT :
!                                  - CHAR_MECA_EPSI_R
!                                  - CHAR_MECA_EPSI_F
!                                  - CHAR_MECA_PESA_R
!                                  - CHAR_MECA_TEMP_R
!                                  - FORC_NODA
!                                  - REFE_FORC_NODA
!                          POUR LES MEMBRANES
!    - ARGUMENTS :
!        DONNEES :      OPTION       -->  OPTION DE CALCUL
!                       NOMTE        -->  NOM DU TYPE ELEMENT
! ----------------------------------------------------------------------
!
    character(len=4) :: fami
    character(len=32) :: phenom
    integer(kind=8) :: nddl, nno, nnos, npg, ndim, ncomp
    integer(kind=8) :: n, kpg
    integer(kind=8) :: ipoids, ivf, idfde, jgano, iret, icompo, itab(1), itemps
    integer(kind=8) :: igeom, icacoq, imate, jvSief, ipesa, iepsin, ivectu
    integer(kind=8) :: icodre1, icodre2
    real(kind=8) :: dff(2, 9), vff(9)
    real(kind=8) :: alpha, beta, h, preten
    aster_logical :: grav
!
!
! -----------------------------------------------------------------
! ---              INITIALISATION DES VARIABLES                 ---
! -----------------------------------------------------------------
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
!
! - PARAMETRES EN ENTREE
! - grav : permet d'utiliser PESANTEUR en STAT_NON_LINE
    grav = (option .eq. 'CHAR_MECA_PESA_R')

    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCACOQU', 'L', icacoq)

    call tecach('N', 'PCOMPOR', 'L', iret, 1, itab)
    icompo = itab(1)
!
    if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', jvSief)
        call jevech('PMATERC', 'L', imate)
!
    else if (option .eq. 'REFE_FORC_NODA') then
        call jevech('PMATERC', 'L', imate)
!
    else if (option .eq. 'CHAR_MECA_EPSI_R') then
        call jevech('PMATERC', 'L', imate)
        call jevech('PEPSINR', 'L', iepsin)
!
    else if (option .eq. 'CHAR_MECA_EPSI_F') then
        call jevech('PMATERC', 'L', imate)
        call jevech('PEPSINF', 'L', iepsin)
        call jevech('PINSTR', 'L', itemps)
!
    else if (grav) then
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
! - DIRECTION DE REFERENCE POUR UN COMPORTEMENT ANISOTROPE
!
    alpha = zr(icacoq+1)*r8dgrd()
    beta = zr(icacoq+2)*r8dgrd()

! - EPAISSEUR ET PRETCONTRAINTES
    h = zr(icacoq)
    preten = zr(icacoq+3)/h
!
! -----------------------------------------------------------------
! ---  VERIFICATION DE LA CORRESPONDANCE MATERIAU / COMPORTMENT ---
! -----------------------------------------------------------------
!
    call rccoma(zi(imate), 'ELAS_MEMBRANE', 0, phenom, icodre1)
    call rccoma(zi(imate), 'ELAS', 0, phenom, icodre2)

    if (icodre1 .eq. 0) then
        if ((icompo .ne. 0) .and. (zk16(icompo+2) (1:5) .ne. 'PETIT')) then
            call utmess('F', 'MEMBRANE_10')
        end if
    elseif (icodre2 .eq. 0) then
        if (((icompo .eq. 0) .or. (zk16(icompo+2) (1:9) .ne. 'GROT_GDEP')) &
            .and. (.not. grav)) then
            call utmess('F', 'MEMBRANE_10')
        end if
    end if
!
! -----------------------------------------------------------------
! ---       DEBUT DE LA BOUCLE SUR LES POINTS DE GAUSS          ---
! -----------------------------------------------------------------
!
    do kpg = 1, npg
!
! --- MISE SOUS FORME DE TABLEAU DES VALEURS ET DES DERIVEES
!     DES FONCTIONS DE FORME
!
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do

        if (icodre1 .eq. 0) then

            call mbxchg(option, fami, nddl, nno, ncomp, kpg, npg, iepsin, itemps, ipoids, igeom, &
                        imate, ipesa, ivectu, jvSief, vff, dff, alpha, beta)

        elseif (icodre2 .eq. 0) then

            if ((option .ne. 'FORC_NODA') .and. (option .ne. 'CHAR_MECA_PESA_R')) then
                call utmess('F', 'MEMBRANE_7')
            end if
            call mbgchg(option, fami, nddl, nno, ncomp, kpg, imate, jvSief, &
                        ipoids, ipesa, igeom, ivectu, vff, dff, h, alpha, beta, preten)

        end if
    end do
!
! - FIN DE LA BOUCLE SUR LES POINTS DE GAUSS
!
end subroutine
