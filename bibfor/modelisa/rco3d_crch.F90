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
!
subroutine rco3d_crch(ligrel, noma, chmlrac, lismaco, &
                      nbmaco, crig, v_epai)

    !
    implicit none
    !
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"

    character(len=19), intent(in) :: ligrel, chmlrac
    character(len=8), intent(in) :: noma
    character(len=24), intent(in) :: lismaco
    real(kind=8), intent(in) :: crig
    integer(kind=8), intent(in) :: nbmaco
    real(kind=8), pointer ::  v_epai(:)

! -------------------------------------------------------
!     CREER LE CHAMP D ENTREE POUR LE TE DU RACCORD
!              COQUE-3D
! -------------------------------------------------------
!  LIGREL        - IN    - K19  - : NOM DU LIGREL
!
!  NOMA          - IN    - K8   - : NOM DU MAILLAGE.
!
!  CHMLRAC       - IN    - K19  - : NOM DE LA SD CHAM_ELEM D ENTREE
!
!  LISMACO       - IN    - K24  - : NOM DE LA SD LISTE_MA CONTENANT
!                                   LES MAILLES DE BORD DE TYPE COQUE.
!
!  NBMACO        - IN    - I    - : NOMBRE DE MAILLES DANS LISMACO.
!
!  CRIG          - IN    - R    - : COEFFICIENT
!
!  V_EPAI        - OUT   - R(*) - : VECTEUR CONTENANT LES ÉPAISSEURS
!                                   ASSOCIÉES AUX MAILLES DE LISMACO.
! -------------------------------------------------------

    integer(kind=8), parameter :: nceld1 = 4
    integer(kind=8), parameter :: nceld2 = 4
    integer(kind=8), parameter :: nceld3 = 4
    integer(kind=8) :: p2
    real(kind=8) :: epais
    character(len=24)  :: chmlrac_celv
    integer(kind=8) :: jv_chmlrac_celv, nb_grel, i_grel
    integer(kind=8) :: decal, i_liel, nb_liel, vale_indx
    integer(kind=8) :: idx, j, k, l
    character(len=24) :: chmlrac_celd
    integer(kind=8), pointer :: v_lmaco(:) => null()
    integer(kind=8), pointer :: v_chmlrac_celd(:) => null()
    integer(kind=8), pointer :: v_ligrel_liel(:) => null()
    integer(kind=8), pointer :: v_list_no_pair(:) => null()
    integer(kind=8) :: num_pair
    integer(kind=8) :: idconnco, nbnoco

    ! retrieve some informations
!
    call jeveuo(lismaco, 'L', vi=v_lmaco)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', p2)
!
    chmlrac_celd = chmlrac//'.CELD'
    chmlrac_celv = chmlrac//'.CELV'
    call jeveuo(chmlrac_celd, 'L', vi=v_chmlrac_celd)
    call jeveuo(chmlrac_celv, 'E', jv_chmlrac_celv)
    nb_grel = v_chmlrac_celd(2)

    do i_grel = 1, nb_grel
        decal = v_chmlrac_celd(nceld1+i_grel)
        nb_liel = v_chmlrac_celd(decal+1)
        call jeveuo(jexnum(ligrel//'.LIEL', i_grel), 'L', vi=v_ligrel_liel)
        do i_liel = 1, nb_liel
            num_pair = -v_ligrel_liel(i_liel)
            call jeveuo(jexnum(ligrel//'.NEMA', num_pair), 'L', vi=v_list_no_pair)

            do k = 1, nbmaco
                call jeveuo(jexnum(noma//'.CONNEX', v_lmaco(k)), 'L', idconnco)
                call jelira(jexnum(noma//'.CONNEX', v_lmaco(k)), 'LONMAX', nbnoco)
                idx = 0
                do l = 1, nbnoco
                    do j = 1, 2
                        if (v_list_no_pair(j) .eq. zi(idconnco+l-1)) then
                            idx = idx+1
                        end if
                    end do
                end do
                if (idx .ge. 2) then
                    exit
                end if
            end do

            ! CHECK NECESSAIRE AVANT DE CONTINUER
            ASSERT(idx .ge. 2)
            ! AFFECTER LA BONNE EPAISSEUR
            epais = v_epai(k)

            vale_indx = jv_chmlrac_celv-1+v_chmlrac_celd(decal+nceld2+nceld3*(i_liel-1)+4)
            zr(vale_indx-1+1) = epais
            zr(vale_indx-1+2) = 0.0d0
            zr(vale_indx-1+3) = 0.0d0
            zr(vale_indx-1+4) = crig
            zr(vale_indx-1+5) = 0.0d0
            zr(vale_indx-1+6) = 0.0d0
        end do
    end do

end subroutine
