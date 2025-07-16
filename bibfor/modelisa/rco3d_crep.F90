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
subroutine rco3d_crep(cara_elem, noma, &
                      lismaco, nbmaco, v_epai)

!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cesexi.h"
#include "asterfort/jedetr.h"
#include "asterfort/detrsd.h"
#include "asterc/indik8.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/utmavo.h"
#include "asterfort/utmess.h"

    character(len=8), intent(in) :: noma
    character(len=8), intent(in) :: cara_elem
    character(len=24), intent(in) :: lismaco
    integer(kind=8), intent(in) :: nbmaco
    real(kind=8), pointer ::  v_epai(:)

! -------------------------------------------------------
!     CREER LE VECTEUR DES EPAISSEURS
!     DES MAILLES DE BORD
! -------------------------------------------------------
!  CARA_ELEM     - IN    - K8   - : NOM DE LA SD CARA_ELEM CONTENANT
!                                   LES CARACTÉRISTIQUES DE L'ÉLÉMENT
!                                   CQOUE
!  NOMA          - IN    - K8   - : NOM DU MAILLAGE
!  LISMACO       - IN    - K24  - : NOM DE LA SD LISTE_MA DES MAILLES
!                                   DE BORD COQUE.
!  NBMACO        - IN    - I    - : NOMBRE DE MAILLES DANS LISMACO.
!  V_EPAI        - IN    - R(*) - : ÉPAISSEUR(S) ASSOCIÉE(S) AUX MAILLES
!                                   DE LISMACO. VECTEUR POUR
!                                   ATTRIBUER UNE ÉPAISSEUR DIFFÉRENTE À
!                                   CHAQUE MAILLE.
! -------------------------------------------------------

    integer(kind=8) :: p3, p4, iret, nb_para_maxi
    integer(kind=8) :: iad1, shell_ep_indx
    real(kind=8) :: shell_ep
    character(len=24)  :: nomavo
    character(len=2)   :: kdim
    integer(kind=8) :: ibid(1), i, j, k, numa
    integer(kind=8), pointer :: v_lmaco(:) => null()
    real(kind=8), pointer :: v_caraelem_cesv(:) => null()
    character(len=8), pointer :: v_caraelem_cesc(:) => null()
    integer(kind=8) :: j_caraelem_cesd, j_caraelem_cesl
    character(len=19) :: cara_elem_s

    ! retrieve some informations
    !
    call jeveuo(lismaco, 'L', vi=v_lmaco)
    !
    ! - Access to elementary characteristics
    !
    cara_elem_s = '&&RACOD3D.CARGEOPO'
    call carces(cara_elem//'.CARCOQUE', 'ELEM', ' ', 'V', cara_elem_s, &
                'A', iret)
    call jeveuo(cara_elem_s//'.CESC', 'L', vk8=v_caraelem_cesc)
    call jeveuo(cara_elem_s//'.CESD', 'L', j_caraelem_cesd)
    call jeveuo(cara_elem_s//'.CESL', 'L', j_caraelem_cesl)
    call jeveuo(cara_elem_s//'.CESV', 'L', vr=v_caraelem_cesv)

    nb_para_maxi = zi(j_caraelem_cesd-1+2)
    shell_ep_indx = indik8(v_caraelem_cesc, 'EP      ', 1, nb_para_maxi)

    !
    nomavo = '&&ORVLMA.MAILLE_VOISINE '
    kdim = '2D'
    call utmavo(noma, kdim, v_lmaco, nbmaco, 'V', &
                nomavo, 0, ibid)
    call jeveuo(jexatr(nomavo, 'LONCUM'), 'L', p4)
    call jeveuo(nomavo, 'L', p3)

    do i = 1, nbmaco
        j = zi(p4+i)-zi(p4-1+i)
        ASSERT(j .gt. 0)
        do k = 1, j
            numa = zi(p3+zi(p4+i-1)-1+k-1)
            call cesexi('C', j_caraelem_cesd, j_caraelem_cesl, numa, 1, &
                        1, shell_ep_indx, iad1)
            if (iad1 .gt. 0) then
                shell_ep = v_caraelem_cesv(iad1)
            else
                call utmess('F', 'CALCULEL3_100', si=numa)
            end if
            ! POUR LE MOMENT ON RECUPERE L'EPAISSEUR DU PREMIER ELEMENT
            v_epai(i) = shell_ep
            exit
        end do
    end do

    ! FIN
    call jedetr(nomavo)
    call detrsd('CHAM_ELEM_S', cara_elem_s)

end subroutine
