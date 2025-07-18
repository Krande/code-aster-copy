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

subroutine verif_affe(modele, sd, non_lin)
    implicit none
!
! person_in_charge: jacques.pellet at edf.fr
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterfort/getvtx.h"
#include "asterfort/exisd.h"
#include "asterfort/assert.h"
#include "asterfort/verif_affe_carte.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"
!
    character(len=*), intent(in) :: modele
    character(len=*), intent(in), optional :: sd
    aster_logical, intent(in), optional ::  non_lin
!
!-----------------------------------------------------------------------
!   But :
!     Emettre des alarmes concernant les affectations douteuses dans les
!     differentes cartes des SD charge et cara_elem
!
!   Entrees:
!     modele     :  sd_modele
!     sd         :  sd_charge ou sd_cara_elem
!
!-----------------------------------------------------------------------
    character(len=8) :: sd_, modele_
    character(len=16) :: verif
    character(len=19) :: carte, ligrmo
    character(len=24) :: typres
    character(len=80) :: comment
    integer(kind=8) :: n1, k, iret
    character(len=5)  :: l_cart_char_meca(24), l_cart_char_ther(11)
    character(len=8)  :: l_cart_cara_elem(14)
    character(len=80) :: l_comm_char_meca(26), l_comm_char_ther(11), l_comm_cara_elem(14)

!-----------------------------------------------------------------------
!
    call jemarq()

    verif = 'OUI'
    call getvtx(' ', 'VERI_AFFE', scal=verif, nbret=n1)
    if (verif .eq. 'NON') goto 999

!   -- cartes des sd_char_meca :
!   -----------------------------
    l_cart_char_meca = (/ &
                       'EFOND', 'EPSIN', 'F1D1D', 'F1D2D', 'F1D3D', 'F2D2D', 'F2D3D', 'F3D3D', &
                       'FCO2D', 'FCO3D', 'FELEC', 'FLUX ', 'FORNO', 'IMPE ', 'ONDE ', 'ONDPL', &
                       'ONDPR', 'PESAN', 'PREFF', 'PRESS', 'ROTAT', 'SIGIN', 'SIINT', 'VFACE'/)

    l_comm_char_meca = ' '
    l_comm_char_meca(22) = 'Chargement provenant du mot cle PRES_REP'

!   -- cartes des sd_char_ther :
!   -----------------------------
    l_cart_char_ther = (/ &
                       'SOURE', 'COEFH', 'FLUNL', 'SOUNL', 'FLUR2', 'FLURE', 'GRAIN', 'HECHP', &
                       'RAYO ', 'T_EXT', 'SOURC'/)

    l_comm_char_ther = ' '
    l_comm_char_ther(1) = 'Chargement provenant du mot cle SOURCE'

!   -- cartes des sd_cara_elem :
!   -----------------------------
    l_cart_cara_elem = (/ &
                       'CARGENBA', 'CARMASSI', 'CARCABLE', 'CARCOQUE', 'CARDISCK', 'CARARCPO', &
                       'CARGENPO', 'CARDISCM', 'CARORIEN', 'CARDISCA', 'CVENTCXF', 'CARPOUFL', &
                       'CARGEOPO', 'CARDINFO'/)

    l_comm_cara_elem = ' '
    l_comm_cara_elem(1) = 'Caracteristiques provenant du mot cle BARRE'
    l_comm_cara_elem(2) = 'Caracteristiques provenant du mot cle MASSIF'
    l_comm_cara_elem(3) = 'Caracteristiques provenant du mot cle CABLE'
    l_comm_cara_elem(4) = 'Caracteristiques provenant du mot cle COQUE'
    l_comm_cara_elem(5) = 'Caracteristiques provenant du mot cle DISCRET de type raideur'
    l_comm_cara_elem(6) = 'Caracteristiques provenant des POUTRES courbes'
    l_comm_cara_elem(7) = 'Caracteristiques provenant du mot cle POUTRE'
    l_comm_cara_elem(8) = 'Caracteristiques provenant du mot cle DISCRET de type masse'
    l_comm_cara_elem(9) = 'Caracteristiques provenant du mot cle ORIENTATION'
    l_comm_cara_elem(10) = 'Caracteristiques provenant du mot cle DISCRET d''amortissement'
    l_comm_cara_elem(11) = 'Caracteristiques provenant du mot cle FCX'
    l_comm_cara_elem(12) = 'Caracteristiques provenant du mot cle POUTRE_FLUI'
    l_comm_cara_elem(13) = 'Caracteristiques provenant du mot cle POUTRE (geometrie)'
    l_comm_cara_elem(14) = 'Caracteristiques provenant du mot cle DISCRET (information)'

    modele_ = modele
    ligrmo = modele_//'.MODELE'

!   -- boucle sur les cartes de la SD :
!   ---------------------------------------
    if (present(sd)) then
        sd_ = sd
        call gettco(sd_, typres)

        if (typres .eq. 'CHAR_MECA') then
            n1 = size(l_cart_char_meca)
            do k = 1, n1
                carte = sd_//'.CHME.'//l_cart_char_meca(k)
                comment = l_comm_char_meca(k)
                call exisd('CARTE', carte, iret)
                if (iret .gt. 0) call verif_affe_carte(ligrmo, carte, comment)
            end do

        elseif (typres .eq. 'CHAR_THER') then
            n1 = size(l_cart_char_ther)
            do k = 1, n1
                carte = sd_//'.CHTH.'//l_cart_char_ther(k)
                comment = l_comm_char_ther(k)
                call exisd('CARTE', carte, iret)
                if (iret .gt. 0) call verif_affe_carte(ligrmo, carte, comment)
            end do

        elseif (typres .eq. 'CARA_ELEM') then
            n1 = size(l_cart_cara_elem)
            do k = 1, n1
                carte = sd_//'.'//l_cart_cara_elem(k)
                comment = l_comm_cara_elem(k)
                call exisd('CARTE', carte, iret)
                if (present(non_lin)) then
                    if (iret .gt. 0) call verif_affe_carte(ligrmo, carte, comment, non_lin=non_lin)
                else
                    if (iret .gt. 0) call verif_affe_carte(ligrmo, carte, comment)
                end if
            end do

        else
            write (6, *) 'A faire ... typres=', sd, typres
            ASSERT(.false.)
        end if

    end if

999 continue

    call jedema()
end subroutine
