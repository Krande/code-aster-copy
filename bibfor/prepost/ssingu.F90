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
subroutine ssingu(nomail, nelem, nbr, ligrmo, alpha, &
                  re, he, chelem)
    implicit none
#include "jeveux.h"
#include "asterfort/cescel.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: nelem, nbr(nelem)
    real(kind=8) :: alpha(nelem), re(nelem), he(nelem)
    character(len=8) :: nomail
    character(len=24) :: ligrmo, chelem
!
!     BUT:
!         STOCKAGE DE LA SINGULARITE (ALPHA), DU RAPPORT DE TAILLE (RE)
!         ET DE LA NOUVELLE TAILLE (TAILLE) DES ELEMENTS DANS CHELEM
!         OPTION : 'SING_ELEM'
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   NOMAIL       : NOM DU MAILLAGE
! IN   NELEM        : NOMBRE D ELEMENTS FINIS DU MAILLAGE
! IN   NBR(NELEM)   : NOMBRE DE COMPOSANTES A STOCKER PAR EF
!      3 SI EF SURFACIQUES EN 2D OU VOLUMIQUES EN 3D
!      0 SINON
! IN   LIGRMO       : NOM DU LIGREL DU MODELE
! IN   ALPHA(NELEM) : DEGRE DE LA SINGULARITE
! IN   RE(NELEM)    : RAPPORT ENTRE ANCIENNE ET NOUVELLE TAILLE
! IN   HE(NELEM)    : NOUVELLE TAILLE
! IN   CHELEM       : CHAM_ELEM QUI VA CONTENIR LE DEGRE ET LA TAILLE
!
!      SORTIE :
!-------------
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: jcesd, jcesl, iad1, iad2, iad3
    integer(kind=8) :: nbcmp, ncmp1, ncmp2, ncmp3
    integer(kind=8) :: icmp, inel, nncp, ibid
    character(len=8) :: nompaz, licmp(3)
    character(len=16) :: opti
    character(len=19) :: chsing
    character(len=8), pointer :: cesc(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
!
    call jemarq()
!
! 1 - CREATION D UN CHAM_ELEM_S CHSING
!
    chsing = '&&SINGUE.SING'
    nompaz = 'SING_R'
    licmp(1) = 'DEGRE'
    licmp(2) = 'RAPPORT'
    licmp(3) = 'TAILLE'
    call cescre('V', chsing, 'ELEM', nomail, nompaz, &
                3, licmp, [-1], [-1], nbr)
!
! 2 - STOCKAGE DANS CHSING DE ALPHA ET RE
!
    call jeveuo(chsing//'.CESC', 'L', vk8=cesc)
    call jelira(chsing//'.CESC', 'LONMAX', nbcmp)
    call jeveuo(chsing//'.CESD', 'L', jcesd)
    call jeveuo(chsing//'.CESL', 'E', jcesl)
    call jeveuo(chsing//'.CESV', 'E', vr=cesv)
!
    do icmp = 1, nbcmp
        if (cesc(icmp) (1:5) .eq. 'DEGRE') ncmp1 = icmp
        if (cesc(icmp) (1:7) .eq. 'RAPPORT') ncmp2 = icmp
        if (cesc(icmp) (1:6) .eq. 'TAILLE') ncmp3 = icmp
    end do
!
    do inel = 1, nelem
        call cesexi('C', jcesd, jcesl, inel, 1, &
                    1, ncmp1, iad1)
        call cesexi('C', jcesd, jcesl, inel, 1, &
                    1, ncmp2, iad2)
        call cesexi('C', jcesd, jcesl, inel, 1, &
                    1, ncmp3, iad3)
        if ((-iad1) .gt. 0) then
            zl(jcesl-iad1-1) = .true.
            cesv(1-iad1-1) = alpha(inel)
        end if
        if ((-iad2) .gt. 0) then
            zl(jcesl-iad2-1) = .true.
            cesv(1-iad2-1) = 1.d0/re(inel)
        end if
        if ((-iad3) .gt. 0) then
            zl(jcesl-iad3-1) = .true.
            cesv(1-iad3-1) = he(inel)
        end if
    end do
!
! 3 - TRANSFORMATION DU CHAM_ELEM_S EN CHAM_ELEM
!
    opti = 'SING_ELEM'
    nompaz = 'PSING_R'
!
    call cescel(chsing, ligrmo(1:19), opti, nompaz, 'NON', &
                nncp, 'G', chelem(1:19), 'F', ibid)
!
    call detrsd('CHAM_ELEM_S', chsing)
!
    call jedema()
!
end subroutine
