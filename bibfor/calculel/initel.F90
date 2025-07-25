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

subroutine initel(ligrel, l_calc_rigi)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/creprn.h"
#include "asterfort/dismoi.h"
#include "asterfort/inigrl.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/typele.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
!
    character(len=19), intent(in) :: ligrel
    aster_logical, optional, intent(out) :: l_calc_rigi
!
! ----------------------------------------------------------------------
!     BUT:
!     INITIALISER LES TYPE_ELEMENTS PRESENTS DANS LE LIGREL (INI00K)
!     CREER (ET REMPLIR) LES OBJETS .PRNM ET/OU .PRNS DU LIGREL.
!
!     IN:
!     LIGREL : NOM DU LIGREL A INITIALISER
!
!     OUT:
!       - INITIALISATION DES ELREFE PRESENTS DANS LE LIGREL
!       - CALCUL DES OBJETS : '.PRNM' ET '.PRNS'
! Out l_calc_rigi : at least one element can support rigidity
!
! ----------------------------------------------------------------------
!
!     VARIABLES LOCALES:
!     ------------------
    integer(kind=8) :: igr, ngr, nmaxob, nbobj, nbprin
    integer(kind=8) :: nbno, jlliel, iconx2
    integer(kind=8) :: nute, nbel, iel, numa, nbnoma, ino, nuno
    parameter(nmaxob=30)
    integer(kind=8) :: adobj(nmaxob)
    character(len=24) :: noobj(nmaxob)
    character(len=1) :: base
    character(len=8) :: exiele, ma, prin, nomail
    character(len=16) :: nomte
    integer(kind=8), pointer :: vprin(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: liel(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
! Initialisation de l'argument optionnel si present
    if (present(l_calc_rigi)) then
        l_calc_rigi = .false.
    end if
!
    call dismoi('EXI_ELEM', ligrel, 'LIGREL', repk=exiele)
    if (exiele(1:3) .eq. 'OUI') then
        call jelira(ligrel//'.LIEL', 'CLAS', cval=base)
    else
!       -- UN LIGREL QUI N'A PAS D'ELEMENTS VIENT FORCEMENT
!          D'UN MODELE QUI DOIT AVOIR DES SOUS-STRUCTURES STATIQUES
        call jelira(ligrel//'.SSSA', 'CLAS', cval=base)
        goto 20
    end if
!
    call jelira(ligrel//'.LIEL', 'NUTIOC', ngr)
    do igr = 1, ngr
        call inigrl(ligrel, igr, nmaxob, adobj, noobj, &
                    nbobj)
    end do
20  continue
!
!
!     -- CALCUL DE .PRNM ET .PRNS :
    call creprn(ligrel, ' ', base, ligrel(1:19)//'.PRNM', ligrel(1:19)//'.PRNS')
!
!
!
!
!     -- ON VERIFIE QUE LES ELEMENTS DE "BORD" SONT COLLES AUX
!        ELEMENTS "PRINCIPAUX" (CEUX QUI CALCULENT LA RIGIDITE):
!     ------------------------------------------------------------
    if ((exiele(1:3) .ne. 'OUI') .or. (.not. present(l_calc_rigi))) goto 90
!
    call jeveuo(ligrel//'.LGRF', 'L', vk8=lgrf)
    call jeveuo(ligrel//'.LIEL', 'L', vi=liel)
    call jeveuo(jexatr(ligrel//'.LIEL', 'LONCUM'), 'L', jlliel)
    ma = lgrf(1)
    call jeveuo(ma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', iconx2)
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbno)
!
!     -- ON COCHE LES NOEUDS PORTES PAR LES ELEMENTS PRINCIPAUX :
    AS_ALLOCATE(vi=vprin, size=nbno)
    do igr = 1, ngr
        nute = typele(ligrel, igr)
        call jenuno(jexnum('&CATA.TE.NOMTE', nute), nomte)
        call dismoi('CALC_RIGI', nomte, 'TYPE_ELEM', repk=prin)
        if (prin .ne. 'OUI') goto 50
        nbel = nbelem(ligrel, igr)
        do iel = 1, nbel
            numa = liel(zi(jlliel+igr-1)+iel-1)
            if (numa .lt. 0) goto 40
            nbnoma = zi(iconx2+numa)-zi(iconx2+numa-1)
            do ino = 1, nbnoma
                nuno = connex(zi(iconx2+numa-1)+ino-1)
                vprin(nuno) = 1
            end do
40          continue
        end do
50      continue
    end do
!
!     -- ON VERIFIE LES NOEUDS DES ELEMENTS NON-PRINCIPAUX (BORD)
    nbprin = 0
    do igr = 1, ngr
        nute = typele(ligrel, igr)
        call jenuno(jexnum('&CATA.TE.NOMTE', nute), nomte)
        call dismoi('CALC_RIGI', nomte, 'TYPE_ELEM', repk=prin)
        nbel = nbelem(ligrel, igr)
        if (prin .eq. 'OUI') then
            if (nbel .gt. 0) nbprin = 1
            goto 80
        end if
        do iel = 1, nbel
            numa = liel(zi(jlliel+igr-1)+iel-1)
            if (numa .lt. 0) goto 70
            nbnoma = zi(iconx2+numa)-zi(iconx2+numa-1)
            do ino = 1, nbnoma
                nuno = connex(zi(iconx2+numa-1)+ino-1)
                if (vprin(nuno) .ne. 1) then
                    nomail = int_to_char8(numa)
                    call utmess('A', 'MODELE1_63', sk=nomail)
                    goto 71
                end if
            end do
71          continue
70          continue
        end do
80      continue
    end do
!
!     -- SI C'EST LE LIGREL DU MODELE, ON VERIFIE QU'IL EXISTE AU MOINS
!        UN ELEMENT PRINCIPAL (QUI CALCULE DE LA RIGIDITE):
    if (present(l_calc_rigi)) then
        l_calc_rigi = .true.
        if (nbprin .eq. 0) then
            l_calc_rigi = .false.
        end if
    end if
!
!
    AS_DEALLOCATE(vi=vprin)
!
!
!
90  continue
    call jedema()
end subroutine
