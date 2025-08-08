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

subroutine veondp(modele, mate, mateco, sddyna, temps, vecelz)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/dbgcal.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/gcncon.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/ndynin.h"
#include "asterfort/ndynkk.h"
#include "asterfort/reajre.h"
#include "asterfort/vemare.h"
#include "asterfort/vrcins.h"
    character(len=*) :: vecelz
    character(len=19) :: sddyna
    character(len=24) :: modele, mate, mateco
    real(kind=8) :: temps
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL)
!
! CALCUL DES VECTEURS ELEMENTAIRES DES ONDES PLANES
!
! ----------------------------------------------------------------------
!
!
! IN  MODELE : NOM DU MODELE
! IN  MATE   : CHAMP DE MATERIAU
! IN  MATECO : CHAMP DE MATERIAU CODE
! IN  SDDYNA : SD DYNAMIQUE
! IN  TEMPS  : INSTANT DE CALCUL
! OUT VECELE : NOM DU VECT_ELEM
!
!
!
!
!
    integer(kind=8) :: nbout, nbin
    parameter(nbout=1, nbin=6)
    character(len=8) :: lpaout(nbout), lpain(nbin), noma
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    integer(kind=8) :: ibid, i, iret, iondp
    character(len=19) :: vecele
    character(len=24) :: chinst
    character(len=24) :: chgeom, ligrmo
    integer(kind=8) :: nchond
    character(len=19) :: chondp
    aster_logical :: debug
    character(len=8) :: newnom
    character(len=16) :: option
    integer(kind=8) :: ifmdbg, nivdbg
    character(len=19) :: chvarc
    character(len=2) :: codret
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('PRE_CALCUL', ifmdbg, nivdbg)
!
! --- INITIALISATIONS
!
    call ndynkk(sddyna, 'CHONDP', chondp)
    nchond = ndynin(sddyna, 'NBRE_ONDE_PLANE')
    vecele = vecelz
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrmo)
    call dismoi('NOM_MAILLA', ligrmo, 'LIGREL', repk=noma)
    chgeom = noma//'.COORDO'
    option = 'ONDE_PLAN'
    chinst = '&&CHINST'
    if (nivdbg .ge. 2) then
        debug = .true.
    else
        debug = .false.
    end if
    call jeveuo(chondp, 'L', iondp)
!
! --- CREATION D'UNE CARTE D'INSTANTS
!
    call mecact('V', chinst, 'MODELE', ligrmo, 'INST_R', &
                ncmp=1, nomcmp='INST', sr=temps)
!
! --- CREATION CHAMP DE VARIABLES DE COMMANDE CORRESPONDANT
!
    chvarc = '&&CHME.ONDPL.CHVARC'
    call vrcins(modele, mate, ' ', temps, chvarc, &
                codret)
!
! --- CHAMPS D'ENTREE
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = mateco(1:19)
    lpain(3) = 'PINSTR'
    lchin(3) = chinst(1:19)
    lpain(4) = 'PONDPLA'
    lpain(5) = 'PONDPLR'
    lpain(6) = 'PVARCPR'
    lchin(6) = chvarc
!
! --- CHAMPS DE SORTIE
!
    lpaout(1) = 'PVECTUR'
    lchout(1) = vecele(1:8)//'.???????'
!
    call detrsd('VECT_ELEM', vecele)
    call vemare('V', vecele, modele(1:8))
!
! -- CALCUL
!
    do i = 1, nchond
        call exisd('CARTE', zk8(iondp+i-1)//'.CHME.ONDPL', iret)
        call exisd('CARTE', zk8(iondp+i-1)//'.CHME.ONDPR', ibid)
!
        if (iret .ne. 0 .and. ibid .ne. 0) then
            lchin(4) = zk8(iondp+i-1)//'.CHME.ONDPL'
            lchin(5) = zk8(iondp+i-1)//'.CHME.ONDPR'
            call gcncon('.', newnom)
            lchout(1) (10:16) = newnom(2:8)
!
            call calcul('S', option, ligrmo, nbin, lchin, &
                        lpain, nbout, lchout, lpaout, 'V', &
                        'OUI')
            call corich('E', lchout(1), ichin_=-1)
!
            if (debug) then
                call dbgcal(option, ifmdbg, nbin, lpain, lchin, &
                            nbout, lpaout, lchout)
            end if
!
            call reajre(vecele, lchout(1), 'V')
        end if
    end do
!
    call detrsd('CHAMP_GD', chvarc)
!
    call jedema()
end subroutine
