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

subroutine ssmage(nomu, option)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/assmam.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/memame.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/sdmpic.h"
#include "asterfort/ssmau2.h"
#include "asterfort/ualfcr.h"
#include "asterfort/utmess.h"
    character(len=8) :: nomu
    character(len=9) :: option
! ----------------------------------------------------------------------
!     BUT: TRAITER LE MOT CLEF "MASS_MECA" (RESP. "AMOR_MECA")
!             DE L'OPERATEUR MACR_ELEM_STAT
!           CALCULER LA MATRICE DE MASSE (OU AMORTISSEMENT)
!           CONDENSEE DU MACR_ELEM_STAT.
!
!     IN: NOMU   : NOM DU MACR_ELEM_STAT
!         OPTION : 'MASS_MECA' OU 'AMOR_MECA'
!
!     OUT: LES OBJETS SUIVANTS DU MACR_ELEM_STAT SONT CALCULES:
!           / NOMU.MAEL_MASS_VALE (SI MASS_MECA)
!           / NOMU.MAEL_AMOR_VALE (SI AMOR_MECA)
!
! ----------------------------------------------------------------------
!
! - Linear case: no behaviour !
    character(len=24), parameter :: compor = ' '
    real(kind=8) :: time
    character(len=1) :: base
    character(len=8) :: nomo, cara, materi, matel, promes
    character(len=14) :: nu
    character(len=19) :: matas
    character(len=24) :: mateco, mater
    integer(kind=8) :: iarefm
    real(kind=8), pointer :: varm(:) => null()
    character(len=24) :: listElemCalc
!-----------------------------------------------------------------------
    call jemarq()
!
! --- ON CREER LES OBJETS DE TRAVAIL SUR LA VOLATILE
    base = 'V'
!
    call jeveuo(nomu//'.REFM', 'E', iarefm)
    nomo = zk8(iarefm-1+1)
    cara = zk8(iarefm-1+4)
    materi = zk8(iarefm-1+3)
    mater = materi
!
    if (materi .eq. '        ') then
        mateco = ' '
    else
        call rcmfmc(materi, mateco, l_ther_=ASTER_FALSE)
    end if
    nu = zk8(iarefm-1+5)
    if (nu(1:8) .ne. nomu) then
        ASSERT(.false.)
    end if
!
    matel = '&&MATEL'
    if (option .eq. 'MASS_MECA') then
        matas = nomu//'.MASSMECA'
    else if (option .eq. 'AMOR_MECA') then
        matas = nomu//'.AMORMECA'
    else
        ASSERT(.false.)
    end if
!
!
    call jeveuo(nomu//'.VARM', 'L', vr=varm)
    time = varm(2)
!
!     -- CALCULS MATRICES ELEMENTAIRES DE MASSE (OU AMORTISSEMENT):
    if (option .eq. 'MASS_MECA') then
        call dismoi('NOM_LIGREL', nomo, 'MODELE', repk=listElemCalc)
        call memame(option, nomo, mater, mateco, &
                    cara, time, compor, matel, &
                    base, listElemCalc)
    else if (option .eq. 'AMOR_MECA') then
        call dismoi('NOM_PROJ_MESU', nomu, 'MACR_ELEM_STAT', repk=promes)
!     --  CAS MODIFICATION STRUCTURALE : CREATION MATRICE PAR SSMAU2
        if (promes .eq. ' ') then
            call utmess('F', 'SOUSTRUC_69')
        end if
    else
        call utmess('F', 'SOUSTRUC_69')
    end if
!
!        -- ASSEMBLAGE:
    if (option .eq. 'MASS_MECA') then
        call assmam('G', matas, 1, matel, [1.d0], &
                    nu, 'ZERO', 1)
!       -- IL FAUT COMPLETER LA MATRICE SI LES CALCULS SONT DISTRIBUES:
        call sdmpic('MATR_ASSE', matas)
        call ualfcr(matas, 'G')
    end if
    call ssmau2(nomu, option)
!
!        -- MISE A JOUR DE .REFM(7) OU REFM(8)
    if (option .eq. 'MASS_MECA') then
        zk8(iarefm-1+7) = 'OUI_MASS'
    else if (option .eq. 'AMOR_MECA') then
        zk8(iarefm-1+8) = 'OUI_AMOR'
    else
        call utmess('F', 'SOUSTRUC_69')
    end if
!
    if (option .eq. 'MASS_MECA') then
        call jedetr(matel)
    end if
    call jedema()
!
end subroutine
