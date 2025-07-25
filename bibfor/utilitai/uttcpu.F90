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

subroutine uttcpu(nommes, action, nomlon)
    implicit none
! person_in_charge: jacques.pellet at edf.fr
#include "jeveux.h"
!
#include "asterfort/assert.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/uttcp0.h"
#include "asterfort/wkvect.h"
    character(len=*) :: nommes, action, nomlon
! ----------------------------------------------------------------------
!  ROUTINE DE MESURE DU TEMPS CPU.
!
! IN  NOMMES (K24)   : NOM (COURT) IDENTIFIANT LA MESURE
!
! IN  ACTION :  ACTION = 'INIT'  LA MESURE DE TEMPS EST (RE)MIS A ZERO
!               ACTION = 'DEBUT' LA MESURE DE TEMPS COMMENCE
!               ACTION = 'FIN'   LA MESURE DE TEMPS S'ARRETE
! IN  NOMLON (K80)   : NOM (LONG) ASSOCIE A LA MESURE (OU ' ').
!              CE NOM EST STOCKE POUR ACTION='INIT'
!              CE NOM SERA EVENTUELLENT IMPRIME PAR UTTCPI
! ----------------------------------------------------------------------
! ON ACCUMULE 7 VALEURS MESUREES POUR CHAQUE MESURE (NOMMES) :
!    TEMPS(1) TEMPS RESTANT EN SECONDES
!    TEMPS(2) NOMBRE D'APPEL A DEBUT/FIN
!    TEMPS(3) TEMPS CPU TOTAL
!    TEMPS(4) TEMPS ELAPSED MOYEN
!    TEMPS(5) TEMPS CPU USER TOTAL
!    TEMPS(6) TEMPS CPU SYSTEME
!    TEMPS(7) TEMPS ELAPSED
! LES VALEURS STOCKEES SONT RECUPERABLES VIA UTTCPR
! ----------------------------------------------------------------------
    integer(kind=8) :: indi, iexi, jvalms, jnoml, jvalmi, k
!
!     -- COMMONS POUR MESURE DE TEMPS :
    integer(kind=8) :: mtpniv, mtpsta, indmax
    parameter(indmax=5)
    character(len=80) :: snolon(indmax)
    real(kind=8) :: valmes(indmax*7), valmei(indmax*7)
    common/mestp1/mtpniv, mtpsta
    common/mestp2/snolon
    common/mestp3/valmes, valmei
! ----------------------------------------------------------------------
!
!     -- POUR CERTAINES MESURES, ON NE PEUT PAS FAIRE DE JEVEUX :
!        ON GARDE ALORS LES INFOS DANS LES COMMON MESTPX
    if (nommes .eq. 'CPU.MEMD.1') then
        indi = 1
    else if (nommes .eq. 'CPU.MEMD.2') then
        indi = 2
    else
        goto 9998
    end if
    ASSERT(indi .le. indmax)
!
    if (action .eq. 'INIT') then
        snolon(indi) = nomlon
!       -- IL FAUT REMETTRE LES COMMON A ZERO :
        do k = 1, 7
            valmes(k) = 0.d0
            valmei(k) = 0.d0
        end do
    end if
!
    call uttcp0(indi, action, 7, valmes(7*(indi-1)+1))
    goto 9999
!
!
!
9998 continue
!     -- INITIALISATION DES OBJETS JEVEUX :
    call jeexin('&&UTTCPU.NOMMES', iexi)
    if (iexi .eq. 0) then
        call jecreo('&&UTTCPU.NOMMES', 'V N K24')
        call jeecra('&&UTTCPU.NOMMES', 'NOMMAX', 100)
        call wkvect('&&UTTCPU.VALMES', 'V V R', 7*100, jvalms)
        call wkvect('&&UTTCPU.VALMEI', 'V V R', 7*100, jvalmi)
        call wkvect('&&UTTCPU.NOMLON', 'V V K80', 100, jnoml)
    else
        call jeveuo('&&UTTCPU.VALMES', 'E', jvalms)
    end if
!
    call jenonu(jexnom('&&UTTCPU.NOMMES', nommes), indi)
    if (indi .eq. 0) then
        if (action .eq. 'INIT') then
            call jecroc(jexnom('&&UTTCPU.NOMMES', nommes))
            call jenonu(jexnom('&&UTTCPU.NOMMES', nommes), indi)
            ASSERT(indi .gt. 0)
            ASSERT(indi .lt. 100)
            call jeveuo('&&UTTCPU.NOMLON', 'E', jnoml)
            zk80(jnoml-1+indi) = nomlon
        else
!         -- LA MESURE N'A PAS ETE INITIALISEE: ON NE LA FAIT PAS
            goto 9999
        end if
    end if
!
    call uttcp0(100+indi, action, 7, zr(jvalms-1+7*(indi-1)+1))
!
9999 continue
end subroutine
