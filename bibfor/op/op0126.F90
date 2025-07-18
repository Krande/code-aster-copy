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
subroutine op0126()
    implicit none
!  P. RICHARD     DATE 13/07/90
!-----------------------------------------------------------------------
!  BUT: TRAITER LA DEFINITION DU MODELE GENERALISE DONNE PAR
!       L'UTILISATEUR ET TRAITER L'ORIENTATION DES MATRICES DE LIAISON
!       PROCEDER AUX VERIFICATIONS SUR LA COHERENCE DE LA DEFINITION
!       DES LIAISONS ET SUR LA COMPATIBILITE DES MACR_ELEM MIS EN JEU
!
!  CONCEPT CREE: MODE_GENE
!
!-----------------------------------------------------------------------
!
!
!
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/arg126.h"
#include "asterfort/callis.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/vecomo.h"
#include "asterfort/versst.h"
    integer(kind=8) :: ival, ibid, nblia, i, iinc, irep11, irep12, irep21, irep22, iopt
    integer(kind=8) :: iret
    character(len=3) :: rep
    character(len=8) :: nomres, sst1, sst2, intf1, intf2, k8bid, option
    character(len=16) :: nomcon, nomope
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
    call infmaj()
    call getres(nomres, nomcon, nomope)
!
!-----TRAITEMENT DES DONNEES UTILISATEUR
!
    call arg126(nomres)
!
!-----VERIFICATION COHERENCE DES SOUS-STRUCTURES ET CREATION DU .DESC
!
    call versst(nomres)
!
!
!-----VERIFICATION DE LA COHERENCE DU MODELE GENERALISE
!
    call getfac('VERIF', ival)
    if (ival .ne. 0) then
        call getvtx('VERIF', 'STOP_ERREUR', iocc=1, scal=rep, nbret=ibid)
        if (rep .eq. 'NON') goto 20
    end if
!
    call getfac('LIAISON', nblia)
!
    do i = 1, nblia
        call getvtx('LIAISON', 'OPTION', iocc=i, scal=option, nbret=iopt)
        call getvtx('LIAISON', 'SOUS_STRUC_1', iocc=i, scal=sst1, nbret=ibid)
        call getvtx('LIAISON', 'SOUS_STRUC_2', iocc=i, scal=sst2, nbret=ibid)
        call getvtx('LIAISON', 'INTERFACE_1', iocc=i, scal=intf1, nbret=ibid)
        call getvtx('LIAISON', 'INTERFACE_2', iocc=i, scal=intf2, nbret=ibid)
        iinc = 0
!     ON TESTE SI LA LIAISON EST INCOMPATIBLE
        call getvtx('LIAISON', 'GROUP_MA_MAIT_1', iocc=i, scal=k8bid, nbret=irep11)
        call getvtx('LIAISON', 'MAILLE_MAIT_1', iocc=i, scal=k8bid, nbret=irep12)
        call getvtx('LIAISON', 'GROUP_MA_MAIT_2', iocc=i, scal=k8bid, nbret=irep21)
        call getvtx('LIAISON', 'MAILLE_MAIT_2', iocc=i, scal=k8bid, nbret=irep22)
        if ((irep11 .ne. 0) .or. (irep12 .ne. 0)) then
            iinc = 1
        else if ((irep21 .ne. 0) .or. (irep22 .ne. 0)) then
            iinc = 2
        end if
!
!       SI ELLE EST COMPATIBLE ON VERIFIE LA COINCIDENCE DES NOEUDS
!       D'INTERFACE, SINON ON FAIT RIEN
        if ((iinc .eq. 0) .and. (option .eq. 'CLASSIQU')) then
            iret = i
            call vecomo(nomres, sst1, sst2, intf1, intf2, &
                        iret, option)
        end if
    end do
20  continue
!
!
!-----ORIENTATION DES MATRICES DE LIAISON
!
    call callis(nomres)
!
!-----VERIFICATION DU MODELE GENERALISE
!
!
end subroutine
