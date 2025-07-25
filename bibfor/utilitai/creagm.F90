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
subroutine creagm(nbmato, nbpart, ma, masd)
!
!    - FONCTION REALISEE:
!       - CREATION DE NOUVEAUX GROUPES DE MAILLES
!              DANS LA STRUCTURE MAILLAGE
!
!    - IN :     NBMATO : NOMBRE DE MAILLES
!               NBPART : NOMBRE DE PARTITION
!               RENUM  : RENUMEROTATION
!               NBMASD : NOMBRE DE MAILLES PAR SOUS DOMAINE
!               MA     : NOM DU MAILLAGE
!               NUMSDM : SOUS DOMAINES DE CHAQUES MAILLE
!
!    - OUT :    MASD   : DONNE LES MAILLES PAR SD
!               IDMASD : INDEX DE MASD
!               NOMSDM : NOM DES GROUP_MA
!
!----------------------------------------------------------------------
! person_in_charge: mathieu.courtois@edf.fr
!
! CORPS DU PROGRAMME
    implicit none
!
!
! DECLARATION VARIABLES D'APPEL
#include "jeveux.h"
#include "asterfort/infniv.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/lxcadr.h"
#include "asterfort/uttcpr.h"
#include "asterfort/uttcpu.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbmato, nbpart, renum, nbmasd, nomsdm, masd, idmasd, numsdm
    character(len=8) :: ma
!
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: id1, isd, nbma, ima, nbre, idma, ifm, niv
    integer(kind=8) :: nbgrsd
    real(kind=8) :: tmps(7)
    character(len=8) :: ktmp
    character(len=24) :: grpema
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    call jemarq()
    call infniv(ifm, niv)
!
    call jeveuo('&&FETSKP.RENUM', 'L', renum)
    call jeveuo('&&FETSKP.NBMASD', 'L', nbmasd)
    call jeveuo('&&FETSKP.NUMSDM', 'E', numsdm)
!
    if (niv .ge. 2) then
        call uttcpu('CPU.CREAGM', 'INIT', ' ')
        call uttcpu('CPU.CREAGM', 'DEBUT', ' ')
    end if
!
    call wkvect('&&FETSKP.ID1', 'V V I', nbpart, id1)
    call wkvect('&&FETSKP.IDMASD', 'V V I', nbpart+1, idmasd)
    call wkvect('&&FETSKP.NOMSDM', 'V V K24', nbpart, nomsdm)
    call wkvect('&&FETSKP.MASD', 'V V I', nbmato, masd)
!
! ------- On RECUPERE LE NOM DES SOUS DOMAINES
!
    do isd = 1, nbpart
        write (ktmp, '(I4)') isd-1
        call lxcadr(ktmp)
        zk24(nomsdm-1+isd) = 'SD'//ktmp
    end do
!
!
! ------- CREATION DU TABLEAU DES PLACES DES GRPMA
!
    zi(idmasd) = 1
    do isd = 2, nbpart+1
        zi(idmasd-1+isd) = zi(idmasd-1+isd-1)+zi(nbmasd-1+isd-1)
    end do
!
! ------- CREATION DU TABLEAU DONNANT LES MAILLES POUR CHAQUE GRMPA
!
    do ima = 1, nbmato
        nbre = zi(numsdm-1+ima)+1
        zi(masd-1+zi(idmasd-1+nbre)+zi(id1-1+nbre)) = zi(renum-1+ima)
        zi(id1-1+nbre) = zi(id1-1+nbre)+1
    end do
!
!
!   -- on cree un "GROUPEMA" par sous-domaine :
    nbgrsd = 0
    do isd = 1, nbpart
        nbma = zi(nbmasd-1+isd)
        if (nbma .gt. 0) nbgrsd = nbgrsd+1
    end do
    call jecrec('&&FETCRF.GROUPEMA', 'V V I', 'NO', 'DISPERSE', 'VARIABLE', &
                nbgrsd)
!
    do isd = 1, nbpart
        grpema = zk24(nomsdm-1+isd)
        nbma = zi(nbmasd-1+isd)
        if (nbma .gt. 0) then
            call jecroc(jexnom('&&FETCRF.GROUPEMA', grpema))
            call jeecra(jexnom('&&FETCRF.GROUPEMA', grpema), 'LONMAX', nbma)
            call jeecra(jexnom('&&FETCRF.GROUPEMA', grpema), 'LONUTI', nbma)
            call jeveuo(jexnom('&&FETCRF.GROUPEMA', grpema), 'E', idma)
            do ima = 0, nbma-1
                zi(idma+ima) = zi(masd+zi(idmasd-1+isd)-1+ima)
            end do
        end if
    end do
!
    call jedetr('&&FETSKP.ID1')
!
    if (niv .ge. 2) then
        call uttcpu('CPU.CREAGM', 'FIN', ' ')
        call uttcpr('CPU.CREAGM', 7, tmps)
        write (ifm, *) '--- CREATION DES GRPMA :', tmps(7)
        write (ifm, *) '  '
    end if
!
    call jedema()
end subroutine
