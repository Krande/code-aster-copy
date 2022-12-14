! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine utims3(comm, sch1, ipos, base)
    implicit none
!     --
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/dbgobj.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelstc.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/uttr24.h"
!
    character(len=*) :: comm, sch1, base
    integer :: ipos
! ----------------------------------------------------------------------
! BUT:
!   IMPRIMER SUR LE FICHIER 'MESSAGE' LE RESUME DES OBJETS JEVEUX (K24)
!   AYANT LA CHAINE SCH1 EN POSITION IPOS DANS LEURS NOMS.
!
!
! IN:
!   COMM   : COMMENTAIRE DE DEBUT DE CHAQUE LIGNE (K8)
!   SCH1   : CHAINE DE CARACTERES CHERCHEE
!   IPOS   : DEBUT DE LA CHAINE DE CARACTERES A CHERCHER
!   BASE   : 'G','V','L',OU ' '(TOUTES)
!
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=24) :: ob1, chain2, ficou
    character(len=8) :: kbid, comm2
    character(len=1) :: bas2
    integer :: long, ifm, nbval, nbobj, i, iret
    character(len=24), pointer :: liste(:) => null()
!
!
    call jemarq()
    bas2 = base
    comm2 = comm
!
!
!     --QUELQUES VERIFICATIONS:
!     -------------------------
    long = len(sch1)
    if (len(sch1) .gt. 24) then
        call utmess('F', 'UTILITAI5_42')
    endif
    if ((ipos.lt.0) .or. (ipos.gt.24)) then
        call utmess('F', 'UTILITAI5_43')
    endif
    if (ipos+len(sch1) .gt. 25) then
        call utmess('F', 'UTILITAI5_44')
    endif
!
!
!     -- DETERMINATION DU NOMBRE DES OBJETS TROUVES :
!    ------------------------------------------------
    if (long .eq. 24) then
        call jeexin(sch1, iret)
        if (iret .gt. 0) then
            nbobj = 1
        else
            nbobj = 0
        endif
    else
        call jelstc(bas2, sch1, ipos, 0, kbid,&
                    nbval)
        nbobj = -nbval
    endif
!
!
!     -- ECRITURE DE L'ENTETE :
!    --------------------------
    ficou = 'MESSAGE'
    ifm = iunifi(ficou)
    chain2 = '????????????????????????'
    chain2(ipos:ipos-1+long) = sch1
    write (ifm,*) ' '
    write (ifm,*) '#AJ1 ====> UTIMS3 DE LA STRUCTURE DE DONNEE : ',&
     &  chain2
    write (ifm,*) '#AJ1 NOMBRE D''OBJETS (OU COLL.) TROUVES :',nbobj
    if (nbobj .eq. 0) goto 20
!
!
!     -- RECHERCHE DES NOMS DES OBJETS VERIFIANT LE CRITERE:
!    -------------------------------------------------------
    AS_ALLOCATE(vk24=liste, size=nbobj)
    if (long .eq. 24) then
        liste(1) = sch1
    else
        call jelstc(bas2, sch1, ipos, nbobj, liste,&
                    nbval)
    endif
!
!     -- ON TRIE PAR ORDRE ALPHABETIQUE:
    call uttr24(liste, nbobj)
!
!
!     -- ECRITURE D'UNE LIGNE D'INFO POUR CHAQUE OBJET
!    -------------------------------------------------------
    write (ifm,1001) comm2,'#AJ2>','NOMOBJ','<','LONUTI',&
     &  'LONMAX','TYPE','IRET','SOMMI','RESUME','SOMMR'
!
    do i = 1, nbobj
        ob1 = liste(i)
        call dbgobj(ob1, 'OUI', ifm, '&&UTIMS3')
    end do
!
!
    AS_DEALLOCATE(vk24=liste)
 20 continue
    call jedema()
!
    1001 format (a8,a4,a5,a24,a1,a8,a8,a5,a5,a10,a10,a15)
end subroutine
