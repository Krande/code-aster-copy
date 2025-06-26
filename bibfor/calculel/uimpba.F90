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

subroutine uimpba(clas, iunmes)
    implicit none
#include "jeveux.h"
!
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/gettco.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelgdq.h"
#include "asterfort/jelira.h"
#include "asterfort/jelstc.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/split_string.h"
    character(len=*) :: clas
    integer(kind=8) :: iunmes
! person_in_charge: jacques.pellet at edf.fr
! ----------------------------------------------------------------------
! BUT:
!   IMPRIMER LA TAILLE DES CONCEPTS STOCKES SUR UNE BASE
!
!  IN    CLAS  : NOM DE LA BASE : 'G', 'V', ..(' ' -> TOUTES LES BASES)
! ----------------------------------------------------------------------
    character(len=8) :: root
    character(len=16) :: typcon
    character(len=19) :: key
    character(len=24) :: kbid, obj
    real(kind=8) :: rlong, mega, taitot
    integer(kind=8) :: i, nbobj, nbval, idx, nbcon, nbsv
    integer(kind=8) ::     nstot
    character(len=24), pointer :: liste_obj(:) => null()
    character(len=24), pointer :: types(:) => null()
    integer(kind=8), pointer :: vnbobj(:) => null()
    integer(kind=8), pointer :: nbsvc(:) => null()
    integer(kind=8), pointer :: nbsvo(:) => null()
    real(kind=8), pointer :: tailcon(:) => null()
    real(kind=8), pointer :: taille(:) => null()
!
!
    mega = 1024*1024
!
!
!     -- 1 : NBOBJ + .LISTE_OBJ :LISTE DES OBJETS :
!     ----------------------------------------------
    nbobj = 0
    call jelstc(clas, ' ', 0, nbobj, kbid, &
                nbval)
    ASSERT(nbval .le. 0)
    if (nbval .eq. 0) goto 999
    nbobj = -nbval
    AS_ALLOCATE(vk24=liste_obj, size=nbobj+1)
    call jelstc(clas, ' ', 0, nbobj, liste_obj, &
                nbval)
!     NBVAL = NBOBJ (+1 EVENTUELLEMENT A CAUSE DE '&&UIMPBA.LISTE_OBJ')
    ASSERT(nbval .eq. nbobj+1 .or. nbval .eq. nbobj)
    nbobj = nbval
!
!     -- 2 : .TAILLE = TAILLE DES OBJETS :
!     --------------------------------------
    AS_ALLOCATE(vr=taille, size=nbobj)
    AS_ALLOCATE(vi=nbsvo, size=nbobj)
    do i = 1, nbobj
        obj = liste_obj(i)
        call jelgdq(obj, rlong, nbsv)
        ASSERT(rlong .gt. 0.d0)
        taille(i) = rlong
        nbsvo(i) = nbsv
    end do
!
!
!     -- 3 : .LSTCON = LISTE DES CONCEPTS (K8) DE .LISTE_OBJ
!     -----------------------------------------------------------
    call jecreo('&&UIMPBA.LSTCON', 'V N K24')
    call jeecra('&&UIMPBA.LSTCON', 'NOMMAX', nbobj)
    do i = 1, nbobj
        obj = liste_obj(i)
        call split_string(obj, ".", root)
        call jenonu(jexnom('&&UIMPBA.LSTCON', root), idx)
        if (idx .eq. 0) then
            call jecroc(jexnom('&&UIMPBA.LSTCON', root))
        end if
    end do
    !
    !
    !     -- 4 : .TAILCON = TAILLE DES CONCEPTS
    !     -----------------------------------------------------------
    call jelira('&&UIMPBA.LSTCON', 'NOMUTI', nbcon)
    AS_ALLOCATE(vr=tailcon, size=nbcon)
    AS_ALLOCATE(vi=nbsvc, size=nbcon)
    AS_ALLOCATE(vi=vnbobj, size=nbcon)
    AS_ALLOCATE(vk24=types, size=nbcon)
    types = " "
    taitot = 0.d0
    nstot = 0
    do i = 1, nbobj
        obj = liste_obj(i)
        call split_string(obj, ".", root)
        call jenonu(jexnom('&&UIMPBA.LSTCON', root), idx)
        ASSERT(idx .gt. 0)
        ASSERT(idx .le. nbcon)
        if (types(idx) .eq. " ") then
            call gettco(root, types(idx), errstop=ASTER_FALSE)
            if (types(idx) .eq. " ") then
                call gettco(obj(1:19), typcon, errstop=ASTER_FALSE)
                if (typcon .ne. " ") then
                    types(idx) = "*"//typcon
                end if
            end if
        end if
        tailcon(idx) = tailcon(idx)+taille(i)
        taitot = taitot+taille(i)
        nbsvc(idx) = nbsvc(idx)+nbsvo(i)
        vnbobj(idx) = vnbobj(idx)+1
        nstot = nstot+nbsvo(i)
    end do
!
!
!     -- 5 : IMPRESSION DU RESULTAT :
!     -----------------------------------------------------------
    write (iunmes, *) '-----------------------------------------------',&
     &                '----------------------------'
    write (iunmes, *) 'Concepts de la base: ', clas
    write (iunmes, *) '   Nom        Type                Taille (Mo)',&
     &                '         Nombre      Nombre de'
    write (iunmes, *) '                                            ',&
     &                '        d''objets       segments'
!
    write (iunmes, 100) 'TOTAL   ', ' ', taitot/mega, nbobj, nstot
    write (iunmes, *) ' '
!
!     -- ON IMPRIME D'ABORD LES CONCEPTS UTILISATEUR :
    do i = 1, nbcon
        if (types(i) .eq. ' ') cycle
        call jenuno(jexnum('&&UIMPBA.LSTCON', i), key)
        write (iunmes, 100) key, types(i), tailcon(i)/mega, vnbobj(i), nbsvc(i)
    end do

!     -- ON IMPRIME ENSUITE LES CONCEPTS CACHES  :
    do i = 1, nbcon
        if (types(i) .ne. ' ') cycle
        call jenuno(jexnum('&&UIMPBA.LSTCON', i), key)
        write (iunmes, 100) key, types(i), tailcon(i)/mega, vnbobj(i), nbsvc(i)
    end do
    write (iunmes, *) '(*) sous-objets non accessible directement'
    write (iunmes, *) '-----------------------------------------------',&
     &                '----------------------------'
!
!
999 continue
    AS_DEALLOCATE(vk24=liste_obj)
    AS_DEALLOCATE(vr=taille)
    call jedetr('&&UIMPBA.LSTCON')
    AS_DEALLOCATE(vr=tailcon)
    AS_DEALLOCATE(vi=nbsvo)
    AS_DEALLOCATE(vi=nbsvc)
    AS_DEALLOCATE(vi=vnbobj)
100 format(4x, a8, 3x, a16, 3x, f12.2, 3x, i12, 3x, i12)
end subroutine
