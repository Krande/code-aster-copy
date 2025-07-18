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

subroutine simult()
    implicit none
!
!     OPERATEUR :   CALC_CHAR_SEISME
!
!     CREE LE VECTEUR SECOND MEMBRE DANS LE CAS D'UN CALCUL SISMIQUE
!     STRUCTURE : MULTI-APPUI
!
!     ------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/compno.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/simul2.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
    real(kind=8) :: xnorm, depl(6)
    character(len=8) :: masse, modsta, mailla, nomnoe
    character(len=16) :: type, nomcmd
    character(len=19) :: resu
    character(len=24) :: magrno
    character(len=8) :: kbid
!     ------------------------------------------------------------------
!
!     --- RECUPERATION DES ARGUMENTS DE LA COMMANDE ---
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idno, ii, in, ldgn
    integer(kind=8) :: nb, nbd, nbdir, nbgr, nbno, nbv
    character(len=24), pointer :: group_no(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    magrno = ' '
    resu = ' '
    call getres(resu, type, nomcmd)
!
!     --- MATRICE DE MASSE ---
!
    call getvid(' ', 'MATR_MASS', scal=masse, nbret=nbv)
    call dismoi('NOM_MAILLA', masse, 'MATR_ASSE', repk=mailla)
!
!     --- QUELLE EST LA DIRECTION ? ---
!
    call getvr8(' ', 'DIRECTION', nbval=0, nbret=nbd)
    nbdir = -nbd
    call getvr8(' ', 'DIRECTION', nbval=nbdir, vect=depl, nbret=nbd)
!
!     --- ON NORMALISE LE VECTEUR ---
    xnorm = 0.d0
    do i = 1, nbdir
        xnorm = xnorm+depl(i)*depl(i)
    end do
    xnorm = sqrt(xnorm)
    if (xnorm .lt. 0.d0) then
        call utmess('F', 'ALGORITH9_81')
    end if
    do i = 1, nbdir
        depl(i) = depl(i)/xnorm
    end do
!
!     --- ON RECUPERE LES MODES STATIQUES ---
!
    call getvid(' ', 'MODE_STAT', scal=modsta, nbret=nbv)
!
!     --- ON RECUPERE LES POINTS D'ANCRAGE ---
!
    call getvem(mailla, 'GROUP_NO', ' ', 'GROUP_NO', 0, &
                0, kbid, nbgr)
    nbgr = -nbgr
    AS_ALLOCATE(vk24=group_no, size=nbgr)
    call getvem(mailla, 'GROUP_NO', ' ', 'GROUP_NO', 0, &
                nbgr, group_no, nbv)
!
!    --- ECLATE LE GROUP_NO EN NOEUD ---
    call compno(mailla, nbgr, group_no, nbno)
    call wkvect('&&SIMULT.NOEUD', 'V V K8', nbno, idno)
    magrno = mailla//'.GROUPENO'
    ii = -1
    do i = 1, nbgr
        call jelira(jexnom(magrno, group_no(i)), 'LONUTI', nb)
        call jeveuo(jexnom(magrno, group_no(i)), 'L', ldgn)
        do in = 0, nb-1
            nomnoe = int_to_char8(zi(ldgn+in))
            ii = ii+1
            zk8(idno+ii) = nomnoe
        end do
    end do
    call simul2(resu, nomcmd, masse, modsta, nbdir, &
                depl, zk8(idno), nbno)
!
! --- MENAGE
    call jedetr('&&SIMULT.NOEUD')
    AS_DEALLOCATE(vk24=group_no)
!
    call jedema()
end subroutine
