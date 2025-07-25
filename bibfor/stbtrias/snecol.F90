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
subroutine snecol(imod, nbnode)
    implicit none
!     =================
!A PRESUPER
!
!     ================================================================
!     !                                                              !
!     !  FONCTION: ECRITURE DES GROUPES DE NOEUDS ASSOCIES           !
!     !            AUX COULEURS                                      !
!     !                                                              !
!     ================================================================
!     !                                                              !
!     !  ROUTINES APPELES : CODENT                                   !
!     !                          : IUNIFI (FONCTION)                 !
!     !                          : CODNOP                            !
!     !                                                              !
!     !  ROUTINE APPELANTE : PRESUP                                  !
!     !                                                              !
!     ================================================================
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/codnop.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=1) :: prfnoe
    character(len=4) :: kbid
    character(len=8) :: chnode, chgrou
    aster_logical :: logiq(256)
    integer(kind=8) :: jpo(256), jnomb(256), jmax(256)
!  ------------ FIN DECLARATION -------------
!
!  -->N  D'UNITE LOGIQUE ASSOCIE AUX FICHIERS
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ic, icmax, icol, imod, inum, ipos
    integer(kind=8) :: j, nbmax, nbno, nbnode, nbtot
    integer(kind=8), pointer :: noeuds(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
!
    prfnoe = 'N'
    icmax = 256
    do i = 1, icmax
        logiq(i) = .false.
        jpo(i) = 0
        jnomb(i) = 0
        jmax(i) = 1000
    end do
!
    nbmax = 1000
    call jeveuo('&&PRESUP.INFO.NOEUDS', 'L', vi=noeuds)
    do i = 1, nbnode
        inum = noeuds((i-1)*3+1)
        call codnop(chnode, prfnoe, 1, 1)
        call codent(inum, 'G', chnode(2:8))
        icol = noeuds((i-1)*3+3)
        ipos = icol+1
        if (ipos .gt. icmax) then
            call utmess('A', 'STBTRIAS_2')
            goto 100
        end if
        if (.not. logiq(ipos)) then
            logiq(ipos) = .true.
            call codent(icol, 'G', kbid)
            call wkvect('&&PRESUP.COUL'//kbid, 'V V K8', nbmax+1, jpo(ipos))
        end if
        nbno = jnomb(ipos)
        nbtot = jmax(ipos)
        if (nbno .ge. nbtot) then
            call codent(icol, 'G', kbid)
            nbtot = nbtot+nbmax
            jmax(ipos) = nbtot
            call juveca('&&PRESUP.COUL'//kbid, nbtot+1)
            call jeveuo('&&PRESUP.COUL'//kbid, 'E', jpo(ipos))
        end if
        jnomb(ipos) = nbno+1
        zk8(jpo(ipos)-1+nbno+1) = chnode
100     continue
    end do
!
! --> ECRITURE DES GROUPES DE NOEUDS PAR COULEUR
!
    do ic = 1, icmax
        if (logiq(ic)) then
            call codent((ic-1), 'G', kbid)
            chgrou = 'COUL_'//kbid
            write (imod, '(A,4X,2A)') 'GROUP_NO', 'NOM=', chgrou
            nbno = jnomb(ic)
            write (imod, '(8(2X,A))') (zk8(jpo(ic)-1+j), j=1, nbno)
            write (imod, '(A)') 'FINSF'
            write (imod, '(A)') '%'
            call jedetr('&&PRESUP.COUL'//kbid)
        end if
    end do
!
    call jedema()
end subroutine
