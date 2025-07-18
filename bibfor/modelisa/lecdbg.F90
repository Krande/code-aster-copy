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

subroutine lecdbg(ifl, icl, iv, rv, cv, &
                  cnl, mcl, nbm, nbg, dim, &
                  nob, irteti)
    implicit none
!       PREMIERE LECTURE DES DONNEES POUR UN MOT CLE DE TYPE GROUPE
!       ----------------------------------------------------------------
!       IN      IFL,ICL,IV,RV,CV,CNL = VOIR LIRITM
!               MCL             = MOTS CLE TYPE DEBUG
!               NBM             = NB DE MOTS CLES TYPE DEBUG
!                               = 1 > ERREUR EN LECTURE
!               DIM             = NB DE NOMS LUS PAR MOT CLE DEBUG
!               NOB             = NOMS LUS
!               NBG             = NIVEAU DEBUG
!               (RETURN)        = MOT CLE SUIVANT (MOT CLE NON RECONNU)
!               (RETURN 1)      = EXIT            (MOT CLE FIN TROUVE)
!               (RETURN 2)      = LIGNE SUIVANTE  (MOT CLE FINSF TROUVE
!                                                  OU ERREUR DETECTE)
!       ----------------------------------------------------------------
!
#include "asterfort/iunifi.h"
#include "asterfort/liritm.h"
#include "asterfort/tesfin.h"
#include "asterfort/tesmcl.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbm
    real(kind=8) :: rv
    character(len=8) :: mcl(nbm)
    integer(kind=8) :: dim(nbm), deblig
    character(len=14) :: cnl
    character(len=24) :: nob(50, nbm), b24, mtc
    character(len=*) :: cv
    save b24
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icl, ifl, ifm, irtet, irteti
    integer(kind=8) :: iv, nbg, numtcl
!-----------------------------------------------------------------------
    data b24/'                        '/
    irteti = 0
!
    ifm = iunifi('MESSAGE')
!
! ----- ITEM = MOT CLE TYPE  DEBUG ?
!
    do i = 1, nbm
        call tesmcl(icl, iv, cv, mcl(i), irtet)
        if (irtet .eq. 1) goto 4
        numtcl = i
        goto 5
4       continue
    end do
    goto 3
!
! ----- LECTURE DES NOMS D OBJETS A DUMPER ?
!
5   continue
    write (ifm, *) ' ----- LECDBG'
6   continue
    call liritm(ifl, icl, iv, rv, cv, &
                cnl, deblig, 1)
    write (ifm, *) '       LIRITM : ICL = ', icl,&
     &  ' IV = ', iv, ' RV = ', rv, ' CV(1:8) = ', cv(1:8), ' DEBLIG =', deblig
    if (deblig .eq. 1) then
        call tesfin(icl, iv, cv, irtet)
        if (irtet .eq. 1) then
            goto 1
        else if (irtet .eq. 2) then
            goto 2
        end if
    end if
!
! - MOT CLE DUMP
!
    if (numtcl .eq. 1) then
        if (icl .ne. 4 .and. icl .ne. 3) then
            call utmess('F', 'MODELISA4_78', sk=mcl(1))
        else if (iv .gt. 24) then
            call utmess('F', 'MODELISA4_79', sk=mcl(1))
        end if
        dim(1) = dim(1)+1
        mtc = b24
        mtc(1:iv) = cv(1:iv)
        nob(dim(1), 1) = mtc
        goto 6
    end if
!
! - MOT CLE DEBUG
!
    if (numtcl .eq. 2) then
        if (icl .ne. 1 .and. icl .ne. 2) then
            call utmess('F', 'MODELISA4_78', sk=mcl(2))
        end if
        if (icl .eq. 1) nbg = iv
        if (icl .eq. 2) nbg = nint(rv)
        if (nbg .gt. 1) nbg = 1
        write (ifm, *) ' ------------ DEBUG NIVEAU ', nbg, ' --------------'
        goto 6
    end if
!
1   continue
    irteti = 1
    goto 999
2   continue
    irteti = 2
    goto 999
3   continue
    irteti = 0
    goto 999
!
999 continue
end subroutine
