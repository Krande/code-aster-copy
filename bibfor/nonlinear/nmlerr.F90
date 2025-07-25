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

subroutine nmlerr(sddisc, action, infz, valr, vali)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=19) :: sddisc
    character(len=1) :: action
    character(len=*) :: infz
    integer(kind=8) :: vali
    real(kind=8) :: valr
!
! ----------------------------------------------------------------------
!
! ROUTINE *_NON_LINE (STRUCTURES DE DONNES - SD DISCRETISATION)
!
! LECTURE/ECRITURE DANS SD STOCKAGE DES INFOS EN COURS DE CALCUL
!
! ----------------------------------------------------------------------
!
! IN  SDDISC : SD DISCRETISATION
! IN  ACTION : 'L' OU 'E'
! IN  INFO   : TYPE D'INFO A STOCKER OU A LIRE
!   1 MXITER               : MAX( ITER_GLOB_MAXI , ITER_GLOB_ELAS )
!   2 MNITER               : MIN( ITER_GLOB_MAXI , ITER_GLOB_ELAS )
!   3 NBITER               : NOMBRE MAX ITERATIONS (Y COMPRIS EXTRAPOL)
!   4 PAS_MINI_ELAS        : PAS_MINI_ELAS
!   5 RESI_GLOB_RELA       : RESI_GLOB_RELA DONNE
!   6 RESI_GLOB_MAXI       : RESI_GLOB_MAXI
!   7 TYPE_RESI            :  =1 ON A DONNE RESI_GLOB_RELA
!                             =2 ON A DONNE RESI_GLOB_MAXI
!                             =3 C'EST (1) ET (2)
!                             =0 ON A RIEN DONNE ==> =1 (DEFAUT)
!   8 INIT_NEWTON_KRYLOV   : RESIDU INITIAL POUR NEWTON KRYLOV
!   9 ITER_NEWTON_KRYLOV   : RESIDU COURANT POUR NEWTON KRYLOV
!  10 ITERSUP              : =3 ON AUTORISE DES ITERATIONS EN PLUS
!
! I/O VALR      : REEL   A ECRIRE OU A LIRE
!     VALI      : ENTIER A ECRIRE OU A LIRE
!
!
!
!
    character(len=24) :: infocv
    integer(kind=8) :: jifcv
    character(len=24) :: info
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    infocv = sddisc(1:19)//'.IFCV'
    call jeveuo(infocv, 'E', jifcv)
    info = infz
!
    ASSERT((action .eq. 'E') .or. (action .eq. 'L'))
!
    if (info .eq. 'MXITER') then
        if (action .eq. 'L') then
            vali = nint(zr(jifcv+1-1))
        else if (action .eq. 'E') then
            zr(jifcv+1-1) = vali
        end if
    else if (info .eq. 'MNITER') then
        if (action .eq. 'L') then
            vali = nint(zr(jifcv+2-1))
        else if (action .eq. 'E') then
            zr(jifcv+2-1) = vali
        end if
    else if (info .eq. 'NBITER') then
        if (action .eq. 'L') then
            vali = nint(zr(jifcv+3-1))
        else if (action .eq. 'E') then
            zr(jifcv+3-1) = vali
        end if
    else if (info .eq. 'PAS_MINI_ELAS') then
        if (action .eq. 'L') then
            valr = zr(jifcv+4-1)
        else if (action .eq. 'E') then
            zr(jifcv+4-1) = valr
        end if
    else if (info .eq. 'RESI_GLOB_RELA') then
        if (action .eq. 'L') then
            valr = zr(jifcv+5-1)
        else if (action .eq. 'E') then
            zr(jifcv+5-1) = valr
        end if
    else if (info .eq. 'RESI_GLOB_MAXI') then
        if (action .eq. 'L') then
            valr = zr(jifcv+6-1)
        else if (action .eq. 'E') then
            zr(jifcv+6-1) = valr
        end if
    else if (info .eq. 'TYPE_RESI') then
        if (action .eq. 'L') then
            vali = nint(zr(jifcv+7-1))
        else if (action .eq. 'E') then
            zr(jifcv+7-1) = vali
        end if
    else if (info .eq. 'INIT_NEWTON_KRYLOV') then
        if (action .eq. 'L') then
            valr = zr(jifcv+8-1)
        else if (action .eq. 'E') then
            zr(jifcv+8-1) = valr
        end if
    else if (info .eq. 'ITER_NEWTON_KRYLOV') then
        if (action .eq. 'L') then
            valr = zr(jifcv+9-1)
        else if (action .eq. 'E') then
            zr(jifcv+9-1) = valr
        end if
    else if (info .eq. 'ITERSUP') then
        if (action .eq. 'L') then
            vali = nint(zr(jifcv+10-1))
        else if (action .eq. 'E') then
            zr(jifcv+10-1) = vali
        end if
    else
        ASSERT(.false.)
    end if
!
    call jedema()
!
end subroutine
