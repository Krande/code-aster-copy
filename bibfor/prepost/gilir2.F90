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
subroutine gilir2(nfic, niv, ndim, nbobo)
! aslint: disable=
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gicnx2.h"
#include "asterfort/gidoma.h"
#include "asterfort/gilig0.h"
#include "asterfort/gilig1.h"
#include "asterfort/gilig2.h"
#include "asterfort/gilig3.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nfic, niv, ndim, nbobo
!
!     BUT: LIRE LE FICHIER DE MAILLAGE GIBI (PROCEDURE SAUVER) :
!
!     IN : NFIC  : UNITE DE LECTURE
!          NIV   : NUMERO DU NIVEAU GIBI
!     OUT: NDIM  : DIMENSION DU PROBLEME (2D OU 3D)
!          NBOBO : NOMBRE D'OBJETS (AU SENS GIBI)
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: r8bid
    integer(kind=8) :: nbobno, ipile, nivo, nberr, nboblu
    character(len=1) :: ityp
    character(len=4) :: k4bid, kbid4
    character(len=6) :: k6bid
    character(len=14) :: kbid14
    aster_logical :: legrno
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iaptin, iret, nbnoto, nbval
!-----------------------------------------------------------------------
    call jemarq()
!
    legrno = .false.
1   continue
    read (nfic, 1001, end=9997) kbid14, kbid4, ityp
!
    if (kbid14 .eq. 'ENREGISTREMENT' .and. kbid4 .eq. 'TYPE') then
!
        if (ityp .eq. '4') then
!
! -- INFORMATIONS GENERALES MAILLAGE ----
!
            read (nfic, 1002) nivo, nberr, ndim
            if (nberr .gt. 0) then
                call utmess('A', 'PREPOST_59')
            end if
            read (nfic, 1003) r8bid
            goto 1
!
        else if (ityp .eq. '7') then
!
! -- INFORMATIONS GENERALES CASTEM 2000 ----
!
            read (nfic, 1004)
            read (nfic, 1004)
            goto 1
!
        else if (ityp .eq. '5') then
!
!          -- ON A TOUT LU ----
!
            goto 9997
!
        else if (ityp .eq. '2') then
!
! -- LECTURE D'UNE PILE  ----
!
            if (niv .le. 6) then
                read (nfic, 1005) k4bid, k6bid, ipile, nbobno, nboblu
            else if (niv .gt. 6) then
                read (nfic, 1006) k4bid, k6bid, ipile, nbobno, nboblu
            end if
!
!
            if (ipile .eq. 0) then
!
!            --- LECTURE DES GROUPES DE NOEUDS NOMMES ---
                legrno = .true.
                call gilig2(nfic, nbobno, niv)
!
!            --- LECTURE DES COORDONNEES ---
                nbval = nboblu*(ndim+1)
                call gilig1(nfic, ndim, nbval, nboblu)
!
                nbnoto = nboblu
!
                call jeexin('&&GILIRE.INDIRECT', iret)
                if (iret .eq. 0) then
                    call wkvect('&&GILIRE.INDIRECT', 'V V I', nboblu, iaptin)
                    do i = 1, nboblu
                        zi(iaptin+i-1) = i
                    end do
                end if
                goto 1
!
            else if (ipile .eq. 32) then
!
!            --- LECTURE DES GROUPES DE NOEUDS NOMMES ---
                if (legrno) goto 1
                call gilig3(nfic, nbobno, niv, nboblu)
                goto 1
!
            else if (ipile .eq. 33) then
!
                read (nfic, 1010) nbval
!
!            --- LECTURE DES COORDONNEES ---
                nboblu = nbval/(ndim+1)
                call gilig1(nfic, ndim, nbval, nboblu)
                nbnoto = nboblu
!
                goto 1
!
            else if (ipile .eq. 1) then
!
!            --- LECTURE DES GROUPES DE MAILLES NOMMEES ---
                call gilig0(nfic, nboblu, nbobno, nbobo, niv)
                goto 1
            end if
!
        end if
        goto 1
!
    else
        goto 1
    end if
!
9997 continue
!
!     -- ON CREE .CONNEX2:
    call gicnx2()
!
!     -- ON CREE .NUMANEW:
    call gidoma(nbnoto)
!
!
!
1001 format(1x, a14, 4x, a4, 3x, a1)
1002 format(7x, i4, 14x, i4, 10x, i4)
1003 format(8x, d12.5)
1004 format(10x)
1005 format(1x, a4, 1x, a6, i4, 18x, i5, 11x, i5)
1006 format(1x, a4, 1x, a6, i4, 18x, i8, 11x, i8)
1010 format(1x, i7)
!
    call jedema()
end subroutine
