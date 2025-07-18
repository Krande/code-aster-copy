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
subroutine dismma(questi, nomobz, repi, repkz, ierd)
    implicit none
!     --     DISMOI(MAILLAGE)
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/gettco.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/ltnotb.h"
#include "asterfort/tbliva.h"
!
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi
    character(len=32) :: repk
    character(len=8) :: nomob
    character(len=*) :: nomobz, repkz
! ----------------------------------------------------------------------
!    IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOB  : NOM D'UN OBJET DE TYPE MAILLAGE
!    OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPK   : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
!
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    complex(kind=8) :: c16b
    character(len=19) :: table
    character(len=16) :: typeco
    character(len=1) :: k1bid
    real(kind=8) :: xmax, xmin, ymax, ymin, zmax, zmin
    integer(kind=8) :: ibid, ier, ilmaco, ism, k, nbma, nbno
    integer(kind=8) :: nbsm, nno, typv
    character(len=8) :: typma
    integer(kind=8), pointer :: dime(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
!
!
!
    call jemarq()
    repk = ' '
    repi = 0
    ierd = 0
!
    nomob = nomobz
    call jeveuo(nomob//'.DIME', 'L', vi=dime)
!
!
    if (questi .eq. 'NB_MA_MAILLA') then
!     ---------------------------------
        repi = dime(3)
!
!
    else if (questi .eq. 'NB_SM_MAILLA') then
!     ---------------------------------
        repi = dime(4)
!
!
    else if (questi .eq. 'NB_NO_MAILLA') then
!     ---------------------------------
        repi = dime(1)
!
!
    else if (questi .eq. 'NB_NL_MAILLA') then
!     ---------------------------------
        repi = dime(2)
!
!
    else if (questi .eq. 'PARALLEL_MESH') then
!     ---------------------------------
        call gettco(nomob, typeco)
        if (typeco .eq. 'MAILLAGE_P') then
            repk = 'OUI'
        else
            repk = 'NON'
        end if
!
!
    else if (questi .eq. 'NB_NO_SS_MAX') then
!     ---------------------------------
        nbsm = dime(4)
        repi = 0
        do ism = 1, nbsm
            call jelira(jexnum(nomob//'.SUPMAIL', ism), 'LONMAX', nno)
            repi = max(repi, nno)
        end do
!
!
    else if (questi .eq. 'Z_CST') then
!     ---------------------------------
        call ltnotb(nomob, 'CARA_GEOM', table)
        call tbliva(table, 0, ' ', [ibid], [0.d0], &
                    [c16b], k1bid, 'ABSO', [0.d0], 'Z_MIN', &
                    k1bid, ibid, zmin, c16b, k1bid, &
                    ier)
        call tbliva(table, 0, ' ', [ibid], [0.d0], &
                    [c16b], k1bid, 'ABSO', [0.d0], 'Z_MAX', &
                    k1bid, ibid, zmax, c16b, k1bid, &
                    ier)
!
        if (zmin .eq. zmax) then
            repk = 'OUI'
        else
            repk = 'NON'
        end if
!
!
    else if (questi .eq. 'Z_ZERO' .or. questi .eq. 'Z_QUASI_ZERO' .or. questi .eq. 'DIM_GEOM') then
!     -----------------------------------------------------------------------------------
        call ltnotb(nomob, 'CARA_GEOM', table)
        call tbliva(table, 0, ' ', [ibid], [0.d0], &
                    [c16b], k1bid, 'ABSO', [0.d0], 'Z_MIN', &
                    k1bid, ibid, zmin, c16b, k1bid, &
                    ier)
        call tbliva(table, 0, ' ', [ibid], [0.d0], &
                    [c16b], k1bid, 'ABSO', [0.d0], 'Z_MAX', &
                    k1bid, ibid, zmax, c16b, k1bid, &
                    ier)
!
        if (zmin .eq. zmax .and. zmin .eq. 0.d0) then
            repk = 'OUI'
        else
            repk = 'NON'
        end if
!
        if (repk == 'NON' .and. questi == 'Z_QUASI_ZERO') then
!       ------------------------------------
            call tbliva(table, 0, ' ', [ibid], [0.d0], &
                        [c16b], k1bid, 'ABSO', [0.d0], 'X_MIN', &
                        k1bid, ibid, xmin, c16b, k1bid, &
                        ier)
            call tbliva(table, 0, ' ', [ibid], [0.d0], &
                        [c16b], k1bid, 'ABSO', [0.d0], 'X_MAX', &
                        k1bid, ibid, xmax, c16b, k1bid, &
                        ier)
            call tbliva(table, 0, ' ', [ibid], [0.d0], &
                        [c16b], k1bid, 'ABSO', [0.d0], 'Y_MIN', &
                        k1bid, ibid, ymin, c16b, k1bid, &
                        ier)
            call tbliva(table, 0, ' ', [ibid], [0.d0], &
                        [c16b], k1bid, 'ABSO', [0.d0], 'Y_MAX', &
                        k1bid, ibid, ymax, c16b, k1bid, &
                        ier)
            if (max(abs(zmax), abs(zmin)) .lt. max(xmax-xmin, ymax-ymin)*1.e-8) then
                repk = 'OUI'
            else
                repk = 'NON'
            end if
            repi = 0
        end if
!
        if (questi .eq. 'DIM_GEOM') then
!       --------------------------------
            repi = dime(6)
!          -- ON RETOURNE 2 SI Z=0. PARTOUT :
            if ((repi .eq. 3) .and. (repk .eq. 'OUI')) then
                repi = 2
            end if
            repk = '???'
        end if
!
!
    else if (questi .eq. 'DIM_GEOM_B') then
!     ----------------------------------------
        repi = dime(6)
        repk = '???'
!
!
    else if (questi .eq. 'NB_NO_MA_MAX') then
!     ----------------------------------------
        nbma = dime(3)
        call jeveuo(jexatr(nomob//'.CONNEX', 'LONCUM'), 'L', ilmaco)
        repi = 0
        do k = 1, nbma
            nbno = zi(ilmaco+k)-zi(ilmaco-1+k)
            repi = max(repi, nbno)
        end do
!
!
    else if ((questi .eq. 'EXI_TRIA3') .or. (questi .eq. 'EXI_TRIA6') &
             .or. (questi .eq. 'EXI_QUAD4') .or. &
             (questi .eq. 'EXI_QUAD8') .or. (questi .eq. 'EXI_QUAD9') &
             .or. (questi .eq. 'EXI_SEG2') .or. (questi .eq. 'EXI_SEG3') .or. &
             (questi .eq. 'EXI_HEXA8') .or. (questi .eq. 'EXI_HEXA20') &
             .or. (questi .eq. 'EXI_HEXA27') &
             .or. (questi .eq. 'EXI_PENTA6') .or. (questi .eq. 'EXI_PENTA15') &
             .or. (questi .eq. 'EXI_TETRA4') .or. (questi .eq. 'EXI_PYRAM13') &
             .or. (questi .eq. 'EXI_POI1')) then
!     ----------------------------------------
        typma = 'XXXX'
        if (questi .eq. 'EXI_TRIA3') typma = 'TRIA3'
        if (questi .eq. 'EXI_TRIA6') typma = 'TRIA6'
        if (questi .eq. 'EXI_QUAD4') typma = 'QUAD4'
        if (questi .eq. 'EXI_QUAD8') typma = 'QUAD8'
        if (questi .eq. 'EXI_QUAD9') typma = 'QUAD9'
        if (questi .eq. 'EXI_SEG2') typma = 'SEG2'
        if (questi .eq. 'EXI_SEG3') typma = 'SEG3'
        if (questi .eq. 'EXI_HEXA8') typma = 'HEXA8'
        if (questi .eq. 'EXI_HEXA20') typma = 'HEXA20'
        if (questi .eq. 'EXI_HEXA27') typma = 'HEXA27'
        if (questi .eq. 'EXI_PENTA6') typma = 'PENTA6'
        if (questi .eq. 'EXI_PENTA15') typma = 'PENTA15'
        if (questi .eq. 'EXI_TETRA4') typma = 'TETRA4'
        if (questi .eq. 'EXI_TETRA10') typma = 'TETRA10'
        if (questi .eq. 'EXI_PYRAM5') typma = 'PYRAM5'
        if (questi .eq. 'EXI_PYRAM13') typma = 'PYRAM13'
        if (questi .eq. 'EXI_POI1') typma = 'POI1'
        ASSERT(typma .ne. 'XXXX')
        call jenonu(jexnom('&CATA.TM.NOMTM', typma), typv)
        ASSERT(typv .gt. 0)
!
        repk = 'NON'
        nbma = dime(3)
        call jeveuo(nomob//'.TYPMAIL', 'L', vi=typmail)
        do k = 1, nbma
            if (typmail(k) .eq. typv) goto 51
        end do
        goto 52
51      continue
        repk = 'OUI'
52      continue
!
!
    else if (questi .eq. 'ONLY_SEG2') then
!     ----------------------------------------
        typma = 'SEG2'
        call jenonu(jexnom('&CATA.TM.NOMTM', typma), typv)
        ASSERT(typv .gt. 0)
!
        repk = 'OUI'
        nbma = dime(3)
        call jeveuo(nomob//'.TYPMAIL', 'L', vi=typmail)
        do k = 1, nbma
            if (typmail(k) .ne. typv) then
                repk = 'NON'
                exit
            end if
        end do
    else
        ASSERT(.false.)
    end if
!
!
!
    repkz = repk
    call jedema()
end subroutine
