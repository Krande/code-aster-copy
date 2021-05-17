! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine irmhpc(idfimd, nomamd, nomast, nbnoeu)
! person_in_charge: nicolas.sellenet at edf.fr
!-----------------------------------------------------------------------
!     ECRITURE DU MAILLAGE -  FORMAT MED - IMPRESSION POUR UN MAILLAGE PARALLELE
!        -  -     -                  -         --
!-----------------------------------------------------------------------
!     ENTREE:
!       IDFIMD  : IDENTIFIANT DU FICHIER MED
!       NOMAMD : NOM DU MAILLAGE MED
!       NOMAST : NOM UTILISATEUR DU MAILLAGE MED
!       NBNOEU : NOMBRE DE NOEUDS DU MAILLAGE
!-----------------------------------------------------------------------
!
implicit none
!
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/as_mmhgnw.h"
#include "asterfort/as_msdcrw.h"
#include "asterfort/as_msdjcr.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"
#include "jeveux.h"
#include "asterf_med.h"
!
! 0.1. ==> ARGUMENTS
!
    med_idt :: idfimd
    integer :: nbnoeu
    character(len=*) :: nomamd, nomast
!
! 0.3. ==> VARIABLES LOCALES
!
!
    integer :: codret, iret
    integer :: jnumno, nbjoin, i_join, nbnoj, jjoinr
    integer :: ifm, nivinf, domdis, rang, nbproc
    integer, pointer :: v_dojoin(:) => null()
    mpi_int :: mrank, msize
!
    character(len=4) :: chrang, chdomdis
    character(len=8) :: k8bid
    character(len=24) :: nonulg, domjoin, nojoin
    character(len=MED_NAME_SIZE) :: nomjoi
    character(len=MED_COMMENT_SIZE) :: descri
!
    aster_logical, pointer :: v_ldomj(:) => null()

!
!====
! 1. PREALABLES
!====
!
    call jemarq()
!
    call infniv(ifm, nivinf)
!
    call asmpi_info(rank = mrank, size = msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
! -- Impression numerotation globale des noeuds
!
    nonulg = nomast//'.NULOGL'
    call jeveuo(nonulg, 'L', jnumno)
!
    call as_mmhgnw(idfimd, nomamd, MED_NODE, MED_NONE, zi(jnumno), nbnoeu, codret)
    ASSERT(codret == 0)
!
! -- Impression des joints
!
    domjoin = nomast//'.DOMJOINTS'
    nbjoin = 0
    call jeexin(domjoin, iret)
    if(iret > 0) then
        call jeveuo(domjoin, 'L', vi=v_dojoin)
        call jelira(domjoin, 'LONMAX', nbjoin, k8bid)
        call wkvect("&&IRMHPC.DOMJOINTS", 'V V L', nbproc, vl=v_ldomj)
        v_ldomj(:) = ASTER_FALSE
        descri = "code_aster"
        call codent(rang, 'G', chrang)
        ASSERT(nbjoin <= 9999)
!
! --- Boucle sur les joints entre les sous-domaines
!
        do i_join = 1, nbjoin
            domdis = v_dojoin(i_join)
            call codent(domdis, 'G', chdomdis)
            nomjoi = " "
            nojoin = " "
            if( .not. v_ldomj(domdis+1) ) then
                if( rang < domdis) then
                    nomjoi(1:8) = chrang // chdomdis
                    nojoin = nomast//'.R'//chdomdis
                else
                    nomjoi(1:8) = chdomdis // chrang
                    nojoin = nomast//'.E'//chdomdis
                end if
            else
                if( rang < domdis) then
                    nomjoi(1:8) = chdomdis // chrang
                    nojoin = nomast//'.E'//chdomdis
                else
                    nomjoi(1:8) = chrang // chdomdis
                    nojoin = nomast//'.R'//chdomdis
                end if
            endif
            v_ldomj(domdis+1) = ASTER_TRUE
!
! --- Creation du joint
!
            call as_msdjcr(idfimd, nomamd, nomjoi, descri, domdis, nomamd, codret)
            ASSERT(codret == 0)
!
! --- Ecriture de la correspondance Noeud, Noeud
!
            call jelira(nojoin, 'LONMAX', nbnoj, k8bid)
            call jeveuo(nojoin, 'L', jjoinr)
!
            call as_msdcrw(idfimd, nomamd, nomjoi, MED_NO_DT, MED_NO_IT, MED_NODE, &
                            MED_NONE, MED_NODE, MED_NONE, nbnoj/2, zi(jjoinr), codret)
            ASSERT(codret == 0)
        end do
        call jedetr("&&IRMHPC.DOMJOINTS")
    endif
!
    call jedema()
!
end subroutine
