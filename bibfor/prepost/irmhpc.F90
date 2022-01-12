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
    character(len=6), parameter :: nompro = 'IRMHPC'
    integer :: codret, iret
    integer :: jnumno, nbjoin, i_join, nbnoj, jjoinr
    integer :: ifm, nivinf, domdis, rang, nbproc, i
    integer, pointer :: v_dojoin(:) => null()
    mpi_int :: mrank, msize
!
    character(len=8) :: chrang, chdomdis
    character(len=8) :: k8bid
    character(len=24) :: nonulg, domjoin, nojoin
    character(len=MED_NAME_SIZE) :: nomjoi
    character(len=MED_COMMENT_SIZE) :: descri

!
!====
! 1. PREALABLES
!====
!
    call jemarq()
!
    call infniv(ifm, nivinf)
    if (nivinf .gt. 1) then
        write (ifm,*) '<',nompro,'> DEBUT IMPRESSION DES JOINTS :'
    endif
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
        descri = "code_aster"
        call codent(rang, 'G', chrang)
        ASSERT(nbjoin <= 9999)
!
! --- Boucle sur les joints entre les sous-domaines
!
        do i_join = 1, nbjoin
            domdis = v_dojoin(i_join)
            call codent(domdis, 'G', chdomdis)
            do i = 1, 2
                nomjoi = " "
                nojoin = " "
                if( i == 1 ) then
                    if( rang < domdis) then
                        nomjoi = chrang // ' ' // chdomdis
                        nojoin = nomast//'.R'//chdomdis
                    else
                        nomjoi = chdomdis // ' ' // chrang
                        nojoin = nomast//'.E'//chdomdis
                    end if
                else
                    if( rang < domdis) then
                        nomjoi = chdomdis // ' ' // chrang
                        nojoin = nomast//'.E'//chdomdis
                    else
                        nomjoi = chrang // ' ' // chdomdis
                        nojoin = nomast//'.R'//chdomdis
                    end if
                endif
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
        end do
    endif
!
    if (nivinf .gt. 1) then
        write (ifm,*) '<',nompro,'> FIN IMPRESSION DES JOINTS.'
    endif
!
    call jedema()
!
end subroutine
