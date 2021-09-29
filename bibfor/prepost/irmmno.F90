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

subroutine irmmno(idfimd, nomamd, ndim, nbnoeu, coordo,&
                  nomnoe, nosdfu)
! person_in_charge: nicolas.sellenet at edf.fr
!-----------------------------------------------------------------------
!     ECRITURE DU MAILLAGE -  FORMAT MED - LES NOEUDS
!        -  -     -                  -         --
!-----------------------------------------------------------------------
!     ENTREE:
!       IDFIMD  : IDENTIFIANT DU FICHIER MED
!       NOMAMD : NOM DU MAILLAGE MED
!       NDIM   : DIMENSION DU PROBLEME (2  OU 3)
!       NBNOEU : NOMBRE DE NOEUDS DU MAILLAGE
!       COORDO : VECTEUR DES COORDONNEES DES NOEUDS
!       NOMNOE : VECTEUR NOMS DES NOEUDS
!-----------------------------------------------------------------------
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "jeveux.h"
#include "asterfort/as_mfrall.h"
#include "asterfort/as_mfrblc.h"
#include "asterfort/as_mfrdea.h"
#include "asterfort/as_mmhcow.h"
#include "asterfort/as_mmhcaw.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    med_idt :: idfimd
    integer :: ndim, nbnoeu
!
    real(kind=8) :: coordo(*)
!
    character(len=*) :: nomnoe(*)
    character(len=*) :: nomamd
    character(len=8) :: nosdfu
!
! 0.2. ==> COMMUNS
!
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter ( nompro = 'IRMMNO' )
!
    integer :: edfuin
    parameter (edfuin=0)
!
    integer :: codret
    integer :: iaux
    integer :: jcoord
    integer :: ifm, nivinf, rang, nbproc, start, filter(1), jaux
    integer :: jnbno, nbnot, jno, nbnol, cmpt
    mpi_int :: mrank, msize
    aster_logical :: lfu
!
    character(len=8) :: saux08
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
    lfu = .false._1
    if(nosdfu.ne.' ') then
        lfu = .true._1
        call jeveuo(nosdfu//'.NBNO', 'L', jnbno)
        call jeveuo(nosdfu//'.NOEU', 'L', jno)
        start = zi(jnbno)
        nbnol = zi(jnbno+1)
        nbnot = zi(jnbno+2)
        call as_mfrall(1, filter, codret)
!
        call as_mfrblc(idfimd, nbnot, 1, ndim, 0,&
                    edfuin, 2, "", start, nbnol,&
                    1, nbnol, 0, filter(1), codret)
        if (codret .ne. 0) then
            saux08='mfrblc'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        endif
    else
        jnbno = 0
        jno = 0
        nbnol = nbnoeu
    endif
!
!====
! 2. ECRITURE DES COORDONNEES DES NOEUDS
!    LA DIMENSION DU PROBLEME PHYSIQUE EST VARIABLE (1,2,3), MAIS
!    ASTER STOCKE TOUJOURS 3 COORDONNEES PAR NOEUDS.
!====
!
!    LE TABLEAU COORDO EST UTILISE AINSI : COORDO(NDIM,NBNOEU)
!    EN FORTRAN, CELA CORRESPOND AU STOCKAGE MEMOIRE SUIVANT :
!    COORDO(1,1), COORDO(2,1), COORDO(3,1), COORDO(1,2), COORDO(2,2),
!    COORDO(3,2), COORDO(1,3), ... , COORDO(1,NBNOEU), COORDO(2,NBNOEU),
!    COORDO(3,NBNOEU)
!    C'EST CE QUE MED APPELLE LE MODE ENTRELACE
!
    call wkvect('&&'//nompro//'.COORDO', 'V V R', nbnol*ndim, jcoord)
!
    if(lfu) then
        cmpt = 0
        do iaux = 0, nbnoeu-1
            if(zi(jno+iaux).gt.0) then
                do jaux = 0, ndim-1
                    zr(jcoord+ndim*cmpt+jaux) = coordo(3*iaux+jaux+1)
                enddo
                cmpt = cmpt + 1
            endif
        enddo
        call as_mmhcaw(idfimd, nomamd, filter(1), zr(jcoord), codret)
    else
        do iaux = 0, nbnoeu-1
            do jaux = 0, ndim-1
                zr(jcoord+ndim*iaux+jaux) = coordo(3*iaux+jaux+1)
            enddo
        enddo
        call as_mmhcow(idfimd, nomamd, zr(jcoord), edfuin, nbnoeu,&
                       codret)
    endif
!
    if (codret .ne. 0) then
        saux08='mmhcow'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    endif
!
    call jedetr('&&'//nompro//'.COORDO')
!
!
!====
! 3. LA FIN
!====
!
    if(lfu) then
        call as_mfrdea(1, filter, codret)
        if (codret .ne. 0) then
            saux08='mfrdea'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        endif
    endif
!
    call jedetr('&&'//nompro//'NOMNOE')
!
    call jedema()
!
end subroutine
