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

subroutine rotchm(profno, vale, tetss, nbss, invsk, &
                  nbnot, nbcmp, iax)
!    P. RICHARD     DATE 10/02/92
!-----------------------------------------------------------------------
!  BUT: EFFECTUER LA ROTATION DE LA PARTIE DES CHAMNO.VALE CORRESPON-
!  -DANT A CHAQUE SOUS-STRUCTURE A PARTIR DU PROFNO GLOBAL, DU
!  MAILLAGE SQUELETTE GLOBAL, ET DU TABLEAU INV-SKELET ET DU VECTEUR
!  DES ANGLES DE ROTATION DE CHAQUE SOUS-STRUCTURE
    implicit none
!
!-----------------------------------------------------------------------
!
! PROFNO   /I/: NOM K19 DU NUME_EQUA GLOBAL
! VALE     /M/: VECTEUR CORRESPONDANT AU .VALE DU CHAMNO COURANT
! TESSS    /I/: VECTEUR DES ANGLE DE ROTATION DES SOUS-STRUCTURES
! NBSS     /I/: NOMBRE DE SOUS-STRUCTURES
! INVSK    /I/: TABLEAU INVERSE-SKELETTE
! NBNOT    /I/: NOMBRE DE NOEUDS GLOBAL
! NBCMP    /I/: NOMBRE DE COMPOSANTE MAX DE LA GRANDEUR SOUS-JACENTE
! IAX      /I/: NUMERO DE L'AXE DE ROTATION
!
!
!
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/intet0.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/nueq_chck.h"
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iax, ibid, icomp, inueq, j
    integer(kind=8) :: k, llnueq, llprno, ltidec, nbcmp, nbcmpm, nbec
    integer(kind=8) :: nbnot, nbss, numsec
    real(kind=8) :: tetac, tetcou
!-----------------------------------------------------------------------
    parameter(nbcmpm=10)
    character(len=6) :: pgc
    character(len=8) :: nomg
    character(len=19) :: profno
    character(len=24) :: prno, nueq
    integer(kind=8) :: invsk(nbnot, 2), ieq(nbcmpm)
    real(kind=8) :: vale(*), tetss(nbss), tet0(nbcmpm, nbcmpm), udep(nbcmpm)
!
!-----------------------------------------------------------------------
    data pgc/'ROTCHM'/
!-----------------------------------------------------------------------
!
!------------------------RECUPERATION DU PRNO DEEQ NUEQ-----------------
!
    call jemarq()
!
!-----RECUPERATION DU NOMBRE DU NOMBRE D'ENTIERS CODES ASSOCIE A DEPL_R
!
    nomg = 'DEPL_R'
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
    if (nbec .gt. 10) then
        call utmess('F', 'MODELISA_94')
    end if
!
    nueq = profno//'.NUEQ'
    prno = profno//'.PRNO'
    call nueq_chck(profno)
!
    call jenonu(jexnom(prno(1:19)//'.LILI', '&MAILLA'), ibid)
    call jeveuo(jexnum(prno, ibid), 'L', llprno)
    call jeveuo(nueq, 'L', llnueq)
!
!----------------------ALLOCATION VECTEUR DECODAGE----------------------
!
    call wkvect('&&'//pgc//'.DECODAGE', 'V V I', nbcmp, ltidec)
!
!---------------------------ROTATION------------------------------------
!
    tetcou = tetss(1)
    call intet0(tetcou, tet0, iax)
!
!
    do i = 1, nbnot
!
        numsec = invsk(i, 1)
        tetac = tetss(numsec)
        if (tetac .ne. tetcou) then
            tetcou = tetac
            call intet0(tetcou, tet0, iax)
        end if
!
        inueq = zi(llprno+(nbec+2)*(i-1))
        call isdeco(zi(llprno+(nbec+2)*(i-1)+2), zi(ltidec), nbcmp)
        icomp = 0
!
        do j = 1, nbcmpm
            if (zi(ltidec+j-1) .gt. 0) then
                icomp = icomp+1
                ieq(j) = zi(llnueq+inueq+icomp-2)
                udep(j) = vale(ieq(j))
            else
                ieq(j) = 0
                udep(j) = 0.d0
            end if
        end do
!
!
        do j = 1, nbcmpm
            if (ieq(j) .gt. 0) then
                vale(ieq(j)) = 0.d0
                do k = 1, nbcmpm
                    vale(ieq(j)) = vale(ieq(j))+tet0(j, k)*udep(k)
                end do
            end if
        end do
!
    end do
!
    call jedetr('&&'//pgc//'.DECODAGE')
!
    call jedema()
end subroutine
