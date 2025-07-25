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
subroutine srlima(mo, mail2d, mail3d, mailto, nbma2d)
    implicit none
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/oriem1.h"
#include "asterfort/reliem.h"
#include "asterfort/utflmd.h"
#include "asterfort/utmamo.h"
#include "asterfort/utmasu.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbma2d
    character(len=8) :: mo
    character(len=24) :: mail2d, mail3d, mailto
!
!  But: construire 3 listes de mailles a partir des donnees
!       utilisateur :
!       liste des mailles 2d
!       liste des mailles 3d sous-jacentes
!       liste de l'ensemble des mailles 2d + mailles 3d sous-jacentes
!
!  in  mo     : nom du modele
!  in/jxout  mail2d : nom objet jeveux contenant la liste des mailles 2d
!  in/jxout  mail3d : nom objet jeveux contenant la liste des mailles 2d
!  in/jxout  mailto : nom objet jeveux contenant la liste des mailles 2d+3d
!  out nbma2d : nombre de mailles 2d trouvees == nb mailles 3d
!
!  Remarque : si une maille 2d n'a pas de maille sous-jacente, le numero
!             de la maille 3d est 0.
!
! ----------------------------------------------------------------------
!
!
    integer(kind=8) :: n1, n2, n3
    integer(kind=8) :: ima, iret, jcesd, jcesl, iad1, numa2d, numa3d
    integer(kind=8) :: nbma, nbmamo, jlima, nbmat, jmato, ndim, ndim2
    character(len=8) :: ma, limocl(3), tymocl(3)
    character(len=24) :: mesmai, limamo
    character(len=19) :: ces, cel
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: listCell2d(:) => null(), listCell3d(:) => null()
!
    data limocl/'TOUT', 'MAILLE', 'GROUP_MA'/
    data tymocl/'TOUT', 'MAILLE', 'GROUP_MA'/
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
!   -- on recupere les mailles de peau :
!   ----------------------------------
    call dismoi('NOM_MAILLA', mo, 'MODELE', repk=ma)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbmat)
    call dismoi('DIM_GEOM', ma, 'MAILLAGE', repi=ndim)
    mesmai = '&&SRLIMA.MAILLU'
!
    call getvtx(' ', 'TOUT', nbret=n1)
    call getvtx(' ', 'MAILLE', nbret=n2)
    call getvtx(' ', 'GROUP_MA', nbret=n3)
!
!   -- on scrute les mot-clé :
!   ------------------------
    if (n1+n2+n3 .ne. 0) then
        call reliem(mo, ma, 'NU_MAILLE', ' ', 0, &
                    3, limocl, tymocl, mesmai, nbma)
    else
!   -- mais s'il n'y en a aucun on fait comme si TOUT='OUI' :
!   -------------------------------------------------------
        call utmamo(mo, nbma, mesmai)
    end if
!
!   -- on ne garde que les mailles surfaciques :
!   ---------------------------------------------
    ndim2 = ndim-1
    call utflmd(ma, mesmai, nbma, ndim2, ' ', &
                nbma2d, mail2d)
    if (nbma2d .gt. 0) then
        call jeveuo(mail2d, 'L', vi=listCell2d)
    else
        call utmess('F', 'CALCULEL5_54')
    end if
!
!   -- on recherche les mailles 3d qui bordent les mailles de peau :
!   -- il faut se limiter aux mailles du modele qui savent calculer SIGM_ELNO :
!   ----------------------------------------------------------------------------
    cel = '&&SRLIMA.CEL_ELNO'
    ces = '&&SRLIMA.CES_ELNO'
    call alchml(mo//'.MODELE', 'SIGM_ELNO', 'PSIEFNOR', 'V', cel, &
                iret, ' ')
    ASSERT(iret .eq. 0)
    call celces(cel, 'V', ces)
    call detrsd('CHAMP', cel)
    call jeveuo(ces//'.CESD', 'L', jcesd)
    call jeveuo(ces//'.CESL', 'L', jcesl)
    nbma = zi(jcesd-1+1)
    limamo = '&&SRLIMA.LIMAIL'
    call wkvect(limamo, 'V V I', nbma, jlima)
    nbmamo = 0
    do ima = 1, nbma
        call cesexi('C', jcesd, jcesl, ima, 1, &
                    1, 1, iad1)
        if (iad1 .gt. 0) then
            nbmamo = nbmamo+1
            zi(jlima-1+nbmamo) = ima
        end if
    end do
    call detrsd('CHAM_ELEM_S', ces)
!
!
!   -- on recherche les mailles 3d associees aux mailles de peau :
!   ---------------------------------------------------------------
    call jeveuo(ma//'.COORDO    .VALE', 'L', vr=vale)
    if (ndim .eq. 3) then
        call utmasu(ma, '3D', nbma2d, listCell2d, mail3d, &
                    vale, nbmamo, zi(jlima), .true._1)
    end if
    if (ndim .eq. 2) then
        call utmasu(ma, '2D', nbma2d, listCell2d, mail3d, &
                    vale, nbmamo, zi(jlima), .true._1)
    end if
    call jeveuo(mail3d, 'L', vi=listCell3d)
!
    call wkvect(mailto, 'V V I', nbma2d*2, jmato)
!
    do ima = 1, nbma2d
        numa2d = listCell2d(ima)
        numa3d = listCell3d(ima)
        if (numa3d .eq. 0) goto 1
!
!       -- si la maille 3d n'est pas du cote "-", on la met a zero (issue22570) :
        if (ndim .eq. 3) then
            call oriem1(ma, '3D', numa2d, numa3d)
        end if
!
        zi(jmato-1+ima) = numa2d
        zi(jmato-1+nbma2d+ima) = numa3d
        listCell3d(ima) = numa3d
1       continue
    end do
!
    call jedetr(mesmai)
    call jedetr(limamo)
!
    call jedema()
!
end subroutine
