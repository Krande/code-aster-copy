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

subroutine crvrc1()
    implicit none
!
!     COMMANDE:  CREA_RESU
!     CREE UNE STRUCTURE DE DONNEE DE TYPE "EVOL_THER"  CONTENANT
!     LA TEMPERATURE SUR LES COUCHES DES COQUES MULTICOUCHE A PARTIR
!     D'UN CHAMP DE FONCTIONS (INST,EPAIS)
!
!
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/cesvar.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    integer :: nbfac, n1, nbinst, kinst, jinst
    integer :: nbin, iret, iexi
    real(kind=8) :: vinst
    character(len=8) :: kbid, resu, carele, paout, lpain(10), tempef
    character(len=8) :: model2, modele
    character(len=16) :: type, oper
    character(len=19) :: ligrmo, chout, chinst
    character(len=24) :: chcara(18), lchin(10)
    real(kind=8), pointer :: linst(:) => null()
    character(len=24), pointer :: celk(:) => null()
!
!----------------------------------------------------------------------
    call jemarq()
!
    call getfac('PREP_VRC1', nbfac)
    if (nbfac .eq. 0) goto 20
    ASSERT(nbfac.eq.1)
!
!
    call getres(resu, type, oper)
    call getvid('PREP_VRC1', 'MODELE', iocc=1, scal=modele, nbret=n1)
    call getvid('PREP_VRC1', 'CARA_ELEM', iocc=1, scal=carele, nbret=n1)
    call getvid('PREP_VRC1', 'CHAM_GD', iocc=1, scal=tempef, nbret=n1)
!
!     -- ON VERIFIE QUE LE CARA_ELEM S'APPUIE BIEN SUR LE MODELE
    call jeexin(carele//'.CANBSP    .CELK', iexi)
    if (iexi .eq. 0) then
        call utmess('F', 'CALCULEL4_14', sk=carele)
    endif
    call jeveuo(carele//'.CANBSP    .CELK', 'L', vk24=celk)
    model2=celk(1)(1:8)
    if (model2 .ne. modele) then
        call utmess('F', 'CALCULEL4_15', sk=carele)
    endif
!
!
!
!     -- INSTANTS DE L'EVOL_THER :
    call getvr8('PREP_VRC1', 'INST', iocc=1, nbval=0, nbret=n1)
    ASSERT(n1.lt.0)
    nbinst = -n1
    AS_ALLOCATE(vr=linst, size=nbinst)
    call getvr8('PREP_VRC1', 'INST', iocc=1, nbval=nbinst, vect=linst,&
                nbret=n1)
!
    call jeexin(resu//'           .DESC', iret)
    if (iret .ne. 0) then
        call utmess('F', 'CALCULEL7_6', sk=resu)
    else
        call rscrsd('G', resu, 'EVOL_THER', nbinst)
    endif
!
    ligrmo = modele//'.MODELE'
    paout = 'PTEMPCR'
    chinst = '&&CRVRC1.CHINST'
    call mecara(carele, chcara)
!
    lpain(1) = 'PNBSP_I'
    lchin(1) = chcara(1) (1:8)//'.CANBSP'
    lpain(2) = 'PTEMPEF'
    lchin(2) = tempef
    lpain(3) = 'PINST_R'
    lchin(3) = chinst
    lpain(4) = 'PCACOQU'
    lchin(4) = chcara(7)
    nbin = 4
!
!     -- BOUCLE SUR LES INSTANTS :
!     --------------------------------
    do kinst = 1,nbinst
        vinst = linst(kinst)
        call mecact('V', chinst, 'MODELE', ligrmo, 'INST_R',&
                    ncmp=1, nomcmp='INST', sr=vinst)
        call rsexch(' ', resu, 'TEMP', kinst, chout,&
                    iret)
!
        call cesvar(carele, ' ', ligrmo, chout)
        call calcul('S', 'PREP_VRC', ligrmo, nbin, lchin,&
                    lpain, 1, chout, paout, 'G',&
                    'OUI')
        call detrsd('CHAM_ELEM_S', chout)
        call rsnoch(resu, 'TEMP', kinst)
        call rsadpa(resu, 'E', 1, 'INST', kinst,&
                    0, sjv=jinst, styp=kbid)
        zr(jinst) = vinst
        call detrsd('CHAMP', chinst)
    end do
!
!
    AS_DEALLOCATE(vr=linst)
!
!
20  continue
    call jedema()
end subroutine
