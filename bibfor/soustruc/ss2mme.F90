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
subroutine ss2mme(modelz, vesstrz, base)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=*), intent(in) :: modelz, vesstrz
    character(len=1), intent(in) :: base
!
! --------------------------------------------------------------------------------------------------
!
! PREPARER LE VECT_ELEM DANS LE CAS DE SOUS-STRUCTURES
!
! --------------------------------------------------------------------------------------------------
!
! IN  NOMO   : NOM DU MODELE
! IN  BASE   : BASE DE CREATION DU VECT_ELEM
! I/O VESSTR : NOM DU VECT_ELEM
!                OUT : VESSTR EST (EVENTUELLEMENT) ENRICHI DE
!                L'OBJET .RELC
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywfact = 'SOUS_STRUC'
    character(len=8) :: mesh, model, k8bid, nosma, nomcas, nomacr
    integer(kind=8) :: nbssa, nbsma, n1, n2, nboc
    integer(kind=8) ::   ialsch, imas
    integer(kind=8) :: ier0, ioc, i, iret
    character(len=19) :: vesstr
    character(len=16) :: valk(2)
    character(len=8), pointer :: lmai(:) => null()
    character(len=8), pointer :: vnomacr(:) => null()
    integer(kind=8), pointer :: sssa(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    vesstr = vesstrz
    model = modelz
    call getfac(keywfact, nboc)
    if (nboc .eq. 0) goto 999

! - Parameters from mesh
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi=nbssa)
    call dismoi('NB_SM_MAILLA', model, 'MODELE', repi=nbsma)
!
    if (nbssa .eq. 0) then
        call utmess('F', 'SOUSTRUC_24')
    end if
!
    call jeveuo(model//'.MODELE    .SSSA', 'L', vi=sssa)
    call jeveuo(mesh//'.NOMACR', 'L', vk8=vnomacr)
!
    call jecrec(vesstr(1:19)//'.RELC', base//' V I', 'NO', 'CONTIG', 'CONSTANT', nboc)
    call jeecra(vesstr(1:19)//'.RELC', 'LONMAX', nbsma)
!
    AS_ALLOCATE(vk8=lmai, size=nbsma)
!
! --- BOUCLE SUR LES CAS_DE_CHARGE
!
    ier0 = 0
    do ioc = 1, nboc
!
        call getvtx(keywfact, 'CAS_CHARGE', iocc=ioc, scal=nomcas, nbret=n1)
        call jecroc(jexnom(vesstr(1:19)//'.RELC', nomcas))
        call jeveuo(jexnom(vesstr(1:19)//'.RELC', nomcas), 'E', ialsch)
!
!       -- CAS : TOUT: 'OUI'
!
        call getvtx(keywfact, 'TOUT', iocc=ioc, scal=k8bid, nbret=n1)
        if (n1 .eq. 1) then
            do i = 1, nbsma
                if (sssa(i) .eq. 1) zi(ialsch-1+i) = 1
            end do
            goto 5
        end if
!
!       -- CAS : MAILLE: L_MAIL
!
        call getvtx(keywfact, 'SUPER_MAILLE', iocc=ioc, nbval=0, nbret=n2)
        if (-n2 .gt. nbsma) then
            call utmess('F', 'SOUSTRUC_25')
        else
            call getvtx(keywfact, 'SUPER_MAILLE', iocc=ioc, nbval=nbsma, vect=lmai, &
                        nbret=n2)
        end if
        do i = 1, n2
            nosma = lmai(i)
            call jenonu(jexnom(mesh//'.SUPMAIL', nosma), imas)
            if (imas .eq. 0) then
                valk(1) = nosma
                valk(2) = mesh
                call utmess('F', 'SOUSTRUC_26', nk=2, valk=valk)
            else
                zi(ialsch-1+imas) = 1
            end if
        end do
!
!       -- ON VERIFIE QUE LES VECTEURS ELEMENTAIRES SONT CALCULES
!
5       continue
        do i = 1, nbsma
            if (zi(ialsch-1+i) .ne. 0) then
                call jenuno(jexnum(mesh//'.SUPMAIL', i), nosma)
                if (sssa(i) .ne. 1) then
                    call utmess('F', 'SOUSTRUC_27', sk=nosma)
                end if
!
                nomacr = vnomacr(i)
                call jeexin(jexnom(nomacr//'.LICA', nomcas), iret)
                if (iret .eq. 0) then
                    ier0 = 1
                    valk(1) = nosma
                    valk(2) = nomcas
                    call utmess('E', 'SOUSTRUC_28', nk=2, valk=valk)
                end if
            end if
        end do
    end do
!
    if (ier0 .eq. 1) then
        call utmess('F', 'SOUSTRUC_29')
    end if
!
    AS_DEALLOCATE(vk8=lmai)
!
999 continue
!
    call jedema()
end subroutine
