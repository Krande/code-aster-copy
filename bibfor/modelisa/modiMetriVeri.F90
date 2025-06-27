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
subroutine modiMetriVeri(noma, ioc, modmai, nutyptu)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
    character(len=8), intent(in) :: noma
    integer(kind=8), intent(in) :: ioc, nutyptu(3)
    character(len=24), intent(in) :: modmai
!     AFFE_CARA_ELEM/POUTRE
!     VERIFICATION QUE MODI_METRIQUE='OUI' N'EST AFFECTE QUE
!     SUR DES TUYAUX
! ----------------------------------------------------------------------
    integer(kind=8) :: ibid, nbma, jma, i, ima, ixma, jdme, itypma, ityp
    character(len=8) :: mmt, typmcl(2)
    character(len=16) :: motfac, motcls(2)
    character(len=24) :: mesmai
    aster_logical :: maok
!     ------------------------------------------------------------------
    call jemarq()
!
    motfac = 'POUTRE'
    mesmai = '&&ACEMMT.MES_MAILLES'
    motcls(1) = 'GROUP_MA'
    motcls(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
    mmt = ' '
!
    call getvtx(motfac, 'MODI_METRIQUE', iocc=ioc, scal=mmt, nbret=ibid)
    if (mmt .eq. 'OUI') then
        call reliem(' ', noma, 'NU_MAILLE', motfac, ioc, &
                    2, motcls, typmcl, mesmai, nbma)
        if (nbma .ne. 0) then
            call jeexin(modmai, ixma)
            if (ixma .ne. 0) then
                call jeveuo(modmai, 'L', jdme)
                call jeveuo(mesmai, 'L', jma)
                do i = 1, nbma
                    ima = zi(jma+i-1)
                    itypma = zi(jdme-1+ima)
                    maok = ASTER_FALSE
                    do ityp = 1, 3
                        if (itypma .eq. nutyptu(ityp)) then
                            maok = ASTER_TRUE
                            exit
                        end if
                    end do
                    if (.not. maok) then
                        call utmess('F', 'MODELISA10_1')
                    end if
                end do
            end if
!
            call jedetr(mesmai)
        end if
    end if
!

    call jedema()
end subroutine
