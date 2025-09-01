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

subroutine pjloin(nbnod, nbnodm, m2, geom2, nbmax, tino2m, tdmin2, lino_loin)
    implicit none

    integer(kind=8), intent(in) :: nbnod, nbnodm, nbmax, tino2m(nbmax), lino_loin(*)
    real(kind=8), intent(in) :: tdmin2(nbmax)
    real(kind=8), intent(in) :: geom2(*)
    character(len=8), intent(in) :: m2

#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterfort/codent.h"
#include "asterfort/crea_maillage.h"
#include "asterfort/detrsd.h"
#include "asterfort/gcncon.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/irmail.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/ulaffe.h"
#include "asterfort/ulnume.h"
#include "asterfort/ulopen.h"
#include "asterfort/utmess.h"
!
!     BUT :
!       Emettre (eventuellement) le message d'alarme de projection sur
!       des mailles lointaines.
!       En INFO=2, on crée le maillage de ces mailles.
!       Cette routine sert a mettre en commun du fortran pour plusieurs routines
! ----------------------------------------------------------------------
!
    character(len=8) ::  madebug, k8bid, kico
    character(len=16) :: alarm, k16bid, nomcmd, formar
    character(len=19) :: k19bid
!
    integer(kind=8) :: vali(2)
    integer(kind=8) :: ii, ino2m, unite, ico, iret, ifm, info, ivers
    real(kind=8) :: valr(4)
    character(len=80) :: valk(2)
    character(len=8) :: nono2
    character(len=80) :: fichier
    save ico
! --- DEB --------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, info)
!
    ico = 0
    if (nbnodm .ne. 0) then
        alarm = 'OUI'
        call getres(k16bid, k16bid, nomcmd)
        if (nomcmd .eq. 'PROJ_CHAMP') then
            call getvtx(' ', 'ALARME', scal=alarm, nbret=iret)
        end if
        if (alarm .eq. 'OUI') then
            ico = ico+1
            call codent(ico, 'D0', kico)
            do ii = 1, nbnod
                ino2m = tino2m(ii)
                nono2 = int_to_char8(ino2m)
                valr(1) = geom2(3*(ino2m-1)+1)
                valr(2) = geom2(3*(ino2m-1)+2)
                valr(3) = geom2(3*(ino2m-1)+3)
                valr(4) = tdmin2(ii)
                call utmess('I', 'CALCULEL5_43', sk=nono2, nr=4, valr=valr)
            end do
            vali(1) = nbnodm
            vali(2) = nbnod
            fichier = 'REPE_OUT/maillage_proj_loin_'//kico//'.med'
            valk(1) = fichier
            call utmess('A', 'CALCULEL5_48', ni=2, vali=vali, nk=1, valk=valk)

!       -- Creation et impression d'un "petit" maillage contenant juste les noeuds
!          lointains. Cela peut aider l'utilisateur a les visualiser.
            if (info .eq. 2) then
                call gcncon('_', madebug)
                call crea_maillage(m2, madebug, 'V', nbno=nbnodm, lino=lino_loin)

                unite = ulnume()
                if (unite .le. 0) call utmess('F', 'UTILITAI5_10')
                call ulaffe(unite, fichier, ' ', 'N', 'O')
                formar = ' '
                k8bid = ' '
                k19bid = ' '
                call irmail('MED', unite, ivers, madebug, ASTER_FALSE, k19bid, 1, formar)
                call ulopen(-unite, k8bid, k8bid, k8bid, k8bid)
                call detrsd('MAILLAGE', madebug)
            end if
        end if
    end if

    call jedema()
end subroutine
