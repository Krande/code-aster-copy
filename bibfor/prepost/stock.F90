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

subroutine stock(resu, chs, nocham, ligrel, tychas, &
                 numord, iouf, numode, masgen, amrge, &
                 numeq)
    implicit none
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/assert.h"
#include "asterfort/cescel.h"
#include "asterfort/cnscno.h"
#include "asterfort/infniv.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsexpa.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: numord, numode
    real(kind=8) :: iouf, masgen, amrge
    character(len=*) :: resu, chs, nocham, ligrel, tychas
!
!     TRANSFERT DES RESULTATS DU CHAMP SIMPLE VERS LE CHAMP VRAI
!     PRESENT DANS LA STRUCTURE DE DONNEES RESULTATS
!
! IN  : CHS    : K19 : NOM DU CHAMP SIMPLE
! IN  : TYCHAS : K4  : TYPE DU CHAMP 'NOEU' 'ELNO' 'ELGA'
! IN  : NOCHAM : K16 : NOM DU CHAMP LU 'DEPL', 'SIEF_ELNO'
! IN  : NUMORD : I   : NUMERO D'ORDRE
! IN  : RESU   : K8  : NOM DE LA STRUCTURE DE DONNEES RESULTATS
! IN  : IOUF   : R   : IOUF DE L'INSTANT OU DE LA FREQUENCE
! IN  : NUMEQ : K19 : PROFIL DE STOCKAGE
!
!----------------------------------------------------------------------
!
!
!
    integer(kind=8) :: iret, jiouf, iad, nncp
    integer(kind=8) :: vali(2), ibid, nivinf, ifm
    real(kind=8) :: depi
    character(len=8) :: k8b, acce
    character(len=24) :: valk(2)
    character(len=16) :: param
    character(len=19) :: nomch, numeq
!
!- RECHERCHE DU NOM DU CHAMP RESULTAT
    depi = r8depi()
!
    call rsexch(' ', resu, nocham, numord, nomch, &
                iret)
!
    if (iret .eq. 100) then
    else if (iret .eq. 0) then
    else if (iret .eq. 110) then
        call rsagsd(resu, 0)
        call rsexch(' ', resu, nocham, numord, nomch, &
                    iret)
    else
        valk(1) = resu
        valk(2) = nomch
        vali(1) = numord
        vali(2) = iret
        call utmess('F', 'PREPOST5_73', nk=2, valk=valk, ni=2, &
                    vali=vali)
    end if
!
! - TRANSFERT DU CHAMP SIMPLE VERS LE CHAMP VRAI
!
!      CALL UTIMSD('MESSAGE',0,.TRUE.,.TRUE.,NUMEQ,1,' ')
!      CALL UTIMSD('MESSAGE',1,.TRUE.,.TRUE.,
!     &             NUMEQ//".PRNO",1,' ')
    if (tychas .eq. 'NOEU') then
        call cnscno(chs, numeq, 'NON', 'G', nomch, &
                    'F', ibid)
    else
        call cescel(chs, ligrel, ' ', ' ', 'OUI', &
                    nncp, 'G', nomch, 'F', ibid)
    end if
!
!
!-    ON NOTE LE CHAMP
!     ---------------------------------
    call infniv(ifm, nivinf)
    if (nivinf .ge. 1) then
        write (ifm, *) '<LRIDEA> LECTURE DU CHAMP  : ', nocham, numord
    end if
    call rsnoch(resu, nocham, numord)
!
!
!-    S: ON STOCKE LES PARAMETRES :
!     ---------------------------------
!
!     S.1 : INST OU FREQ :
!     --------------------
    acce = 'INST'
    call rsexpa(resu, 0, 'FREQ', iret)
    if (iret .gt. 0) acce = 'FREQ'
!
    call rsexpa(resu, 0, acce, iret)
    ASSERT(iret .gt. 0)
    call rsadpa(resu, 'E', 1, acce, numord, &
                0, sjv=jiouf, styp=k8b)
    zr(jiouf) = iouf
!
!     S.2 : NUME_MODE, MASS_GENE et AMOR_GENE :
!     ------------------------------
    param = 'NUME_MODE'
    call rsexpa(resu, 2, param, iret)
    if (iret .gt. 0) then
        call rsadpa(resu, 'E', 1, param, numord, &
                    0, sjv=iad, styp=k8b)
        zi(iad) = numode
    end if
!
    param = 'MASS_GENE'
    call rsexpa(resu, 2, param, iret)
    if (iret .gt. 0) then
        call rsadpa(resu, 'E', 1, param, numord, &
                    0, sjv=iad, styp=k8b)
        zr(iad) = masgen
    end if
!
    param = 'AMOR_REDUIT'
    call rsexpa(resu, 2, param, iret)
    if (iret .gt. 0) then
        call rsadpa(resu, 'E', 1, param, numord, &
                    0, sjv=iad, styp=k8b)
        zr(iad) = amrge
    end if
!
    param = 'AMOR_GENE'
    call rsexpa(resu, 2, param, iret)
    if (iret .gt. 0) then
        call rsadpa(resu, 'E', 1, param, numord, &
                    0, sjv=iad, styp=k8b)
        if (amrge .lt. 1.d92) then
            zr(iad) = 2*amrge*masgen*depi*iouf
        else
            zr(iad) = 0.d0
        end if
    end if
!
    param = 'RIGI_GENE'
    call rsexpa(resu, 2, param, iret)
    if (iret .gt. 0) then
        call rsadpa(resu, 'E', 1, param, numord, &
                    0, sjv=iad, styp=k8b)
        if (masgen .lt. 1.d92) then
            zr(iad) = masgen*(depi*iouf)**2
        else
            zr(iad) = 0.d0
        end if
    end if
!
    param = 'OMEGA2'
    call rsexpa(resu, 2, param, iret)
    if (iret .gt. 0) then
        call rsadpa(resu, 'E', 1, param, numord, &
                    0, sjv=iad, styp=k8b)
        zr(iad) = (depi*iouf)**2
    end if
!
!
!
end subroutine
