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

subroutine ssvalm(statut, option, mo, ma, isma, &
                  jresl, nbvel)
!
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/matrot.h"
#include "asterfort/ssrone.h"
#include "asterfort/ssvaro.h"
#include "asterfort/wkvect.h"
    character(len=8) :: mo, ma
    character(len=*) :: option, statut
    integer(kind=8) :: isma, jresl
! ----------------------------------------------------------------------
!     BUT:
!
!         DONNER L'ADRESSE DE L'OBJET JEVEUX CONTENANT LA MATRICE
!         ELEMENTAIRE (CONDENSEE) CORRESPONDANT A L'OPTION : OPTION
!         (CET OBJET N'EST PAS FORCEMENT LA MATRICE CONDENSEE : XP_EE
!          CAR IL PEUT Y AVOIR ROTATION/SYMETRIE DE LA SOUS-STRUCTURE)
!
!     IN:    STATUT : 'DEBUT','FIN', OU ' '(COURANT)
!     ---
!
!            LES STATUTS : 'DEBUT' ET 'FIN' SERVENT A PREPARER LA BOUCLE
!                       SUR LES MATRICES ELEMENTAIRES. (OBLIGATOIRES !)
!            EXEMPLE:
!
!              DO 1  ISMA= 1, NB_MAILLES
!
!           1  CONTINUE
!
!
!
!            OPTION(K16):  'RIGI_MECA', 'MASS_MECA', 'AMOR_MECA',
!            MO(K8) : NOM DU MODELE
!            MA(K8) : NOM DU MAILLAGE
!            ISMA   : NUMERO DE LA (SUPER)MAILLE
!
!     OUT:   JRESL : ADRESSE DE LA MATRICE CONDENSEE.
!     ----   NBVEL   : NOMBRE DE VALEURS DE CETTE MATRICE.
!
!
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=8) :: rota
    character(len=16) :: optio2
    character(len=8) :: nomacr
    real(kind=8) :: lambda(6, 6), angl(3), pgl(3, 3)
    character(len=24) :: nomob
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadesm, ianmcr, iavmat
    integer(kind=8) :: iret, j, jsma, nbsma, nbssa, nbvel, nddle
    integer(kind=8) :: ndim, nmxval
    integer(kind=8), pointer :: sssa(:) => null()
    real(kind=8), pointer :: para_r(:) => null()
!-----------------------------------------------------------------------
    optio2 = option
!
!     -- SI APPEL INITIAL : ON ALLOUE UN OBJET SUFFISANT :
!     ----------------------------------------------------
    if (statut(1:5) .eq. 'DEBUT') then
        call dismoi('NB_SM_MAILLA', mo, 'MODELE', repi=nbsma)
        call dismoi('NB_SS_ACTI', mo, 'MODELE', repi=nbssa)
        if (nbssa .gt. 0) then
            call jeveuo(mo//'.MODELE    .SSSA', 'L', vi=sssa)
            call jeveuo(ma//'.NOMACR', 'L', ianmcr)
            nmxval = 0
            do jsma = 1, nbsma
                if (sssa(jsma) .eq. 1) then
                    nomacr = zk8(ianmcr-1+jsma)
                    call jeveuo(nomacr//'.DESM', 'L', iadesm)
                    nddle = zi(iadesm-1+4)
!             --LA DIMENSION DES MATRICES CONDENSEES EST DONNEE PAR
!               LA RIGIDITE:
                    ndim = nddle*(nddle+1)/2
                    nmxval = max(nmxval, ndim)
                end if
            end do
            if (nmxval .gt. 0) then
                call wkvect('&&SSVALM.VALEURS', 'V V R', nmxval, jresl)
            end if
        end if
    end if
!
!
!     -- SI APPEL FINAL : ON DETRUIT LES OBJETS DE TRAVAIL :
!     ------------------------------------------------------
    if (statut(1:3) .eq. 'FIN') then
        call jedetr('&&SSVALM.VALEURS')
        call jeexin('&&SSVARO.IINO', iret)
        if (iret .gt. 0) call jedetr('&&SSVARO.IINO')
    end if
!
!
!     -- SI APPEL COURANT : ON RECOPIE OU ON TOURNE :
!     -----------------------------------------------
    if (statut(1:1) .eq. ' ') then
        call jeveuo(ma//'.NOMACR', 'L', ianmcr)
        nomacr = zk8(ianmcr-1+isma)
        call jeveuo(nomacr//'.DESM', 'L', iadesm)
        nddle = zi(iadesm-1+4)
!
        if (optio2(1:4) .eq. 'RIGI' .or. optio2 == "MECA_DDLM_R") then
!          NOMOB=NOMACR//'.KP_EE'
            nomob = nomacr//'.MAEL_RAID_VALE'
        else if (optio2(1:4) .eq. 'MASS') then
!          NOMOB=NOMACR//'.MP_EE'
            nomob = nomacr//'.MAEL_MASS_VALE'
        else if (optio2(1:4) .eq. 'AMOR') then
!          NOMOB=NOMACR//'.AP_EE'
            nomob = nomacr//'.MAEL_AMOR_VALE'
        else
            ASSERT(.false.)
        end if
!
!       call jeveuo(nomob, 'L', iavmat)
        call jeveuo(jexnum(nomob, 1), 'L', iavmat)
        nbvel = nddle*(nddle+1)/2
!
!
!       -- RECOPIE (OU ROTATION):
!       -------------------------
        call jeveuo('&&SSVALM.VALEURS', 'E', jresl)
        call ssrone(ma, isma, rota)
!
        if (rota(1:3) .eq. 'NON') then
!         RECOPIE:
            do i = 1, nbvel
                zr(jresl-1+i) = zr(iavmat-1+i)
            end do
        else if (rota(1:3) .eq. 'OUI') then
!         ROTATION:
            call jeveuo(ma//'.PARA_R', 'L', vr=para_r)
            angl(1) = para_r(14*(isma-1)+4)
            angl(2) = para_r(14*(isma-1)+5)
            angl(3) = para_r(14*(isma-1)+6)
            call matrot(angl, pgl)
            do i = 1, 3
                do j = 1, 3
                    lambda(i, j) = pgl(i, j)
                    lambda(i, j+3) = 0.d0
                    lambda(i+3, j) = 0.d0
                    lambda(i+3, j+3) = pgl(i, j)
                end do
            end do
            call ssvaro(lambda, 'LG', .true._1, 'EXTE', nomacr, &
                        iavmat, jresl)
        else
            ASSERT(.false.)
        end if
!
        call jelibe(nomob)
    end if
!
!
!
!
end subroutine
