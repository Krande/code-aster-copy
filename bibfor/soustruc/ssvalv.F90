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

subroutine ssvalv(statut, nomcas, mo, ma, isma, &
                  idresl, long, instap)
!
! INSPI  SSVALM
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/matrot.h"
#include "asterfort/ssrone.h"
#include "asterfort/ssvaro.h"
#include "asterfort/ssvau1.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: mo, ma
    character(len=*) :: statut
    character(len=8) :: nomcas
    integer(kind=8) :: isma, idresl
    real(kind=8), optional :: instap
! ----------------------------------------------------------------------
!     BUT:
!
!         DONNER L'ADRESSE DE L'OBJET JEVEUX CONTENANT LE VECTEUR
!         ELEMENTAIRE (CONDENSE) CORRESPONDANT A NOMCHAR
!         (CET OBJET N'EST PAS FORCEMENT LE VECTEUR CONDENSE : XP_EE
!          CAR IL PEUT Y AVOIR ROTATION/SYMETRIE DE LA SOUS-STRUCTURE)
!
!     IN:    STATUT : 'DEBUT','FIN', OU ' '(COURANT)
!     ---
!
!            LES STATUTS : 'DEBUT' ET 'FIN' SERVENT A PREPARER LA BOUCLE
!                       SUR LES VECTEURS ELEMENTAIRES. (OBLIGATOIRES !)
!            EXEMPLE:
!
!              DO 1  ISMA= 1, NB_MAILLES
!
!           1  CONTINUE
!
!
!            NOMCAS    :   NOM DU CAS DE CHARGE
!            MO(K8) : NOM DU MODELE
!            MA(K8) : NOM DU MAILLAGE
!            ISMA   : NUMERO DE LA (SUPER)MAILLE
!
!     OUT:   IDRESL : ADRESSE DU VECTEUR CONDENSE.
!     ----   LONG   : NOMBRE DE VALEURS DE CE VECTEUR.
!
!
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=8) :: nomacr, rota
    real(kind=8) :: lambda(6, 6), angl(3), pgl(3, 3)
!
!
!
!     1- SI APPEL INITIAL : ON ALLOUE UN OBJET SUFFISANT :
!     ----------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadesm, ialica, ialich, ianmcr
    integer(kind=8) :: idres2, iret, j, jsma, long, nbsma
    integer(kind=8) :: nbssa, nddle, nddli, nddlt, nmxval
    integer(kind=8) :: nbinst, jresu, jinst
    character(len=19) :: resuge
    real(kind=8) :: inst1, inst2, coecor
    real(kind=8), pointer :: para_r(:) => null()
    integer(kind=8), pointer :: sssa(:) => null()
!-----------------------------------------------------------------------
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
                    nddli = zi(iadesm-1+5)
                    nddlt = nddli+nddle
                    nmxval = max(nmxval, nddlt)
                end if
            end do
            if (nmxval .gt. 0) then
                call wkvect('&&SSVALV.VALEURS', 'V V R', nmxval, idresl)
!           --          '&&SSVALV.VALTEMP' EST UN VECTEUR DE TRAVAIL :
                call wkvect('&&SSVALV.VALTEMP', 'V V R', nmxval, idres2)
            end if
        end if
    end if
!
!
!     2- SI APPEL FINAL : ON DETRUIT L OBJET DE TRAVAIL :
!     ---------------------------------------------------
    if (statut(1:3) .eq. 'FIN') then
        call jedetr('&&SSVALV.VALEURS')
        call jedetr('&&SSVALV.VALTEMP')
        call jeexin('&&SSVARO.IINO', iret)
        if (iret .gt. 0) call jedetr('&&SSVARO.IINO')
    end if
!
!
!     3- SI APPEL COURANT :
!     ---------------------
    if (statut(1:1) .eq. ' ') then
        call jeveuo(ma//'.NOMACR', 'L', ianmcr)
        nomacr = zk8(ianmcr-1+isma)
        call jeveuo(nomacr//'.DESM', 'L', iadesm)
        nddle = zi(iadesm-1+4)
        nddli = zi(iadesm-1+5)
        nddlt = nddli+nddle
        long = nddle
        call jeveuo(jexnom(nomacr//'.LICH', nomcas), 'L', ialich)
!
!
!       3.1- ON DETERMINE SI ON DOIT FAIRE LA ROTATION:
!       -----------------------------------------------
        call ssrone(ma, isma, rota)
!
!
!       3.2- RECOPIE (OU ROTATION) DE .LICA DANS .VALEURS :
!       ---------------------------------------------------
        call jeveuo('&&SSVALV.VALEURS', 'E', idresl)
        call jeveuo(jexnom(nomacr//'.LICA', nomcas), 'L', ialica)
!
        if (rota(1:3) .eq. 'NON') then
            if (zk8(ialich-1+2) .eq. ' ' .and. present(instap)) then
                resuge = zk8(ialich-1+3)
                call jelira(resuge//'.DISC', 'LONMAX', nbinst)
                call jeveuo(resuge//'.DISC', 'L', jinst)
                call jeveuo(resuge//'.DEPL', 'L', jresu)
                do j = 1, nbinst
                    if ((instap-1.e-12) .le. zr(jinst+j-1)) goto 10
                end do
10              continue
                if (j .eq. 1) then
                    do i = nddli+1, nddlt
                        zr(idresl-1+i) = zr(jresu-1+i)
                    end do
                elseif (j .le. nbinst) then
                    inst1 = zr(jinst+j-2)
                    inst2 = zr(jinst+j-1)
                    coecor = (instap-inst1)/(inst2-inst1)
                    do i = nddli+1, nddlt
                        zr(idresl-1+i) = coecor*zr(jresu-1+(j-1)*nddle+i) &
                                         +(1.d0-coecor)*zr(jresu-1+(j-2)*nddle+i)
                    end do
                else
                    do i = nddli+1, nddlt
                        zr(idresl-1+i) = zr(jresu-1+(nbinst-1)*nddle+i)
                    end do
                end if
            else
!         RECOPIE DU VECTEUR DEJA CONDENSE :
                do i = nddli+1, nddlt
                    zr(idresl-1+i) = zr(ialica-1+nddlt+i)
                end do
            end if
!
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
!
            if (zk8(ialich-1+1) (1:3) .eq. 'NON') then
!
!           -- LE CHARGEMENT N'EST PAS "SUIVEUR" :
                call ssvaro(lambda, 'GL', .false._1, 'TOUS', nomacr, &
                            ialica, idresl)
                call jeveuo('&&SSVALV.VALTEMP', 'E', idres2)
                call ssvau1(nomacr, idresl, idres2)
                call ssvaro(lambda, 'LG', .false._1, 'EXTE', nomacr, &
                            idres2, idresl)
!
            else if (zk8(ialich-1+1) (1:3) .eq. 'OUI') then
!
!           -- LE CHARGEMENT EST "SUIVEUR" :
                call ssvaro(lambda, 'LG', .false._1, 'EXTE', nomacr, &
                            ialica+nddlt, idresl)
            else
                call utmess('F', 'SOUSTRUC_47')
            end if
!
        else
            ASSERT(.false.)
        end if
!
!
        idresl = idresl+nddli
!
    end if
!
end subroutine
