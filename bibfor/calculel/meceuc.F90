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

! aslint: disable=W1306
!
subroutine meceuc(stop, option, caraez, ligrel, &
                  nin, lchin, lpain, nou, lchou, &
                  lpaou, base)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assach.h"
#include "asterfort/assert.h"
#include "asterfort/barych.h"
#include "asterfort/calcul.h"
#include "asterfort/chlici.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/sepach.h"
!
    integer(kind=8) :: nin, nou
    character(len=1) :: stop
    character(len=8) :: carael
    character(len=*) :: base, option
    character(len=*) :: lchin(*), lchou(*), lpain(*), lpaou(*), ligrel, caraez
!
! --------------------------------------------------------------------------------------------------
!
!  BUT : CETTE ROUTINE EST UN INTERMEDIAIRE VERS LA ROUTINE CALCUL.F
!        POUR CERTAINES OPTIONS DONT LE RESULTAT EST COMPLEXE, ON
!        APPELLE 2 FOIS CALCUL EN AYANT SEPARE LES CHAMPS COMPLEXES"IN"
!        EN 2 : PARTIE REELLE ET PARTIE IMAGINAIRE
!
! --------------------------------------------------------------------------------------------------
!
!     ENTREES:
!        STOP   :  /'S' : ON S'ARRETE SI AUCUN ELEMENT FINI DU LIGREL
!                         NE SAIT CALCULER L'OPTION.
!                  /'C' : ON CONTINUE SI AUCUN ELEMENT FINI DU LIGREL
!                         NE SAIT CALCULER L'OPTION. IL N'EXISTE PAS DE
!                         CHAMP "OUT" DANS CE CAS.
!        OPTIO  :  NOM D'1 OPTION
!        LIGREL :  NOM DU LIGREL SUR LEQUEL ON DOIT FAIRE LE CALCUL
!        NIN    :  NOMBRE DE CHAMPS PARAMETRES "IN"
!        NOU    :  NOMBRE DE CHAMPS PARAMETRES "OUT"
!        LCHIN  :  LISTE DES NOMS DES CHAMPS "IN"
!        LCHOU  :  LISTE DES NOMS DES CHAMPS "OUT"
!        LPAIN  :  LISTE DES NOMS DES PARAMETRES "IN"
!        LPAOU  :  LISTE DES NOMS DES PARAMETRES "OUT"
!        BASE   :  'G' , 'V' OU 'L'
!
!     SORTIES:
!       ALLOCATION ET CALCUL DES OBJETS CORRESPONDANT AUX CHAMPS "OUT"
!     CETTE ROUTINE MET EN FORME LES CHAMPS EVENTUELLEMENT COMPLEXE
!     POUR L APPEL A CALCUL
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: chdecr(nin), chdeci(nin), ch19, chr, chi, ch1, ch2
    character(len=19) :: lchinr(nin), lchini(nin)
    character(len=16) :: optio2
    character(len=8) :: nomgd
    integer(kind=8) :: k, iexi, iexi1, iexi2
    integer(kind=8) :: inddec(nin)
    aster_logical :: lcmplx, lsspt, ldbg, lopdec
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    optio2 = option
    carael = caraez
    ch1 = '&&MECEUC.CH1'
    ch2 = '&&MECEUC.CH2'
    ldbg = .false.
!
!
!     -- 0. LA ROUTINE N'EST UTILE QUE POUR CERTAINES OPTIONS :
!     ---------------------------------------------------------
    lopdec = .false.
    if (option .eq. 'ECIN_ELEM') lopdec = .true.
    if (option .eq. 'EFGE_ELGA') lopdec = .true.
    if (option .eq. 'EFGE_ELNO') lopdec = .true.
    if (option .eq. 'ENEL_ELGA') lopdec = .true.
    if (option .eq. 'ENEL_ELNO') lopdec = .true.
    if (option .eq. 'EPOT_ELEM') lopdec = .true.
    if (option .eq. 'EPSI_ELGA') lopdec = .true.
    if (option .eq. 'EPSI_ELNO') lopdec = .true.
    if (option .eq. 'SIEF_ELGA') lopdec = .true.
    if (option .eq. 'SIEF_ELNO') lopdec = .true.
    if (option .eq. 'SIGM_ELGA') lopdec = .true.
    if (option .eq. 'SIGM_ELNO') lopdec = .true.
    if (option .eq. 'SIPM_ELNO') lopdec = .true.
    if (option .eq. 'SIPO_ELNO') lopdec = .true.
    if (.not. lopdec) then
        call calcul(stop, optio2, ligrel, nin, lchin, &
                    lpain, nou, lchou, lpaou, base, &
                    'OUI')
        goto 999
    end if
!
!     -- 1. Y-A-T-IL DES CHAMPS "IN" COMPLEXES ?
!           SI OUI, IL FAUT LES DECOUPER
!     -----------------------------------------------------------
    lcmplx = .false.
    do k = 1, nin
        inddec(k) = 0
        if (lpain(k) .eq. ' ') cycle
        ch19 = lchin(k)
        if (ch19 .eq. ' ') cycle
        if (ldbg) call chlici(ch19, 19)
        call exisd('CHAMP', ch19, iexi)
        if (iexi .eq. 0) cycle
        call dismoi('NOM_GD', ch19, 'CHAMP', repk=nomgd)
!        -- MECHPO CREE PARFOIS UN CHAMP DE FORC_C
!           IL M'EMBETE ! COMMENT SAVOIR S'IL EST PERTINENT ?
        if (nomgd .eq. 'FORC_C') cycle
        if (nomgd(5:6) .eq. '_C') then
            lcmplx = .true.
            inddec(k) = 1
            chr = '&&MECEUC.CHXX.R'
            chi = '&&MECEUC.CHXX.I'
            call codent(k, 'D0', chr(12:13))
            call codent(k, 'D0', chi(12:13))
            call sepach(carael, ch19, 'V', chr, chi)
            chdecr(k) = chr
            chdeci(k) = chi
        end if
    end do
!
!
!     -- 2. S'IL N'Y A AUCUN CHAMP COMPLEXE, C'EST FACILE :
!     -------------------------------------------------------
    if (.not. lcmplx) then
        call calcul(stop, optio2, ligrel, nin, lchin, &
                    lpain, nou, lchou, lpaou, base, &
                    'OUI')
        goto 999
    end if
!
!
!     -- 3. LE CHAMP "OUT" EST-IL A SOUS-POINTS ?
!     -------------------------------------------
    ASSERT(nou .eq. 1)
    call exisd('CHAM_ELEM_S', lchou(1), iexi)
    lsspt = (iexi .ne. 0)
!
!
!     -- 4.0 ON PREPARE LCHINR ET LCHINI :
!     ------------------------------------
    do k = 1, nin
        if (inddec(k) .eq. 0) then
            lchinr(k) = lchin(k)
            lchini(k) = lchin(k)
        else
            lchinr(k) = chdecr(k)
            lchini(k) = chdeci(k)
        end if
    end do
!
!
!     -- 4.1 APPEL A CALCUL AVEC LES PARTIES REELLES :
!     ------------------------------------------------
    if (lsspt) call copisd('CHAM_ELEM_S', 'V', lchou(1), ch1)
    call calcul(stop, optio2, ligrel, nin, lchinr, &
                lpain, nou, ch1, lpaou, 'V', &
                'OUI')
!
!
!     -- 4.2 APPEL A CALCUL AVEC LES PARTIES IMAGINAIRES :
!     ----------------------------------------------------
    if (lsspt) call copisd('CHAM_ELEM_S', 'V', lchou(1), ch2)
    call calcul(stop, optio2, ligrel, nin, lchini, &
                lpain, nou, ch2, lpaou, 'V', &
                'OUI')
!
!
!     -- 4.3 SI STOP='C' ET QUE CH1 ET CH2 N'EXISTENT PAS :
!     -----------------------------------------------------
    if (stop .eq. 'C') then
        call exisd('CHAM_ELEM', ch1, iexi1)
        call exisd('CHAM_ELEM', ch2, iexi2)
        if (iexi1 .eq. 0) then
            ASSERT(iexi2 .eq. 0)
            goto 888
        end if
    end if
!
!
!     -- 6.  ASSEMBLAGE (R,I) OU CUMUL (R+I) :
!     -----------------------------------------
    if ((optio2 .eq. 'SIEF_ELNO') .or. (optio2 .eq. 'SIGM_ELGA') .or. (optio2 .eq. 'SIGM_ELNO') &
        .or. (optio2 .eq. 'EFGE_ELGA') .or. (optio2 .eq. 'EFGE_ELNO') .or. &
        (optio2 .eq. 'SIPM_ELNO') .or. (optio2 .eq. 'SIPO_ELNO') .or. (optio2 .eq. 'EPSI_ELNO') &
        .or. (optio2 .eq. 'EPSI_ELGA') .or. (optio2 .eq. 'STRX_ELGA') .or. &
        (optio2 .eq. 'SIEF_ELGA')) then
        call assach(ch1, ch2, base, lchou(1))
    elseif ((optio2 .eq. 'EPOT_ELEM') .or. (optio2 .eq. 'ENEL_ELGA') .or. &
            (optio2 .eq. 'ENEL_ELNO') .or. (optio2 .eq. 'ECIN_ELEM')) then
        call barych(ch1, ch2, 1.d0, 1.d0, lchou(1), &
                    'G')
    else
        ASSERT(.false.)
    end if
!
!     -- 7. MENAGE :
!     --------------
888 continue
    do k = 1, nin
        if (inddec(k) .ne. 0) then
            call detrsd('CHAMP', chdecr(k))
            call detrsd('CHAMP', chdeci(k))
        end if
    end do
!
    call detrsd('CHAMP', ch1)
    call detrsd('CHAMP', ch2)
    call detrsd('CHAM_ELEM_S', ch1)
    call detrsd('CHAM_ELEM_S', ch2)
!
999 continue
    call jedema()
end subroutine
