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
subroutine irsspt(cesz, unite, nbmat, nummai, nbcmp, &
                  nomcmp, lsup, linf, lmax, lmin, &
                  borinf, borsup)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cesexi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*) :: cesz, nomcmp(*)
    integer(kind=8) :: unite, nbmat, nummai(*), nbcmp
    real(kind=8) :: borinf, borsup
    aster_logical :: lsup, linf, lmax, lmin
! ---------------------------------------------------------------------
! BUT: IMPRIMER LES VALEURS MIN/MAX DES COMPOSANTES D'UN CHAM_ELEM_S
!      A DES SOUS-POINTS
! ---------------------------------------------------------------------
!     ARGUMENTS:
! CESZ   IN/JXIN  K19 : SD CHAM_ELEM_S A IMPRIMER
! UNITE  IN       I   : NUMERO DE L'UNITE LOGIQUE D'IMPRESSION
! NBMAT  IN       I   : /0 : ON IMPRIME TOUTES LES MAILLES
! NBMAT  IN       I   : >0 : ON N'IMPRIME QUE LES MAILLES DE NUMMAI(*)
!                            DE 1 A NBMAT
! NUMMAI IN      V(I) : NUMEROS DES MAILLES A IMPRIMER (SI NBMAT >0)
! NBCMP  IN       I   : NOMBRE DE COMPOSANTES A IMPRIMER
! NOMCMP IN       K8  : NOMS DES COMPOSANTES A IMPRIMER
! LSUP   IN       L   : =.TRUE.  INDIQUE PRESENCE D'UNE BORNE SUPERIEURE
! BORSUP IN       R8  : VALEUR DE LA BORNE SUPERIEURE
! LINF   IN       L   : =.TRUE.  INDIQUE PRESENCE D'UNE BORNE INFERIEURE
! BORINF IN       R8  : VALEUR DE LA BORNE INFERIEURE
! LMAX   IN       L   : =.TRUE.  INDIQUE IMPRESSION VALEUR MAXIMALE
! LMIN   IN       L   : =.TRUE.  INDIQUE IMPRESSION VALEUR MINIMALE
! ---------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: jcesd, jcesl, ncmpc, icmp, i
    integer(kind=8) :: ncmp, nbmac, nbma, nbpt, nbsp, j, ipt, isp, iad
    integer(kind=8) :: ispmin, ispmax, ispmi2, ispma2, ispmi3, ispma3
    integer(kind=8) :: iptmin, iptmax, iptmi2, iptma2
    integer(kind=8) :: imamin, imamax, ima
    real(kind=8) :: vspmi3, vspma3, valr
    real(kind=8) :: vptmi2, vptma2, vmamin, vmamax
    character(len=8) :: ma, noma
    character(len=19) :: ces
    aster_logical :: lmamin, lmamax, lptmin, lptmax, lspmin, lspmax
    integer(kind=8), pointer :: num_cmp_cham(:) => null()
    integer(kind=8), pointer :: num_mail_cham(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    character(len=8), pointer :: cesk(:) => null()
    character(len=8), pointer :: cesc(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!     ---------------
    ces = cesz
!
    call jeveuo(ces//'.CESK', 'L', vk8=cesk)
    call jeveuo(ces//'.CESD', 'L', jcesd)
    call jeveuo(ces//'.CESC', 'L', vk8=cesc)
    call jeveuo(ces//'.CESV', 'L', vr=cesv)
    call jeveuo(ces//'.CESL', 'L', jcesl)
    call jelira(ces//'.CESC', 'LONMAX', ncmpc)
    ma = cesk(1)
    nbmac = zi(jcesd-1+1)
!
    AS_ALLOCATE(vi=num_cmp_cham, size=ncmpc)
    AS_ALLOCATE(vi=num_mail_cham, size=nbmac)
!
!
! --- ON RECUPERE LES COMPOSANTES AD-HOC:
!     ------------------------------------
!     SI L'UTILISATEUR A RENSEIGNE NOM_CMP
    if (nbcmp .ne. 0) then
        ncmp = 0
        do i = 1, nbcmp
            icmp = indik8(cesc, nomcmp(i), 1, ncmpc)
            if (icmp .ne. 0) then
                num_cmp_cham(ncmp+1) = icmp
                ncmp = ncmp+1
            end if
        end do
    else
!       SINON TOUT_CMP='OUI'
        do i = 1, ncmpc
            num_cmp_cham(i) = i
        end do
        ncmp = ncmpc
    end if
!
! --- ON RECUPERE LES MAILLES AD-HOC:
!     -------------------------------
!     SI L'UTILISATEUR A RENSEIGNE MAILLE/GROUP_MA
    if (nbmat .ne. 0) then
        do i = 1, nbmat
            num_mail_cham(i) = nummai(i)
        end do
        nbma = nbmat
    else
!        SINON
        do i = 1, nbmac
            num_mail_cham(i) = i
        end do
        nbma = nbmac
    end if
!
!
! --- RECUPERATION DES VALEURS MIN/MAX
!     -------------------------------
!
!     BOUCLE SUR LES COMPOSANTES
    do i = 1, ncmp
        icmp = num_cmp_cham(i)
!
!       LMAxxx : BOOLEEN INDIQUANT LE PREMIER PASSAGE
!       LORS DU DES MAILLES POUR STOCKER LES VALEURS
        lmamin = .true.
        lmamax = .true.
!
!       BOUCLE SUR LES MAILLES
        do j = 1, nbma
            ima = num_mail_cham(j)
            nbpt = zi(jcesd-1+5+4*(ima-1)+1)
            nbsp = zi(jcesd-1+5+4*(ima-1)+2)
!
!         LPTxxx : BOOLEEN INDIQUANT LE PREMIER PASSAGE
!         LORS DU PARCOURT DES POINTS POUR STOCKER LES VALEURS MIN/MAX
            lptmin = .true.
            lptmax = .true.
!
!         BOUCLE SUR LES POINTS
            do ipt = 1, nbpt
!
!           VSPMA3: VALEUR MAX SUR TOUS LES SOUS-POINTS
!           VSPMI3: VALEUR MIN SUR TOUS LES SOUS-POINTS
!           ISPMA3: NUMERO DU SOUS_POINT ASSOCIE A VSPMA3
!           ISPMI3: NUMERO DU SOUS_POINT ASSOCIE A VSPMI3
!
!           LSPxxx : BOOLEEN INDIQUANT LE PREMIER PASSAGE
!           LORS DU PARCOURT DES SOUS-POINTS POUR STOCKER LES VALEURS
                lspmin = .true.
                lspmax = .true.
!
!           BOUCLE SUR LES SOUS-POINTS:
                do isp = 1, nbsp
                    call cesexi('C', jcesd, jcesl, ima, ipt, &
                                isp, icmp, iad)
                    if (iad .gt. 0) then
!
                        valr = cesv(iad)
!
!                SI VALE_MAX
                        if (lmax) then
!                  SI BORNE_SUP
                            if (lsup) then
                                if ((valr-borsup) .gt. 0.d0) goto 80
                            end if
!                  SI BORNE_INF
                            if (linf) then
                                if ((valr-borinf) .lt. 0.d0) goto 80
                            end if
!                  PREMIER PASSAGE
                            if (lspmax) then
                                vspma3 = valr
                                ispma3 = isp
                                lspmax = .false.
                            else
                                if (valr .gt. vspma3) then
                                    vspma3 = valr
                                    ispma3 = isp
                                end if
                            end if
                        end if
!
!                SI VALE_MIN
                        if (lmin) then
!                  SI BORNE_SUP
                            if (lsup) then
                                if ((valr-borsup) .gt. 0.d0) goto 80
                            end if
!                  SI BORNE_INF
                            if (linf) then
                                if ((valr-borinf) .lt. 0.d0) goto 80
                            end if
!                  PREMIER PASSAGE
                            if (lspmin) then
                                vspmi3 = valr
                                ispmi3 = isp
                                lspmin = .false.
                            else
                                if (valr .lt. vspmi3) then
                                    vspmi3 = valr
                                    ispmi3 = isp
                                end if
                            end if
                        end if
!
                    end if
!
!           FIN BOUCLE SUR LES SOUS-POINTS
80                  continue
                end do
!
!           VPTMA2: VALEUR MAX SUR TOUS LES POINTS
!           VPTMI2: VALEUR MIN SUR TOUS LES POINTS
!           IPTMA2: NUMERO DU POINT ASSOCIE A VPTMA2
!           IPTMI2: NUMERO DU POINT ASSOCIE A VPTMI2
!           ISPMA2: NUMERO DU SOUS_POINT ASSOCIE A IPTMA2
!           ISPMI2: NUMERO DU SOUS_POINT ASSOCIE A IPTMI2
!
!           SI VALE_MAX
                if (lmax .and. .not. lspmax) then
!             PREMIER PASSAGE
                    if (lptmax) then
                        iptma2 = ipt
                        vptma2 = vspma3
                        ispma2 = ispma3
                        lptmax = .false.
                    else
!               ON REACTUALISE LA VALEUR MAX
                        if (vptma2 .lt. vspma3) then
                            vptma2 = vspma3
                            iptma2 = ipt
                            ispma2 = ispma3
                        end if
                    end if
                end if
!
!           SI VALE_MIN
                if (lmin .and. .not. lspmin) then
!             PREMIER PASSAGE
                    if (lptmin) then
                        iptmi2 = ipt
                        vptmi2 = vspmi3
                        ispmi2 = ispmi3
                        lptmin = .false.
                    else
!               ON REACTUALISE LA VALEUR MIN
                        if (vptmi2 .gt. vspmi3) then
                            vptmi2 = vspmi3
                            iptmi2 = ipt
                            ispmi2 = ispmi3
                        end if
                    end if
                end if
!
!         FIN BOUCLE SUR LES POINTS
            end do
!
!         VMAMAX: VALEUR MAX SUR TOUTES LES MAILLES
!         VMAMIN: VALEUR MIN SUR TOUTES LES MAILLES
!         IMAMAX: NUMERO DE LA MAILLE ASSOCIEE A VMAMAX
!         IMAMIN: NUMERO DE LA MAILLE ASSOCIEE A VMAMIN
!         IPTMAX: NUMERO DU POINT ASSOCIE A IMAMAX
!         IPTMIN: NUMERO DU POINT ASSOCIE A IMAMIN
!         ISPMAX: NUMERO DU SOUS_POINT ASSOCIE A IPTMAX
!         ISPMIN: NUMERO DU SOUS_POINT ASSOCIE A IPTMIN
!
!         SI VALE_MAX
            if (lmax .and. .not. lptmax) then
!           PREMIER PASSAGE
                if (lmamax) then
                    imamax = ima
                    vmamax = vptma2
                    iptmax = iptma2
                    ispmax = ispma2
                    lmamax = .false.
                else
!             ON REACTUALISE LA VALEUR MAX
                    if (vmamax .lt. vptma2) then
                        vmamax = vptma2
                        iptmax = iptma2
                        ispmax = ispma2
                        imamax = ima
                    end if
                end if
            end if
!
!         SI VALE_MIN
            if (lmin .and. .not. lptmin) then
!           PREMIER PASSAGE
                if (lmamin) then
                    imamin = ima
                    vmamin = vptmi2
                    iptmin = iptmi2
                    ispmin = ispmi2
                    lmamin = .false.
                else
!             ON REACTUALISE LA VALEUR MIN
                    if (vmamin .gt. vptmi2) then
                        vmamin = vptmi2
                        iptmin = iptmi2
                        ispmin = ispmi2
                        imamin = ima
                    end if
                end if
            end if
!
!     FIN BOUCLE SUR LES MAILLES
        end do
!
!
!     IMPRESSIONS
!     -----------
!
        if (lmax .and. .not. lmamax) then
            noma = int_to_char8(imamax)
            write (unite, *) ' '
            write (unite, 2000) 'LA VALEUR MAXIMALE DE ', cesc(1+icmp- &
                                                               1), 'EST: ', vmamax
            write (unite, 2001) 'OBTENUE DANS LA MAILLE ', noma,&
     &    'AU SOUS_POINT ', ispmax, ' DU POINT ', iptmax
        end if
!
        if (lmin .and. .not. lmamin) then
            noma = int_to_char8(imamin)
            write (unite, *) ' '
            write (unite, 2000) 'LA VALEUR MINIMALE DE ', cesc(1+icmp- &
                                                               1), 'EST: ', vmamin
            write (unite, 2001) 'OBTENUE DANS LA MAILLE ', noma,&
     &    'AU SOUS_POINT ', ispmin, ' DU POINT ', iptmin
        end if
!
    end do
!
2000 format(3(a), e12.5)
2001 format(3(a), i3, a, i3)
!
    AS_DEALLOCATE(vi=num_cmp_cham)
    AS_DEALLOCATE(vi=num_mail_cham)
!
    call jedema()
!
end subroutine
