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

subroutine cnsprj(cns1z, correz, basez, cns2z, iret)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*) :: cns1z, correz, basez, cns2z
    integer(kind=8) :: iret
! ------------------------------------------------------------------
! BUT : PROJETER UN CHAM_NO_S  SUR UN AUTRE MAILLAGE
! ------------------------------------------------------------------
!     ARGUMENTS:
! CNS1Z  IN/JXIN  K19 : CHAM_NO_S A PROJETER
! CORREZ IN/JXIN  K16 : NOM DE LA SD CORRESP_2_MAILLA
! BASEZ  IN       K1  : BASE DE CREATION POUR CNS2Z : G/V
! CNS2Z  IN/JXOUT K19 : CHAM_NO_S RESULTAT DE LA PROJECTION
! IRET   OUT      I   : CODE RETOUR :
!                       0 -> OK
!                       1 -> ECHEC DE LA PROJECTION
! ------------------------------------------------------------------
!    ON NE TRAITE QUE LES CHAMPS REELS (R8) OU COMPLEXES (C16)
!
!
! REMARQUE :
!   LA PROJECTION EST APPROCHEE DANS LE CAS OU TOUS LES NOEUDS
!   DE LA MAILLE M1 NE PORTENT PAS LES MEMES DDLS.
!   LE TRAITEMENT EST EXPLIQUE DANS LA REPONSE A LA FICHE 9259
!
    character(len=24) :: valk(3)
!     ------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=1) :: base
    character(len=3) :: tsca
    character(len=8) :: ma1, ma2, nomgd, nomcmp, nomno2
    character(len=16) :: corres
    character(len=19) :: cns1, cns2
    integer(kind=8) :: jcns1l, jcns1v
    integer(kind=8) :: jcns2c, jcns2l, jcns2v, jcns2k
    integer(kind=8) :: nbno1, ncmp, gd, nbno2
    integer(kind=8) :: idecal, ino2, icmp, ico1, ico2, ino1, nuno1, kalarm
    real(kind=8) :: v1, v2, coef1, coetot, vrmoy
    complex(kind=8) :: v1c, v2c, vcmoy
    aster_logical :: lexact
    integer(kind=8), pointer :: pjef_nu(:) => null()
    character(len=8), pointer :: cns1k(:) => null()
    integer(kind=8), pointer :: pjef_nb(:) => null()
    character(len=8), pointer :: cns1c(:) => null()
    character(len=24), pointer :: pjxx_k1(:) => null()
    real(kind=8), pointer :: pjef_cf(:) => null()
    integer(kind=8), pointer :: cns1d(:) => null()
    integer(kind=8), pointer :: cns2d(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
    cns1 = cns1z
    cns2 = cns2z
    base = basez
    corres = correz
    iret = 0
    kalarm = 0
!
!
!------------------------------------------------------------------
!     1- RECUPERATION DES OBJETS ET INFORMATIONS DE CNS1 :
!     ----------------------------------------------------
!
    call jeveuo(cns1//'.CNSK', 'L', vk8=cns1k)
    call jeveuo(cns1//'.CNSD', 'L', vi=cns1d)
    call jeveuo(cns1//'.CNSC', 'L', vk8=cns1c)
    call jeveuo(cns1//'.CNSV', 'L', jcns1v)
    call jeveuo(cns1//'.CNSL', 'L', jcns1l)
!
    ma1 = cns1k(1)
    nomgd = cns1k(2)
    nbno1 = cns1d(1)
    ncmp = cns1d(2)
!
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
!
!
!------------------------------------------------------------------
!     2- RECUPERATION DES OBJETS ET INFORMATIONS DE CORRES :
!     ----------------------------------------------------
    call jeveuo(corres//'.PJXX_K1', 'L', vk24=pjxx_k1)
    call jeveuo(corres//'.PJEF_NB', 'L', vi=pjef_nb)
    call jeveuo(corres//'.PJEF_NU', 'L', vi=pjef_nu)
    call jeveuo(corres//'.PJEF_CF', 'L', vr=pjef_cf)
!
    ma2 = pjxx_k1(2)
!
!
!------------------------------------------------------------------
!     3- QUELQUES VERIFS :
!     ------------------------
    if (tsca .ne. 'R' .and. tsca .ne. 'C') then
!        -- ON NE TRAITE QUE LES CHAMPS R/C :
        iret = 1
        goto 60
!
    end if
!     TEST SUR IDENTITE DES 2 MAILLAGES
    ASSERT(pjxx_k1(1) .eq. ma1)
!
    call jenonu(jexnom('&CATA.GD.NOMGD', nomgd), gd)
    if (gd .eq. 0) then
        call utmess('F', 'CALCULEL_67', sk=nomgd)
    end if
!
!
!------------------------------------------------------------------
!     4- ALLOCATION DE CNS2 :
!     ------------------------
    call detrsd('CHAM_NO_S', cns2)
    call cnscre(ma2, nomgd, ncmp, cns1c, base, &
                cns2)
    call jeveuo(cns2//'.CNSK', 'L', jcns2k)
    call jeveuo(cns2//'.CNSD', 'L', vi=cns2d)
    call jeveuo(cns2//'.CNSC', 'L', jcns2c)
    call jeveuo(cns2//'.CNSV', 'E', jcns2v)
    call jeveuo(cns2//'.CNSL', 'E', jcns2l)
!
    nbno2 = cns2d(1)
!
!------------------------------------------------------------------
!     5- CALCUL DES VALEURS DE CNS2 :
!     -------------------------------
    idecal = 0
    do ino2 = 1, nbno2
        nbno1 = pjef_nb(ino2)
        if (nbno1 .eq. 0) goto 50
        do icmp = 1, ncmp
!
!          -- ON COMPTE (ICO1) LES NOEUDS PORTANT LE DDL :
!             ON COMPTE AUSSI (ICO2) CEUX DONT LE COEF EST > 0
!             ON CALCULE LA VALEUR MOYENNE SUR LA MAILLE (VXMOY)
!             ON CALCULE LA SOMME DES COEF > 0 (COETOT)
            ico1 = 0
            ico2 = 0
            vrmoy = 0.d0
            vcmoy = dcmplx(0.d0, 0.d0)
            coetot = 0.d0
            do ino1 = 1, nbno1
                nuno1 = pjef_nu(1+idecal-1+ino1)
                coef1 = pjef_cf(1+idecal-1+ino1)
                if (zl(jcns1l-1+(nuno1-1)*ncmp+icmp)) then
                    ico1 = ico1+1
                    if (coef1 .gt. 0.d0) then
                        ico2 = ico2+1
                        coetot = coetot+coef1
                    end if
                    if (tsca .eq. 'R') then
                        vrmoy = vrmoy+zr(jcns1v-1+(nuno1-1)*ncmp+icmp)
!
                    else
                        vcmoy = vcmoy+zc(jcns1v-1+(nuno1-1)*ncmp+icmp)
                    end if
                end if
            end do
            if (ico1 .eq. 0) goto 40
            zl(jcns2l-1+(ino2-1)*ncmp+icmp) = .true.
!
!
!         -- SI COETOT EST FAIBLE, LA PROJECTION N'EST PAS PRECISE :
!            L'EMISSION DE L'ALARME EST COUTEUSE, ON LA LIMITE :
            if (coetot .lt. 1.d-3 .and. kalarm .le. 6) then
                kalarm = kalarm+1
                nomno2 = int_to_char8(ino2)
                nomcmp = cns1c(icmp)
                valk(1) = nomgd
                valk(2) = nomno2
                valk(3) = nomcmp
                call utmess('A', 'CALCULEL4_9', nk=3, valk=valk)
            end if
!
!
!          -- 3 CAS DE FIGURE POUR L'INTERPOLATION :
!          ----------------------------------------
            if (ico1 .eq. nbno1) then
!            1 : NORMAL ON PREND TOUS LES NOEUDS N1
                lexact = .true.
                coetot = 1.d0
!
            else if (ico2 .gt. 0) then
!            2 : ON PREND LES NOEUDS N1 DE COEF > 0
                lexact = .false.
!
            else
!            3 : ON FAIT UNE MOYENNE ARITHMETIQUE
                if (tsca .eq. 'R') then
                    zr(jcns2v-1+(ino2-1)*ncmp+icmp) = vrmoy/ico1
!
                else
                    zc(jcns2v-1+(ino2-1)*ncmp+icmp) = vcmoy/ico1
                end if
                goto 40
!
            end if
!
!
            if (tsca .eq. 'R') then
                v2 = 0.d0
                do ino1 = 1, nbno1
                    nuno1 = pjef_nu(1+idecal-1+ino1)
                    coef1 = pjef_cf(1+idecal-1+ino1)
                    if (zl(jcns1l-1+(nuno1-1)*ncmp+icmp)) then
                        if (lexact .or. coef1 .gt. 0) then
                            v1 = zr(jcns1v-1+(nuno1-1)*ncmp+icmp)
                            v2 = v2+coef1*v1
                        end if
                    end if
                end do
                zr(jcns2v-1+(ino2-1)*ncmp+icmp) = v2/coetot
!
            else if (tsca .eq. 'C') then
                v2c = dcmplx(0.d0, 0.d0)
                do ino1 = 1, nbno1
                    nuno1 = pjef_nu(1+idecal-1+ino1)
                    coef1 = pjef_cf(1+idecal-1+ino1)
                    if (zl(jcns1l-1+(nuno1-1)*ncmp+icmp)) then
                        if (lexact .or. coef1 .gt. 0) then
                            v1c = zc(jcns1v-1+(nuno1-1)*ncmp+icmp)
                            v2c = v2c+coef1*v1c
                        end if
                    end if
                end do
                zc(jcns2v-1+(ino2-1)*ncmp+icmp) = v2c/coetot
            end if
40          continue
        end do
        idecal = idecal+nbno1
50      continue
    end do
!
60  continue
    call jedema()
end subroutine
