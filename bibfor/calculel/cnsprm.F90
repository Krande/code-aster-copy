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

subroutine cnsprm(cns1z, basez, cns2z, iret)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: iret
    character(len=*) :: cns1z, basez, cns2z
! ------------------------------------------------------------------
! BUT : PROJETER UN CHAM_NO_S  SUR UN MAILLAGE MESURE
! ------------------------------------------------------------------
!     ARGUMENTS:
! CNS1Z  IN/JXIN  K19 : CHAM_NO_S A PROJETER
! BASEZ  IN       K1  : BASE DE CREATION POUR CNS2Z : G/V/L
! CNS2Z  IN/JXOUT K19 : CHAM_NO_S RESULTAT DE LA PROJECTION
! IRET   OUT      I   : IRET = 0 : OK / IRET = 1 : PB
! ------------------------------------------------------------------
!    ON NE TRAITE QUE LES CHAMPS REELS (R8) OU COMPLEXES (C16)
!
!     ------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    aster_logical :: newk, axe
    character(len=1) :: base
    character(len=3) :: tsca, dir
    character(len=8) :: ma1, ma2, ma3, nomgd, promes, macrel
    character(len=8) :: model3, kcmp, nono, kcmp2
    character(len=8) :: basemo, mail, numddl, licmp, ksto
    character(len=16) :: k16bid
    character(len=19) :: cns1, cns2, trav
    character(len=24) :: vnoeud, vrange, vnoeum, vrangm, vmes, vsu, k24bid
    character(len=24) :: valk(2), vorien, vref, vrefpm
    integer(kind=8) :: jcns1l, jcns1v, icmp1, icmp2, ier
    integer(kind=8) :: jcns2l, jcns2v, jcns2k, jcns2d, lvsu, lcmp
    integer(kind=8) :: ncmp, ibid, gd, ncmp2, ino2, icmp, ino1, icmpd
    integer(kind=8) :: isma, lori, lref, lrefms
    integer(kind=8) :: iddl, jddl, imod, ipos, iposi, iposj, lnoeud, lrange
    integer(kind=8) :: lnoeum, lrangm, nbmesu, nbord, nddle, lmesu, ltrav
    real(kind=8) :: v1, v2, coef1, valx, valy, valz, eps
    complex(kind=8) :: v1c, v2c
    character(len=8), pointer :: cns1k(:) => null()
    character(len=8), pointer :: nomacr(:) => null()
    character(len=8), pointer :: cns1c(:) => null()
    character(len=8), pointer :: cns2c(:) => null()
    integer(kind=8), pointer :: cns1d(:) => null()
    aster_logical :: lcolle, lcolle2
!     ------------------------------------------------------------------
!
    call jemarq()
!
    cns1 = cns1z
    cns2 = cns2z
    base = basez
!
! RECUPERATION DES OBJETS ET INFORMATIONS DE CNS1 :
!
    call jeveuo(cns1//'.CNSK', 'L', vk8=cns1k)
    call jeveuo(cns1//'.CNSD', 'L', vi=cns1d)
    call jeveuo(cns1//'.CNSC', 'L', vk8=cns1c)
    call jeveuo(cns1//'.CNSV', 'L', jcns1v)
    call jeveuo(cns1//'.CNSL', 'L', jcns1l)
!
! MA1 : MAILLAGE DE LA MODIFICATION
    ma1 = cns1k(1)
    nomgd = cns1k(2)
    ncmp = cns1d(2)
!
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
!
    call jeveuo(ma1//'.NOMACR', 'L', vk8=nomacr)
    lcolle = .false.
    call jeexin(ma1//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
!
    call getvtx(' ', 'SUPER_MAILLE', scal=mail, nbret=ibid)
!
    call jenonu(jexnom(ma1//'.SUPMAIL', mail), isma)
    if (isma .le. 0) then
        valk(1) = mail
        valk(2) = ma1
        call utmess('F', 'CALCULEL5_53', nk=2, valk=valk)
    end if
    macrel = nomacr(isma)
!
    call dismoi('NOM_PROJ_MESU', macrel, 'MACR_ELEM_STAT', repk=promes)
!
! RECUPERATION DES ELEMENTS RELATIFS A LA MESURE
    vnoeud = promes//'.PROJM    .PJMNO'
    vrange = promes//'.PROJM    .PJMRG'
    vorien = promes//'.PROJM    .PJMOR'
    vrefpm = promes//'.PROJM    .PJMRF'
!
    call jeveuo(vnoeud, 'L', lnoeud)
    call jelira(vnoeud, 'LONUTI', nbmesu)
!
! MODEL3 : MODELE MESURE
    call jeveuo(vrange, 'L', lrange)
    call jeveuo(vrefpm, 'L', lrefms)
    k16bid = zk16(lrefms-1+1)
    model3 = k16bid(1:8)
!
    vref = macrel//'.PROJM    .PJMRF'
    call jeveuo(vref, 'L', lref)
    k16bid = zk16(lref-1+3)
    basemo = k16bid(1:8)
!
! POUR LES ORIENTATIONS DES CAPTEURS
    call jeveuo(vorien, 'L', lori)
!
! BASEMO : POUR LA RECUPERATION DU MAILLAGE DU MODELE SUPPORT (MA2)
!
    call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=k24bid)
    numddl = k24bid(1:8)
    call dismoi('NOM_MAILLA', numddl, 'NUME_DDL', repk=ma2)
    lcolle2 = .false.
    call jeexin(ma2//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle2 = .true.
    end if
!
    call dismoi('NOM_MAILLA', model3, 'MODELE', repk=ma3)
!
!  QUELQUES VERIFS :
    if (tsca .ne. 'R' .and. tsca .ne. 'C') then
!        -- ON NE TRAITE QUE LES CHAMPS R/C :
        iret = 1
        goto 999
    end if
!
    call jenonu(jexnom('&CATA.GD.NOMGD', nomgd), gd)
    if (gd .eq. 0) then
        call utmess('F', 'CALCULEL_67', sk=nomgd)
    end if
!
! ALLOCATION DE CNS2 :
    call detrsd('CHAM_NO_S', cns2)
!
! FAIRE APPEL A VRANGE POUR LA LISTE DES CMP MESURE
! ON FAIT L UNION DES CMP DE CNS1 ET VRANGE
!
    licmp = '&&LICMP'
    call wkvect(licmp, 'V V K8', 3*ncmp, lcmp)
    do icmp = 1, ncmp
        zk8(lcmp-1+icmp) = cns1c(icmp)
    end do
    ncmp2 = ncmp
    do iddl = 1, nbmesu
        kcmp = zk8(lrange-1+iddl)
        newk = .true.
        do icmp = 1, ncmp2
            ksto = zk8(lcmp-1+icmp)
            if (kcmp .eq. ksto) newk = .false.
        end do
        if (newk) then
            ncmp2 = ncmp2+1
            zk8(lcmp-1+ncmp2) = kcmp
        end if
    end do
!
    call cnscre(ma3, nomgd, ncmp2, zk8(lcmp), base, &
                cns2)
    call jeveuo(cns2//'.CNSK', 'L', jcns2k)
    call jeveuo(cns2//'.CNSD', 'L', jcns2d)
    call jeveuo(cns2//'.CNSC', 'L', vk8=cns2c)
    call jeveuo(cns2//'.CNSV', 'E', jcns2v)
    call jeveuo(cns2//'.CNSL', 'E', jcns2l)
!
! LISTE DES NOEUDS DU MACRO ELEMENT
    vnoeum = macrel//'.PROJM    .PJMNO'
    vrangm = macrel//'.PROJM    .PJMRG'
    call jeveuo(vnoeum, 'L', lnoeum)
    call jeveuo(vrangm, 'L', lrangm)
    call jelira(vnoeum, 'LONUTI', nddle)
!
! INVERSE DE LA MATRICE DE PASSAGE : VSU = (TIT*PHI)-1
    vsu = macrel//'.PROJM    .PJMIG'
    call jeveuo(vsu, 'L', lvsu)
    call jelira(vsu, 'LONUTI', nbord)
    nbord = nbord/nddle
! NBORD : NOMBRE DE NUMERO D'ORDRE (MODE MESURE)
!
! RECUPERATION DES MODES MESURES
    vmes = macrel//'.PROJM    .PJMMM'
    call jeveuo(vmes, 'L', lmesu)
!
    trav = '&TRAV'
    call wkvect(trav, 'V V R', nbmesu*nddle, ltrav)
! CALCUL DU PRODUIT : PHI*VSU
    do iddl = 1, nbmesu
        do jddl = 1, nddle
            ipos = (jddl-1)*nbmesu+iddl
            zr(ltrav-1+ipos) = 0.d0
            do imod = 1, nbord
                iposi = (imod-1)*nbmesu+iddl
                iposj = (jddl-1)*nbord+imod
                zr(ltrav-1+ipos) = zr(ltrav-1+ipos)+zr(lmesu-1+iposi)*zr(lvsu-1+iposj)
            end do
        end do
    end do
!
!
! INITIALISATION A ZERO
    v2 = 0.d0
    v2c = dcmplx(0.d0, 0.d0)
!
    do iddl = 1, nbmesu
        ino2 = zi(lnoeud-1+iddl)
        do icmp = 1, ncmp2
            zl(jcns2l-1+(ino2-1)*ncmp2+icmp) = .true.
            if (tsca .eq. 'R') then
                zr(jcns2v-1+(ino2-1)*ncmp2+icmp) = v2
            else
                zc(jcns2v-1+(ino2-1)*ncmp2+icmp) = v2c
            end if
        end do
    end do
!
!
! PROJECTION DU CHAMP SUIVANT LA DIRECTION DE MESURE
!
    do iddl = 1, nbmesu
        ino2 = zi(lnoeud-1+iddl)
        kcmp2 = zk8(lrange-1+iddl)
!
        do icmp = 1, ncmp2
            if (cns2c(icmp) .eq. kcmp2) then
                icmp2 = icmp
                goto 60
            end if
        end do
60      continue
!
        v2 = 0.d0
        v2c = dcmplx(0.d0, 0.d0)
!
        do jddl = 1, nddle
            ino1 = zi(lnoeum-1+jddl)
            kcmp = zk8(lrangm-1+jddl)
! ICI ON SUPPOSE QUE LES NOEUDS INTERFACES ONT LE MEME NOM
            nono = int_to_char8(ino1, lcolle2, ma2, 'NOEUD')
            ino1 = char8_to_int(nono, lcolle, ma1, 'NOEUD')
!
            do icmp = 1, ncmp
                if (cns1c(icmp) .eq. kcmp) then
                    icmp1 = icmp
                    goto 160
                end if
            end do
160         continue
!
            coef1 = zr(ltrav-1+(jddl-1)*nbmesu+iddl)
!
            if (tsca .eq. 'R') then
                v1 = zr(jcns1v-1+(ino1-1)*ncmp+icmp1)
                v2 = v2+coef1*v1
            else
                v1c = zc(jcns1v-1+(ino1-1)*ncmp+icmp1)
                v2c = v2c+coef1*v1c
            end if
!
        end do
!
        zl(jcns2l-1+(ino2-1)*ncmp2+icmp2) = .true.
        if (tsca .eq. 'R') then
            zr(jcns2v-1+(ino2-1)*ncmp2+icmp2) = v2
        else
            zc(jcns2v-1+(ino2-1)*ncmp2+icmp2) = v2c
        end if
!
! VERIFICATION SI LA MESURE EST SUR UN DES AXES DE COORDONNEES
! CERTAINS UTILISATEURS SONT HABITUES AUX CMP DX, DY, DZ
        if ((kcmp2 .eq. 'D1') .or. (kcmp2 .eq. 'D2') .or. (kcmp2 .eq. 'D3')) then
            valx = zr(lori-1+(iddl-1)*3+1)
            valy = zr(lori-1+(iddl-1)*3+2)
            valz = zr(lori-1+(iddl-1)*3+3)
!
            valx = abs(valx)
            valy = abs(valy)
            valz = abs(valz)
!
            eps = 1.d2*r8prem()
            axe = .false.
            if ((valy .lt. eps) .and. (valz .lt. eps)) then
                dir = 'DX'
                axe = .true.
            end if
            if ((valx .lt. eps) .and. (valz .lt. eps)) then
                dir = 'DY'
                axe = .true.
            end if
            if ((valx .lt. eps) .and. (valy .lt. eps)) then
                dir = 'DZ'
                axe = .true.
            end if
!
            if (axe) then
                do icmp = 1, ncmp2
                    if (cns2c(icmp) .eq. dir) then
                        icmpd = icmp
                        goto 260
                    end if
                end do
260             continue
!
                zl(jcns2l-1+(ino2-1)*ncmp2+icmpd) = .true.
                if (tsca .eq. 'R') then
                    zr(jcns2v-1+(ino2-1)*ncmp2+icmpd) = v2
                else
                    zc(jcns2v-1+(ino2-1)*ncmp2+icmpd) = v2c
                end if
            end if
        end if
!
    end do
!
    call jedetr(trav)
    call jedetr(licmp)
!
    iret = 0
!
999 continue
    call jedema()
end subroutine
