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

subroutine nocart(carte, code, ncmp, groupma, mode, &
                  nma, limano, limanu, ligrel, rapide, &
                  jdesc, jnoma, jncmp, jnoli, jvale, &
                  jvalv, jnocmp, ncmpmx, nec, ctype, &
                  jlima0, jlimac, lontav)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
!
#include "asterfort/agcart.h"
#include "asterfort/assert.h"
#include "asterfort/editgd.h"
#include "asterfort/jeagco.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/jecroc.h"
#include "asterfort/nbec.h"
#include "asterfort/char8_to_int.h"
!
    character(len=*), intent(in) :: carte
    integer(kind=8), intent(in) :: code
    integer(kind=8), intent(in) :: ncmp
    character(len=*), intent(in), optional :: groupma
    character(len=*), intent(in), optional :: mode
    integer(kind=8), intent(in), optional :: nma
    character(len=*), intent(in), optional :: limano(*)
    integer(kind=8), intent(in), optional :: limanu(*)
    character(len=*), intent(in), optional :: ligrel
!
!
!   -- arguments optionnels pour gagner du CPU :
    character(len=3), intent(in), optional :: rapide
    integer(kind=8), intent(inout), optional :: jdesc
    integer(kind=8), intent(inout), optional :: jnoma
    integer(kind=8), intent(inout), optional :: jncmp
    integer(kind=8), intent(inout), optional :: jnoli
    integer(kind=8), intent(inout), optional :: jvale
    integer(kind=8), intent(inout), optional :: jvalv
    integer(kind=8), intent(in), optional :: jnocmp
    integer(kind=8), intent(in), optional :: ncmpmx
    integer(kind=8), intent(in), optional :: nec
    character(len=8), intent(in), optional :: ctype
    integer(kind=8), intent(inout), optional :: jlima0
    integer(kind=8), intent(inout), optional :: jlimac
    integer(kind=8), intent(inout), optional :: lontav
! ----------------------------------------------------------------------
!     entrees:
!     --------
!     carte : nom de la carte a enrichir
!
!     code : 1: 'tout(es)' les mailles du maillage.
!           -1: 'tout(es)' les mailles suppl. d'1 ligrel.
!            2: 1 groupe_ma du maillage.
!            3: 1 liste de mailles du maillage.
!           -3: 1 liste de mailles tardives d'1 ligrel.
!
!     ncmp : nombre des composantes decrites
!            dans chin.ncmp et chin.valv
!            -- remarque : on peut sur-dimensionner ncmp a condition
!                          que la liste des noms de cmps contienne des
!                          "blancs". les cmps reellement notees sont
!                          les composantes non-blanches.
!
!     groupma : nom d' 1 group_ma du maillage
!              ( uniquement si code= 2)
!
!     mode : 'nom' ou 'num' :
!             si 'nom' on utilise la liste limano (noms des mailles)
!                 ( uniquement si code= 3)
!             si 'num' on utilise la liste limanu (numero des mailles)
!                 ( uniquement si code= 3 ou -3)
!
!     nma  : nombre de mailles dans limano ou limanu
!                 ( uniquement si code= 3 ou -3)
!
!     limano : noms des mailles du groupe_tardif (code=3)
!
!     limanu : numeros des mailles du groupe_tardif (code=3 ou -3)
!
!     ligrel : nom du ligrel ou sont eventuellement definies les mailles
!         tardives que l'on veut affecter.
!         ligrel est utilise uniquement si code=-3 ou code=-1
!
!     Tous les arguments suivants sont facultatifs.
!     ---------------------------------------------
!     Ils ne sont pas documentes. Il ne doivent etre renseignes que dans le
!     cas ou la routine nocart est appelee de (trop) nombreuses fois.
!     Exemple d'utilisation : aflrch.F90
!
!     rapide: 'OUI' / 'NON'
!     jdesc : adresse de l'objet chin.DESC
!     jnoma : adresse de l'objet chin.NOMA
!     jnoli : adresse de l'objet chin.NOLI
!     jvale : adresse de l'objet chin.VALE
!     jvalv : adresse de l'objet chin.VALV
!     jlima0, jlimac : adresses pour l'objet l'objet chin.LIMA
!     + jncmp, jnocmp, ncmpmx, nec, ctype, lontav
!
!     Attention : si rapide='OUI', il faut que l'appelant fasse appel a
!                 jeecra / NUTIOC apres le denier appel a nocart
!                 (objet .LIMA)
!
!     sorties:
!        on enrichit le contenu de la carte chin
!
! ----------------------------------------------------------------------
    character(len=24) :: groupe
    character(len=19) :: chin, nomlig
    character(len=8) :: nomail
    integer(kind=8) :: nedit, ngdmx, iaddg, gr, dim, i, numero, gd
    integer(kind=8) :: jlima02, jlimac2, jlima
    integer(kind=8) :: ldim, jdesc2, jnoma2, jncmp2, jnoli2, jvale2, jvalv2
    integer(kind=8) :: jnocmp2, ncmpmx2, nec2
    character(len=8) :: ma, base, mode2, ctype2
    character(len=24) :: lima, trav
    aster_logical :: laggr, lrapid
    integer(kind=8) :: lontap, lontav2, jbid
!-----------------------------------------------------------------------
    chin = carte
    lima = chin//'.LIMA'
!
    laggr = .false.
    nomlig = ' '
    groupe = ' '
!
    lrapid = .false.
    if (present(rapide)) then
        if (rapide .eq. 'OUI') lrapid = .true.
    end if
!
!
!   -- gestion des derniers arguments facultatifs (pour gains CPU) :
!   ----------------------------------------------------------------
    if (.not. lrapid) then
        call jeveuo(chin//'.DESC', 'E', jdesc2)
        call jeveuo(chin//'.NOMA', 'E', jnoma2)
        call jeveuo(chin//'.NCMP', 'E', jncmp2)
        call jeveuo(chin//'.NOLI', 'E', jnoli2)
        call jeveuo(chin//'.VALE', 'E', jvale2)
        call jeveuo(chin//'.VALV', 'E', jvalv2)
        call jelira(chin//'.VALV', 'TYPELONG', cval=ctype2)
        call jeveuo(chin//'.LIMA', 'E', jlima02)
        call jelira(chin//'.LIMA', 'LONT', lontav2)
        call jeveuo(jexatr(chin//'.LIMA', 'LONCUM'), 'E', jlimac2)
        gd = zi(jdesc2-1+1)
        nec2 = nbec(gd)
        call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', jnocmp2)
        call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx2)
    else
        ASSERT(present(jdesc))
        ASSERT(present(jnoma))
        ASSERT(present(jncmp))
        ASSERT(present(jnoli))
        ASSERT(present(jvale))
        ASSERT(present(jvalv))
        ASSERT(present(ctype))
        ASSERT(present(lontav))
        ASSERT(present(jnocmp))
        ASSERT(present(ncmpmx))
        ASSERT(present(nec))
        ASSERT(present(jlima0))
        ASSERT(present(jlimac))
        jdesc2 = jdesc
        jnoma2 = jnoma
        jncmp2 = jncmp
        jnoli2 = jnoli
        jvale2 = jvale
        jvalv2 = jvalv
        ctype2 = ctype
        jlima02 = jlima0
        jlimac2 = jlimac
        lontav2 = lontav
        jnocmp2 = jnocmp
        ncmpmx2 = ncmpmx
        nec2 = nec
    end if
!
!
!   -- verification des arguments :
!   --------------------------------------------
    if (code .eq. 1) then
    else if (code .eq. -1) then
        ASSERT(present(ligrel))
        nomlig = ligrel
    else if (code .eq. 2) then
        ASSERT(present(groupma))
        groupe = groupma
    else if (code .eq. 3) then
        ASSERT(present(mode))
        ASSERT(mode .eq. 'NUM' .or. mode .eq. 'NOM')
        mode2 = mode
        ASSERT(present(nma))
        ASSERT(nma .gt. 0)
        if (mode2 .eq. 'NUM') then
            ASSERT(present(limanu))
        else
            ASSERT(present(limano))
        end if
    else if (code .eq. -3) then
        ASSERT(present(ligrel))
        nomlig = ligrel
        ASSERT(present(nma))
        ASSERT(nma .gt. 0)
        ASSERT(present(limanu))
        mode2 = 'NUM'
    else
        ASSERT(.false.)
    end if
!
!
    ma = zk8(jnoma2-1+1)
    ngdmx = zi(jdesc2-1+2)
    nedit = zi(jdesc2-1+3)+1
!
    if (lrapid .and. nedit .eq. 1) zi(jlimac2) = 1
!
!
!   -- faut-il agrandir la carte ?
!   -------------------------------
    if (nedit .gt. ngdmx) then
        laggr = .true.
        ngdmx = 2*ngdmx
        call agcart(ngdmx, chin)
! agcart a deplace les objets .DESC, .NOLI et .VALE :
        call jeveuo(chin//'.DESC', 'E', jdesc2)
        call jeveuo(chin//'.NOLI', 'E', jnoli2)
        call jeveuo(chin//'.VALE', 'E', jvale2)
        if (present(jdesc)) jdesc = jdesc2
        if (present(jnoli)) jnoli = jnoli2
        if (present(jvale)) jvale = jvale2
    end if
!
    zi(jdesc2-1+3) = nedit
!
    if ((code .eq. -1) .or. (code .eq. -3)) then
        ASSERT(nomlig(1:8) .ne. ' ')
        zk24(jnoli2-1+nedit) = nomlig
    end if
!
!   -- appel a editgd qui remplit .vale et met a jour  le desc_grandeur
    iaddg = 3+2*ngdmx+(nedit-1)*nec2+1
    call editgd(ncmp, nedit, zi(jdesc2-1+iaddg), ncmpmx2, ctype2, &
                jnocmp2, jncmp2, jvalv2, jvale2)
!
!
!   -- mise a jour de .desc :
!   --------------------------
    zi(jdesc2-1+3+2*nedit-1) = code
    dim = 0
    if (abs(code) .eq. 1) then
!        -- on note le numero d'entite conventionnel relatif a "TOUT":   9999
        zi(jdesc2-1+3+2*nedit) = 9999
    else if (code .eq. 2) then
        call jenonu(jexnom(ma//'.GROUPEMA', groupma), gr)
        zi(jdesc2-1+3+2*nedit) = gr
    else if (abs(code) .eq. 3) then
        zi(jdesc2-1+3+2*nedit) = nedit
        dim = nma
    else
        ASSERT(.false.)
    end if
!
!
!   -- mise a jour de lima :
!   ------------------------
!
!   -- faut-il agrandir .lima ?
    lontap = zi(jlimac2-1+nedit)+max(dim, 1)
    if (lontap .gt. lontav2) then
        laggr = .true.
        lontap = max(2*lontav2, lontap)
    end if
    if (laggr) then
        trav = chin//'.TRAV'
        call jeveuo(jexatr(lima, 'LONCUM'), 'E', jbid)
        call jedupo(lima, 'V', trav, .false._1)
        call jeveuo(jexatr(trav, 'LONCUM'), 'E', jbid)
        call jelira(lima, 'CLAS', cval=base)
        call jedetr(lima)
        call jeagco(trav, lima, ngdmx, lontap, base)
        call jeveuo(jexatr(lima, 'LONCUM'), 'E', jbid)
        call jedetr(trav)
        call jeveuo(lima, 'E', jlima02)
        if (present(jlima0)) jlima0 = jlima02
        call jeveuo(jexatr(lima, 'LONCUM'), 'E', jlimac2)
        if (present(jlimac)) jlimac = jlimac2
        call jelira(lima, 'LONT', lontav2)
        if (present(lontav)) lontav = lontav2
    end if
!
!
! -- mise a jour des longueurs cumulees pour lima :
    ldim = max(dim, 1)
    if (.not. lrapid) then
        call jeecra(jexnum(lima, nedit), 'LONMAX', ival=ldim)
    else
        zi(jlimac2-1+nedit+1) = zi(jlimac2-1+nedit)+ldim
    end if
!
! -- on "cree" un nouvel objet dans la collection lima :
    if (.not. lrapid) then
        call jecroc(jexnum(lima, nedit))
    else
! le cout du jeecra est trop important, il ne faut le
! faire qu'une seule fois quand la carte est terminee
! call jeecra(lima,'NUTIOC',ival= nedit)
    end if
!
!
! -- on remplit lima :
    jlima = jlima02-1+zi(jlimac2-1+nedit)
    do i = 1, dim
        if (mode2 .eq. 'NUM') then
!           --  mailles numerotees ( du maillage (>0) ou tardives(<0) )
            zi(jlima-1+i) = limanu(i)
        else if (mode2 .eq. 'NOM') then
!           --  mailles nommees du maillage
            nomail = limano(i)
            numero = char8_to_int(nomail)
            zi(jlima-1+i) = numero
        else
            ASSERT(.false.)
        end if
    end do
!
end subroutine
