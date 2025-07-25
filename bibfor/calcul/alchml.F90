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

subroutine alchml(ligrel_, option_, nompar_, base_, cel_, &
                  iret, dcel_)
    implicit none
! person_in_charge: jacques.pellet at edf.fr

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/digdel.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/modat2.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/nopar2.h"
#include "asterfort/scalai.h"
#include "asterfort/typele.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/nucalc.h"

    character(len=*) :: ligrel_, dcel_, cel_, base_, option_, nompar_
    character(len=19) :: ligrel, cel, dcel
    character(len=16) :: option
    character(len=8) :: nompar
    character(len=1) :: base
    integer(kind=8) :: iret
!----------------------------------------------------------------------
!  But : Creer un cham_elem "vierge"

!  Arguments :
!  -----------
!    ligrel_ in/jxin  k19 : sd ligrel sur lequel on alloue le cham_elem
!    option_ in       k16 : nom de l'option servant a decrire le cham_elem
!    nompar_ in       k8  : nom du parametre (in ou out) de l'option
!                           servant a decrire le cham_elem
!                           - si nompar_ n'est pas fourni (' '), on se sert
!                             du parametre "out" s'il est unique.
!                           sinon : erreur <F>
!    base_   in       k1  : 'G','V'
!    cel_    in/jxout k19 : sd cham_elem a creer
!    iret    out      i   : code retour :
!                          0 -> le champ a ete cree
!                          1 -> le champ n'a pas ete cree car
!                               aucun type_elem du ligrel ne connait
!                               le parametre de l'option

!  Arguments supplementaires pour allouer un cham_elem "etendu" :
!  --------------------------------------------------------------
!    dcel_   in/jxin  k19 :
!      sd cham_elem_s permettant de creer un cham_elem "etendu".
!      la grandeur associee a dcel_ doit etre "DCEL_I"
!      et les cmps doivent etre "NPG_DYN" et "NCMP_DYN" (dans cet ordre)

!    si dcel_ = ' '
!       le cham_elem n'est pas etendu.

!    si dcel_ /= ' '   :  le cham_elem est etendu :
!       Les mailles "tardives" sont alors interdites dans ligrel_

!----------------------------------------------------------------------
    character(len=8) :: scal, nomgd
    character(len=16) :: nomte1, nomte, ma, ma2, kbid
    character(len=24) :: nomolo
    character(len=24) :: valk(3)
    integer(kind=8) :: ngrel, igrel, te, te1, mode, long, jceld, ncmpv, debgrl
    integer(kind=8) :: gd, jcelk, iopt, iprem, nel, iel, lgcata, nbspt
    integer(kind=8) :: ncdyn, lgchel, numc
    integer(kind=8) :: ibid, modmx, iamolo, itych, itych1
    integer(kind=8) :: illiel, jdcesd, jdcesl
    integer(kind=8) :: ima, ncmpv2, kk, ityplo
    aster_logical :: letendu
    integer(kind=8), pointer :: liel(:) => null()
    integer(kind=8), pointer :: cesv(:) => null()
    character(len=8), pointer :: cesk(:) => null()
    character(len=8), pointer :: cesc(:) => null()
!   ------------------------------------------------------------------

    call jemarq()
    cel = cel_
    dcel = dcel_
    ligrel = ligrel_
    option = option_
    nompar = nompar_
    if (nompar .eq. ' ') call nopar2(option, ' ', 'OUT', nompar)
    base = base_

    ngrel = nbgrel(ligrel)
    call jenonu(jexnom('&CATA.OP.NOMOPT', option), iopt)
    call jeveuo(ligrel//'.LIEL', 'L', vi=liel)
    call jeveuo(jexatr(ligrel//'.LIEL', 'LONCUM'), 'L', illiel)
    call dismoi('NOM_MAILLA', ligrel, 'LIGREL', repk=ma)

!   1- le cham_elem doit-il etre cree ?
!   ----------------------------------
    modmx = 0
    do igrel = 1, ngrel
        te = typele(ligrel, igrel)
        mode = modat2(iopt, te, nompar)
        modmx = max(modmx, mode)
    end do
    if (modmx .eq. 0) then
        iret = 1
        goto 60
    else
        iret = 0
    end if

!   2- quelle est la grandeur associee au cham_elem ?
!   -----------------------------------------------------------
    call jeveuo(jexnum('&CATA.TE.MODELOC', modmx), 'L', iamolo)
    gd = zi(iamolo-1+2)
    call jenuno(jexnum('&CATA.GD.NOMGD', gd), nomgd)
    scal = scalai(gd)

!   3- doit-on creer un champ etendu ? (letendu.eq..true.)
!   ------------------------------------------------------
    if (dcel_ .eq. ' ') then
!       -- champ non-etendu :
        letendu = .false.
    else
!       -- champ etendu : dcel est fourni par l'appelant
        letendu = .true.
        dcel = dcel_
    end if

!   3.1 si champ etendu : on recupere quelques adresses :
!   --------------------------------------------------------
    if (letendu) then
        call jeveuo(dcel//'.CESK', 'L', vk8=cesk)
        call jeveuo(dcel//'.CESC', 'L', vk8=cesc)
        call jeveuo(dcel//'.CESD', 'L', jdcesd)
        call jeveuo(dcel//'.CESL', 'L', jdcesl)
        call jeveuo(dcel//'.CESV', 'L', vi=cesv)

!       -- quelques verifications :
        ma2 = cesk(1)
        if (ma2 .ne. ma) then
            valk(1) = ma2
            valk(2) = ma
            call utmess('F', 'CALCUL_39', nk=2, valk=valk)
        end if
        ASSERT(zi(jdcesd-1+2) .eq. 2)
        ASSERT(zi(jdcesd-1+3) .eq. 1)
        ASSERT(zi(jdcesd-1+4) .eq. 1)

        kbid = cesk(2)
        ASSERT(kbid .eq. 'DCEL_I')

        kbid = cesc(1)
        ASSERT(kbid .eq. 'NPG_DYN')
        kbid = cesc(2)
        ASSERT(kbid .eq. 'NCMP_DYN')
    end if

!   4- objet .CELD :
!   -------------------
    long = 4+ngrel+4*ngrel
    do igrel = 1, ngrel
        nel = nbelem(ligrel, igrel)
        te = typele(ligrel, igrel)
        mode = modat2(iopt, te, nompar)
        if (mode .gt. 0) then
            long = long+4*nel
        end if
    end do

    call wkvect(cel//'.CELD', base//' V I', long, jceld)
    call jeecra(cel//'.CELD', 'DOCU', cval='CHML')

    zi(jceld-1+1) = gd
    zi(jceld-1+2) = ngrel
    zi(jceld-1+3) = 1
    zi(jceld-1+4) = 0

!   ncmpv: longueur de .CELV  (+1)
    ncmpv = 1

!   debgrl: debut de description du grel dans .CELD
    debgrl = 4+ngrel

    iprem = 0
    do igrel = 1, ngrel
        nel = nbelem(ligrel, igrel)
        te = typele(ligrel, igrel)
        mode = modat2(iopt, te, nompar)
        zi(jceld-1+4+igrel) = debgrl
        zi(jceld-1+debgrl+1) = nel
        zi(jceld-1+debgrl+2) = mode
        ncmpv2 = ncmpv

!       -- faut-il verifier que numc != -1 ?
        if (igrel .eq. 1) then
            numc = nucalc(iopt, te, 1)
        else
            numc = nucalc(iopt, te, 0)
        end if
        if (numc .eq. -1 .or. numc .eq. -2) then
            call jenuno(jexnum('&CATA.TE.NOMTE', te), nomte)
            valk(1) = nomte
            valk(2) = option
            if (numc .eq. -1) then
                call utmess('F', 'CALCUL_37', nk=2, valk=valk)
            else
                call utmess('A', 'CALCUL_41', nk=2, valk=valk)
            end if
        end if

        if (mode .gt. 0) then
            iprem = iprem+1
            call jeveuo(jexnum('&CATA.TE.MODELOC', mode), 'L', iamolo)
            ityplo = zi(iamolo-1+1)
            if (ityplo .gt. 3) then
                call jenuno(jexnum('&CATA.TE.NOMMOLOC', mode), nomolo)
                call utmess('F', 'CALCUL_33', sk=nomolo)
            end if

            itych = zi(iamolo-1+1)

            if (iprem .eq. 1) then
                itych1 = itych
                te1 = te
            else
                if (itych*itych1 .lt. 0) goto 50
            end if

            lgcata = digdel(mode)
            zi(jceld-1+debgrl+3) = lgcata

            do iel = 1, nel

!               4.1 calcul de nbspt et ncdyn pour chaque element :
!               --------------------------------------------------

!               -- cas d'un cham_elem etendu :
                if (letendu) then
                    ima = liel(zi(illiel+igrel-1)+iel-1)

!                   -- si la maille appartient au maillage, on se sert
!                      de dcel, sinon on arrete le code:
                    if (ima .gt. 0) then

                        call cesexi('C', jdcesd, jdcesl, ima, 1, &
                                    1, 1, kk)
                        if (kk .gt. 0) then
                            nbspt = max(cesv(kk), 1)
                        else
                            nbspt = 1
                        end if

                        ncdyn = 0
                        if (nomgd .eq. 'VARI_R') then
                            call cesexi('C', jdcesd, jdcesl, ima, 1, &
                                        1, 2, kk)
                            if (kk .gt. 0) ncdyn = cesv(kk)
                        end if

                    else
!                       -- cas des mailles tardives :
                        nbspt = 1
                        ncdyn = 0
                        if (nomgd .eq. 'VARI_R') ncdyn = 1
                    end if

                else
!                   -- cas d'un cham_elem non-etendu :
                    nbspt = 1
                    ncdyn = 0
                    if (nomgd .eq. 'VARI_R') ncdyn = 1
                end if

!               4.2 affectation des valeurs dans CELD :
!               ---------------------------------------
                zi(jceld-1+debgrl+4+(iel-1)*4+1) = nbspt
                zi(jceld-1+debgrl+4+(iel-1)*4+2) = ncdyn
                zi(jceld-1+3) = max(zi(jceld-1+3), nbspt)
                zi(jceld-1+4) = max(zi(jceld-1+4), ncdyn)

                lgchel = lgcata*nbspt*max(1, ncdyn)
                zi(jceld-1+debgrl+4+(iel-1)*4+3) = lgchel
                zi(jceld-1+debgrl+4+(iel-1)*4+4) = ncmpv
                ncmpv = ncmpv+lgchel
            end do
        end if

        zi(jceld-1+debgrl+4) = ncmpv-ncmpv2
        if (mode .gt. 0) then
            debgrl = debgrl+4+4*nel
        else
            debgrl = debgrl+4
        end if
    end do

!   5- objet .CELV:
!   ----------------
    call wkvect(cel//'.CELV', base//' V '//scal(1:4), ncmpv-1, ibid)

!   6- objet .CELK:
!   ----------------
    call wkvect(cel//'.CELK', base//' V K24', 7, jcelk)
    zk24(jcelk-1+1) = ligrel
    zk24(jcelk-1+2) = option
    zk24(jcelk-1+6) = nompar
    if (itych1 .eq. 1) then
        zk24(jcelk-1+3) = 'ELEM'
    else if (itych1 .eq. 2) then
        zk24(jcelk-1+3) = 'ELNO'
    else if (itych1 .eq. 3) then
        zk24(jcelk-1+3) = 'ELGA'
    else
        ASSERT(.false.)
    end if
    zk24(jcelk-1+4) = ' '
    zk24(jcelk-1+5) = ' '
    zk24(jcelk-1+7) = 'MPI_COMPLET'

    goto 60

!   7- section erreur:
!   ------------------------
50  continue

!   ce cas de figure ne devrait plus exister apres la verif dans
!   caver1.F90  (coherence des type_elem avec l'option):
    call jenuno(jexnum('&CATA.TE.NOMTE', te), nomte)
    call jenuno(jexnum('&CATA.TE.NOMTE', te1), nomte1)
    valk(1) = option
    valk(2) = nomte1
    valk(3) = nomte
    call utmess('F', 'CALCUL_36', nk=3, valk=valk)

!   8- fin normale:
!   ----------------
60  continue

    call jedema()

end subroutine
