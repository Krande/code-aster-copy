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
subroutine char_eval_fonc(load, mesh, geomDime, param)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cescar.h"
#include "asterfort/cesexi.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/nocart.h"
#include "asterfort/rccome.h"
!

    character(len=8), intent(in) :: load, mesh
    integer(kind=8), intent(in) :: geomDime
    character(len=5), intent(in) :: param
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load 'FORCE_COQUE_FO'
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
! In  mesh             : mesh
! In  geomDime         : dimension of space
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ibid, jcesdf, jcesdr, nbcmpf, nbcmpr, nbmail, igeom
    integer(kind=8) :: jconne, jtabco, jceslf, jceslr, ii, iadr1, iadr2, nbno, adrm
    integer(kind=8) :: icompo, inoeu, nunoe, kk, iret, iad, jj, i, jvalv, n1, nbpara
    integer(kind=8) :: jcesdcoq, nbcmpcoq, jceslcoq, jcesdmat, jceslmat, icodn
    integer(kind=8) :: jvalk, jvalr, nbcste, icste
    real(kind=8) :: valr(5), fresu
    character(len=8) :: nomval(5), nomfct, nmcmpf, carele, chmat
    character(len=19) :: carte, cartefonc, celsreel, celsfonc, connex, cartco
    character(len=19) :: celscoq, cartmat, celsmat
    character(len=24) :: k24bid, rcvalk, rcvalr
    character(len=11) :: k11
    character(len=8), pointer :: vncmp(:) => null()
    character(len=8), pointer :: cesvf(:) => null()
    real(kind=8), pointer :: cesvcoq(:) => null()
    character(len=8), pointer :: cesvmat(:) => null()
    real(kind=8), pointer :: cesvr(:) => null()
    character(len=8), pointer :: cescf(:) => null()
    character(len=8), pointer :: cescr(:) => null()
    character(len=8), pointer :: cesccoq(:) => null()
    aster_logical :: lcoor, l_cara, l_mater

    data nomval/'X', 'Y', 'Z', 'EP', 'RHO'/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

    ASSERT(geomDime .eq. 3)
!
!   copie de la carte de fonction et destruction
    carte = load(1:8)//'.CHME.'//param(1:5)
    call exisd('CARTE', carte, iret)
    ASSERT(iret .ne. 0)
    cartefonc = '&&CHAREVALFONC.CART'
    call copisd('CHAMP', 'V', carte, cartefonc)
    call detrsd('CARTE', carte)

!   création de la carte à valeurs réelles et initialisation à 0
    call alcart('V', carte, mesh, 'FORC_R')
!
    call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
!
    vncmp(1) = 'FX'
    vncmp(2) = 'FY'
    vncmp(3) = 'FZ'
    vncmp(4) = 'MX'
    vncmp(5) = 'MY'
    vncmp(6) = 'MZ'
    vncmp(7) = 'REP'
    vncmp(8) = 'PLAN'
    vncmp(9) = 'MGX'
    vncmp(10) = 'MGY'
    vncmp(11) = 'MGZ'

    do i = 1, 11
        zr(jvalv-1+i) = 0.d0
    end do
    call nocart(carte, 1, 11)

!   transformation des cartes en cham_elem_s
    celsreel = '&&CHAREVALFONC.CHSR'
    celsfonc = '&&CHAREVALFONC.CHSF'
    call carces(carte, 'ELEM', ' ', 'V', celsreel, &
                'A', ibid)
    call carces(cartefonc, 'ELEM', ' ', 'V', celsfonc, &
                'A', ibid)
    call detrsd('CARTE', carte)
    call detrsd('CARTE', cartefonc)
    call jeveuo(celsfonc//'.CESD', 'L', jcesdf)
    call jeveuo(celsreel//'.CESD', 'L', jcesdr)
    nbcmpf = zi(jcesdf+1)
    nbcmpr = zi(jcesdr+1)
!     ADRESSE DES NOMS DES COMPOSANTES
    call jeveuo(celsfonc//'.CESC', 'L', vk8=cescf)
    call jeveuo(celsreel//'.CESC', 'L', vk8=cescr)
!
! --- information sur le maillage
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbmail)
!
    k24bid = mesh//'.COORDO    .VALE'
    call jeveuo(k24bid, 'L', igeom)
    connex = mesh//'.CONNEX'
    call jeveuo(connex, 'L', jconne)
    call jeveuo(jexatr(connex, 'LONCUM'), 'L', jtabco)
!
    call jeveuo(celsfonc//'.CESL', 'L', jceslf)
    call jeveuo(celsreel//'.CESL', 'L', jceslr)
    call jeveuo(celsfonc//'.CESV', 'L', vk8=cesvf)
    call jeveuo(celsreel//'.CESV', 'E', vr=cesvr)

!   récupération de l'épaisseur si CARA_ELEM présent
    l_cara = ASTER_FALSE
    call getvid(' ', 'CARA_ELEM', scal=carele, nbret=n1)
    if (n1 .ne. 0) then
        l_cara = ASTER_TRUE
        cartco = carele//'.CARCOQUE'
        call exisd('CARTE', cartco, iret)
        ASSERT(iret .ne. 0)
        celscoq = '&&CHAREVALFONC.COQU'
        call carces(cartco, 'ELEM', ' ', 'V', celscoq, &
                    'A', ibid)
        call jeveuo(celscoq//'.CESD', 'L', jcesdcoq)
        nbcmpcoq = zi(jcesdcoq+1)
!       ADRESSE DES NOMS DES COMPOSANTES
        call jeveuo(celscoq//'.CESC', 'L', vk8=cesccoq)
        call jeveuo(celscoq//'.CESL', 'L', jceslcoq)
        call jeveuo(celscoq//'.CESV', 'L', vr=cesvcoq)
    end if
!   récupération de RHO si CHAM_MATER présent
    l_mater = ASTER_FALSE
    call getvid(' ', 'CHAM_MATER', scal=chmat, nbret=n1)
    if (n1 .ne. 0) then
        l_mater = ASTER_TRUE
        cartmat = chmat//'.CHAMP_MAT '
        call exisd('CARTE', cartmat, iret)
        ASSERT(iret .ne. 0)
!        call etenca(carte, ligrmo, iret)
        celsmat = '&&CHAREVALFONC.MATE'
        call carces(cartmat, 'ELEM', ' ', 'V', celsmat, &
                    'A', ibid)
        call jeveuo(celsmat//'.CESD', 'L', jcesdmat)
!        nbcmpcoq = zi(jcesdcoq+1)
!       ADRESSE DES NOMS DES COMPOSANTES
!        call jeveuo(celscoq//'.CESC', 'L', vk8=cesccoq)
        call jeveuo(celsmat//'.CESL', 'L', jceslmat)
        call jeveuo(celsmat//'.CESV', 'L', vk8=cesvmat)
    end if
! --- TRAITEMENT
!        BOUCLE SUR TOUTES LES MAILLES
!           BOUCLE SUR LES COMPOSANTES AVEC FONCTIONS DE CELSFONC
!              SI LA FONCTION EXISTE
!                 CALCUL DE LA POSITION DU CDG DE LA MAILLE
!                 CALCUL DE LA FONCTION
!                 AFFECTE LE RESULTAT A LA COMPOSANTE DE CELSREEL
    do ii = 1, nbmail
        lcoor = ASTER_FALSE
        nbpara = 3
        do jj = 1, nbcmpf
            nmcmpf = cescf(jj)
            call cesexi('C', jcesdf, jceslf, ii, 1, &
                        1, jj, iad)
            if (iad .gt. 0) then
                nomfct = cesvf(iad)
                if (nmcmpf .eq. 'REP') then
                    if (nomfct .eq. 'GLOBAL') then
                        fresu = 0.d0
                    elseif (nomfct .eq. 'LOCAL') then
                        fresu = 1.d0
                    elseif (nomfct .eq. 'VENT') then
                        fresu = 2.d0
                    elseif (nomfct .eq. 'LOCAL_PR') then
                        fresu = 3.d0
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                else if (nomfct(1:2) .ne. '&&') then
                    if (.not. lcoor) then
                        lcoor = .true.
                        iadr1 = zi(jtabco-1+ii)
                        iadr2 = zi(jtabco-1+ii+1)
                        nbno = iadr2-iadr1
                        adrm = jconne-1+iadr1
!                    CENTRE DE GRAVITE DE LA MAILLE
                        do icompo = 1, 3
                            valr(icompo) = 0.0d0
                            do inoeu = 1, nbno
                                nunoe = zi(adrm-1+inoeu)
                                valr(icompo) = valr(icompo)+zr(igeom+3*(nunoe-1)+icompo-1)
                            end do
                            valr(icompo) = valr(icompo)/nbno
                        end do
!                       epaisseur de la maille
                        if (l_cara) then
                            kk = indik8(cesccoq, 'EP', 1, nbcmpcoq)
                            call cesexi('C', jcesdcoq, jceslcoq, ii, 1, &
                                        1, kk, iad)
                            if (iad .gt. 0) then
                                nbpara = nbpara+1
                                valr(nbpara) = cesvcoq(iad)
                            end if
                        end if
!                       rho
                        if (l_mater) then
                            call cesexi('C', jcesdmat, jceslmat, ii, 1, &
                                        1, 1, iad)
                            if (iad .gt. 0) then
                                call rccome(cesvmat(iad), 'ELAS', icodn, &
                                            k11_ind_nomrc=k11)
                                if (icodn .eq. 0) then
                                    rcvalk = cesvmat(iad)//k11//'.VALK'
                                    rcvalr = cesvmat(iad)//k11//'.VALR'
                                    call jeveuo(rcvalk, 'L', jvalk)
                                    call jeveuo(rcvalr, 'L', jvalr)
                                    call jelira(rcvalr, 'LONMAX', nbcste)
                                    do icste = 1, nbcste
                                        if (zk16(jvalk+icste-1) .eq. 'RHO') then
                                            nbpara = nbpara+1
                                            valr(nbpara) = zr(jvalr+icste-1)
                                        end if
                                    end do
                                end if
                            end if
                        end if
                    end if
                    call fointe('F', nomfct, nbpara, nomval, valr, &
                                fresu, iret)
                else
                    cycle
                end if
!               ON RECHERCHE DANS CELSREEL LA COMPOSANTE CORRESPONDANTE
                kk = indik8(cescr, nmcmpf, 1, nbcmpr)
                call cesexi('S', jcesdr, jceslr, ii, 1, &
                            1, kk, iad)
                cesvr(iad) = fresu
            end if
        end do
    end do
!     CONSTRUCTION DE LA CARTE DES REELS A PARTIR DE celsreel
    call cescar(celsreel, carte, 'G')
!     DESTRUCTION DES CHELEM_S
    call detrsd('CHAM_ELEM_S', celsreel)
    call detrsd('CHAM_ELEM_S', celsfonc)
    if (l_cara) call detrsd('CHAM_ELEM_S', celscoq)
    if (l_mater) call detrsd('CHAM_ELEM_S', celsmat)

    call jedema()
end subroutine
