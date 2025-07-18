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

subroutine coqucf(nomu)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cescar.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/fointe.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
!
    character(len=8) :: nomu
!
!                          AFFE_CARA_ELEM
!
!     POUR LES COQUES ET LES GRILLES IL PEUT EXISTER UNE CARTE DE
!     FONCTIONS DE (X,Y,Z) DONNANT :
!        EPAISSEUR
!        SECTION
!        EXCENTREMENT
!
!     IL FAUT EVALUER LES FONCTIONS AU CENTRE DE GRAVITE DE L'ELEMENT
!     ET REMPLIR LA CARTE DE REELS
! ----------------------------------------------------------------------
!
!
    integer(kind=8) :: ifm, niv, iret, ibid, ii, jj, kk
    integer(kind=8) :: jcesdf, jcesdo, nbmail, adrm, iad
    integer(kind=8) :: nbcmpf, nbcmpo, icompo, inoeu, nbno, nunoe, igeom
    integer(kind=8) :: jceslf, jceslo
    integer(kind=8) :: jconne, iadr1, iadr2, jtabco
    real(kind=8) :: valr(3), fresu
    character(len=8) :: nomma, nmcmpf, nomval(3), nomfct
    character(len=19) :: cartco, cartcf, celsco, celscf, connex
    character(len=24) :: k24bid
    aster_logical :: lcoor
    character(len=8), pointer :: cesvf(:) => null()
    real(kind=8), pointer :: cesvo(:) => null()
    character(len=8), pointer :: cescf(:) => null()
    character(len=8), pointer :: cesco(:) => null()
!
    data nomval/'X', 'Y', 'Z'/
! ----------------------------------------------------------------------
    call jemarq()
!
! --- CARTE POUR LES FONCTIONS
    cartcf = nomu//'.CARCOQUF'
    call exisd('CARTE', cartcf, iret)
!     SI LA CARTE DE FONCTIONS N'EXISTE PAS : RIEN A FAIRE
    if (iret .eq. 0) goto 999
!
! --- CARTE POUR LES VALEURS REELLES
    cartco = nomu//'.CARCOQUE'
    call exisd('CARTE', cartco, iret)
!     SI LA CARTE DE REELS N'EXISTE PAS : BIZARRE
    ASSERT(iret .ne. 0)
!     POUR TRACE DES INFORMATIONS DANS LE FICHIER DE MESSAGES
    call infniv(ifm, niv)
!
    celsco = '&&COQUCF.CHSCO'
    celscf = '&&COQUCF.CHSCF'
! --- TRANSFORMATION DE LA CARTE DE REELS EN ELEM_S
    call carces(cartco, 'ELEM', ' ', 'V', celsco, &
                'A', ibid)
! --- TRANSFORMATION DE LA CARTE DE FONCTIONS EN ELEM_S
    call carces(cartcf, 'ELEM', ' ', 'V', celscf, &
                'A', ibid)
    if (niv .ge. 2) then
        write (ifm, '(A)') ' '
        write (ifm, '(A)') 'TRANSFORMATION DES CARTES EN CHELEM_S'
        write (ifm, '(A)') '  <'//cartco//'>  <'//cartcf//'>'
    end if
!
! --- VERIFICATION QUE LES NOMS DANS LES DEUX ELEM_S EXISTENT
!        RECUPERATION DES NOMS DES COMPOSANTES DANS CELSCF
!        VERIFICATION DE L'EXISTANCE DANS CELSCO
    call jeveuo(celscf//'.CESD', 'L', jcesdf)
    call jeveuo(celsco//'.CESD', 'L', jcesdo)
    nbcmpf = zi(jcesdf+1)
    nbcmpo = zi(jcesdo+1)
!     ADRESSE DES NOMS DES COMPOSANTES
    call jeveuo(celscf//'.CESC', 'L', vk8=cescf)
    call jeveuo(celsco//'.CESC', 'L', vk8=cesco)
85  format('(', i3, 'A9)')
    if (niv .ge. 2) then
        write (ifm, '(A,I3)') 'COMPOSANTES DE <'//cartco//'> ', nbcmpo
        write (k24bid, 85) nbcmpo
        write (ifm, k24bid) (cesco(jj+1), jj=0, nbcmpo-1)
        write (ifm, '(A,I3)') 'COMPOSANTES DE <'//cartcf//'> ', nbcmpf
        write (k24bid, 85) nbcmpf
        write (ifm, k24bid) (cescf(jj+1), jj=0, nbcmpf-1)
    end if
!
    iret = 0
    do ii = 1, nbcmpf
        nmcmpf = cescf(ii)
        jj = indik8(cesco, nmcmpf, 1, nbcmpo)
        if (jj .ne. 0) then
            if (niv .ge. 2) then
                write (ifm, '(A)') 'COMPOSANTE <'//nmcmpf//'>'
                write (ifm, '(A,I5)') '   INDEX DANS CARCOQUF ', ii
                write (ifm, '(A,I5)') '   INDEX DANS CARCOQUE ', jj
            end if
        else
            iret = iret+1
        end if
    end do
    ASSERT(iret .eq. 0)
!
! --- INFORMATIONS SUR LE MAILLAGE
    call dismoi('NOM_MAILLA', cartco, 'CARTE', repk=nomma)
    call dismoi('NB_MA_MAILLA', nomma, 'MAILLAGE', repi=nbmail)
!
    k24bid = nomma//'.COORDO    .VALE'
    call jeveuo(k24bid, 'L', igeom)
    connex = nomma//'.CONNEX'
    call jeveuo(connex, 'L', jconne)
    call jeveuo(jexatr(connex, 'LONCUM'), 'L', jtabco)
!
    call jeveuo(celscf//'.CESL', 'L', jceslf)
    call jeveuo(celsco//'.CESL', 'L', jceslo)
    call jeveuo(celscf//'.CESV', 'L', vk8=cesvf)
    call jeveuo(celsco//'.CESV', 'E', vr=cesvo)
! --- TRAITEMENT
!        BOUCLE SUR TOUTES LES MAILLES
!           BOUCLE SUR LES COMPOSANTES AVEC FONCTIONS DE CELSCF
!              SI LA FONCTION EXISTE
!                 CALCUL DE LA POSITION DU CDG DE LA MAILLE
!                 CALCUL DE LA FONCTION
!                 AFFECTE LE RESULTAT A LA COMPOSANTE DE CELSCO
    if (niv .ge. 2) then
        write (ifm, '(A)') 'VALEURS DES FONCTIONS'
        write (ifm, 90)
    end if
    do ii = 1, nbmail
        lcoor = .false.
        do jj = 1, nbcmpf
            nmcmpf = cescf(jj)
            call cesexi('C', jcesdf, jceslf, ii, 1, &
                        1, jj, iad)
            if (iad .gt. 0) then
                nomfct = cesvf(iad)
                if (nomfct(1:2) .ne. '&&') then
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
                    end if
                    call fointe('F', nomfct, 3, nomval, valr, &
                                fresu, iret)
                    if (niv .ge. 2) then
                        write (ifm, 91) ii, jj, (valr(icompo), icompo=1, 3), &
                            nomfct, fresu
                    end if
!                 ON RECHERCHE DANS CELSCO LA COMPOSANTE CORRESPONDANTE
                    kk = indik8(cesco, nmcmpf, 1, nbcmpo)
                    call cesexi('S', jcesdo, jceslo, ii, 1, &
                                1, kk, iad)
                    cesvo(iad) = fresu
                end if
            end if
        end do
    end do
!
!     DESTRUCTION DE LA CARTE DES REELS, DES FONCTIONS
    call detrsd('CARTE', cartco)
    call detrsd('CARTE', cartcf)
!     CONSTRUCTION DE LA CARTE DES REELS A PARTIR DE CELSCO
    call cescar(celsco, cartco, 'G')
!     DESTRUCTION DES CHELEM_S
    call detrsd('CHAM_ELEM_S', celsco)
    call detrsd('CHAM_ELEM_S', celscf)
!
90  format(' MAILLE_NB  : [ [  CENTRE DE GRAVITE ],', &
           '  FONCTION  ,  VALEUR       ]')
91  format("'", i7, "_", i1, "' : [ [", 3(e18.10, ","), "], '", &
           a, "' ,", e18.10, "],")
!
999 continue
    call jedema()
end subroutine
