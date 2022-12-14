! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine utmavo(mail, kdim, lima, nlima, base,&
                  nomz, nbmavo, mailvo)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cncinv.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    integer :: lima(*), nlima, nbmavo, mailvo(*)
    character(len=1) :: base
    character(len=2) :: kdim
    character(len=8) :: mail
    character(len=*) :: nomz
!
!     DETERMINE LES MAILLES VOISINES D'UNE LISTE DE MAILLES
!     MAILLE --> LISTE DES MAILLES VOISINES
!
!   ARGUMENT EN ENTREE
!   ------------------
!     MAIL   : NOM DE L'OJB REPRESENTANT LE MAILLAGE
!     KDIM   : '3D' RECHERCHE LES MAILLES 3D VOISINES
!              '2D' RECHERCHE LES MAILLES 2D VOISINES
!              '1D' RECHERCHE LES MAILLES 1D VOISINES
!              '  ' RECHERCHE TOUTES LES MAILLES VOISINES
!     LIMA   : LISTE DES NUMEROS DE MAILLES
!     NLIMA  : NOMBRE DE MAILLES
!     BASE   : BASE DE CREATION
!     NOMZ   : NOM DE L' OJB A CREER
!     MAILVO : SI ORIE_PEAU_3D ("GROUP_MA_VOLU"):
!                  = LISTE DES MAILLES VOLUMIQUES
!                    UTILES A LA REORIENTATION
!              SI ORIE_PEAU_2D ("GROUP_MA_SURF"):
!                  = LISTE DES MAILLES SURFACIQUES
!                    UTILES A LA REORIENTATION
!              SINON: MAILVO N'EST PAS UTILISE
!     NBMAVO : NB DE MAILLES DE MAILVO
!
!   ORGANISATION
!   ------------
!     TYPE : XC V I ACCES(NUMEROTE) LONG(VARIABLE)
!-----------------------------------------------------------------------
!
    integer :: ibid, nare, numa, nbno, nbmat, ino, nuno, p2
    integer :: i, j, k, jmail, nbman, adrvlc, acncin, ima, ii
    integer :: adra, iad, jtr1(1000), nutyma, iexinv
    character(len=8) :: type
    character(len=24) :: nom, ncninv
    integer, pointer :: trav2(:) => null()
    integer, pointer :: typmail(:) => null()
    integer, pointer :: connex(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    nom = nomz
!
! --- APPEL A LA CONNECTIVITE INVERSE
!
    ibid = 0
    ncninv = '&&UTMAVO.CONNEC_INV'
!     EST-CE QUE LA CONNECTIVITE INVERSE A DEJA ETE CALCULEE ?
    call jeexin(ncninv, iexinv)
    if (nbmavo .eq. 0) then
        if (iexinv .eq. 0) call cncinv(mail, [ibid], 0, 'V', ncninv)
    else
!        ON FORCE LE CALCUL DE LA CONNECTIVITE INVERSE
        call jedetr(ncninv)
        call cncinv(mail, mailvo, nbmavo, 'V', ncninv)
    endif
    call jeveuo(jexatr(ncninv, 'LONCUM'), 'L', adrvlc)
    call jeveuo(ncninv, 'L', acncin)
!
! --- APPEL A LA CONNECTIVITE
!
    call jeveuo(jexatr(mail//'.CONNEX', 'LONCUM'), 'L', p2)
    call jeveuo(mail//'.CONNEX', 'L', vi=connex)
!
    call jeveuo(mail//'.TYPMAIL', 'L', vi=typmail)
!
! --- DIMENSIONNEMENT DE LA SD
!
    AS_ALLOCATE(vi=trav2, size=nlima)
    nare = 0
    do i = 1, nlima
        numa = lima(i)
        nbno = zi(p2+numa+1-1) - zi(p2+numa-1)
        iad = zi(p2+numa-1)
        nbmat = 0
        do ino = 1, nbno
            nuno = connex(1+iad-1+ino-1)
            nbman = zi(adrvlc+nuno+1-1) - zi(adrvlc+nuno-1)
            adra = zi(adrvlc+nuno-1)
            do j = 1, nbman
                ii = zi(acncin+adra-1+j-1)
!              -- SI UN NOEUD EST ORPHELIN : II=0
!                 (PAS D'OBJET JEVEUX DE LONG=0)
                if (ii .eq. 0) then
                    ASSERT(nbman.eq.1)
                    goto 120
                endif
                if (nbmavo .eq. 0) then
                    ima=ii
                else
                    ima=mailvo(ii)
                endif
                if (ima .eq. numa) goto 120
                nutyma = typmail(ima)
                call jenuno(jexnum('&CATA.TM.NOMTM', nutyma), type)
                if (type(1:4) .eq. 'HEXA') then
                    if (kdim .eq. '2D' .or. kdim .eq. '1D') goto 120
                else if (type(1:4).eq.'PENT') then
                    if (kdim .eq. '2D' .or. kdim .eq. '1D') goto 120
                else if (type(1:4).eq.'PYRA') then
                    if (kdim .eq. '2D' .or. kdim .eq. '1D') goto 120
                else if (type(1:4).eq.'TETR') then
                    if (kdim .eq. '2D' .or. kdim .eq. '1D') goto 120
                else if (type(1:4).eq.'QUAD') then
                    if (kdim .eq. '3D' .or. kdim .eq. '1D') goto 120
                else if (type(1:4).eq.'TRIA') then
                    if (kdim .eq. '3D' .or. kdim .eq. '1D') goto 120
                else if (type(1:3).eq.'SEG') then
                    if (kdim .eq. '3D' .or. kdim .eq. '2D') goto 120
                else if (type(1:3).eq.'POI') then
                    if (kdim .ne. '  ') goto 120
                else
                    call utmess('F', 'PREPOST4_89', sk=type)
                endif
                do k = 1, nbmat
                    if (jtr1(k) .eq. ima) goto 120
                end do
                nbmat = nbmat + 1
                jtr1(nbmat) = ima
120             continue
            end do
        end do
        trav2(i) = nbmat
        nare = nare + max(nbmat,1)
    end do
!
! --- CREATION DE LA SD
!
    call jecrec(nom, base//' V I', 'NU', 'CONTIG', 'VARIABLE',&
                nlima)
    call jeecra(nom, 'LONT', nare)
!
! --- ON REMPLIT LA SD
!
    do i = 1, nlima
        numa = lima(i)
        nbno = zi(p2+numa+1-1) - zi(p2+numa-1)
        iad = zi(p2+numa-1)
        call jecroc(jexnum(nom, i))
        if (trav2(i) .eq. 0) then
            call jeecra(jexnum(nom, i), 'LONMAX', 1)
            call jeecra(jexnum(nom, i), 'LONUTI', 0)
            goto 200
        else
            call jeecra(jexnum(nom, i), 'LONMAX', trav2(i))
            call jeecra(jexnum(nom, i), 'LONUTI', trav2(i))
            call jeveuo(jexnum(nom, i), 'E', jmail)
        endif
!
        nbmat = 0
        do ino = 1, nbno
            nuno = connex(1+iad-1+ino-1)
            nbman = zi(adrvlc+nuno+1-1) - zi(adrvlc+nuno-1)
            adra = zi(adrvlc+nuno-1)
            do j = 1, nbman
                ii = zi(acncin+adra-1+j-1)
                if (ii .eq. 0) goto 220
!
                if (nbmavo .eq. 0) then
                    ima=ii
                else
                    ima=mailvo(ii)
                endif
                if (ima .eq. numa) goto 220
                nutyma = typmail(ima)
                call jenuno(jexnum('&CATA.TM.NOMTM', nutyma), type)
                if (type(1:4) .eq. 'HEXA') then
                    if (kdim .eq. '2D') goto 220
                    if (kdim .eq. '1D') goto 220
                else if (type(1:4).eq.'PENT') then
                    if (kdim .eq. '2D') goto 220
                    if (kdim .eq. '1D') goto 220
                else if (type(1:4).eq.'PYRA') then
                    if (kdim .eq. '2D') goto 220
                    if (kdim .eq. '1D') goto 220
                else if (type(1:4).eq.'TETR') then
                    if (kdim .eq. '2D') goto 220
                    if (kdim .eq. '1D') goto 220
                else if (type(1:4).eq.'QUAD') then
                    if (kdim .eq. '3D') goto 220
                    if (kdim .eq. '1D') goto 220
                else if (type(1:4).eq.'TRIA') then
                    if (kdim .eq. '3D') goto 220
                    if (kdim .eq. '1D') goto 220
                else if (type(1:3).eq.'SEG') then
                    if (kdim .eq. '3D') goto 220
                    if (kdim .eq. '2D') goto 220
                else if (type(1:3).eq.'POI') then
                    if (kdim .ne. '  ') goto 220
                else
                    call utmess('F', 'PREPOST4_89', sk=type)
                endif
                do k = 1, nbmat
                    if (zi(jmail-1+k) .eq. ima) goto 220
                end do
                nbmat = nbmat + 1
                zi(jmail-1+nbmat) = ima
220             continue
            end do
        end do
200     continue
    end do
!
    call jedetr('&&UTMAVO.TRAV1')
    AS_DEALLOCATE(vi=trav2)
    if (nbmavo .ne. 0) call jedetr(ncninv)
    call jedema()
!
end subroutine
