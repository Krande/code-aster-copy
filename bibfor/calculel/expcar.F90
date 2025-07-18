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

subroutine expcar(carte)
    implicit none
!
! EXTENSION COMPLETE DE LA CARTE( EN VUE D'1 COMPRESSION ULTERIEURE )
! ( LORSQUE LES CMPS D'1 GRANDEUR N'ONT PAS ETE DONNEES SIMULTANEMENT)
!
!-----------------------------------------------------------------------
!
!     ARGUMENTS:
!     ----------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mecumu.h"
#include "asterfort/melima.h"
#include "asterfort/nbec.h"
#include "asterfort/scalai.h"
!
    character(len=19) :: carte
! ----------------------------------------------------------------------
!     ENTREES:
!       CARTE : NOM D'1 CARTE A ETENDRE PROVISOIREMENT
!               ( AVANT COMPRESSION)
!     SORTIES:
!      ON A CREE QUELQUES OBJETS SUR LA VOLATILE ...
! ----------------------------------------------------------------------
    aster_logical :: dejavu
    character(len=8) :: scal, noma
    character(len=24) :: noli
!
!
!     -- RECUPERATION DES OBJETS JEVEUX DE LA CARTE:
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i1, i2, i3, i4, i5, iadesc, iadgp
    integer(kind=8) :: ialima, iavale, iavalp, icode
    integer(kind=8) :: iedit, ient, igd, ima, iret, j
    integer(kind=8) :: nbedit, nbgdmx, nbma, nbmato, ncmpmx, nec, num1
    integer(kind=8) :: num2, numat
    integer(kind=8), pointer :: numt(:) => null()
    character(len=24), pointer :: vnoli(:) => null()
    character(len=8), pointer :: vnoma(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call jeveuo(carte//'.DESC', 'L', iadesc)
    call jeveuo(carte//'.LIMA', 'L', ialima)
    call jeveuo(carte//'.VALE', 'L', iavale)
    call jeveuo(carte//'.NOMA', 'L', vk8=vnoma)
!
    noma = vnoma(1)
!
    igd = zi(iadesc-1+1)
    nec = nbec(igd)
!     -- SCAL = I,R,C,K8,...
    scal = scalai(igd)
!
!     -- NCMPMX : NOMBRE MAXIMAL DE CMP POUR LA GRANDEUR.
!     ----------------------------------------------------
    call jelira(jexnum('&CATA.GD.NOMCMP', igd), 'LONMAX', ncmpmx)
!
!     -- NBEDIT : NOMBRE DE VALEURS EDITEES DANS LA CARTE:
!     ----------------------------------------------------
    nbedit = zi(iadesc-1+3)
!
!     -- CREATION DE L'OBJET .NUMT:
!     -----------------------------
!     NUMT(IEDIT) := V(I)  (DIM=3)
!     V(1) : NUMERO GLOBAL DANS .VALP ET .DGP
!            DE LA 1ERE MAILLE AFFECTEE PAR LA GRANDEUR IEDIT.
!     V(2) : NUMERO DE LA DER MAILLE AFFECTEE.
!     V(3) : 1 --> CETTE GRANDEUR EST ASSOCIEE A 1 PAQUET DE MAILLES
!                 NON TRAITE POUR LES IEDIT PRECEDENTS.
!            0 --> SINON
!     SI CODE(IEDIT)= 1:  TTES LES MAILLES DU MAILLAGE'
!     SI CODE(IEDIT)=-1:  TTES LES MAILLES SUPPL. DU LIGREL NOLI(IEDIT)
!     SI CODE(IEDIT)= 2:  LES MAILLES D'1 GROUPE NOMME.
!     SI CODE(IEDIT)= 3:  LES MAILLES DU MAILLAGE D'1 LISTE TARDIVE.
!     SI CODE(IEDIT)=-3:  LES MAILLES SUPPL. D'1 LISTE TARDIVE.
!
    call jeexin(carte//'.NOLI', iret)
    ASSERT(iret .ne. 0)
    call jeveuo(carte//'.NOLI', 'L', vk24=vnoli)
    call jecreo(carte//'.NUMT', 'V V I')
    call jeecra(carte//'.NUMT', 'LONMAX', 3*nbedit)
    call jeveuo(carte//'.NUMT', 'E', vi=numt)
    nbmato = 0
    do iedit = 1, nbedit
        icode = zi(iadesc-1+3+2*(iedit-1)+1)
        noli = vnoli(iedit)
        dejavu = .false.
        do j = iedit-1, 1, -1
            if (noli .eq. vnoli(j)) then
                dejavu = .true.
                goto 3
            end if
        end do
3       continue
        if (dejavu) then
            numt(3*(iedit-1)+1) = numt(3*(j-1)+1)
            numt(3*(iedit-1)+2) = numt(3*(j-1)+2)
            numt(3*(iedit-1)+3) = 0
        else
            numt(3*(iedit-1)+1) = nbmato+1
            if (noli(1:8) .eq. '        ') then
!              MAILLES DU MAILLAGE:
                call jelira(noma//'.TYPMAIL', 'LONMAX', nbma)
                nbmato = nbmato+nbma
            else
!              MAILLES SUPPLEMENTAIRES D'1 LIGREL:
!              TEST : ON DOIT AVOIR ICODE=3 POUR DES MAILLES TARDIVES
                ASSERT(icode .eq. -3)
                call dismoi('NB_MA_SUP', noli, 'LIGREL', repi=nbma)
                nbmato = nbmato+nbma
            end if
            numt(3*(iedit-1)+2) = nbmato
            numt(3*(iedit-1)+3) = 1
        end if
    end do
!
!     -- ALOCATION DES OBJETS DE TRAVAIL : .VALP ET .DGP
!
    call jecreo(carte//'.VALP', 'V V '//scal(1:4))
    call jecreo(carte//'.DGP ', 'V V I')
    call jeecra(carte//'.VALP', 'LONMAX', ncmpmx*nbmato)
    call jeecra(carte//'.DGP ', 'LONMAX', nec*nbmato)
    call jeveuo(carte//'.VALP', 'E', iavalp)
    call jeveuo(carte//'.DGP ', 'E', iadgp)
!
!     --REMPLISSAGE DE .VALP ET .DGP:
!     -------------------------------
!
    nbgdmx = zi(iadesc-1+2)
    do iedit = 1, nbedit
        i1 = iavale+(iedit-1)*ncmpmx
        i3 = iadesc-1+3+2*nbgdmx+(iedit-1)*nec+1
        num1 = numt((iedit-1)*3+1)
        num2 = numt((iedit-1)*3+2)
        icode = zi(iadesc-1+3+2*(iedit-1)+1)
        ient = zi(iadesc-1+3+2*(iedit-1)+2)
        noli = vnoli(iedit)
        if (abs(icode) .eq. 1) then
!            -- GROUPE 'TOUT':
            nbma = num2-num1+1
        else
            call melima(carte, noma, icode, ient, i5, &
                        nbma)
!           --I5 : ADRESSE DANS ZI DE LA LISTE DES MAILLES A TRAITER.
!           --NBMA: NOMBRE DE MAILLES A TRAITER.
        end if
        do ima = 1, nbma
            if (abs(icode) .eq. 1) then
                numat = ima
            else
                numat = num1-1+abs(zi(i5-1+ima))
            end if
            i2 = iavalp+(numat-1)*ncmpmx
            i4 = iadgp-1+(numat-1)*nec+1
            call mecumu(scal, ncmpmx, i1, i2, nec, &
                        zi(i3), zi(i4))
        end do
    end do
!
!
    call jedema()
end subroutine
