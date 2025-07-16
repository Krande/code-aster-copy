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
subroutine orvlse(noma, listCellNume, nbCell, norien, vect, &
                  noeud)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/indiis.h"
#include "asterfort/infniv.h"
#include "asterfort/iorim1.h"
#include "asterfort/ioriv3.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmavo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nbCell, noeud, norien
    integer(kind=8), pointer :: listCellNume(:)
    character(len=8) :: noma
    real(kind=8) :: vect(*)
!.======================================================================
!
!   ORVLSE  --  LE BUT EST QUE TOUTES LES MAILLES DE TYPE SEG
!               DE LA LISTE SOIENT
!               ORIENTEES SUIVANT LE VECTEUR DIRECTEUR.
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NOMA           IN    K8      NOM DU MAILLAGE
!    LISTMA         IN    I       LISTE DES MAILLES A REORIENTER
!    NBMAIL         IN    I       NB DE MAILLES DE LA LISTE
!    NORIEN        VAR            NOMBRE DE MAILLES REORIENTEES
!    VECT           IN    R       VECTEUR DIRECTEUR
!    NOEUD          IN    I       NOEUD D'ORIENTATION
!.========================= DEBUT DES DECLARATIONS ====================
! -----  VARIABLES LOCALES
    integer(kind=8) :: cellTypeNume, lori, jori, nori, kori, iliste
    integer(kind=8) :: iCell, numail, cellNume, norieg, lliste, zero, ibid(1)
    integer(kind=8) :: im1, im2, ico
    integer(kind=8) :: p1, p2, ifm, niv, p3, p4
    integer(kind=8) :: jdesm1, jdesm2
    integer(kind=8) :: nbmavo, indi, im3, jcoor
    integer(kind=8) :: nbmaor, ii, kdeb
    aster_logical :: reorie
    character(len=2) :: kdim
    character(len=8) :: cellTypeName, nomail
    character(len=24) :: nomavo
    character(len=24) :: valk(2)
    integer(kind=8), pointer :: typmail(:) => null()
!
#define pasori(iCell) zi(lori-1+iCell).eq.0
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
!
    call infniv(ifm, niv)
    reorie = .true.
    zero = 0
!
! --- VECTEUR DU TYPE DES MAILLES DU MAILLAGE :
!     ---------------------------------------
    call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
!
! --- COORDONNEES DES NOEUDS DU MAILLAGE :
!     ----------------------------------
    call jeveuo(noma//'.COORDO    .VALE', 'L', jcoor)
!
! --- APPEL A LA CONNECTIVITE :
!     -----------------------
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', p2)
    call jeveuo(noma//'.CONNEX', 'E', p1)
!
!     ALLOCATIONS :
!     -----------
    call wkvect('&&ORVLMA.ORI1', 'V V I', nbCell, lori)
    call wkvect('&&ORVLMA.ORI2', 'V V I', nbCell, jori)
    call wkvect('&&ORVLMA.ORI3', 'V V I', nbCell, nori)
    call wkvect('&&ORVLMA.ORI4', 'V V I', nbCell, kori)
    call wkvect('&&ORVLMA.ORI5', 'V V I', nbCell, kdeb)
!
! --- VERIFICATION DU TYPE DES MAILLES
! --- (ON DOIT AVOIR DES MAILLES DE PEAU) :
!     -----------------------------------
    do iCell = 1, nbCell
        zi(lori-1+iCell) = 0
        cellNume = listCellNume(iCell)
        zi(nori-1+iCell) = zi(p2+cellNume)-zi(p2-1+cellNume)
        zi(kori-1+iCell) = zi(p2+cellNume-1)
!
! ---   TYPE DE LA MAILLE COURANTE :
!       --------------------------
        cellTypeNume = typmail(cellNume)
        call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
!
        if (cellTypeName(1:3) .ne. 'SEG') then
            nomail = int_to_char8(cellNume)
            valk(1) = nomail
            valk(2) = cellTypeName
            call utmess('F', 'MODELISA5_93', nk=2, valk=valk)
        end if
    end do
!ok --- on teste le type de maille de la liste
!       si ce n'est pas un seg => erreur
!
! --- RECUPERATION DES MAILLES VOISINES DU GROUP_MA :
!     ---------------------------------------------
    kdim = '1D'
    nomavo = '&&ORVLMA.MAILLE_VOISINE '
    call utmavo(noma, kdim, listCellNume, nbCell, 'V', &
                nomavo, zero, ibid)
    call jeveuo(jexatr(nomavo, 'LONCUM'), 'L', p4)
    call jeveuo(nomavo, 'L', p3)
!
!
! --- ON TESTE SI LA POUTRE CONTIENT UN EMBRANCHEMENT
    do iCell = 1, nbCell
        numail = zi(p4+iCell)-zi(p4+iCell-1)
        if (numail .ge. 3) then
            norieg = 0
            do ii = 1, numail
                cellNume = zi(p3+zi(p4+iCell-1)-1+ii-1)
                do ico = 1, nbCell
                    if (cellNume .eq. listCellNume(ico)) then
                        norieg = norieg+1
                        goto 12
                    end if
                end do
12              continue
            end do
            if (norieg .ge. 3) then
                call utmess('F', 'MODELISA4_84')
            end if
        end if
    end do
!
!
! --- PREMIER PASSAGE: METTRE LES MAILLES AYANT LE NOEUD DANS
!     LA BONNE ORIENTATION
!
    norieg = 0
!
    nbmaor = 0
    do iCell = 1, nbCell
        cellNume = listCellNume(iCell)
        jdesm1 = zi(kori-1+iCell)
!
! ------ VERIFICATION QUE LE NOEUD EST DANS LA MAILLE
        ico = ioriv3(zi(p1+jdesm1-1), noeud, vect, zr(jcoor))
!
! ------ LA MAILLE NE CONTIENT PAS LE NOEUD
        if (ico .eq. 0) then
!
! ------ LA MAILLE A ETE REORIENTEE
        else if (ico .lt. 0) then
            nbmaor = nbmaor+1
            zi(kdeb+nbmaor-1) = iCell
            zi(lori-1+iCell) = 1
            if (niv .eq. 2) then
                nomail = int_to_char8(cellNume)
                write (ifm, *) 'LA MAILLE '//nomail//&
     &                       ' A ETE ORIENTEE PAR RAPPORT AU VECTEUR'
            end if
            norieg = norieg+1
!
! ------ LA MAILLE A LA BONNE ORIENTATION
        else
            nbmaor = nbmaor+1
            zi(kdeb+nbmaor-1) = iCell
            zi(lori-1+iCell) = 1
            if (niv .eq. 2) then
                nomail = int_to_char8(cellNume)
                write (ifm, *) 'LA MAILLE '//nomail//&
     &                       ' EST ORIENTEE PAR RAPPORT AU VECTEUR'
            end if
        end if
!
    end do
    if (nbmaor .eq. 0) then
        call utmess('F', 'MODELISA6_1')
    end if
!
    do ii = 1, nbmaor
        lliste = 0
        iliste = 0
        zi(jori+lliste) = zi(kdeb+ii-1)
!
! --- ON ORIENTE TOUTES LES MAILLES DU CONNEXE
!
200     continue
!
        im1 = zi(jori+iliste)
        jdesm1 = zi(kori-1+im1)
! --- ON ESSAYE D'ORIENTER LES MAILLES VOISINES
        nbmavo = zi(p4+im1)-zi(p4-1+im1)
        do im3 = 1, nbmavo
            indi = zi(p3+zi(p4+im1-1)-1+im3-1)
            im2 = indiis(listCellNume, indi, 1, nbCell)
            if (im2 .eq. 0) goto 210
            numail = listCellNume(im2)
            if (pasori(im2)) then
                jdesm2 = zi(kori-1+im2)
!           VERIFICATION DE LA CONNEXITE ET REORIENTATION EVENTUELLE
                ico = iorim1(zi(p1+jdesm1-1), zi(p1+jdesm2-1), reorie)
!           SI MAILLES CONNEXES
                if (ico .ne. 0) then
                    zi(lori-1+im2) = 1
                    lliste = lliste+1
                    zi(jori+lliste) = im2
                    if (reorie .and. niv .eq. 2) then
                        nomail = int_to_char8(numail)
                        if (ico .lt. 0) then
                            write (ifm, *) 'LA MAILLE ', nomail, ' A ETE REORIENTEE'
                        else
                            write (ifm, *) 'LA MAILLE ', nomail, ' EST ORIENTEE'
                        end if
                    end if
                end if
!
!           SI ORIENTATIONS CONTRAIRES
                if (ico .lt. 0) norieg = norieg+1
!
            end if
210         continue
        end do
        iliste = iliste+1
        if (iliste .le. lliste) goto 200
    end do
!
! --- ON VERIFIE QU'ON A BIEN TRAITE TOUTES LES MAILLES
!
    do iCell = 1, nbCell
        if (pasori(iCell)) then
            call utmess('F', 'MODELISA6_2')
        end if
    end do
!
    norien = norien+norieg
!
    call jedetr('&&ORVLMA.ORI1')
    call jedetr('&&ORVLMA.ORI2')
    call jedetr('&&ORVLMA.ORI3')
    call jedetr('&&ORVLMA.ORI4')
    call jedetr('&&ORVLMA.ORI5')
    call jedetr(nomavo)
!
    call jedema()
end subroutine
