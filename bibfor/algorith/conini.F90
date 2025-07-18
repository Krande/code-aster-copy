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
subroutine conini(ma, noecon, maicon, marcon, nbmar, &
                  nbnoe, nbmarc, nommar, jmicor, mbcor, &
                  nomtyr, nbgco, io8gco)
!
!  ROUTINE CONINI
!    ROUTINE DE PREPARATION DE TABLEAUX PERMETTANT D'OPTIMISER
!    LA ROUTINE CONORI
!  DECLARATIONS
!    NBGCO  : NOMBRE DE GROUPE DE FISSURE
!    IC     : NB D'OCCURENCE DE GROUP_MA
!    ICOC   : INDICE COURANT DES CONNEX     POUR UNE MAILLE FISSURE
!    ICOR   : INDICE COURANT DES CONNEX     POUR UNE MAILLE REFERENCE
!    IDUM   : ENTIER DE TRAVAIL
!    IGCO   : INDICE COURANT SUR LES GROUPE DE FISSURE
!    IGMA   : INDICE COURANT SUR LES GROUP_MA
!    IMAC   : INDICE COURANT DES MAILLES    POUR UNE MAILLE FISSURE
!    IMAG   : INDICE COURANT SUR LES MAILLES D UN GROUPE
!    IMAR   : INDICE COURANT DES MAILLES    POUR UNE MAILLE REFERENCE
!    IMICOC : INDICE DE  CONNEX    DANS ZK8 POUR UNE MAILLE FISSURE
!    IMICOR : INDICE DE  CONNEX    DANS ZK8 POUR UNE MAILLE REFERENCE
!    IMIGMA : INDICE DE  GROUPEMA  DANS ZI
!    IMITYC : INDICE DE  TYPMAIL   DANS ZI  POUR UNE MAILLE FISSURE
!    IMITYR : INDICE DE  TYPMAIL   DANS ZI  POUR UNE MAILLE REFERENCE
!    INOC   : NUMERO D UN NOEUD             POUR UNE MAILLE FISSURE
!    INOR   : NUMERO D UN NOEUD             POUR UNE MAILLE REFERENCE
!    IO8GCO : INDICE DE OP0154 NOGCO DANS ZK8
!    ITYC   : INDICE COURANT DU TYPE        POUR UNE MAILLE FISSURE
!    ITYR   : INDICE COURANT DU TYPE        POUR UNE MAILLE REFERENCE
!    JEXNOM : FUNCTION D ASTER
!    JEXNUM : FUNCTION D ASTER
!    KBID   : CHARACTER DE TRAVAIL
!    KMAC   : NOM D UNE MAILLE              POUR UNE MAILLE FISSURE
!    KMAR   : NOM D UNE MAILLE              POUR UNE MAILLE REFERENCE
!    KNOC   : NOM D UN NOEUD                POUR UNE MAILLE FISSURE
!    KNOR   : NOM D UN NOEUD                POUR UNE MAILLE REFERENCE
!    KTYC   : NOM DU TYPE                   POUR UNE MAILLE FISSURE
!    KTYR   : NOM DU TYPE                   POUR UNE MAILLE REFERENCE
!    LOCONT : LOGICAL PRECISANT SI LA MAILLE EST UNE MAILLE FISSURE
!    LOMODI : LOGICAL PRECISANT SI LA MAILLE EST UNE MAILLE MODIFIE
!    LOREOR : LOGICAL PRECISANT SI LA MAILLE EST UNE MAILLE REORIENTEE
!    MA     : L OBJET DU MAILLAGE
!    MACOC  : TABLEAU DES NOMS DES NOEUDS   POUR UNE MAILLE FISSURE
!    MACOR  : TABLEAU DES NOMS DES NOEUDS   POUR UNE MAILLE REFERENCE
!    NBCOC  : NOMBRE DE CONNEX              POUR UNE MAILLE FISSURE
!    NBCOR  : NOMBRE DE CONNEX              POUR UNE MAILLE REFERENCE
!    NBMAG  : NOMBRE DE MAILLE DANS UN GROUP_MA
!    NBMAR  : NOMBRE DE MAILLE              POUR UNE MAILLE REFERENCE
!
!  MOT_CLEF : ORIE_FISSURE
!
!
    implicit none
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/infniv.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: io8gco, nbgco, igco
    integer(kind=8) :: imigma, igma
    integer(kind=8) :: nbmag, imag
    integer(kind=8) :: imac
    integer(kind=8) :: imityc, ityc
    integer(kind=8) :: imicoc, nbcoc, icoc
    integer(kind=8) :: inoc
    integer(kind=8) :: nbmar, imar
    integer(kind=8) :: imityr, ityr
    integer(kind=8) :: imicor, nbcor, icor
    integer(kind=8) :: iatyma
!
    character(len=8) :: kmac, ktyc, kmar, ktyr
    character(len=24) :: valk(2)
    character(len=8) :: ma
!
    aster_logical :: inval
    aster_logical :: cas2d, cas3d
!
    integer(kind=8) :: nbnoe, noecon(nbnoe), maicon(nbmar), marcon(nbmar)
    integer(kind=8) :: mbcor(nbmar), jmicor(nbmar)
    character(len=8) :: nommar(nbmar), nomtyr(nbmar)
    integer(kind=8) :: ierr, ifm, imai, inoe, inor, itest, nbcom
    integer(kind=8) :: nbmarc, niv
!-----------------------------------------------------------------------
!     TYPES VALIDES POUR LES MAILLES DE REFERENCE
#define valid() (cas2d .and. (ktyr(:4).eq.'TRIA'.or.ktyc(: \
4   ) .eq. 'QUAD')) .or. (cas3d .and. \
    (ktyr(:5) .eq. 'PENTA' .or. ktyr(:4) .eq. 'HEXA' .or. ktyr(:\
5   ) .eq. 'PYRAM' .or. ktyr(:5) .eq. 'TETRA'))
!
!
    inval = .false.
    cas2d = .false.
    cas3d = .false.
!CC     ON COMMENTE JEMARQ CAR ADRESSES PASSEES EN ARGUMENT
!CC      CALL JEMARQ()
    call infniv(ifm, niv)
!
!     ==================================================================
!
    do inoe = 1, nbnoe
        noecon(inoe) = 0
    end do
!
    do imai = 1, nbmar
        maicon(imai) = 0
    end do
!
    nbmarc = 0
    ierr = 0
30  continue
!     ------------------------------------------------------------------
!     BOUCLE SUR LES GROUPE_MA_FISSURE
!     ------------------------------------------------------------------
    do igco = 1, nbgco
!     ------------------------------------------------------------------
!     RECHERCHE D EXISTENCE DU GROUP_MA_FISSURE CONSIDERE
!     ------------------------------------------------------------------
        call jenonu(jexnom(ma//'.GROUPEMA', zk24(io8gco+igco-1)), igma)
!
        if (igma .eq. 0) then
!     ------------------------------------------------------------------
!     TRAITEMENT DU CAS DE NON-EXISTENCE
!     ------------------------------------------------------------------
            call utmess('I', 'ALGORITH2_26', sk=zk24(io8gco+igco-1))
!
        else
!     ------------------------------------------------------------------
!     TRAITEMENT DU CAS D EXISTENCE
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
!     RECHERCHE DE L'ADRESSE DU GROUP_MA DANS ZI
!     ------------------------------------------------------------------
            call jeveuo(jexnum(ma//'.GROUPEMA', igma), 'L', imigma)
!     ------------------------------------------------------------------
!     RECHERCHE DU NOMBRE DE MAILLE DU GROUP_MA
!     ------------------------------------------------------------------
            call jelira(jexnum(ma//'.GROUPEMA', igma), 'LONMAX', nbmag)
!     ------------------------------------------------------------------
!     BOUCLE SUR LES MAILLES DU GROUP_MA
!     ------------------------------------------------------------------
            do imag = 1, nbmag
                imac = zi(imigma+imag-1)
                maicon(imac) = maicon(imac)+1
!     ------------------------------------------------------------------
!     RECHERCHE DU NOM DE LA MAILLE
!     ------------------------------------------------------------------
                kmac = int_to_char8(imac)
!     ------------------------------------------------------------------
!     RECHERCHE DE L'ADRESSE DU TYPE DE LA MAILLE DANS ZI
!     ------------------------------------------------------------------
                call jeveuo(ma//'.TYPMAIL', 'L', iatyma)
                imityc = iatyma-1+imac
                ityc = zi(imityc)
!     ------------------------------------------------------------------
!     RECHERCHE DU TYPE DE LA MAILLE DANS CATA.TM.NOMTM
!     ------------------------------------------------------------------
                call jenuno(jexnum('&CATA.TM.NOMTM', ityc), ktyc)
!
                if (ktyc(:5) .eq. 'QUAD4' .or. ktyc(:5) .eq. 'QUAD8') then
                    cas2d = .true.
                    if (ierr .ne. 0) write (ifm, *) 'MAILLE 2D : ', kmac, ' DE TYPE ', ktyc
                elseif (ktyc(:5) .eq. 'PENTA' .or. ktyc(:4) .eq. 'HEXA') &
                    then
                    cas3d = .true.
                    if (ierr .ne. 0) write (ifm, *) 'MAILLE 3D : ', kmac, ' DE TYPE ', ktyc
                else
                    inval = .true.
                    valk(1) = kmac
                    valk(2) = ktyc
                    call utmess('E', 'ALGORITH2_27', nk=2, valk=valk)
                end if
!
!     ------------------------------------------------------------------
!     RECHERCHE DE L ADRESSE DES CONNEXIONS DE LA MAILLE
!     ------------------------------------------------------------------
                call jeveuo(jexnum(ma//'.CONNEX', imac), 'E', imicoc)
!     ------------------------------------------------------------------
!     RECHERCHE DU NOMBRE DE CONNEXIONS DE LA MAILLE
!     ------------------------------------------------------------------
                call jelira(jexnum(ma//'.CONNEX', imac), 'LONMAX', nbcoc)
!
!     ------------------------------------------------------------------
!     BOUCLE SUR LES CONNEXIONS DE LA MAILLE
!     ------------------------------------------------------------------
                do icoc = 1, nbcoc
                    inoc = zi(imicoc+icoc-1)
                    noecon(inoc) = noecon(inoc)+1
                end do
            end do
!     ------------------------------------------------------------------
        end if
!     ------------------------------------------------------------------
    end do
    if (inval) then
        call utmess('F', 'ALGORITH2_28')
    end if
!
    if (cas2d .and. cas3d) then
        if (ierr .eq. 0) then
!       ON RETOURNE DANS LA BOUCLE AVEC DEMANDE DE MESSAGES
            ierr = 1
            goto 30
!
        else
            call utmess('F', 'ALGORITH2_29')
        end if
    end if
    if (cas2d) itest = 2
    if (cas3d) itest = 3
!
!     ------------------------------------------------------------------
!     BOUCLE SUR LES MAILLES DU MAILLAGE
!     ------------------------------------------------------------------
    do imar = 1, nbmar
        if (maicon(imar) .ne. 0) goto 80
!     ------------------------------------------------------------------
!     RECHERCHE DU NOM DE LA MAILLE
!     ------------------------------------------------------------------
        kmar = int_to_char8(imar)
        nommar(imar) = kmar
!
!     ------------------------------------------------------------------
!     RECHERCHE DE L ADRESSE DES CONNEXIONS DE LA MAILLE
!     ------------------------------------------------------------------
        call jeveuo(jexnum(ma//'.CONNEX', imar), 'L', imicor)
        jmicor(imar) = imicor
!     ------------------------------------------------------------------
!     RECHERCHE DU NOMBRE DE CONNEXIONS DE LA MAILLE
!     ------------------------------------------------------------------
        call jelira(jexnum(ma//'.CONNEX', imar), 'LONMAX', nbcor)
        mbcor(imar) = nbcor
!     ------------------------------------------------------------------
!     BOUCLE SUR LES CONNEXIONS DE LA MAILLE
!     ------------------------------------------------------------------
        nbcom = 0
        do icor = 1, nbcor
            inor = zi(imicor+icor-1)
            if (noecon(inor) .ne. 0) nbcom = nbcom+1
!
        end do
        if (nbcom .ge. itest) then
!     ------------------------------------------------------------------
!     RECHERCHE DE L'ADRESSE DU TYPE DE LA MAILLE DANS ZI
!     ------------------------------------------------------------------
            call jeveuo(ma//'.TYPMAIL', 'L', iatyma)
            imityr = iatyma-1+imar
            ityr = zi(imityr)
!     ------------------------------------------------------------------
!     RECHERCHE DU TYPE DE LA MAILLE DANS CATA.TM.NOMTM
!     ------------------------------------------------------------------
            call jenuno(jexnum('&CATA.TM.NOMTM', ityr), ktyr)
            nomtyr(imar) = ktyr
!
            if (valid()) then
                nbmarc = nbmarc+1
                marcon(nbmarc) = imar
            end if
        end if
80      continue
    end do
!     ==================================================================
!CC      ON COMMENTE JEMARQ CAR ADRESSES PASSEES EN ARGUMENT
!CC      CALL JEDEMA()
end subroutine
