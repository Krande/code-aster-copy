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
subroutine conori(ma)
!
!  ROUTINE CONORI
!    ROUTINE RACINE D'ORIENTATION DES MAILLES DE FISSURE
!  DECLARATIONS
!    NBGCO  : NOMBRE DE GROUPE DE FISSURE
!    IC     : NB D'OCCURENCE DE GROUP_MA
!    ICOC   : INDICE COURANT DES CONNEX     POUR UNE MAILLE FISSURE
!    ICOR   : INDICE COURANT DES CONNEX     POUR UNE MAILLE REFERENCE
!    IDUM   : ENTIER DE TRAVAIL
!    IFM    : IUNIT DU FICHIER MESSAGE
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
!    NBGMA  : NOMBRE DE GROUP_MA
!    NBMAG  : NOMBRE DE MAILLE DANS UN GROUP_MA
!    NBMAR  : NOMBRE DE MAILLE              POUR UNE MAILLE REFERENCE
!    NBNOMX : NOMBRE DE NOEUD MAXIMUM POUR UNE MAILLE ( 100 )
!    NIV    : NIVEAU D'IMPRESSION (OPTION INFO)
!
!
!
    implicit none
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/chkmsg.h"
#include "asterfort/conini.h"
#include "asterfort/contac.h"
#include "asterfort/getvem.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: idum, ic, ifm, niv, ichk
    integer(kind=8) :: io8gco, nbgco, igco
    integer(kind=8) :: imigma, nbgma, igma
    integer(kind=8) :: nbmag, imag
    integer(kind=8) :: imac
    integer(kind=8) :: imityc, ityc
    integer(kind=8) :: imicoc, nbcoc, icoc
    integer(kind=8) :: inoc
    integer(kind=8) :: nbmar, imar
    integer(kind=8) :: imicor, nbcor, icor
    integer(kind=8) :: iatyma
!
    character(len=8) :: kmac, ktyc, knoc, kmar, ktyr, knor
    character(len=8) :: ma, kbid
!
    aster_logical :: lomodi, loreo0, loreor, lomod0, locor0, lface, lface0
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ikmar, iktyr, imai, imarc, imaz, inoe
    integer(kind=8) :: inor, jmb, jmic, nbmac, nbmarc, nbnoe, nbnomx
!
!-----------------------------------------------------------------------
    parameter(nbnomx=100)
    character(len=8) :: macor(nbnomx+2), macoc(nbnomx+2), macos(nbnomx+2)
!
    call jemarq()
!
!     ==================================================================
!     ------------------------------------------------------------------
!     FORMAT ET UNIT D ECRITURE ET NIVEAU D'IMPRESSION
!     ------------------------------------------------------------------
    call infniv(ifm, niv)
!     ------------------------------------------------------------------
!     RECHERCHE DU NOMBRE DE GROUP_MA_FISSURE DANS .COMM
!     ------------------------------------------------------------------
    ic = 1
    call getvem(ma, 'GROUP_MA', 'ORIE_FISSURE', 'GROUP_MA', ic, &
                0, kbid, nbgco)
    nbgco = -nbgco
!
!     ==================================================================
    if (nbgco .ne. 0) then
!     ------------------------------------------------------------------
!     RECHERCHE DU NOMBRE DE GROUP_MA DANS .MAIL
!     ------------------------------------------------------------------
        call jelira(ma//'.GROUPEMA', 'NUTIOC', nbgma)
        if (niv .eq. 2) then
            write (ifm, *) ' '
            write (ifm, *) ' LA LISTE DES GROUP_MA '
            write (ifm, *) ' '
        end if
!     ------------------------------------------------------------------
!     RECHERCHE DES NOMS DES GROUP_MA DANS .MAIL
!     ------------------------------------------------------------------
        do igma = 1, nbgma
            call jenuno(jexnum(ma//'.GROUPEMA', igma), kbid)
            if (niv .eq. 2) then
                write (ifm, *) '   GROUP_MA     : ', kbid
            end if
        end do
        write (ifm, *) ' '
!     ------------------------------------------------------------------
!     CREATION D UN TABLEAU DE TRAVAIL
!     ------------------------------------------------------------------
        call wkvect('&&OP0154.NOGCO', 'V V K24', nbgco, io8gco)
!     ------------------------------------------------------------------
!     RECHERCHE DES NOMS DES GROUP_MA_FISSURE DANS .COMM
!     ------------------------------------------------------------------
        call getvem(ma, 'GROUP_MA', 'ORIE_FISSURE', 'GROUP_MA', ic, &
                    nbgco, zk24(io8gco), idum)
        if (niv .eq. 2) then
            write (ifm, *) ' '
            write (ifm, *) ' LA LISTE DES ORIE_FISSURE'
            write (ifm, *) ' '
            do igco = 1, nbgco
                write (ifm, *) '   ORIE_FISSURE: ', zk24(io8gco+igco-1)
            end do
            write (ifm, *) ' '
        end if
!     ------------------------------------------------------------------
        call jelira(ma//'.TYPMAIL', 'LONMAX', nbmar)
        call jelira(ma//'.COORDO    .VALE', 'LONMAX', nbnoe)
        nbnoe = nbnoe/3
!
        call wkvect('&&OP0154.NOE', 'V V I', nbnoe, inoe)
        call wkvect('&&OP0154.MAI', 'V V I', nbmar, imai)
        call wkvect('&&OP0154.MAR', 'V V I', nbmar, imaz)
        call wkvect('&&OP0154.KMR', 'V V K8', nbmar, ikmar)
        call wkvect('&&OP0154.KTR', 'V V K8', nbmar, iktyr)
        call wkvect('&&OP0154.IMI', 'V V I', nbmar, jmic)
        call wkvect('&&OP0154.MBL', 'V V I', nbmar, jmb)
        call conini(ma, zi(inoe), zi(imai), zi(imaz), nbmar, &
                    nbnoe, nbmarc, zk8(ikmar), zi(jmic), zi(jmb), &
                    zk8(iktyr), nbgco, io8gco)
        write (ifm, *) 'NOMBRE DE MAILLES DE REFERENCE TESTEES : ', &
            nbmarc
!
!     ==================================================================
!     ------------------------------------------------------------------
!     BOUCLE SUR LES GROUPE_MA_FISSURE
!     ------------------------------------------------------------------
        do igco = 1, nbgco
!     ------------------------------------------------------------------
!     RECHERCHE D EXISTENCE DU GROUP_MA_FISSURE CONSIDERE
!     ------------------------------------------------------------------
            call jenonu(jexnom(ma//'.GROUPEMA', zk24(io8gco+igco-1)), igma)
!
            if (niv .eq. 2) then
                write (ifm, *) ' '
                write (ifm, *) ' TRAITEMENT DE ', zk24(io8gco+igco-1)
                write (ifm, *) ' '
            end if
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
                if (niv .eq. 2) then
                    write (ifm, *) '   LA LISTE DES MAILLES DU GROUPE '
                    write (ifm, *) ' '
                end if
!
!     ------------------------------------------------------------------
!     BOUCLE SUR LES MAILLES DU GROUP_MA
!     ------------------------------------------------------------------
                do imag = 1, nbmag
                    imac = zi(imigma+imag-1)
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
                    if (niv .eq. 2) then
                        write (ifm, *) '     MAILLE NU : ', imag, ' NOM : ', kmac,&
     &            ' ORDRE : ', imac, ' TYPE : ', ityc, ' TYPE : ', ktyc
                    end if
                    macoc(1) = kmac
                    macoc(2) = ktyc
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
!     ------------------------------------------------------------------
!     RECHERCHE DU NOM DU NOEUD
!     ------------------------------------------------------------------
                        knoc = int_to_char8(inoc)
                        macoc(icoc+2) = knoc
                    end do
!
!     ------------------------------------------------------------------
!     SAUVEGARDE DE LA MAILLE DE FISSURE
!     ------------------------------------------------------------------
                    do idum = 1, nbcoc+2
                        macos(idum) = macoc(idum)
                    end do
!     ==================================================================
!     ------------------------------------------------------------------
!     BOUCLE SUR LES MAILLES DU MAILLAGE
!     ------------------------------------------------------------------
                    lomodi = .false.
                    loreor = .false.
                    nbmac = 0
                    do imarc = 1, nbmarc
                        imar = zi(imaz-1+imarc)
!     ------------------------------------------------------------------
!     RECHERCHE DU NOM DE LA MAILLE
!     ------------------------------------------------------------------
                        kmar = zk8(ikmar-1+imar)
!     ------------------------------------------------------------------
!     RECHERCHE DU TYPE DE LA MAILLE
!     ------------------------------------------------------------------
                        ktyr = zk8(iktyr-1+imar)
!
                        macor(1) = kmar
                        macor(2) = ktyr
!
!     ------------------------------------------------------------------
!     RECHERCHE DE L ADRESSE DES CONNEXIONS DE LA MAILLE
!     ------------------------------------------------------------------
                        imicor = zi(jmic-1+imar)
!     ------------------------------------------------------------------
!     RECHERCHE DU NOMBRE DE CONNEXIONS DE LA MAILLE
!     ------------------------------------------------------------------
                        nbcor = zi(jmb-1+imar)
!     ------------------------------------------------------------------
!     BOUCLE SUR LES CONNEXIONS DE LA MAILLE
!     ------------------------------------------------------------------
                        do icor = 1, nbcor
                            inor = zi(imicor+icor-1)
!     ------------------------------------------------------------------
!     RECHERCHE DU NOM DU NOEUD
!     ------------------------------------------------------------------
                            knor = int_to_char8(inor)
                            macor(icor+2) = knor
                        end do
!     ==================================================================
!     ------------------------------------------------------------------
!     APPEL DE CONTAC
!                      ORIENTATION DE LA MAILLE FISSURE SI NECESSAIRE
!                      LOMODI = .TRUE.  SI MODIFICATION
!                      LOMODI = .FALSE. SINON
!     ------------------------------------------------------------------
                        lomod0 = .false.
                        locor0 = .false.
                        loreo0 = .false.
                        call contac(macor, nbcor, macoc, nbcoc, lface0, &
                                    lomod0, locor0, loreo0, ma)
                        if (loreo0) lface0 = .not. lface0
                        if (locor0 .or. lomod0) then
                            nbmac = nbmac+1
                            if (niv .eq. 2) then
                                write (ifm, *) 'LA MAILLE DE FISSURE   ', macoc(1),&
     &                ' DE TYPE ', macoc(2)
                                write (ifm, *) (macoc(i+2), i=1, nbcoc)
                                write (ifm, *) 'S''APPUIE SUR LA MAILLE ', macor(1),&
     &                ' DE TYPE ', macor(2)
                                write (ifm, *) (macor(i+2), i=1, nbcor)
                                if (lface0) then
                                    write (ifm, *) 'PAR SA FACE INFERIEURE'
                                else
                                    write (ifm, *) 'PAR SA FACE SUPERIEURE'
                                end if
                                if (lomod0) then
                                    write (ifm, *)&
     &                  'UNE REORIENTATION POUR L''APPUI A EU LIEU'
                                end if
                                if (loreo0) then
                                    write (ifm, *)&
     &                  'UNE REORIENTATION POUR LA NORMALE A EU LIEU'
                                end if
                                write (ifm, *)
                            end if
                            if (nbmac .eq. 3) then
                                call utmess('F', 'ALGORITH2_30')
                            end if
                            if (nbmac .eq. 2 .and. (lface0 .eqv. lface)) then
                                call utmess('F', 'ALGORITH2_30')
                            end if
                            lface = lface0
                            if (lomod0) lomodi = .true.
                            if (loreo0) loreor = .true.
                            if (nbmac .eq. 2 .and. (lomod0 .or. loreo0)) then
                                call utmess('F', 'ALGORITH2_30')
                            end if
                        end if
!
!     ==================================================================
                    end do
                    if (nbmac .eq. 0) then
                        call utmess('F', 'ALGORITH2_30')
                    end if
!
                    if (lomodi .or. loreor) then
!     ------------------------------------------------------------------
!     ECRITURE DES MAILLES MODIFIEES
!     ------------------------------------------------------------------
                        if (niv .eq. 2) then
                            write (ifm, *) ' '
                            write (ifm, *) '       MODIFICATION DE LA MAILLE'
                            write (ifm, *) ' '
                            write (ifm, *) '       AVANT'
                            write (ifm, 9000) (macos(idum), idum=1, nbcoc+ &
                                               2)
                            write (ifm, *) '       APRES'
                            write (ifm, 9000) (macoc(idum), idum=1, nbcoc+ &
                                               2)
                            write (ifm, *) ' '
                        end if
!
                        do icoc = 1, nbcoc
                            knoc = macoc(icoc+2)
!     ------------------------------------------------------------------
!     RECHERCHE DE L'ORDRE DU NOEUD
!     ------------------------------------------------------------------
                            inoc = char8_to_int(knoc)
!     ------------------------------------------------------------------
!     MODIFICATION DE L ORIENTATION DE LA MAILLE
!     ------------------------------------------------------------------
                            zi(imicoc+icoc-1) = inoc
                        end do
                    end if
!     ==================================================================
!
                end do
!     ------------------------------------------------------------------
            end if
!     ------------------------------------------------------------------
        end do
!
    end if
!     ------------------------------------------------------------------
!     ==================================================================
!     EMISSION D'UNE ERREUR <F> SI UNE ERREUR <E> S'EST PRODUITE
    call chkmsg(0, ichk)
    call jedema()
!
9000 format(6x, 6(2x, a8), (/, 26x, 4(2x, a8)))
end subroutine
