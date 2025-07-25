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

subroutine rvcohe(xdicmp, xdncmp, vcheff, i, ier)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/gettco.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    character(len=24) :: xdicmp, xdncmp, vcheff
    integer(kind=8) :: i, ier
!     VERIFICATION DE COHERENCE DES ARGUMENTS D' APPEL DE LA COMMANDE
!     DE CALCUL DU POST TRAITEMENT (OP0051)
!       1. EXISTENCE DES CHAM_GD MIS EN JEU
!       2. LEGALITE DES COMPOSANTES MISES EN JEU POUR CES CHAM_GD
!       3. CONCORDANCE ENTRE LES MAILLAGES DES CHAMPS ET DES COURBES
!          OU DES NOEUDS
!     ------------------------------------------------------------------
! IN  XDNCMP : K : NOM DE OJB XD V K8 DES NOMS DE CMP MISES EN JEU
! IN  XDICMP : K : NOM DE OJB XD V I  DES NUMS DE CMP (0 <=> ILLEGALE)
! IN  VCHEFF : K : NOM DU OJB S  V K24 DES NOMS DE CHAMPS EFFECTIF
! IN  I      : I : NUMERO DE L' OCCURENCE A TRAITER
! OUT IER    : I : CODE RETOUR, 1 --> OK, 0 --> KO
!     ------------------------------------------------------------------
!
!
    character(len=24) :: ncheff, ndesc, valk(7), nomgrn
    character(len=19) :: nchp19
    character(len=16) :: nchsym, tresu
    character(len=15) :: nrepnd
    character(len=8) :: nresu, nomcmp, nmaich, nomnd
    character(len=4) :: docu
    integer(kind=8) :: acheff, alneud, anumcp, anomcp, nbcmp
    integer(kind=8) :: nbgrpn, nbneud, grel, nbgrel, jceld, amod, mod
    integer(kind=8) :: j, k, n1, ngrn, iexi, ier2
    aster_logical :: chelok, parMesh, lnomnoe
    character(len=24), pointer :: grpn(:) => null()
!
!=====================================================================
!
    call jemarq()
    ier = 1
    call jeveuo(vcheff, 'L', acheff)
    ncheff = zk24(acheff+i-1)
!
    if (ncheff(1:1) .eq. '&') then
!
!        CAS D'UN CHAMP SYMBOLIQUE ILLEGAL OU DE NON EXISTENCE
!                  D' UN CHAMP EFFECTIF ASSOCIE
!
        ier = 0
        call getvid('ACTION', 'RESULTAT', iocc=i, scal=nresu, nbret=n1)
        call getvtx('ACTION', 'NOM_CHAM', iocc=i, scal=nchsym, nbret=n1)
        call gettco(nresu, tresu)
!
        valk(1) = nchsym
        valk(2) = nresu
        valk(3) = tresu
        call utmess('F', 'POSTRELE_46', nk=3, valk=valk, si=i)
!
    else
!
!        LE CHAMP SYMBOLIQUE EXISTE ET UN CHAMP EFFECTIF LUI
!        CORRESPOND OU UN CHAMP EFFECTIF EST DIRECTEMENT ARGUMENT
!        VERIFICATION POUR LES CHAM_ELEM DU CARACTERE "AUX NOEUDS"
!
        nchp19 = ncheff(1:19)
        call dismoi("DOCU", nchp19, "CHAMP", repk=docu)
!
        if (docu .eq. 'CHML') then
            ndesc = nchp19//'.CELD'
            call jeveuo(ndesc, 'L', jceld)
            nbgrel = zi(jceld+2-1)
            chelok = .true.
            grel = 0
10          continue
            if ((chelok) .and. (grel .lt. nbgrel)) then
                grel = grel+1
                mod = zi(jceld-1+zi(jceld-1+4+grel)+2)
                if (mod .ne. 0) then
                    call jeveuo(jexnum('&CATA.TE.MODELOC', mod), 'L', amod)
                    chelok = (zi(amod-1+1) .eq. 2)
                end if
                goto 10
            end if
            if (.not. chelok) then
                call utmess('F', 'POSTRELE_47', si=i)
            end if
        end if
!
!        --- VERIFICATION SUR LES CMPS ---
        call jelira(jexnum(xdicmp, i), 'LONMAX', nbcmp)
        call jeveuo(jexnum(xdicmp, i), 'L', anumcp)
        do j = 1, nbcmp, 1
            if (zi(anumcp+j-1) .eq. 0) then
                call jeveuo(jexnum(xdncmp, i), 'L', anomcp)
                nomcmp = zk8(anomcp+j-1)
                call utmess('F', 'POSTRELE_48', sk=nomcmp, si=i)
            end if
        end do
!
!        --- VERIFICATION DE CONCORDANCE DES MAILLAGES ---
        call dismoi('NOM_MAILLA', ncheff, 'CHAMP', repk=nmaich)
!       /* LE LIEU DU POST TRAITEMENT EST UN ENSEMBLE DE NOEUDS */
!       VERIFICATION D' EXISTENCE DES NOEUDS DANS LE MAILLAGE DU CHP
        call getvtx('ACTION', 'GROUP_NO', iocc=i, nbval=0, nbret=nbgrpn)
        call getvtx('ACTION', 'NOEUD', iocc=i, nbval=0, nbret=nbneud)
        nbgrpn = -nbgrpn
        nbneud = -nbneud
        parMesh = isParallelMesh(nmaich)
        if (parMesh) then
            if (nbneud .ne. 0) then
                ASSERT(.false.)
            end if
        end if
        if (nbgrpn .ne. 0) then
            call jecreo('&&OP0051.NOM.GRPN', 'V V K24')
            call jeecra('&&OP0051.NOM.GRPN', 'LONMAX', nbgrpn)
            call jeveuo('&&OP0051.NOM.GRPN', 'E', vk24=grpn)
            call getvtx('ACTION', 'GROUP_NO', iocc=i, nbval=nbgrpn, vect=grpn, &
                        nbret=n1)
            do k = 1, nbgrpn, 1
                nomgrn = grpn(k)
                iexi = 0
                ngrn = 0
                call jeexin(nmaich//'.GROUPENO', iexi)
                if (iexi .eq. 0) then
                    ier = 0
                else
                    call jelira(nmaich//'.GROUPENO', 'NUTIOC', ngrn)
                    if (ngrn .eq. 0) then
                        ier = 0
                    else
                        call jenonu(jexnom(nmaich//'.GROUPENO', nomgrn), n1)
                        if (n1 .eq. 0) then
                            ier = 0
                        end if
                    end if
                end if
                if (parMesh) then
                    call asmpi_comm_vect('MPI_MAX', 'I', sci=ier)
                    if (ier .eq. 0) then
                        exit
                    end if
                else
                    if (ier .eq. 0) then
                        exit
                    end if
                end if
            end do
            call jedetr('&&OP0051.NOM.GRPN')
        end if
        if (nbneud .ne. 0) then
            call wkvect('&&OP0051.NOM.NEUD', 'V V K8', nbneud, alneud)
            call getvtx('ACTION', 'NOEUD', iocc=i, nbval=nbneud, vect=zk8(alneud), &
                        nbret=n1)
            nrepnd = nmaich//'.NOMNOE'
            call jeexin(nrepnd, ier2)
            lnomnoe = .false.
            if (ier2 .ne. 0) then
                lnomnoe = .true.
            end if
            do k = 1, nbneud, 1
                nomnd = zk8(alneud+k-1)
                if (lnomnoe) then
                    call jenonu(jexnom(nrepnd, nomnd), n1)
                else
                    n1 = char8_to_int(nomnd)
                end if
                if (n1 .eq. 0) then
                    call utmess('A', 'POSTRELE_51', sk=nomnd, si=i)
                    ier = 0
                    exit
                end if
            end do
            call jedetr('&&OP0051.NOM.NEUD')
        end if
    end if
    call jedema()
end subroutine
