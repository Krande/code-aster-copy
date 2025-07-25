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
subroutine dismcm(questi, nomobz, repi, repkz, ierd)
    implicit none
!     --     DISMOI(CHAM_MATER)
!     ARGUMENTS:
!     ----------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismca.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jelstc.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi
    character(len=*) :: nomobz, repkz
    character(len=6) :: k6
    character(len=32) :: repk
    character(len=8) :: nomob
! ----------------------------------------------------------------------
!     IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE TYPE NUM_DDL
!     OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
!
! ----------------------------------------------------------------------
!
!
    character(len=8) :: mater, nomf, nomp, novarc
    character(len=32) :: nomrc
    character(len=16) :: ktyp
    character(len=24) :: quest2, nomobj(100)
    character(len=19) :: nomcar2
    aster_logical :: trouve
    integer(kind=8) :: nbmax, izone, i
!-----------------------------------------------------------------------
    integer(kind=8) :: ianorc, iaobj, iaprol, iavale, iavalk, if, ii
    integer(kind=8) :: iii, imax, irc, iret, jdesc, lonobj, n
    integer(kind=8) :: n1, nbrc, nbzone, nc, nf, nmat, nr
    integer(kind=8) :: n2, nbvarc, nedit, kvarc, kedit
    character(len=8), pointer :: cvrcvarc(:) => null()
    character(len=16), pointer :: vale(:) => null()
!
!-----------------------------------------------------------------------
    parameter(nbmax=30)
!
    call jemarq()
    repk = ' '
    repi = 0
    ierd = 0
!
    nomob = nomobz
!
!
    if (questi .eq. 'NOM_MAILLA') then
!     -----------------------------------
        call dismca(questi, nomob//'.CHAMP_MAT', repi, repk, ierd)
!
!
    elseif (questi .eq. 'EXI_AMOR_ALPHA' .or. questi .eq. 'EXI_AMOR_NOR' &
            .or. questi .eq. 'EXI_AMOR_TAN' .or. questi .eq. 'EXI_AMOR_HYST' &
            .or. questi .eq. 'EXI_AMOR_BETA') then
!     ---------------------------------------------------------------
        call jeveuo(nomob//'.CHAMP_MAT .VALE', 'L', iavale)
        call jelira(nomob//'.CHAMP_MAT .VALE', 'LONMAX', nmat)
        call jeveuo(nomob//'.CHAMP_MAT .DESC', 'L', jdesc)
        trouve = .false.
        quest2 = questi
        nbzone = zi(jdesc+2)
!
        do izone = 1, nbzone
            do imax = 1, nbmax
                i = (izone-1)*nbmax+imax
                mater = zk8(iavale-1+i)
                if (mater .eq. ' ') goto 40
                if (mater .eq. 'TREF=>') goto 40
!
                call jelstc('G', mater, 1, 100, nomobj, &
                            n)
                if (n .lt. 0) then
                    call utmess('F', 'UTILITAI_54')
                end if
                do ii = 1, n
                    if (nomobj(ii) (20:24) .eq. '.VALK') then
                        call jeveuo(nomobj(ii), 'L', iaobj)
                        call jelira(nomobj(ii), 'LONMAX', lonobj)
                        do iii = 1, lonobj
                            if (zk16(iaobj-1+iii) .eq. quest2(5:14)) then
                                trouve = .true.
                                goto 50
!
                            end if
                        end do
                    end if
                end do
            end do
40          continue
        end do
50      continue
        repk = 'NON'
        if (trouve) repk = 'OUI'
!
!
    else if (questi .eq. 'EXI_ANISO') then
!     -----------------------------------
        repk = 'NON'
        call jeveuo(nomob//'.CHAMP_MAT .VALE', 'L', iavale)
        call jelira(nomob//'.CHAMP_MAT .VALE', 'LONMAX', nmat)
        call jeveuo(nomob//'.CHAMP_MAT .DESC', 'L', jdesc)
        nbzone = zi(jdesc+2)
!
        do izone = 1, nbzone
            do imax = 1, nbmax
                i = (izone-1)*nbmax+imax
                mater = zk8(iavale-1+i)
                if (mater .eq. ' ') goto 80
                if (mater .eq. 'TREF=>') goto 80
!
                call jeveuo(mater//'.MATERIAU.NOMRC', 'L', ianorc)
                call jelira(mater//'.MATERIAU.NOMRC', 'LONUTI', nbrc)
                do irc = 1, nbrc
                    nomrc = zk32(ianorc-1+irc) (1:10)
                    if (nomrc .eq. 'ELAS_COQUE') then
                        repk = 'OUI'
                        goto 90
!
                    else if (nomrc .eq. 'THER_COQUE') then
                        repk = 'OUI'
                        goto 90
!
                    else if (nomrc .eq. 'ELAS_ORTH') then
                        repk = 'OUI'
                        goto 90
!
                    else if (nomrc .eq. 'THER_ORTH') then
                        repk = 'OUI'
                        goto 90
!
                    else if (nomrc .eq. 'ELAS_COQMU') then
                        repk = 'OUI'
                        goto 90
!
                    else if (nomrc .eq. 'THER_COQMU') then
                        repk = 'OUI'
                        goto 90
!
                    end if
                end do
            end do
80          continue
        end do
90      continue
!
!
    else if (questi .eq. 'THER_F_INST') then
!     --------------------------------------
        repk = 'NON'
        call jeveuo(nomob//'.CHAMP_MAT .VALE', 'L', iavale)
        call jelira(nomob//'.CHAMP_MAT .VALE', 'LONMAX', nmat)
        call jeveuo(nomob//'.CHAMP_MAT .DESC', 'L', jdesc)
        nbzone = zi(jdesc+2)
!
        do izone = 1, nbzone
            do imax = 1, nbmax
                i = (izone-1)*nbmax+imax
                mater = zk8(iavale-1+i)
                if (mater .eq. ' ') goto 140
                if (mater .eq. 'TREF=>') goto 140
!
                call jeveuo(mater//'.MATERIAU.NOMRC', 'L', ianorc)
                call jelira(mater//'.MATERIAU.NOMRC', 'LONUTI', nbrc)
                do irc = 1, nbrc
                    nomrc = zk32(ianorc-1+irc)
!
!            -- SI LE MATERIAU EST ISSU DE LA COMMANDE DEFI_COQU_MULT :
                    if (nomrc(5:10) .eq. '_COQMU') goto 120
!
                    if (nomrc(1:4) .ne. 'THER') goto 110
                    call codent(irc, 'D0', k6)
                    call jeveuo(mater//'.CPT.'//k6//'.VALK', 'L', iavalk)
                    call jelira(mater//'.CPT.'//k6//'.VALK', 'LONUTI', n1)
                    call jelira(mater//'.CPT.'//k6//'.VALR', 'LONUTI', nr)
                    call jelira(mater//'.CPT.'//k6//'.VALC', 'LONUTI', nc)
                    nf = (n1-nr-nc)/2
                    do if = 1, nf
                        nomf = zk16(iavalk-1+nr+nc+nf+if) (1:8)
                        call jeveuo(nomf//'           .PROL', 'L', iaprol)
                        if (zk24(iaprol-1+1) .eq. 'NAPPE') then
!              -- CAS D'UNE FONCTION A 2 VARIABLES :
                            if (zk24(iaprol-1+3) .eq. 'INST') repk = 'OUI'
                            if (zk24(iaprol-1+7) .eq. 'INST') repk = 'OUI'
                        else
!              -- CAS D'UNE FONCTION A 1 VARIABLE :
                            if (zk24(iaprol-1+3) .eq. 'INST') repk = 'OUI'
                        end if
                    end do
110                 continue
                end do
120             continue
            end do
140         continue
        end do
!
!
!     -- CETTE QUESTION N'EXISTE PLUS. IL NE FAUT PAS L'UTILISER.
!        JE LA CONSERVE JUSTE LE TEMPS DE FAIRE MA RESTIT LA MEME
!        SEMAINE QUE SEBASTIEN MEUNIER QUI MODIFIE VECTME.F
    else if (questi .eq. 'ELAS_F_TEMP') then
!     --------------------------------------
        repk = '???'
        repi = -99999
!
!
    else if (questi == 'ELAS_FO' .or. questi == 'NU_FO') then
!     --------------------------------------
        repk = 'NON'
        call jeexin(nomob//'.CHAMP_MAT .VALE', iret)
        if (iret == 0) goto 200
        call jeveuo(nomob//'.CHAMP_MAT .VALE', 'L', iavale)
        call jelira(nomob//'.CHAMP_MAT .VALE', 'LONMAX', nmat)
        call jeveuo(nomob//'.CHAMP_MAT .DESC', 'L', jdesc)
        nbzone = zi(jdesc+2)
!
        do izone = 1, nbzone
            do imax = 1, nbmax
                i = (izone-1)*nbmax+imax
                mater = zk8(iavale-1+i)
                if (mater .eq. ' ') goto 190
                if (mater .eq. 'TREF=>') goto 190
!
                call jeveuo(mater//'.MATERIAU.NOMRC', 'L', ianorc)
                call jelira(mater//'.MATERIAU.NOMRC', 'LONUTI', nbrc)
                do irc = 1, nbrc
                    nomrc = zk32(ianorc-1+irc)
!
!                   -- Si le materiau est issu de la commande DEFI_COQU_MULT :
                    if (nomrc(5:10) .eq. '_COQMU') goto 170
!
                    if (nomrc(1:4) .ne. 'ELAS') goto 160
                    if (nomrc .eq. 'ELAS_DHRC') goto 160
!
                    call codent(irc, 'D0', k6)
                    call jeveuo(mater//'.CPT.'//k6//'.VALK', 'L', iavalk)
                    call jelira(mater//'.CPT.'//k6//'.VALK', 'LONUTI', n1)
                    call jelira(mater//'.CPT.'//k6//'.VALR', 'LONUTI', nr)
                    call jelira(mater//'.CPT.'//k6//'.VALC', 'LONUTI', nc)
                    nf = (n1-nr-nc)/2
!
! loop over parameters for which function values were given
                    do if = 1, nf
!
! name of the parameter
                        nomp = zk16(iavalk-1+nr+nc+if) (1:8)
!
! name of the function
                        nomf = zk16(iavalk-1+nr+nc+nf+if) (1:8)
                        call jeveuo(nomf//'           .PROL', 'L', iaprol)
                        if (zk24(iaprol-1+1) .eq. 'CONSTANT') then
!                           -- cas d'une fonction constante :
                        else
!                           -- cas d'une fonction variable :
                            if (questi == 'ELAS_FO') repk = 'OUI'
                            if (questi == 'NU_FO' .and. nomp(1:2) == 'NU') then
                                repk = 'OUI'
                            end if
                        end if
                    end do
160                 continue
                end do
170             continue
            end do
190         continue
        end do
200     continue
!
!
    else if (questi .eq. 'EXI_VARC') then
!     --------------------------------------
        repk = 'NON'
        call jeexin(nomob//'.CVRCVARC', iret)
        if (iret .ne. 0) then
            repk = 'OUI'
        end if
!
!
    else if (questi .eq. 'VARC_F_INST') then
!     --------------------------------------
        repk = 'NON'
        call jeexin(nomob//'.CVRCVARC', iret)
        if (iret .ne. 0) then
            call jeveuo(nomob//'.CVRCVARC', 'L', vk8=cvrcvarc)
            call jelira(nomob//'.CVRCVARC', 'LONMAX', nbvarc)
            do kvarc = 1, nbvarc
                novarc = cvrcvarc(kvarc)
                nomcar2 = nomob//'.'//novarc
                nomcar2 = nomcar2(1:17)//'.2'
                call jeveuo(nomcar2//'.DESC', 'L', jdesc)
                nedit = zi(jdesc-1+3)
                call jeveuo(nomcar2//'.VALE', 'L', vk16=vale)
                call jelira(nomcar2//'.VALE', 'LONMAX', n1)
                n2 = n1/nedit
                ASSERT(n1 .eq. nedit*n2)
                do kedit = 1, nedit
                    ktyp = vale((kedit-1)*n2+2)
                    ASSERT(ktyp .eq. 'EVOL' .or. ktyp .eq. 'CHAMP')
                    if (ktyp .eq. 'EVOL') then
                        repk = 'OUI'
                        goto 210
                    end if
                end do
            end do
        end if
210     continue
!
!
    else
!     --------------------------------------
        ierd = 1
    end if
!
    repkz = repk
    call jedema()
end subroutine
