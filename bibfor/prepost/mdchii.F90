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
subroutine mdchii(idfimd, nochmd, typent, typgeo, prefix, &
                  nbtv, codret)
!_____________________________________________________________________
! person_in_charge: nicolas.sellenet at edf.fr
! ======================================================================
!     FORMAT MED - CHAMP - INFORMATIONS - FICHIER CONNU PAR IDENTIFIANT
!            --    --      -                                -
!     DONNE LE NOMBRE DE TABLEAUX DE VALEURS ET LEURS CARACTERISTIQUES
!     TEMPORELLES POUR UN CHAMP ET UN SUPPORT GEOMETRIQUE
!-----------------------------------------------------------------------
!      ENTREES:
!        IDFIMD : IDENTIFIANT DU FICHIER MED
!        NOCHMD : NOM MED DU CHAMP A LIRE
!        TYPENT : TYPE D'ENTITE AU SENS MED
!        TYPGEO : TYPE DE SUPPORT AU SENS MED
!      ENTREES/SORTIES:
!        PREFIX : BASE DU NOM DES STRUCTURES
!                 POUR LE TABLEAU NUMERO I
!                 PREFIX//'.NUME' : T(2I-1) = NUMERO DE PAS DE TEMPS
!                                   T(2I)   = NUMERO D'ORDRE
!                 PREFIX//'.INST' : T(I) = INSTANT S'IL EXISTE
!      SORTIES:
!        NBTV   : NOMBRE DE TABLEAUX DE VALEURS DU CHAMP
!        CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!_____________________________________________________________________
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_mfdcsi.h"
#include "asterfort/as_mfdfdi.h"
#include "asterfort/as_mfdfin.h"
#include "asterfort/as_mfdncn.h"
#include "asterfort/as_mfdnfc.h"
#include "asterfort/as_mfdnfd.h"
#include "asterfort/as_mfdonp.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    med_idt :: idfimd
    integer(kind=8) :: nbtv
    integer(kind=8) :: typent, typgeo
    integer(kind=8) :: codret
!
    character(len=19) :: prefix
    character(len=*) :: nochmd
!
! 0.2. ==> COMMUNS
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='MDCHII')
!
    integer(kind=8) :: ednopt
    integer(kind=8) :: vali(2)
    parameter(ednopt=-1)
    integer(kind=8) :: ednono
    parameter(ednono=-1)
    integer(kind=8) :: mfloat
    parameter(mfloat=6)
    integer(kind=8) :: iterma
    parameter(iterma=1)
!
    integer(kind=8) :: ifm, nivinf, nbtv2, ncmp
    integer(kind=8) :: iaux, nbcham, lnochm, nbcmfi
    integer(kind=8) :: adncmp, aducmp, jaux, existc
    integer(kind=8) :: finupt, finuno, numsau
    integer(kind=8) :: adnume, adinst, nchmed, jnptno, jpasdt
    integer(kind=8) :: typech, nbcmp, jcmp, junit, nseqca, npro
    aster_logical :: ilocal
!
    character(len=8) :: saux08
    character(len=21) :: saux21
    character(len=64) :: saux64, nomcha, nomprf, nomloc, nomam2, nomamd
!
    real(kind=8) :: fiinst
!
    call jemarq()
    call infniv(ifm, nivinf)
!
!====
! 1. LE CHAMP EST-IL PRESENT ?
!====
!
! 1.1. ==> NBCHAM : NOMBRE DE CHAMPS DANS LE FICHIER
!
    call as_mfdnfd(idfimd, nbcham, codret)
    if (codret .ne. 0) then
        saux08 = 'mfdnfd'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
! 1.2. ==> RECHERCHE DU CHAMP VOULU
!
    lnochm = lxlgut(nochmd)
    existc = 0
    ilocal = .true.
    numsau = 0
!
    do iaux = 1, nbcham
!
! 1.2.1. ==> NBCMFI : NOMBRE DE COMPOSANTES DANS LE FICHIER POUR
!                     LE CHAMP NUMERO IAUX
!
        call as_mfdnfc(idfimd, iaux, nbcmfi, codret)
        if (codret .ne. 0) then
            saux08 = 'mfdnfc'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
! 1.2.2. ==> POUR LE CHAMP NUMERO IAUX, ON RECUPERE :
!            SAUX64 : NOM DU CHAMP
!            ZK16(ADNCMP) : NOM DE SES NBCMFI COMPOSANTES
!            ZK16(ADUCMP) : UNITES DE SES NBCMFI COMPOSANTES
!
        call codent(iaux, 'G', saux08)
        call wkvect('&&'//nompro//'N'//saux08, 'V V K16', nbcmfi, adncmp)
        call wkvect('&&'//nompro//'U'//saux08, 'V V K16', nbcmfi, aducmp)
        saux64 = ' '
        call as_mfdfdi(idfimd, iaux, saux64, jaux, zk16(adncmp), &
                       zk16(aducmp), nseqca, codret)
        if (codret .ne. 0 .or. jaux .ne. mfloat) then
            if (codret .ne. 0) then
                saux08 = 'mfdfdi'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
!         TYPE INCORRECT
            if (saux64 .eq. nochmd .and. jaux .ne. mfloat) then
                vali(1) = jaux
                call utmess('A+', 'MED_84', si=vali(1))
                call utmess('F', 'MED_75')
            end if
        end if
!
! 1.2.3. ==> COMPARAISON DU NOM DU CHAMP
!
        jaux = lxlgut(saux64)
!
        if (jaux .eq. lnochm) then
            if (saux64(1:jaux) .eq. nochmd(1:lnochm)) then
                numsau = iaux
                existc = 1
            end if
        end if
!
    end do
    if (numsau .eq. 0) then
        call utmess('F', 'MED_43', sk=nochmd(1:lnochm))
    end if
!
    if (existc .ne. 1) then
        call as_mfdnfd(idfimd, nchmed, codret)
        call utmess('F+', 'MED_57', sk=nochmd(1:lnochm))
        do iaux = 1, nchmed
            call as_mfdnfc(idfimd, iaux, nbcmp, codret)
            call wkvect('&&MDCHII.NOMCMP_K16', 'V V K16', nbcmp, jcmp)
            call wkvect('&&MDCHII.UNITCMP', 'V V K16', nbcmp, junit)
            call as_mfdfdi(idfimd, iaux, nomcha, typech, zk16(jcmp), &
                           zk16(junit), nseqca, codret)
            call utmess('F+', 'MED2_2', sk=nomcha)
            call jedetr('&&MDCHII.NOMCMP_K16')
            call jedetr('&&MDCHII.UNITCMP')
        end do
        call utmess('F', 'VIDE_1')
    end if
!
!====
! 2. NOMBRE DE TABLEAUX DE VALEURS ASSOCIES AU CHAMP
!====
!
    call as_mfdncn(idfimd, nochmd, ncmp, codret)
    call wkvect('&&MDCHII.CNAME', 'V V K16', ncmp, jcmp)
    call wkvect('&&MDCHII.CUNIT', 'V V K16', ncmp, junit)
!
    call as_mfdfin(idfimd, nochmd, nomam2, nbtv, zk16(junit), &
                   zk16(jcmp), codret)
!
    if (codret .ne. 0) then
        saux08 = 'mfdfin'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
    call jedetr('&&MDCHII.CNAME')
    call jedetr('&&MDCHII.CUNIT')
    if (nivinf .gt. 1) then
        write (ifm, 20001) nochmd, nbtv
    end if
20001 format&
& ('LE CHAMP MED ', a, ' CONTIENT ', i6, ' TABLEAUX SUR LES NOEUDS :')
!
!====
! 3. ALLOCATION DES TABLEAUX
!====
!
    if (nbtv .gt. 0) then
!
        call wkvect(prefix//'.NUME', 'V V I', 2*nbtv, adnume)
        call wkvect(prefix//'.INST', 'V V R', nbtv, adinst)
!
    end if
!
!====
! 4. POUR CHAQUE TABLEAU :
!====
!
    call codent(numsau, 'G', saux08)
    call jeexin('&&MDCHI2.'//saux08, codret)
    if (codret .eq. 0) then
        call wkvect('&&MDCHI2.'//saux08, 'V V I', 2*nbtv, jnptno)
        call wkvect('&&MDCHI3.'//saux08, 'V V R', nbtv, jpasdt)
        do iaux = 1, nbtv
!
! 4.1. ==> LECTURE
!    . NUMERO, UNITE ET VALEUR DU PAS DE TEMPS : FINUPT, FIINST
!    . NUMERO D'ORDRE : FINUNO
!
            call as_mfdcsi(idfimd, nochmd, iaux, finupt, finuno, &
                           fiinst, codret)
!
            if (codret .ne. 0) then
                saux08 = 'mfdcsi'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
            if (finupt .ne. finuno) then
                saux21 = 'readMedFileToResults'
                call utmess('F', 'MED_12', sk=saux21)
            end if
            zi(jnptno+iaux*2-2) = finupt
            zi(jnptno+iaux*2-1) = finuno
            zr(jpasdt+iaux-1) = fiinst
        end do
    else
        call jeveuo('&&MDCHI2.'//saux08, 'L', jnptno)
        call jeveuo('&&MDCHI3.'//saux08, 'L', jpasdt)
    end if
!
    nbtv2 = nbtv
    do 40, iaux = 1, nbtv
!
        finupt = zi(jnptno+iaux*2-2)
        finuno = zi(jnptno+iaux*2-1)
        fiinst = zr(jpasdt+iaux-1)
!
        if (.not. ilocal) then
            call utmess('F', 'MED_60')
        end if
!
!       ON REGARDE DANS LE FICHIER MED SI DES VALEURS SONT REFERENCEES
!       POUR TYPENT ET TYPGEO, SI CE N'EST PAS LE CAS, ON RETIRE
!       CET INSTANT DE LA LISTE
        call as_mfdonp(idfimd, nochmd, finupt, finuno, typent, &
                       typgeo, iterma, nomamd, nomprf, nomloc, &
                       npro, codret)
!
        if (npro .eq. 0) then
            nbtv2 = nbtv2-1
            goto 40
        end if
!
        if (nivinf .gt. 1) then
            if (finupt .eq. ednopt) then
                write (ifm, 40012)
            else
                write (ifm, 40022) finupt, fiinst
            end if
            if (finuno .eq. ednono) then
                write (ifm, 40013)
            else
                write (ifm, 40023) finuno
            end if
        end if
!
! 4.2. ==> ARCHIVAGE DANS LES TABLEAUX EXPLOITES PAR LE PROGRAMME
!          APPELANT
!
        zi(adnume+2*iaux-2) = finupt
        zi(adnume+2*iaux-1) = finuno
        if (finupt .ne. ednopt) then
            zr(adinst+iaux-1) = fiinst
        end if
!
40  end do
!
!     SOIT ON A TROUVE POUR TOUS LES INSTANTS DES VALEURS DANS
!     LE FICHIER MED, SOIT ON EN A TROUVE AUCUNE
!     LE CAS INTERMEDIAIRE EST PROBLEMATIQUE ET NON TRAITE
    ASSERT((nbtv .eq. nbtv2) .or. (nbtv2 .eq. 0))
    if (nbtv2 .eq. 0) then
        call jedetr(prefix//'.NUME')
        call jedetr(prefix//'.INST')
        nbtv = 0
    end if
!
40012 format(2x, '. AUCUNE INDICATION DE PAS DE TEMPS')
40022 format(2x, '. PAS DE TEMPS NUMERO ', i5,&
&       /, 2x, '. INSTANT : ', g13.5)
40013 format(2x, '. AUCUNE INDICATION DE NUMERO D''ORDRE',/)
40023 format(2x, '. NUMERO D''ORDRE : ', i5,/)
!
! --- MENAGE
    do 70, iaux = 1, nbcham
        call codent(iaux, 'G', saux08)
        call jedetr('&&'//nompro//'N'//saux08)
        call jedetr('&&'//nompro//'U'//saux08)
70  end do
!
    call jedema()
end subroutine
