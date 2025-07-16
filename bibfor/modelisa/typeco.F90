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

subroutine typeco(char, noma)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/cfnumm.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mmelin.h"
#include "asterfort/mmelty.h"
#include "asterfort/mminfi.h"
#include "asterfort/mmssfr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: noma, char
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES MAILLEES - LECTURE DONNEES)
!
! CONSTRUCTION DU TABLEAU POUR TYPE DE NOEUD/MAILLES
!
! ----------------------------------------------------------------------
!
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  CHAR   : NOM UTILISATEUR DU CONCEPT DE CHARGE
!
! ----------------------------------------------------------------------
!
    character(len=24) :: typeno, typema, maescl
    integer(kind=8) :: jtypno, jtypma, jmaesc
    integer(kind=8) :: ztypm, ztypn, zmaes
    integer(kind=8) :: nzoco, nnoco, nmaco, ntmae, iform
    integer(kind=8) :: izone
    integer(kind=8) :: jdecnm, jdecne, inoe, inom, nbnoe, nbnom
    integer(kind=8) :: jdecmm, jdecme, imae, imam, nbmae, nbmam
    integer(kind=8) :: posnom, posnoe
    integer(kind=8) :: posmam, posmae, nummae, nummam
    integer(kind=8) :: indmae, indmam
    integer(kind=8) :: ino, ima, posno, posma
    integer(kind=8) :: ndexfr(1), typint, nptm
    character(len=24) :: defico
    character(len=8) :: alias, nommae, nommam
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    defico = char(1:8)//'.CONTACT'
    nzoco = cfdisi(defico, 'NZOCO')
    nnoco = cfdisi(defico, 'NNOCO')
    nmaco = cfdisi(defico, 'NMACO')
    ntmae = cfdisi(defico, 'NTMAE')
    iform = cfdisi(defico, 'FORMULATION')
!
! --- ACCES OBJETS JEVEUX
!
    typeno = defico(1:16)//'.TYPENO'
    typema = defico(1:16)//'.TYPEMA'
    maescl = defico(1:16)//'.MAESCL'
    ztypm = cfmmvd('ZTYPM')
    ztypn = cfmmvd('ZTYPN')
    zmaes = cfmmvd('ZMAES')
!
! --- CREATION DES TABLEAUX
!
    call wkvect(typeno, 'G V I', ztypn*nnoco, jtypno)
    call wkvect(typema, 'G V I', ztypm*nmaco, jtypma)
    call wkvect(maescl, 'G V I', zmaes*ntmae, jmaesc)
!
! --- REMPLISSAGE DU TABLEAU TYPE_NOEUD
!
    do izone = 1, nzoco
        nbnoe = mminfi(defico, 'NBNOE', izone)
        nbnom = mminfi(defico, 'NBNOM', izone)
        jdecne = mminfi(defico, 'JDECNE', izone)
        jdecnm = mminfi(defico, 'JDECNM', izone)
!
        do inom = 1, nbnom
            posnom = jdecnm+inom
            zi(jtypno+ztypn*(posnom-1)+1-1) = 1
            zi(jtypno+ztypn*(posnom-1)+2-1) = izone
        end do
!
        do inoe = 1, nbnoe
            posnoe = jdecne+inoe
            zi(jtypno+ztypn*(posnoe-1)+1-1) = -1
            zi(jtypno+ztypn*(posnoe-1)+2-1) = izone
        end do
!
    end do
!
! --- REMPLISSAGE DU TABLEAU TYPE_MAILLE - TYPE MAITRE OU ESCLAVE
!
    indmae = 0
    indmam = 0
    do izone = 1, nzoco
        nbmae = mminfi(defico, 'NBMAE', izone)
        nbmam = mminfi(defico, 'NBMAM', izone)
        jdecme = mminfi(defico, 'JDECME', izone)
        jdecmm = mminfi(defico, 'JDECMM', izone)
!
        do imam = 1, nbmam
            posmam = jdecmm+imam
            indmam = indmam+1
            zi(jtypma+ztypm*(posmam-1)+1-1) = 1
            zi(jtypma+ztypm*(posmam-1)+2-1) = indmam
            if (iform .eq. 2) then
                call cfnumm(defico, posmam, nummam)
                call mmelty(noma, nummam, alias)
                if (alias .eq. 'PO1') then
                    nommam = int_to_char8(nummam)
                    call utmess('F', 'CONTACT3_2', sk=nommam)
                end if
            end if
        end do
!
        do imae = 1, nbmae
            posmae = jdecme+imae
            indmae = indmae+1
            zi(jtypma+ztypm*(posmae-1)+1-1) = -1
            zi(jtypma+ztypm*(posmae-1)+2-1) = indmae
            if (iform .eq. 2) then
                call cfnumm(defico, posmae, nummae)
                call mmelty(noma, nummae, alias)
                if (alias .eq. 'PO1') then
                    nommae = int_to_char8(nummae)
                    call utmess('F', 'CONTACT3_2', sk=nommae)
                end if
            end if
        end do
    end do
!
! --- REMPLISSAGE DU TABLEAU DES MAILLES ESCLAVES MAESC
!
    indmae = 0
    do izone = 1, nzoco
        nbmae = mminfi(defico, 'NBMAE', izone)
        jdecme = mminfi(defico, 'JDECME', izone)
!
        do imae = 1, nbmae
            posmae = jdecme+imae
            indmae = indmae+1
            call cfnumm(defico, posmae, nummae)
            if (iform .eq. 2) then
                typint = mminfi(defico, 'INTEGRATION', izone)
                call mmelin(noma, nummae, typint, nptm)
                call mmssfr(defico, izone, posmae, ndexfr(1))
            else
                nptm = 0
                ndexfr(1) = 0
            end if
            zi(jmaesc+zmaes*(indmae-1)+1-1) = posmae
            zi(jmaesc+zmaes*(indmae-1)+2-1) = izone
            zi(jmaesc+zmaes*(indmae-1)+3-1) = nptm
            zi(jmaesc+zmaes*(indmae-1)+4-1) = ndexfr(1)
        end do
    end do
!
! --- VERIFS: TYPENO ET TYPEMA SANS TROUS !
!
    do ino = 1, nnoco
        posno = ino
        if (zi(jtypno+ztypn*(posno-1)+1-1) .eq. 0) then
            ASSERT(.false.)
        end if
    end do
!
    do ima = 1, nmaco
        posma = ima
        if (zi(jtypma+ztypm*(posma-1)+1-1) .eq. 0) then
            ASSERT(.false.)
        end if
    end do
!
    call jedema()
end subroutine
