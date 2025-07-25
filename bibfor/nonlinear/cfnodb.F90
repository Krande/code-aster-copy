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
subroutine cfnodb(sdcont)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/utlisi.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=8), intent(in) :: sdcont
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES MAILLEES - LECTURE DONNEES - ELIMINATION)
!
! DETECTION DE NOEUDS APPARTENANT A DEUX SURFACES DE CONTACT
!
! DEUX CAS SONT DETECTES : AU SEIN D'UNE MEME ZONE
!                          ENTRE 2 SURFACES ESCLAVES
!
! ----------------------------------------------------------------------
!
! IN  CHAR   : NOM UTILISATEUR DU CONCEPT DE CHARGE
!
!
!
!
    character(len=24) :: defico
    aster_logical :: lcalc, lliss
    integer(kind=8) :: nzoco, nnoco, iform
    character(len=24) :: nodbl, nodbl2
    integer(kind=8) :: jnodbl, jnodb2
    character(len=24) :: contno, sansno, psans
    integer(kind=8) :: jnoco, jsans, jpsans
    integer(kind=8) :: izone, ibid(1), ndoubl, nvdbl
    integer(kind=8) :: izonea, izoneb, nvdba, nvdbb
    integer(kind=8) :: nbnoe, jdecne, nbnoea, nbnoeb
    integer(kind=8) :: nbnom, jdecnm, jdecea, jdeceb
    integer(kind=8) :: nsans, jdecs, nsansa, jdecsa, nsansb, jdecsb
    integer(kind=8) :: vali(3)
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    defico = sdcont(1:8)//'.CONTACT'
    nzoco = cfdisi(defico, 'NZOCO')
    nnoco = cfdisi(defico, 'NNOCO')
    iform = cfdisi(defico, 'FORMULATION')
!
! --- OBJETS TEMPORAIRES
!
    nodbl = '&&CFNODB.NODBL'
    call wkvect(nodbl, 'V V I', nnoco, jnodbl)
    nodbl2 = '&&CFNODB.NODBL2'
    call wkvect(nodbl2, 'V V I', nnoco, jnodb2)
!
! --- ACCES AU TABLEAU DES NOEUDS DE CONTACT
!
!
    contno = defico(1:16)//'.NOEUCO'
    call jeveuo(contno, 'L', jnoco)
    if (iform .ne. 5) then
        sansno = defico(1:16)//'.SSNOCO'
        call jeveuo(sansno, 'L', jsans)
        psans = defico(1:16)//'.PSSNOCO'
        call jeveuo(psans, 'L', jpsans)
    end if
!
! ----------------------------------------------------------------------
!
! --- PREMIER CAS : NOEUDS COMMUNS DANS UNE MEME ZONE DE CONTACT
!
    if (iform .ne. 5) then
        do izone = 1, nzoco
            nbnoe = mminfi(defico, 'NBNOE', izone)
            nbnom = mminfi(defico, 'NBNOM', izone)
            jdecne = mminfi(defico, 'JDECNE', izone)
            jdecnm = mminfi(defico, 'JDECNM', izone)
            lcalc = mminfl(defico, 'CALCUL', izone)
            if (.not. lcalc) then
                goto 100
            end if
            call utlisi('INTER', zi(jnoco+jdecne), nbnoe, zi(jnoco+jdecnm), nbnom, &
                        zi(jnodbl), nnoco, ndoubl)
            if (ndoubl .ne. 0) then
                if (ndoubl .gt. 0) then
! --------- LES NOEUDS COMMUNS SONT-ILS EXCLUS PAR SANS_NOEUD ?
                    nsans = zi(jpsans+izone)-zi(jpsans+izone-1)
                    jdecs = zi(jpsans+izone-1)
                    call utlisi('DIFFE', zi(jnodbl), ndoubl, zi(jsans+jdecs), nsans, &
                                ibid, 1, nvdbl)
! --------- NON !
                    if (nvdbl .ne. 0) then
                        vali(1) = izone
                        vali(2) = abs(nvdbl)
                        call utmess('F', 'CONTACT2_13', ni=2, vali=vali)
                    end if
                else
                    ASSERT(.false.)
                end if
            end if
100         continue
        end do
    end if
!
! ----------------------------------------------------------------------
!
! --- SECOND CAS : NOEUDS COMMUNS A DEUX SURFACES ESCLAVES
!
    do izonea = 1, nzoco
        lcalc = mminfl(defico, 'CALCUL', izonea)
        if (.not. lcalc) then
            goto 200
        end if
        do izoneb = izonea+1, nzoco
            lcalc = mminfl(defico, 'CALCUL', izoneb)
            if (.not. lcalc) then
                goto 201
            end if
            nbnoea = mminfi(defico, 'NBNOE', izonea)
            nbnoeb = mminfi(defico, 'NBNOE', izoneb)
            jdecea = mminfi(defico, 'JDECNE', izonea)
            jdeceb = mminfi(defico, 'JDECNE', izoneb)
            call utlisi('INTER', zi(jnoco+jdecea), nbnoea, zi(jnoco+jdeceb), nbnoeb, &
                        zi(jnodbl), nnoco, ndoubl)
            if (ndoubl .ne. 0) then
                if (ndoubl .gt. 0) then
                    if (iform .eq. 1) then
! ------------- LES NOEUDS COMMUNS SONT-ILS EXCLUS PAR LA ZONE A ?
                        nsansa = zi(jpsans+izonea)-zi(jpsans+izonea-1)
                        jdecsa = zi(jpsans+izonea-1)
                        call utlisi('DIFFE', zi(jnodbl), ndoubl, zi(jsans+jdecsa), nsansa, &
                                    zi(jnodb2), nnoco, nvdba)
                        if (nvdba .ne. 0) then
                            if (nvdba .gt. 0) then
! ----------------- LES NOEUDS RESTANTS SONT-ILS EXCLUS PAR LA ZONE B ?
                                nsansb = zi(jpsans+izoneb)-zi(jpsans+izoneb-1)
                                jdecsb = zi(jpsans+izoneb-1)
                                call utlisi('DIFFE', zi(jnodb2), nvdba, zi(jsans+jdecsb), nsansb, &
                                            ibid, 1, nvdbb)
                                if (nvdbb .ne. 0) then
                                    vali(1) = izonea
                                    vali(2) = izoneb
                                    vali(3) = abs(nvdbb)
                                    call utmess('A', 'CONTACT2_15', ni=3, vali=vali)
                                end if
                            else
                                ASSERT(.false.)
                            end if
                        end if
                    else if (iform .eq. 2 .or. iform .eq. 5) then
                        lliss = cfdisl(defico, 'LISSAGE')
                        vali(1) = izonea
                        vali(2) = izoneb
                        vali(3) = abs(ndoubl)
                        if ((iform .eq. 5 .and. lliss) .or. (iform .eq. 2)) then
                            call utmess('F', 'CONTACT2_16', ni=3, vali=vali)
                        end if
                    else
                        ASSERT(.false.)
                    end if
                else
                    ASSERT(.false.)
                end if
            end if
201         continue
        end do
200     continue
    end do
!
! --- MENAGE
!
    call jedetr(nodbl)
    call jedetr(nodbl2)
!
    call jedema()
end subroutine
