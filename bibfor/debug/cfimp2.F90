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

subroutine cfimp2(sdcont_defi, sdcont_solv, mesh, iliai, typeou)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/apinfi.h"
#include "asterfort/apnomp.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisr.h"
#include "asterfort/cfnoap.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8) :: mesh
    character(len=24) :: sdcont_defi, sdcont_solv
    integer(kind=8) :: iliai
    character(len=3) :: typeou
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE DISCRETE - APPARIEMENT - UTILITAIRE)
!
! IMPRESSION DE L'ACTIVATION/DESACTIVATION DE LA LIAISON ESCLAVE/MAITRE
!
! ----------------------------------------------------------------------
!
!
! IN  DEFICO : SD DE DEFINITION DU CONTACT
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  NOMA   : NOM DU MAILLAGE
! IN  ILIAI  : NUMERO DE LA LIAISON (INDICE DANS LE TABLEAU GLOBAL DE
!              TOUTE LES LIAISONS POSSIBLES -APPARIEES-)
! IN  TYPEOU : LIEU OU L'OPERATION A ETE FAITE
!                'ACT' : ACTIVATION LIAISON DE CONTACT
!                'LIB' : DESACTIVATION LIAISON DE CONTACT
!                'NEG' : SUPPRESSION D'UNE LIAISON A PRESSION NEGATIVE
!                'GLI' : SUPPRESSION D'UNE LIAISON GLISSANTE
!                'ADH' : AJOUT D'UNE LIAISON ADHERENTE
!                'PIV' : SUPPRESSION D'UNE LIAISON A PIVOT NUL
!                'ALJ' : ALARME LORSQU'UN JEU DEPASSE LA VALEUR SEUIL
!                         DANS LE CAS DU CONTACT GLISSIERE
!                'SIN' : AFFICHAGE DE LA LIAISON PROVOQUANT UNE MATRICE
!                         DE CONTACT SINGULIERE
!                'AGC' : ALARME LORSQUE LE NOMBRE D'ITERATIONS MAX
!                         DU GCP EST DEPASSE
!
!
!
!
    integer(kind=8) :: ifm, niv
    character(len=24) :: numlia
    integer(kind=8) :: jnumli
    integer(kind=8) :: ip, izone
    integer(kind=8) :: entapp
    character(len=8) :: nomapp
    character(len=16) :: nomnoe
    integer(kind=8) :: typapp
    character(len=16) :: etalia
    character(len=4) :: type2
    character(len=19) :: sdappa
    character(len=38) :: nomlia
    real(kind=8) :: aljeu
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
!
! --- ACCES SD CONTACT
!
    numlia = sdcont_solv(1:14)//'.NUMLIA'
    call jeveuo(numlia, 'L', jnumli)
!
! --- INITIALISATIONS
!
    aljeu = cfdisr(sdcont_defi, 'ALARME_JEU')
    nomlia = ' '
    etalia = ' '
!
! --- SD APPARIEMENT
!
    sdappa = sdcont_solv(1:14)//'.APPA'
!
    if (niv .ge. 2) then
!
! ----- POINT DE CONTACT
!
        ip = zi(jnumli+4*(iliai-1)+1-1)
!
! ----- NOM DU NOEUD ESCLAVE
!
        call apnomp(sdappa, ip, nomnoe)
!
! ----- INFOS APPARIEMENT
!
        call apinfi(sdappa, 'APPARI_TYPE', ip, typapp)
        call apinfi(sdappa, 'APPARI_ENTITE', ip, entapp)
        call apinfi(sdappa, 'APPARI_ZONE', ip, izone)
!
! ----- NOM ET TYPE DU MAITRE
!
        call cfnoap(mesh, sdcont_defi, typapp, entapp, nomapp, &
                    type2)
!
! ----- NOM DE LA LIAISON
!
        write (nomlia, 100) iliai, '( ', nomnoe, type2, nomapp, '): '
!
! ----- AFFICHAGE LIAISON
!
        if (typeou .eq. 'ACT') then
            etalia = ' CONT - AJOUTE  '
            write (ifm, 202) nomlia, etalia
!
        else if (typeou .eq. 'LIB') then
            etalia = ' CONT - SUPPRIME'
            write (ifm, 202) nomlia, etalia
!
        else if (typeou .eq. 'NEG') then
            etalia = ' PRES. NEGATIVE '
            write (ifm, 200) nomlia, etalia, ' - TYPE CONT.'
!
        else if (typeou .eq. 'PIV') then
            etalia = ' PIVOT NUL      '
            write (ifm, 200) nomlia, etalia, ' - TYPE CONT.'
!
        else if (typeou .eq. 'GLI') then
            etalia = ' GLIS - SUPPRIME'
            write (ifm, 200) nomlia, etalia, ' - TYPE CONT.'
!
        else if (typeou .eq. 'ADH') then
            etalia = ' ADHE - AJOUTE  '
            write (ifm, 200) nomlia, etalia, ' - TYPE CONT.'
!
        else if (typeou .eq. 'ALJ') then
            etalia = ' DECOLLE DU JEU '
            write (ifm, 201) nomlia, etalia, aljeu
!
        else if (typeou .eq. 'AGC') then
            etalia = ' INTERPENETRE   '
            write (ifm, 202) nomlia, etalia
!
        else
            ASSERT(.false.)
        end if
    end if
!
100 format(i5, a2, a16, a4, a8, a3)
200 format(' <CONTACT><CALC> LIAISON ', a38, a16, a20)
201 format(' <CONTACT><CALC> LIAISON ', a38, a16, e12.3)
202 format(' <CONTACT><CALC> LIAISON ', a38, a16)
!
    call jedema()
!
end subroutine
