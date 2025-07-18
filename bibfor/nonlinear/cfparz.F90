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

subroutine cfparz(ds_contact, iliai, coefff, coefpn, coefpt, &
                  coefte, dissup, izone, ip, numnoe, &
                  posnoe)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_Contact), intent(in) :: ds_contact
    real(kind=8) :: coefff, coefpn, coefpt, coefte, dissup
    integer(kind=8) :: iliai, ip, izone, numnoe, posnoe
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES DISCRETES - APPARIEMENT)
!
! CARACTERISTIQUES DES LIAISONS POUR LA ZONE
!
! ----------------------------------------------------------------------
!
! IN  ILIAI  : INDICE DE LA LIAISON COURANTE
! In  ds_contact       : datastructure for contact management
! IN  DISSUP : JEU FICTIF DE LA ZONE
! IN  COEFPN : COEFFICIENT DE PENALISATION DE CONTACT
! IN  COEFPT : COEFFICIENT DE PENALISATION DE FROTTEMENT
! IN  COEFFF : COEFFICIENT DE FROTTEMENT
! IN  COEFTE : COEFFICIENT THETA POUR FROTTEMENT
! IN  POSNOE : INDICES DANS CONTNO DU NOEUD ESCLAVE
! IN  NUMNOE : NUMERO ABSOLU DU NOEUD ESCLAVE
! IN  IP     : INDICE DU POINT DANS LA SD APPARIEMENT
! IN  IZONE  : NUMERO DE LA ZONE DE CONTACT
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: ztacf
    character(len=24) :: tacfin, jeusup
    integer(kind=8) :: jtacf, jjsup
    character(len=24) :: jeuite, numlia
    integer(kind=8) :: jjeuit, jnumli
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    jeuite = ds_contact%sdcont_solv(1:14)//'.JEUITE'
    tacfin = ds_contact%sdcont_solv(1:14)//'.TACFIN'
    jeusup = ds_contact%sdcont_solv(1:14)//'.JSUPCO'
    numlia = ds_contact%sdcont_solv(1:14)//'.NUMLIA'
!
    call jeveuo(jeuite, 'E', jjeuit)
    call jeveuo(tacfin, 'E', jtacf)
    call jeveuo(jeusup, 'E', jjsup)
    call jeveuo(numlia, 'E', jnumli)
!
    ztacf = cfmmvd('ZTACF')
!
! --- JEU FICTIF DE LA ZONE
!
    zr(jjsup+iliai-1) = dissup
!
! --- ADDITION DU JEU FICTIF DE LA ZONE
!
    zr(jjeuit+3*(iliai-1)+1-1) = zr(jjeuit+3*(iliai-1)+1-1)-dissup
!
! --- RECOPIE CARACTERISTIQUES ZONE -> NOEUDS ESCLAVES
!
    zr(jtacf+ztacf*(iliai-1)+0) = coefff
    zr(jtacf+ztacf*(iliai-1)+1) = coefpn
    zr(jtacf+ztacf*(iliai-1)+2) = coefpt
    zr(jtacf+ztacf*(iliai-1)+3) = coefte
!
! --- SAUVEGARDE DANS LA SD APPARIEMENT
!
    zi(jnumli+4*(iliai-1)+1-1) = ip
    zi(jnumli+4*(iliai-1)+2-1) = posnoe
    zi(jnumli+4*(iliai-1)+3-1) = numnoe
    zi(jnumli+4*(iliai-1)+4-1) = izone
!
    call jedema()
end subroutine
