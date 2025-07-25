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
subroutine iremed(fileUnit, dsNameZ, lResu, &
                  fieldListNb, fieldListType, fieldMedListType, &
                  storeListNb, storeListIndx, &
                  paraListNb, paraListName, &
                  cmpListNb, cmpListName, &
                  cellUserNb, cellUserNume, &
                  nodeUserNb, nodeUserNume, &
                  cplxFormat, lVariName, caraElem, &
                  lfichUniq, lNomCas)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/asmpi_barrier_wrap.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_split_comm.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/irchme.h"
#include "asterfort/irmpar.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mdnoch.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutrg.h"
#include "asterfort/utcmp3.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/bool_to_int.h"
!
    integer(kind=8), intent(in) :: fileUnit
    character(len=19), intent(in) :: dsNameZ
    aster_logical, intent(in) :: lResu
    integer(kind=8), intent(in) :: fieldListNb
    character(len=16), pointer :: fieldListType(:)
    character(len=80), pointer :: fieldMedListType(:)
    integer(kind=8), intent(in) :: storeListNb
    integer(kind=8), pointer :: storeListIndx(:)
    integer(kind=8), intent(in) :: paraListNb
    character(len=16), pointer :: paraListName(:)
    integer(kind=8), intent(in) :: cmpListNb
    character(len=8), pointer :: cmpListName(:)
    integer(kind=8), intent(in) :: cellUserNb
    integer(kind=8), pointer :: cellUserNume(:)
    integer(kind=8), intent(in) :: nodeUserNb
    integer(kind=8), pointer :: nodeUserNume(:)
    character(len=*), intent(in) ::  cplxFormat
    aster_logical, intent(in) :: lVariName
    character(len=8), intent(in) :: caraElem
    aster_logical, intent(in) :: lfichUniq, lNomCas
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE D'UN CONCEPT SUR FICHIER MED
!
! --------------------------------------------------------------------------------------------------
!
! IN  NOMCON : K8  : NOM DU CONCEPT A IMPRIMER
! IN  IFICHI : IS  : UNITE LOGIQUE D'ECRITURE
! IN  NOCHAM : K*  : NOM DES CHAMPS A IMPRIMER ( EX 'DEPL', .... )
! IN  NOVCMP : K*  : NOM JEVEUX DE L'OBJET CONTENANT LES NOMS MED
! IN  PARTIE : K4  : IMPRESSION DE LA PARTIE COMPLEXE OU REELLE DU CHAMP
! IN  LIORDR : K*  : LISTE DES NUMEROS D'ORDRE A IMPRIMER
! IN  LRESU  : L   : INDIQUE SI NOMCON EST UN CHAMP OU UN RESULTAT
! IN  NBNOEC : I   : NOMBRE DE NOEUDS A IMPRIMER
! IN  LINOEC : I   : NUMEROS DES NOEUDS A IMPRIMER
! IN  NBMAEC : I   : NOMBRE DE MAILLES A IMPRIMER
! IN  LIMAEC : I   : NUMEROS DES MAILLES A IMPRIMER
! IN  NOMCMP : K*  : NOMS DES COMPOSANTES A IMPRIMER
! IN  CARAEL : K*  : NOM DU CARA_ELEM
! IN  NOPARA : K16 : NOM D'UN PARAMATRE A AJOUTER
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6) :: chnumo, tsca
    character(len=8) :: typech, nomgd, dsName, sdcarm, carel2, valk2(2), saux08
    character(len=16) :: fieldType, fieldResult
    character(len=19) :: fieldName, cesnsp, cescoq, cesfib, cesori, cestuy
    character(len=24) :: valk(2), resultType
    character(len=64) :: fieldNameMed
    integer(kind=8) :: storeIndx, iField, iStore, iret, ibid, codret, iaux
    integer(kind=8) :: lnochm, ierd, nbCmpDyna, minval, nbNomCas

    aster_logical :: lfirst, l_mult_model, l_vari_name, l_meta_name
    integer(kind=8), pointer :: cmpListNume(:) => null()
    character(len=24), pointer :: celk(:) => null()
    mpi_int :: world, newcom, color, key, ierror
    parameter(cesnsp='&&IREMED.CANBSP')
    parameter(cescoq='&&IREMED.CARCOQUE')
    parameter(cesfib='&&IREMED.CAFIBR')
    parameter(cesori='&&IREMED.CARORIEN')
    parameter(cestuy='&&IREMED.CARGEOPO')
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    dsName = dsNameZ
    lfirst = .true.
    nbCmpDyna = 0
!
    sdcarm = ' '
    if (caraElem .ne. ' ' .and. caraElem .ne. '-') then
        if (lResu) then
            call dismoi('CARA_ELEM', dsName, 'RESULTAT', repk=carel2, arret='C', &
                        ier=ierd)
            if (carel2 .ne. caraElem) then
                valk2(1) = caraElem
                valk2(2) = carel2
                call utmess('F', 'MED_64', sk=valk2(1))
            end if
        end if
        sdcarm = '&&IREMED'
        call exisd('CHAMP', caraElem//'.CANBSP', iret)
        if (iret .ne. 0) then
            call celces(caraElem//'.CANBSP', 'V', cesnsp)
        end if
        call exisd('CARTE', caraElem//'.CARORIEN', iret)
        if (iret .ne. 0) then
            call carces(caraElem//'.CARORIEN', 'ELEM', ' ', 'V', cesori, &
                        'A', iret)
        end if
        call exisd('CARTE', caraElem//'.CARCOQUE', iret)
        if (iret .ne. 0) then
            call carces(caraElem//'.CARCOQUE', 'ELEM', ' ', 'V', cescoq, &
                        'A', iret)
        end if
        call exisd('CARTE', caraElem//'.CARGEOPO', iret)
        if (iret .ne. 0) then
            call carces(caraElem//'.CARGEOPO', 'ELEM', ' ', 'V', cestuy, &
                        'A', iret)
        end if
        call exisd('CHAM_ELEM', caraElem//'.CAFIBR', iret)
        if (iret .ne. 0) then
            call celces(caraElem//'.CAFIBR', 'V', cesfib)
        end if
    end if
!
! - Create parameters in MED file
!
    if (paraListNb .gt. 0) then
        call irmpar(dsName, fileUnit, paraListNb, paraListName)
    end if
!
! - Loop on fields
!
    do iField = 1, fieldListNb
        fieldType = fieldListType(iField)
        l_mult_model = ASTER_FALSE
        l_vari_name = ASTER_FALSE
        l_meta_name = ASTER_FALSE
! ----- Loop on storing slots
        do iStore = 1, storeListNb
            storeIndx = storeListIndx(iStore)
!
!         --- SI VARIABLE DE TYPE RESULTAT = RESULTAT COMPOSE :
!             VERIFICATION CORRESPONDANCE ENTRE NUMERO D'ORDRE
!             UTILISATEUR ORDR(IORDR) ET NUMERO DE RANGEMENT IRET
! AU CAS OU ON NE PASSE PAS EN DESSOUS ON INITIALISE LORDR A FALSE
            if (lResu) then
                call rsutrg(dsName, storeIndx, iret, ibid)
                if (iret .eq. 0) then
!             - MESSAGE NUMERO D'ORDRE NON LICITE
                    if (lfirst) then
                        call codent(storeIndx, 'G', chnumo)
                        call utmess('A', 'PREPOST2_46', sk=chnumo)
                    end if
                    goto 22
                end if
            end if
! --------- Get name of field
            if (lResu) then
                call rsexch(' ', dsName, fieldType, storeIndx, fieldName, iret)
                if (lfichUniq) then
                    minval = iret
                    call asmpi_comm_vect('MPI_MIN', 'I', sci=minval)
                    if (minval .eq. 0) then
                        call asmpi_comm('GET', world)
                        if (iret .ne. 0) then
                            color = 1
                            key = 1
                            call asmpi_split_comm(world, color, key, "TMPCOMM", newcom)
                            call asmpi_comm('SET', newcom)
                            goto 21
                        else
                            color = 0
                            key = 0
                            call asmpi_split_comm(world, color, key, "TMPCOMM", newcom)
                            call asmpi_comm('SET', newcom)
                        end if
                    else
                        newcom = 0
                        goto 21
                    end if
                else
                    if (iret .ne. 0) goto 21
                end if
            else
                fieldName = dsNameZ
                newcom = 0
            end if
!
!         * IMPRESSION DU CHAMP (CHAM_NO OU CHAM_ELEM)
!             LE CHAMP EST UN CHAM_GD SIMPLE SI LRESU=.FALSE. OU
!             LE CHAMP EST LE CHAM_GD CHAM(ISY) DE NUMERO D'ORDRE
!             ORDR(IORDR) ISSU DE LA SD_RESULTAT NOMCON
!
!       --- TYPE DU CHAMP A IMPRIMER (CHAM_NO OU CHAM_ELEM)
            call dismoi('TYPE_CHAMP', fieldName, 'CHAMP', repk=typech)
            call dismoi('TYPE_RESU', dsName, 'RESULTAT', repk=resultType)
!
            if ((typech(1:4) .eq. 'NOEU') .or. (typech(1:2) .eq. 'EL')) then
            else if (typech(1:4) .eq. 'CART') then
!            GOTO 9999
            else
                valk(1) = typech
                valk(2) = fieldName
                if (resultType(1:9) .eq. 'MODE_GENE' .or. resultType(1:9) .eq. 'HARM_GENE') then
                    call utmess('A+', 'PREPOST_87', nk=2, valk=valk)
                    call utmess('A', 'PREPOST6_36')
                else
                    call utmess('A', 'PREPOST_87', nk=2, valk=valk)
                end if
            end if
!
!         --- NOM DE LA GRANDEUR ASSOCIEE AU CHAMP CHAM19
            call dismoi('NOM_GD', fieldName, 'CHAMP', repk=nomgd)
            call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)

            if ((typech(1:4) .eq. 'CART' .and. tsca .ne. 'R')) then
                if (.not. lResu) then
                    valk(1) = tsca
                    call utmess('A+', 'PREPOST_91', nk=1, valk=valk)
                end if
                goto 999
            end if
!
            if (cmpListNb .ne. 0) then
                if ((nomgd .eq. 'VARI_R') .and. (typech(1:2) .eq. 'EL')) then
                    AS_ALLOCATE(vi=cmpListNume, size=cmpListNb)
                    call utcmp3(cmpListNb, cmpListName, cmpListNume)
                    AS_DEALLOCATE(vi=cmpListNume)
                end if
            end if
! --------- Get name of field in MED file
            fieldNameMed = ' '
            if (.not. associated(fieldMedListType)) then
                call mdnoch(fieldNameMed, lnochm, bool_to_int(lResu), dsName, fieldType, codret)
            else
                fieldNameMed = fieldMedListType(iField) (1:64)
            end if
!
! --------- Change name and nume_ordre for mult_elas sd
            if (lNomCas) then
                call rsadpa(dsName, 'L', 1, "NOM_CAS", storeIndx, 0, styp=saux08, sjv=iaux, istop=0)

                fieldResult = zk16(iaux)
                nbNomCas = len_trim(fieldResult(1:8))
                fieldNameMed(1:8) = repeat('_', 8)
                fieldNameMed(1:nbNomCas) = fieldResult(1:nbNomCas)
                storeIndx = 1
            end if

!
!         -- TRAITEMENT SPECIFIQUE POUR LES CHAMPS ISSUE DE PROJ_CHAMP
!            METHODE='SOUS_POINT'
            if (typech .eq. 'ELGA') then
                call jeveuo(fieldName//'.CELK', 'L', vk24=celk)
                if (celk(2) .eq. 'INI_SP_MATER') then
                    call utmess('A', 'MED2_9', sk=fieldType)
                    codret = 0
                    goto 999
                end if
            end if
!
!         -- ON LANCE L'IMPRESSION:
!         -------------------------
!
            if (.not. lResu) dsName = ' '
            call irchme(fileUnit, fieldName, cplxFormat, fieldNameMed, dsName, &
                        fieldType, typech, storeIndx, cmpListNb, cmpListName, &
                        nodeUserNb, nodeUserNume, &
                        cellUserNb, cellUserNume, lVariName, &
                        sdcarm, caraElem, paraListNb, paraListName, &
                        nbCmpDyna, lfichUniq, codret)
!
999         continue
!
            if (codret .ne. 0 .and. codret .ne. 100 .and. &
                codret .ne. 200 .and. codret .ne. 300) then
                valk(1) = fieldType
                valk(2) = 'MED'
                if (fieldType .ne. 'COMPORTEMENT') then
                    call utmess('A', 'PREPOST_90', nk=2, valk=valk)
                end if
            end if
            if (codret .eq. 100) then
                l_mult_model = ASTER_TRUE
            end if
            if (codret .eq. 200) then
                l_vari_name = ASTER_TRUE
            end if
            if (codret .eq. 400) then
                l_meta_name = ASTER_TRUE
            end if
!
22          continue
!
21          continue
            if (lfichUniq .and. newcom .ne. 0) then
                call asmpi_barrier_wrap(world, ierror)
                call asmpi_comm('SET', world)
                if (newcom .ne. 0) then
                    call asmpi_comm('FREE', newcom)
                end if
            end if
        end do
        if (l_mult_model) then
            valk(1) = fieldType
            call utmess('I', 'MED_30', sk=valk(1))
        end if
        if (l_vari_name) then
            call utmess('A', 'MED2_7')
        end if
        if (l_meta_name) then
            call utmess('A', 'MED2_10')
        end if
        lfirst = .false.
!
    end do
!
    call detrsd('CHAM_ELEM_S', cesnsp)
    call detrsd('CHAM_ELEM_S', cescoq)
    call detrsd('CHAM_ELEM_S', cesfib)
    call detrsd('CHAM_ELEM_S', cesori)
    call detrsd('CHAM_ELEM_S', cestuy)
!
    call jedema()
end subroutine
