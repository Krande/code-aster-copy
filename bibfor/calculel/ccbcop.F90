! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine ccbcop(resultIn    , resultOut,&
                  listStoreJv , nbStore,&
                  listOptionJv, nbOption)
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/calcop.h"
#include "asterfort/ccfnrn.h"
#include "asterfort/ccvrch.h"
#include "asterfort/dismoi.h"
#include "asterfort/indk16.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/medom1.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rslesd.h"
#include "asterfort/rsnopa.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/deprecated_option.h"
!
integer, intent(in) :: nbStore, nbOption
character(len=8), intent(in) :: resultOut, resultIn
character(len=19), intent(in) :: listStoreJv, listOptionJv
!
! --------------------------------------------------------------------------------------------------
!
!  CALC_CHAMP - BOUCLE SUR LA LISTE D'OPTION ET APPEL A CALCOP
!
! --------------------------------------------------------------------------------------------------
!
!  ROUTINE PREPARANT L'APPEL A CALCOP
!
! IN  :
!   RESUIN K8   NOM DE LA SD IN
!   RESUOU K8   NOM DE LA SD OUT
!   LISORD K19  NOM DE LA LISTE DES NUMEROS D'ORDRE
!   NBORDR I    NOMBRE DE NUMEROS D'ORDRE
!   LISOOP K19  NOM DE LA LISTE DES OPTIONS A CALCULER
!   NBROPT I    LONGUEUR DE LA LISTE D'OPTIONS

!
! --------------------------------------------------------------------------------------------------
!
    integer :: iret, posopt
    integer :: ifm, niv, numeStore0
    integer :: nbac, nbpa, nbpara, jvPara
    integer :: iStore, iPara, iadou, iadin, numeStore, jopt, iOption
    integer, pointer :: listStore(:) => null()
    character(len=8) :: paraType, model, caraElem
    character(len=8) :: answer
    character(len=16) :: option, resultType
    character(len=24) :: nompar
    aster_logical :: exipla, newResult, lforc_noda
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
!
    call gettco(resultIn, resultType)
!
    call jeexin(resultOut//'           .DESC', iret)
    newResult = iret .eq. 0
!
    if ((resultIn.ne.resultOut) .and. (.not.newResult)) then
        call utmess('F', 'CALCULEL_18')
    endif
!
    call jeveuo(listStoreJv, 'L', vi = listStore)

! - Get parameters at initial state
    numeStore0 = listStore(1)
    call rslesd(resultIn, numeStore0, model_ = model, cara_elem_ = caraElem)
    if (model .eq. ' ') then
        call utmess('F', 'CALCULEL2_44')
    endif

! - New output result
    if (newResult) then
        call rscrsd('G', resultOut, resultType, nbStore)
        call titre()
    endif
!
!     ON VERIFIE QUE CARA_ELEM EST RENSEIGNES POUR LES COQUES
    exipla=.false.
    call dismoi('EXI_COQ3D', model, 'MODELE', repk=answer)
    if (answer(1:3) .eq. 'OUI') exipla=.true.
    call dismoi('EXI_PLAQUE', model, 'MODELE', repk=answer)
    if (answer(1:3) .eq. 'OUI') exipla=.true.
    if (exipla .and. caraElem .eq. ' ') then
        call utmess('A', 'CALCULEL2_94')
        goto 30
    endif
!
!     RECOPIE DES PARAMETRES DANS LA NOUVELLE SD RESULTAT
    if (newResult) then
        nompar='&&CCBCOP.NOMS_PARA '
        call rsnopa(resultIn, 2, nompar, nbac, nbpa)
        nbpara=nbac+nbpa
        call jeveuo(nompar, 'L', jvPara)
        do iStore = 1, nbStore
            numeStore = listStore(iStore)
            do iPara = 1, nbpara
                call rsadpa(resultIn, 'L', 1, zk16(jvPara+iPara-1), numeStore,&
                            1, sjv=iadin, styp=paraType, istop=0)
                call rsadpa(resultOut, 'E', 1, zk16(jvPara+iPara-1), numeStore,&
                            1, sjv=iadou, styp=paraType)
                if (paraType(1:1) .eq. 'I') then
                    zi(iadou)=zi(iadin)
                else if (paraType(1:1).eq.'R') then
                    zr(iadou)=zr(iadin)
                else if (paraType(1:1).eq.'C') then
                    zc(iadou)=zc(iadin)
                else if (paraType(1:3).eq.'K80') then
                    zk80(iadou)=zk80(iadin)
                else if (paraType(1:3).eq.'K32') then
                    zk32(iadou)=zk32(iadin)
                else if (paraType(1:3).eq.'K24') then
                    zk24(iadou)=zk24(iadin)
                else if (paraType(1:3).eq.'K16') then
                    zk16(iadou)=zk16(iadin)
                else if (paraType(1:2).eq.'K8') then
                    zk8(iadou)=zk8(iadin)
                endif
            end do
        end do
        call jedetr(nompar)
    endif
!
    call jeexin(listOptionJv, iret)
    if(iret .ne. 0) then
        call jeveuo(listOptionJv, 'L', jopt)
    end if
!
!     VERIFICATION DE LA PRESENCE D'UN EXCIT DANS LE FICHIER
!     DE COMMANDE OU DES CHARGES DANS LA SD RESULTAT
    posopt = indk16(zk16(jopt), 'FORC_NODA', 1, nbOption)
    lforc_noda = .true.
    if (posopt.eq.0) lforc_noda = .false.
    call ccvrch(resultIn, numeStore0, lforc_noda)
!
!     BOUCLE SUR LES OPTIONS DEMANDEES PAR L'UTILISATEUR
    do iOption = 1, nbOption
!
        option = zk16(jopt+iOption-1)
        call deprecated_option(option)

        if (option .eq. ' ') cycle
!
        if ((option.eq.'FORC_NODA') .or. (option.eq.'REAC_NODA')) then
            call ccfnrn(option, resultIn, resultOut, listStoreJv, nbStore,&
                        resultType)
        else
            call calcop(option, listOptionJv, resultIn, resultOut, listStoreJv,&
                        nbStore, resultType, iret, tldist=.True._1)
            ASSERT(iret .eq. 0)

        endif

    end do
!
 30 continue
!
    call jedema()
!
end subroutine
