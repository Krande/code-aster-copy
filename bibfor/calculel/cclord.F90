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
subroutine cclord(numeOptEff, nbStore, listStore, jvBaseName, isOptionFromUser, &
                  numeStoreMin, numeStoreMax, resultIn, resultOut)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: numeOptEff
    integer(kind=8), intent(in) :: nbStore
    integer(kind=8), pointer :: listStore(:)
    character(len=8), intent(in) :: jvBaseName
    aster_logical, intent(in) :: isOptionFromUser
    integer(kind=8), intent(in) ::  numeStoreMin, numeStoreMax
    character(len=8), intent(in) :: resultIn, resultOut
!
! --------------------------------------------------------------------------------------------------
!
!  CALC_CHAMP
!
!  Create list of storing index to compute given option
!
! --------------------------------------------------------------------------------------------------
!
!  CREATION D'UNE LISTE DE NUMEROS D'ORDRE POUR L'OPTION EN ARGUMENT
!  EN PRENANT EN COMPTE SA LISTE DE DEPENDANCE
!
! In  numeOptEff        : index of option in list of option
! In  nbStore           : number of storing indexes
! In  listStore         : list of storing indexes
! In  jvBaseName        : string base for names of JEVEUX objects
! In  isOptionFromUser  : option is required by user (not in dependencies)
! In  numeStoreMin      : minimum storing index in list of storing indexes to compute option
! In  numeStoreMax      : maximum storing index in list of storing indexes to compute option
! In  resultIn          : name of datastructure for input results
! In  resultOut         : name of datastructure for output results
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jlisop, jliori, jlidep, jordop, ierd, inddeb, indfin
    integer(kind=8) :: nbStoreToCompute, iStore, curmax, curmin, iter, decal, numeStore, jlnoin
    integer(kind=8) :: jordo2, jlisde
    character(len=24) :: listStoreOptJv
    character(len=1) :: isodep
    character(len=5) :: numeOptStr
    character(len=16) :: option
    character(len=19) :: fieldName
    character(len=24) :: noliop, nolori, noldep, noliin, nolisd
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Generate name of object for list of storing index to compute for this option
    call codent(numeOptEff, 'D0', numeOptStr)
    listStoreOptJv = jvBaseName//'.OP'//numeOptStr

! - Access to objects
    noliop = jvBaseName//'.LISOPT'
    nolori = jvBaseName//'.LISORI'
    noldep = jvBaseName//'.LISDEP'
    noliin = jvBaseName//'.LNOINS'
    nolisd = jvBaseName//'.ISODEP'
    call jeveuo(noliop, 'L', jlisop)
    call jeveuo(nolori, 'L', jliori)
    call jeveuo(noldep, 'L', jlidep)
    call jeveuo(noliin, 'E', jlnoin)
    call jeveuo(nolisd, 'L', jlisde)

! - Parameters for current option
    option = zk16(jlisop+numeOptEff-1)
    inddeb = zi(jliori+2*numeOptEff-2)
    indfin = zi(jliori+2*numeOptEff-1)
    isodep = zk8(jlisde+numeOptEff-1) (1:1)

! - Create list of storing index to compute given option
    call wkvect(listStoreOptJv, 'V V I ', nbStore+3, jordop)
!
    if (inddeb .ne. 0) then
!       CAS 1 : CETTE OPTION DEPEND D'AUTRES OPTIONS A CALCULER
!               AUQUEL CAS, IL FAUT REGARDER COMMENT ELLE EN DEPEND
!               ET LA LISTE DES NUMEROS D'ORDRE DE SES PARENTS
        call jeveuo(noliin, 'E', jlnoin)
        curmax = numeStoreMax
        curmin = numeStoreMin
        do iter = inddeb, indfin
            call jeveuo(zk24(jlnoin+iter-1), 'L', jordo2)
            if (zk8(jlidep+iter-1) .eq. 'NP1') then
                decal = -1
            else if (zk8(jlidep+iter-1) .eq. 'NM1') then
                decal = +1
            else
                decal = 0
            end if
!         LA LISTE DE NUMEROS D'ORDRE PROVIENT DE OP0058
!         ELLE EST DONC CROISSANTE
            curmax = min(curmax, zi(jordo2+2)+decal)
            curmin = max(curmin, zi(jordo2+1)+decal)
        end do
!
        nbStoreToCompute = 0
        do iStore = 1, nbStore
            numeStore = listStore(iStore)
            if ((isodep .eq. '-') .and. (numeStore .eq. numeStoreMin)) then
                goto 30
            else if ((isodep .eq. '+') .and. (numeStore .eq. numeStoreMax)) then
                goto 30
            end if
            if (numeStore .ge. curmin) then
                if (numeStore .gt. curmax) goto 40
                fieldName = ' '
                ierd = 1
                if (.not. isOptionFromUser) then
                    call rsexch(' ', resultIn, option, numeStore, fieldName, ierd)
                end if
                if (ierd .ne. 0) then
                    call rsexch(' ', resultOut, option, numeStore, fieldName, ierd)
                end if
                if (ierd .ne. 0) then
                    nbStoreToCompute = nbStoreToCompute+1
                    zi(jordop+nbStoreToCompute+2) = numeStore
                end if
            end if
30          continue
        end do
40      continue
        zi(jordop+1) = curmin
        zi(jordop+2) = curmax
        zi(jordop) = nbStoreToCompute
    else
!       CAS 2 : AUCUNE DEPENDANCE, ON PEUT DONC RECOPIER LA LISTE
!               DES NUMEROS D'ORDRE EN VERIFIANT QUE L'OPTION
!               N'EXISTE NI DANS RESUIN NI DANS RESUOU
        nbStoreToCompute = 0
        do iStore = 1, nbStore
            numeStore = listStore(iStore)
            if ((isodep .eq. '-') .and. (numeStore .eq. numeStoreMin)) then
                goto 20
            else if ((isodep .eq. '+') .and. (numeStore .eq. numeStoreMax)) then
                goto 20
            end if
            fieldName = ' '
            ierd = 1
            if (.not. isOptionFromUser) then
                call rsexch(' ', resultIn, option, numeStore, fieldName, ierd)
            end if
            if (ierd .ne. 0) then
                call rsexch(' ', resultOut, option, numeStore, fieldName, ierd)
            end if
!
            if (ierd .ne. 0) then
                nbStoreToCompute = nbStoreToCompute+1
                zi(jordop+nbStoreToCompute+2) = numeStore
            end if
20          continue
        end do
        zi(jordop+1) = listStore(1)
        zi(jordop+2) = listStore(nbStore)
        zi(jordop) = nbStoreToCompute
    end if
!
    zk24(jlnoin+numeOptEff-1) = listStoreOptJv
!
    call jedema()
!
end subroutine
