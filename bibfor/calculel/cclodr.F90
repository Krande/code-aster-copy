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
subroutine cclodr(numeOptEff, nbStore, listStore, jvBaseName, numeStoreMin, &
                  numeStoreMax, resultIn, resultOut, lacalc)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
!
    integer(kind=8), intent(in) :: numeOptEff
    integer(kind=8), intent(in) :: nbStore
    integer(kind=8), pointer :: listStore(:)
    character(len=8), intent(in) :: jvBaseName
    integer(kind=8), intent(in) ::  numeStoreMin, numeStoreMax
    character(len=8), intent(in) :: resultIn, resultOut
    integer(kind=8), pointer :: lacalc(:)
!
! --------------------------------------------------------------------------------------------------
!
!  CALC_CHAMP
!
!  DETERMINATION LISTE OPTIONS AVEC DEPENDANCE REDUITE
!
! --------------------------------------------------------------------------------------------------
!
!  MODIFICATION DE LACALC EN METTANT DES 0 LORSQUE L'OPTION NE DOIT
!   PAS ETRE CALCULEE
!
! In  numeOptEff        : index of option in list of option
! In  nbStore           : number of storing indexes
! In  listStore         : list of storing indexes
! In  jvBaseName        : string base for names of JEVEUX objects
! In  numeStoreMin      : minimum storing index in list of storing indexes to compute option
! In  numeStoreMax      : maximum storing index in list of storing indexes to compute option
! In  resultIn          : name of datastructure for input results
! In  resultOut         : name of datastructure for output results
! Ptr lacalc            : pointer to list of options to have to compute
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jlisop, jliori, jlidep, ierd, inddeb, indfin
    integer(kind=8) :: iStore, curmax, curmin, iter, decal, numeStore, jlnoin
    integer(kind=8) :: jordo2, jlisde
    character(len=1) :: isodep
    character(len=16) :: option
    character(len=24) :: fieldName, noliop, nolori, noldep, noliin, nolisd
    aster_logical :: exitor
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Access to objects
    noliop = jvBaseName//'.LISOPT'
    nolori = jvBaseName//'.LISORI'
    noldep = jvBaseName//'.LISDEP'
    noliin = jvBaseName//'.LNOINS'
    nolisd = jvBaseName//'.ISODEP'
    call jeveuo(noliop, 'L', jlisop)
    call jeveuo(nolori, 'L', jliori)
    call jeveuo(noldep, 'L', jlidep)
    call jeveuo(noliin, 'L', jlnoin)
    call jeveuo(nolisd, 'L', jlisde)

! - Parameters for current option
    option = zk16(jlisop+numeOptEff-1)
    inddeb = zi(jliori+2*numeOptEff-2)
    indfin = zi(jliori+2*numeOptEff-1)
    isodep = zk8(jlisde+numeOptEff-1) (1:1)
!
    if (inddeb .ne. 0) then
!       CAS 1 : CETTE OPTION DEPEND D'AUTRES OPTIONS A CALCULER
!               AUQUEL CAS, IL FAUT REGARDER COMMENT ELLE EN DEPEND
!               ET LA LISTE DES NUMEROS D'ORDRE DE SES PARENTS
        call jeveuo(noliin, 'L', jlnoin)
        curmax = numeStoreMax
        curmin = numeStoreMin
        do iter = inddeb, indfin
            call jeveuo(zk24(jlnoin+iter-1), 'L', jordo2)
!
            if (zk8(jlidep+iter-1) .eq. 'NP1') then
                decal = -1
            else if (zk8(jlidep+iter-1) .eq. 'NM1') then
                decal = +1
            else
                decal = 0
            end if
!
!         LA LISTE DE NUMEROS D'ORDRE PROVIENT DE OP0058
!         ELLE EST DONC CROISSANTE
            curmax = min(curmax, zi(jordo2+2)+decal)
            curmin = max(curmin, zi(jordo2+1)+decal)
        end do
!
        exitor = .true.
        if (lacalc(numeOptEff) .eq. 1) then
            do iStore = 1, nbStore
                numeStore = listStore(iStore)
                if ((isodep .eq. '-') .and. (numeStore .eq. numeStoreMin)) then
                    goto 30
                elseif ((isodep .eq. '+') .and. (numeStore .eq. numeStoreMax)) &
                    then
                    goto 30
                end if
                if (numeStore .ge. curmin) then
                    if (numeStore .gt. curmax) goto 40
                    fieldName = ' '
                    call rsexch(' ', resultIn, option, numeStore, fieldName, ierd)
                    if (ierd .ne. 0) then
                        call rsexch(' ', resultOut, option, numeStore, fieldName, ierd)
                    end if
                    if (ierd .ne. 0) then
                        exitor = .false.
                    end if
                end if
30              continue
            end do
        end if
40      continue
        if (exitor) then
            lacalc(inddeb:indfin) = 0
        end if
    end if
!
    call jedema()
!
end subroutine
