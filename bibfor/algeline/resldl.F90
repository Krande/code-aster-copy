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

subroutine resldl(solveu, nommat, vcine, nsecm, rsolu, &
                  csolu, prepos)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/csmbgg.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mcconl.h"
#include "asterfort/mrconl.h"
#include "asterfort/mtdscr.h"
#include "asterfort/rldlg3.h"
#include "asterfort/utmess.h"
#include "asterfort/isParallelMatrix.h"
!
    character(len=*) :: nommat, vcine
    integer(kind=8) :: nsecm
    real(kind=8) :: rsolu(*)
    complex(kind=8) :: csolu(*)
    aster_logical :: prepos
!
! BUT : RESOUDRE UN SYSTEME LINEAIRE D'EQUATIONS (REEL OU COMPLEXE)
!       SOLVEUR = 'LDLT' OU 'MULT_FRONT'
!-----------------------------------------------------------------------
! IN/JXIN  K19 SOLVEU : SD_SOLVEUR
! IN/JXIN  K19 NOMMAT : MATR_ASSE PREMIER MEMBRE DU SYSTEME LINEAIRE
! IN/JXIN  K*  VCINE  : CHAMP ASSOCIE AUX CHARGES CINEMATIQUES (OU ' ')
! IN       I   NSECM  :  N : NOMBRE DE SECONDS MEMBRES
! IN/OUT   R   RSOLU(*,NSECM)  :
!        EN ENTREE : VECTEUR DE REELS CONTENANT LES SECONDS MEMBRES
!        EN SORTIE : VECTEUR DE REELS CONTENANT LES SOLUTIONS
! IN/OUT   C   CSOLU(*,NSECM)  : IDEM RSOLU POUR LES COMPLEXES.
! IN      LOG  PREPOS : SI .TRUE. ON FAIT LES PRE ET POSTTRAITEMENTS DE
!           MISE A L'ECHELLE DU RHS ET DE LA SOLUTION (MRCONL) ET DE LA
!           PRISE EN COMPTE DES AFFE_CHAR_CINE (CSMBGG).
!           SI .FALSE. ON NE LES FAIT PAS (PAR EXEMPLE EN MODAL).
!-----------------------------------------------------------------------
!
!
    character(len=19) :: nomma2
    character(len=8) :: type
    character(len=16) :: metres
    character(len=19) :: vci19, solveu
    complex(kind=8) :: cbid
    integer(kind=8) :: k, kdeb, idvalc, lmat, neq, nimpo
    character(len=24), pointer :: slvk(:) => null()
    aster_logical :: l_parallel_matrix
    cbid = dcmplx(0.d0, 0.d0)
!     ------------------------------------------------------------------
!
    call jemarq()
    vci19 = vcine
    nomma2 = nommat
!
    call jeveuo(solveu//'.SLVK', 'L', vk24=slvk)
    metres = slvk(1) (1:16)
!
    l_parallel_matrix = isParallelMatrix(nomma2)
    if (l_parallel_matrix) then
        call utmess('F', 'FACTOR_93')
    end if
!
    call mtdscr(nomma2)
    call jeveuo(nomma2(1:19)//'.&INT', 'E', lmat)
    if (lmat .eq. 0) then
        call utmess('F', 'ALGELINE3_40')
    end if
!
    neq = zi(lmat+2)
    nimpo = zi(lmat+7)
    if (vci19 .eq. ' ') then
! --- SI ON NE FAIT PAS LES PREPOS, ON NE SE PREOCCUPE PAS DES
!     AFFE_CHAR_CINE. DONC C'EST NORMAL QUE L'INFO SOIT INCOHERENTE
!     A CE NIVEAU
        if ((nimpo .ne. 0) .and. prepos) then
            call utmess('F', 'ALGELINE3_41')
        end if
        idvalc = 0
    else
        call jeveuo(vci19//'.VALE', 'L', idvalc)
        call jelira(vci19//'.VALE', 'TYPE', cval=type)
        if (((type .eq. 'R') .and. (zi(lmat+3) .ne. 1)) .or. &
            ((type .eq. 'C') .and. (zi(lmat+3) .ne. 2))) then
            call utmess('F', 'ALGELINE3_42')
        end if
    end if
!
    if (zi(lmat+3) .eq. 1) then
        type = 'R'
    else if (zi(lmat+3) .eq. 2) then
        type = 'C'
    else
        ASSERT(.false.)
    end if
!
!
!
!
    if (type .eq. 'R') then
!     ----------------------------------------
        if (prepos) then
!         MISE A L'ECHELLE DES LAGRANGES DANS LE SECOND MEMBRE
            call mrconl('MULT', lmat, 0, 'R', rsolu, &
                        nsecm)
            if (idvalc .ne. 0) then
                do k = 1, nsecm
                    kdeb = (k-1)*neq+1
                    call csmbgg(lmat, rsolu(kdeb), zr(idvalc), [cbid], [cbid], &
                                'R')
                end do
            end if
        end if
        call rldlg3(metres, lmat, rsolu, [cbid], nsecm)
        if (prepos) then
!         MISE A L'ECHELLE DES LAGRANGES DANS LA SOLUTION
            call mrconl('MULT', lmat, 0, 'R', rsolu, &
                        nsecm)
        end if
!
!
    else if (type .eq. 'C') then
!     ----------------------------------------
        if (prepos) then
!         MISE A L'ECHELLE DES LAGRANGES DANS LE SECOND MEMBRE
            call mcconl('MULT', lmat, 0, 'C', csolu, &
                        nsecm)
            if (idvalc .ne. 0) then
                do k = 1, nsecm
                    kdeb = (k-1)*neq+1
                    call csmbgg(lmat, [0.d0], [0.d0], csolu(kdeb), zc(idvalc), &
                                'C')
                end do
            end if
        end if
        call rldlg3(metres, lmat, [0.d0], csolu, nsecm)
        if (prepos) then
!         MISE A L'ECHELLE DES LAGRANGES DANS LA SOLUTION
            call mcconl('MULT', lmat, 0, 'C', csolu, &
                        nsecm)
        end if
    end if
!
!
    call jedema()
end subroutine
