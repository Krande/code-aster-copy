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

subroutine cbval2(nbcomb, typcst, const, lmat, typres, &
                  lres, ddlexc)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/pteddl.h"
#include "asterfort/rrssm2.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbcomb, lmat(*), lres
    character(len=*) :: typres, ddlexc, typcst(*)
    real(kind=8) :: const(*)
!     COMBINAISON LINEAIRE DES .VALM DES MATRICES
!       *  LES MATRICES SONT SUPPOSEES AVOIR LE MEME TYPE DE STOCKAGE
!          (MORSE) MAIS ELLES ONT DES PROFILS DIFFERENTS
!       *  POUR L'INSTANT ON NE TRAITE QUE LE CAS DE MATRICES
!          SYMETRIQUES, REELLES.
!       *  LES SCALAIRES SONT PRIS REELS POUR L'INSTANT
!     -----------------------------------------------------------------
! IN  I  NBCOMB = NOMBRE DE MATRICES A COMBINER
! IN  R  CONST  = TABLEAU DE R*8    DES COEFICIENTS
! IN  I  LMAT = TABLEAU DES POINTEURS DES MATRICES
! IN  K* TYPRES = TYPE DES MATRICES   (R)
! IN  I  LRES = POINTEUR DE MATRICE RESULTAT
! IN  K* DDLEXC = NOM DES DDLS A EXCLURE (CONCRETEMENT IL S'AGIT
!                                         DES LAGRANGE)
!
!     -----------------------------------------------------------------
!
    aster_logical :: symr, symi, symrl
!     -----------------------------------------------------------------
    integer(kind=8) :: lgbloc
    character(len=1) :: clas, typmat
    character(len=8) :: nomddl
    character(len=14) :: numr, numi
    character(len=19) :: matres, mati
    character(len=24) :: valmi, valmr
    integer(kind=8) :: neq, mxddl, lddl, jsmhcr
    integer(kind=8) :: iconst, imat, jsmhci, jvlmi1, jvlmr1, k
    integer(kind=8) :: jvlmi2, jvlmr2
    real(kind=8) :: zero
    character(len=24), pointer :: refai(:) => null()
    character(len=24), pointer :: refar(:) => null()
    integer(kind=8), pointer :: smdii(:) => null()
    integer(kind=8), pointer :: smdir(:) => null()
!     -----------------------------------------------------------------
    call jemarq()
    zero = 0.d0
!
    nomddl = ddlexc
    matres = zk24(zi(lres+1)) (1:19)
    ASSERT(typres .eq. 'R' .or. typres .eq. 'C')
    neq = zi(lres+2)
    valmr = matres//'.VALM'
    lgbloc = zi(lres+14)
    call jelira(matres//'.REFA', 'CLAS', cval=clas)
    call jeveuo(matres//'.REFA', 'L', vk24=refar)
    ASSERT(refar(9) (1:1) .eq. 'M')
    symr = refar(9) .eq. 'MS'
!
    mxddl = 1
    call dismoi('NOM_NUME_DDL', matres, 'MATR_ASSE', repk=numr)
    call wkvect('&&CBVAL2', 'V V I', neq*mxddl, lddl)
    call pteddl('NUME_DDL', numr, mxddl, nomddl, neq, &
                tabl_equa=zi(lddl))
!
!
    call jeveuo(numr//'.SMOS.SMHC', 'L', jsmhcr)
    call jeveuo(numr//'.SMOS.SMDI', 'L', vi=smdir)
    call jeveuo(jexnum(valmr, 1), 'E', jvlmr1)
    if (typres(1:1) .eq. 'R') then
        do k = 1, lgbloc
            zr(jvlmr1-1+k) = zero
        end do
        if (.not. symr) then
            call jeveuo(jexnum(valmr, 2), 'E', jvlmr2)
            do k = 1, lgbloc
                zr(jvlmr2-1+k) = zero
            end do
        end if
!
    else
        call utmess('F', 'ALGELINE_5')
    end if
!
!
! --- BOUCLE SUR LES MATRICES A COMBINER :
!     ----------------------------------
    iconst = 1
    do imat = 1, nbcomb
        ASSERT(typcst(imat) .eq. 'R')
        mati = zk24(zi(lmat(imat)+1)) (1:19)
        call dismoi('NOM_NUME_DDL', mati, 'MATR_ASSE', repk=numi)
        call jeveuo(numi//'.SMOS.SMHC', 'L', jsmhci)
        call jeveuo(numi//'.SMOS.SMDI', 'L', vi=smdii)
        call jeveuo(mati//'.REFA', 'L', vk24=refai)
        valmi = mati//'.VALM'
        symi = refai(9) .eq. 'MS'
        symrl = symr
        if (.not. symi) then
            ASSERT(.not. symrl)
        end if
        call jelira(valmi, 'TYPE', cval=typmat)
        call jeveuo(jexnum(valmi, 1), 'L', jvlmi1)
        ASSERT(typmat .eq. 'R')
        if (.not. symi) call jeveuo(jexnum(valmi, 2), 'L', jvlmi2)
!
        if (typres(1:1) .eq. 'R') then
            if (typmat .eq. 'R') then
                call rrssm2(neq, zi4(jsmhcr), zi4(jsmhci), smdir, smdii, &
                            zi(lddl), const(iconst), zr(jvlmi1), zr(jvlmr1))
                if (.not. symr) then
                    if (.not. symi) then
                        call rrssm2(neq, zi4(jsmhcr), zi4(jsmhci), smdir, smdii, &
                                    zi(lddl), const(iconst), zr(jvlmi2), zr(jvlmr2))
!
                    else
                        call rrssm2(neq, zi4(jsmhcr), zi4(jsmhci), smdir, smdii, &
                                    zi(lddl), const(iconst), zr(jvlmi1), zr(jvlmr2))
                    end if
                end if
            end if
        end if
        iconst = iconst+1
        call jelibe(jexnum(valmi, 1))
        if (.not. symi) call jelibe(jexnum(valmi, 2))
    end do
!
    call jedetr('&&CBVAL2')
!
    call jedema()
!
end subroutine
