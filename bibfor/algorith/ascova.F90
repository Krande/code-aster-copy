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
subroutine ascova(detr, vachar, fomulz, npara, vpara, &
                  typres, cnchar, basez)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/fointc.h"
#include "asterfort/fointe.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/vtcmbl.h"
#include "asterfort/wkvect.h"
    character(len=*) :: fomulz, npara, typres, cnchar, detr
    character(len=24) :: vachar
    real(kind=8) :: vpara, tval(1)
    character(len=1), intent(in), optional :: basez
! ----------------------------------------------------------------------
!  BUT :   AJOUTER / COMBINER DES VECTEURS ASSEMBLES (CHAM_NO)
!
! IN  DETR    : / 'D' : ON DETRUIT  LE VACHAR  (CONSEILLE)
!               / 'G' : ON NE DETRUIT PAS LE VACHAR
! IN/JXVAR  VACHAR  : LISTE DES VECTEURS ASSEMBLES
! IN  FOMULT  : LISTE DES FONCTIONS MULTIPLICATIVES
! IN  NPARA   : NOM DU PARAMETRE
! IN  VPARA   : VALEUR DU PARAMETRE
! IN  TYPRES  : TYPE DES VECTEURS ET DU CHAM_NO RESULTANT 'R' OU 'C'
! VAR/JXOUT  CNCHAR  : CHAM_NO RESULTAT
!
! REMARQUES :
! --------------
! CETTE ROUTINE SERT A COMBINER LES CHAM_NO D'UN "VACHAR"
! UN VACHAR EST UN VECTEUR DE K24 CONTENANT UNE LISTE DE CHAM_NO.
! ON PEUT OBTENIR UN VACHAR PAR ASASVE PAR EXEMPLE.
!
! ATTENTION: SI DETR='D', CETTE ROUTINE DETRUIT LES CHAM_NO DU VACHAR
! =========  APRES LES AVOIR COMBINES. ELLE DETRUIT EGALEMENT LE VACHAR.
!
! POUR TENIR EVENTUELLEMENT COMPTE D'UNE FONC_MULT, ASCOVA
! UTILISE SYSTEMATIQUEMENT LA ROUTINE CORICH SUR LES CHAM_NO
! A COMBINER. IL FAUT QUE LES CHAM_NO DU VACHAR AIENT ETE
! RENSEIGNES PAR CORICH.
!
!
! SI CNCHAR=' ' LE NOM DU CHAMP RESULTAT SERA :
!   CNCHAR=VACHAR(1:8)//'.ASCOVA'
!
!
!
!
    integer(kind=8) :: kk, nbvec, nchar, iret, jvec, jfonct, jcoef, jtype, k
    integer(kind=8) :: icha, ier, n1, npuis, n2, ltyp
    real(kind=8) :: valres, valre, valim, dgrd, omega, phase
    aster_logical :: fct
    character(len=19) :: chamno
    character(len=24) :: fomult
    complex(kind=8) :: calpha
    character(len=1) :: base
!
    call jemarq()
    fomult = fomulz
    dgrd = r8dgrd()
!
    if (present(basez)) then
        base = basez
    else
        base = 'V'
    end if
!
!
!     -- ON VERIFIE QUE LE VACHAR A LES BONNES PROPRIETES:
!     ----------------------------------------------------
    call jeexin(vachar, iret)
    ASSERT(iret .ne. 0)
    call jelira(vachar, 'LONMAX', nbvec)
    ASSERT(nbvec .ne. 0)
    call jeveuo(vachar, 'L', jvec)
!
!
    call jeexin(fomult, iret)
    if (iret .eq. 0) then
        fct = .false.
        nchar = 0
    else
        fct = .true.
        call jelira(fomult, 'LONMAX', nchar)
        ASSERT(nchar .ne. 0)
        call jeveuo(fomult, 'L', jfonct)
        call jelira(fomult, 'LTYP', ltyp)
    end if
!
!
!
!     -- CAS DES CHAM_NO REELS :
!     ----------------------------------------------------
    if (typres(1:1) .eq. 'R') then
        call wkvect('&&ASCOVA.COEF', 'V V R8', nbvec, jcoef)
        call wkvect('&&ASCOVA.TYPE', 'V V K8', nbvec, jtype)
        do k = 1, nbvec
!
            chamno = zk24(jvec+k-1) (1:19)
            call corich('L', chamno, ichout_=icha)
!
            ASSERT((icha .ne. 0) .and. (icha .ge. -2))
!
            if (icha .eq. -1) then
                valres = 1.d0
            else if (icha .eq. -2) then
                valres = 0.d0
            else
                ASSERT(icha .le. nchar)
                valres = 1.d0
                if (fct) then
                    if (ltyp .eq. 8) then
                        if (zk8(jfonct+icha-1) .ne. ' ') then
                            call fointe('F ', zk8(jfonct+icha-1), 1, [npara], &
                                        [vpara], valres, ier)
                        end if
                    else
                        if (zk24(jfonct+icha-1) .ne. ' ') then
                            call fointe('F ', zk24(jfonct+icha-1), 1, [npara], &
                                        [vpara], valres, ier)
                        end if
                    end if
                end if
            end if
!
            zr(jcoef+k-1) = valres
            zk8(jtype+k-1) = 'R'
        end do
!
!
!     -- CAS DES CHAM_NO COMPLEXES :
!     ----------------------------------------------------
    else
        omega = r8depi()*vpara
        kk = 0
        call wkvect('&&ASCOVA.COEF', 'V V R8', 2*nbvec, jcoef)
        call wkvect('&&ASCOVA.TYPE', 'V V K8', nbvec, jtype)
        do k = 1, nbvec
!
            phase = 0.d0
            npuis = 0
            calpha = cmplx(1.d0, 0.d0)

            call getvr8('EXCIT', 'PHAS_DEG', iocc=k, scal=phase, nbret=n1)
            call getvis('EXCIT', 'PUIS_PULS', iocc=k, scal=npuis, nbret=n2)
            if (n1 .ne. 0) calpha = exp(dcmplx(0.d0, phase*dgrd))
            if (n2 .ne. 0 .and. npuis .ne. 0) calpha = calpha*omega**npuis
!
            chamno = zk24(jvec+k-1) (1:19)
            call corich('L', chamno, ichout_=icha)
!
            ASSERT((icha .ne. 0) .and. (icha .ge. -2))
!
            if (icha .eq. -1) then
                valre = 1.d0
                valim = 0.d0
            else if (icha .eq. -2) then
                valre = 0.d0
                valim = 0.d0
            else
                ASSERT(icha .le. nchar)
                valre = 1.d0
                valim = 0.d0
                tval(1) = vpara
                if (fct .and. zk24(jfonct+icha-1) .ne. ' ' .and. &
                    zk24(jfonct+icha-1) .ne. '&&ACDOCH') then
                    call fointc('F', zk24(jfonct+icha-1) (1:8), 1, npara, tval, &
                                valre, valim, ier)
                end if
            end if
!
            zk8(jtype+k-1) = 'C'
            kk = kk+1
            zr(jcoef+kk-1) = valre*dble(calpha)-valim*dimag(calpha)
            kk = kk+1
            zr(jcoef+kk-1) = valim*dble(calpha)+valre*dimag(calpha)
        end do
    end if
!
!
!     COMBINAISON LINEAIRES DES CHAM_NO :
!     -----------------------------------
    if (cnchar .eq. ' ') cnchar = vachar(1:8)//'.ASCOVA'
    call vtcmbl(nbvec, zk8(jtype), zr(jcoef), zk8(jtype), zk24(jvec), &
                zk8(jtype), cnchar, base)
    call jedetr('&&ASCOVA.COEF')
    call jedetr('&&ASCOVA.TYPE')
!
!
!     DESTRUCTION DU VACHAR :
!     -----------------------------------
    if (detr .eq. 'D') then
        do k = 1, nbvec
            call corich('S', zk24(jvec+k-1))
            call detrsd('CHAMP_GD', zk24(jvec+k-1))
        end do
        call jedetr(vachar)
    else if (detr .eq. 'G') then
!       -- EN PRINCIPE UTILISE PAR DYNA_VIBRA//HARM
    else
        ASSERT(.false.)
    end if
!
!
    call jedema()
end subroutine
