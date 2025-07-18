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

subroutine celcel(transf, cel1, base, cel2)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=*) :: transf, cel1, base, cel2
! person_in_charge: jacques.pellet at edf.fr
! ----------------------------------------------------------------------
!  BUT : TRANSFORMER UN CHAM_ELEM (CEL1) EN UN CHAM_ELEM (CEL2) QUI
!        A DE NOUVELLES PROPRIETES. PAR EXEMPLE POUR POUVOIR ETRE
!        IMPRIME, POST-TRAITE, ...
!
!  IN        TRANSF K*  : NOM DE LA TRANSFORMATION :
!
!  'NBVARI_CST' : SI LES ELEMENTS DU CHAM_ELEM (VARI_R) N'ONT
!                   PAS LE MEME NOMBRE DE DE CMPS, ON LES ALIGNE
!                   TOUS SUR LE NOMBRE MAX (SITUATION AVANT 10/99)
!
!  'PAS_DE_SP'  : SI DES ELEMENTS DU CHAM_ELEM ONT DES SOUS-POINTS
!                 ON MET LE MODE LOCAL DE LEUR GREL A 0.
!                 C'EST COMME S'ILS AVAIENT ETE EFFACES.
!
!
!  IN/JXIN   CEL1   K19 : CHAM_ELEM A TRANSFORMER
!  IN        BASE   K1  : /'V'/'G'/ : BASE DE CREATION DE CEL2
!  IN/JXOUT  CEL2   K19 : CHAM_ELEM APRES TRANSFORMATION
! ----------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    integer(kind=8) :: ima, ipt, ispt, icmp, nbma, ibid
    integer(kind=8) :: jcesd1, jcesl1, jcesc1
    integer(kind=8) :: jcesd2, jcesl2, nbgrel, igrel, debugr, nbel
    integer(kind=8) :: nbpt, nbspt, ncmp2, iad1, iad2, imolo, nbspmx, iel, nbsp
    integer(kind=8) :: ncmpg, nbvamx, nncp, ico
    character(len=19) :: ces1, ces2, ligrel, cel11, cel22
    character(len=16) :: optini, nompar
    character(len=8) :: ma, nomgd, typces, kbid
    character(len=24) :: valk(3)
    integer(kind=8), pointer :: celd(:) => null()
    character(len=24), pointer :: celk(:) => null()
    character(len=8), pointer :: cesk(:) => null()
    real(kind=8), pointer :: cesv1(:) => null()
    real(kind=8), pointer :: cesv2(:) => null()
    integer(kind=8), pointer :: vnbpt(:) => null()
    integer(kind=8), pointer :: vnbspt(:) => null()
! -DEB------------------------------------------------------------------
    call jemarq()
!
    call dismoi('NOM_GD', cel1, 'CHAMP', repk=nomgd)
!
!
    if (transf .eq. 'NBVARI_CST') then
!     =================================
!
        if (nomgd .ne. 'VARI_R') then
!         -- IL N'Y A RIEN A FAIRE : NBVARI EST CST !
            call copisd('CHAMP_GD', base, cel1, cel2)
            goto 80
        end if
!
!       1- ON TRANSFORME CEL1 EN CHAM_ELEM_S : CES1
!       -------------------------------------------
        ces1 = '&&CELCEL.CES1'
        call celces(cel1, 'V', ces1)
        call jeveuo(ces1//'.CESD', 'L', jcesd1)
        call jeveuo(ces1//'.CESL', 'L', jcesl1)
        call jeveuo(ces1//'.CESV', 'L', vr=cesv1)
        call jeveuo(ces1//'.CESC', 'L', jcesc1)
        call jeveuo(ces1//'.CESK', 'L', vk8=cesk)
!
!       2- ON ALLOUE UN CHAM_ELEM_S PLUS GROS: CES2
!       -------------------------------------------
        ces2 = '&&CELCEL.CES2'
        ma = cesk(1)
        typces = cesk(3)
        nbma = zi(jcesd1-1+1)
        ncmpg = zi(jcesd1-1+2)
        nbvamx = zi(jcesd1-1+5)
        ASSERT(ncmpg .eq. nbvamx)
!
!       2.1 : CALCUL DE 2 VECTEURS CONTENANT LE NOMBRE DE
!             POINTS DE SOUS-POINTS DES MAILLES
!       ---------------------------------------------------
        AS_ALLOCATE(vi=vnbpt, size=nbma)
        AS_ALLOCATE(vi=vnbspt, size=nbma)
        do ima = 1, nbma
            vnbpt(ima) = zi(jcesd1-1+5+4*(ima-1)+1)
            vnbspt(ima) = zi(jcesd1-1+5+4*(ima-1)+2)
        end do
!
!       2.2 : ALLOCATION DE CES2 :
!       ---------------------------------------------------
        call cescre('V', ces2, typces, ma, nomgd, &
                    -nbvamx, kbid, vnbpt, vnbspt, [-nbvamx])
        call jeveuo(ces2//'.CESD', 'L', jcesd2)
        call jeveuo(ces2//'.CESL', 'E', jcesl2)
        call jeveuo(ces2//'.CESV', 'E', vr=cesv2)
!
!
!
!
!       3- ON RECOPIE LES VALEURS DE CES1 DANS CES2 :
!       ---------------------------------------------
        do ima = 1, nbma
            nbpt = zi(jcesd1-1+5+4*(ima-1)+1)
            nbspt = zi(jcesd1-1+5+4*(ima-1)+2)
!
            ncmp2 = zi(jcesd2-1+5+4*(ima-1)+3)
!
            do ipt = 1, nbpt
                do ispt = 1, nbspt
                    do icmp = 1, ncmp2
                        call cesexi('C', jcesd1, jcesl1, ima, ipt, &
                                    ispt, icmp, iad1)
                        call cesexi('C', jcesd2, jcesl2, ima, ipt, &
                                    ispt, icmp, iad2)
                        ASSERT(iad2 .lt. 0)
                        zl(jcesl2-1-iad2) = .true.
                        if (iad1 .gt. 0) then
                            cesv2(1-1-iad2) = cesv1(iad1)
!
                        else
                            cesv2(1-1-iad2) = 0.d0
                        end if
                    end do
                end do
            end do
        end do
!
!
!
!       4- ON TRANSFORME CES2 EN CHAM_ELEM : CEL2
!       -------------------------------------------
        cel11 = cel1
        call jeveuo(cel11//'.CELK', 'L', vk24=celk)
        ligrel = celk(1)
        optini = celk(2)
        nompar = celk(6)
        call cescel(ces2, ligrel, optini, nompar, 'NON', &
                    nncp, base, cel2, 'F', ibid)
!
!
!       5- MENAGE :
!       -------------------------------------------
        AS_DEALLOCATE(vi=vnbpt)
        AS_DEALLOCATE(vi=vnbspt)
        call detrsd('CHAM_ELEM_S', ces1)
        call detrsd('CHAM_ELEM_S', ces2)
!
!
    else if (transf .eq. 'PAS_DE_SP') then
!     =====================================
!
        call copisd('CHAMP_GD', base, cel1, cel2)
        cel22 = cel2
        call jeveuo(cel22//'.CELD', 'E', vi=celd)
        nbgrel = celd(2)
!
!       -- ON MET A ZERO LE MODE LOCAL DES GRELS QUI ONT DES
!          SOUS-POINTS :
        ico = 0
        do igrel = 1, nbgrel
            debugr = celd(4+igrel)
            nbel = celd(debugr+1)
            imolo = celd(debugr+2)
            if (imolo .gt. 0) then
                nbspmx = 0
                do iel = 1, nbel
                    nbsp = celd(debugr+4+4*(iel-1)+1)
                    nbspmx = max(nbspmx, nbsp)
                end do
                if (nbspmx .gt. 1) then
                    celd(debugr+2) = 0
                else
                    ico = ico+1
                end if
            end if
        end do
        if (ico .eq. 0) then
            valk(1) = cel1
            valk(2) = nomgd
            call utmess('F', 'CALCULEL2_40', nk=2, valk=valk)
        end if
!
!
    else
!       CAS RESTANT A PROGRAMMER ...
        ASSERT(.false.)
    end if
!
80  continue
    call jedema()
end subroutine
