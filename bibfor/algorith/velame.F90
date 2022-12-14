! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine velame(modele, charge, infcha, depmoz, vecelz)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/exisd.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/vtgpld.h"
    character(len=*) :: vecelz, depmoz
    character(len=24) :: modele, charge, infcha
! ----------------------------------------------------------------------
!     CALCUL DES VECTEURS ELEMENTAIRES DES FORCES DE LAPLACE
!     PRODUIT UN VECT_ELEM DEVANT ETRE ASSEMBLE PAR LA ROUTINE ASASVE
!
! IN  MODELE  : NOM DU MODELE
! IN  CHARGE  : LISTE DES CHARGES
! IN  INFCHA  : INFORMATIONS SUR LES CHARGEMENTS
! IN  DEPMOI  : DEPLACEMENT A L'INSTANT TEMMOI
! OUT/JXOUT  VECELZ  : VECT_ELEM RESULTAT.
!
!   ATTENTION :
!   -----------
!   LE VECT_ELEM (VECELZ) RESULTAT A 2 PARTICULARITES :
!   1) LE NOM DES RESUELEM COMMENCE PAR '&&ASASVE.'
!   2) LES RESUELEM VERIFIENT LA  PROPRIETE :
!      AU CHARGEMENT ELEMENTAIRE
!      (ICHA=0 SI IL N'Y A PAS DE CHARGE)
!
!
!
    character(len=8) :: nomcha, lpain(3), paout
    character(len=8) :: lcmp(2), newnom
    character(len=16) :: option
    character(len=19) :: resuel, vecele, depmoi
    character(len=24) :: chgeom, chlapl, chgeo2
    character(len=24) :: ligrmo, ligrch, lchin(3), kcmp(2)
    integer :: iret, nchar
    aster_logical :: bidon
    integer :: icha, ifla, j, jchar, jinf, lonlis
!-----------------------------------------------------------------------
    call jemarq()
    newnom = '.0000000'
!
    vecele = vecelz
    if (vecele .eq. ' ') then
        vecele = '&&VELAME'
    endif
    resuel = '&&VELAME.???????'
    depmoi = depmoz
!
    bidon = .true.
    lonlis = 0
    call jeexin(charge, iret)
    if (iret .ne. 0) then
        call jelira(charge, 'LONMAX', nchar)
        if (nchar .ne. 0) then
            bidon = .false.
            call jeveuo(charge, 'L', jchar)
            call jeveuo(infcha, 'L', jinf)
            lonlis = (zi(jinf+2*nchar+2))*nchar
            if (lonlis .eq. 0) bidon = .true.
        endif
    endif
!
!
!     -- ALLOCATION DU VECT_ELEM RESULTAT :
!     -------------------------------------
    call detrsd('VECT_ELEM', vecele)
    call memare('V', vecele, modele(1:8), ' ', ' ',&
                'CHAR_MECA')
    call reajre(vecele, ' ', 'V')
    if (bidon) goto 40
!
    ligrmo = modele(1:8)//'.MODELE'
    call megeom(modele(1:8), chgeom)
!
!     REACTUALISATION DE LA GEOMETRIE SI DEPMOI EXISTE
    if (depmoi .ne. ' ') then
        chgeo2 = '&&VELAME.CH_GEOMER'
        call vtgpld('CUMU', 1.d0, chgeom, depmoi, 'V',&
                    chgeo2)
    else
        chgeo2 = chgeom
    endif
!
!
    option = 'CHAR_MECA_FRLAPL'
    lpain(1) = 'PFLAPLA'
    lpain(2) = 'PGEOMER'
    lchin(2) = chgeo2
    lpain(3) = 'PLISTMA'
    paout = 'PVECTUR'
!
    ifla = 0
    do icha = 1, nchar
        nomcha = zk24(jchar+icha-1) (1:8)
        ligrch = nomcha//'.CHME.LIGRE'
        lchin(3) (1:17) = ligrch(1:13)//'.FL1'
        do j = 1, 99
            call codent(j, 'D0', lchin(3) (18:19))
            lchin(3) = lchin(3) (1:19)//'.DESC'
            call exisd('CHAMP_GD', lchin(3), iret)
            if (iret .ne. 0) then
                if (ifla .eq. 0) then
                    chlapl = '&&VELAME.CH_FLAPLA'
                    lcmp(1) = 'NOMAIL'
                    lcmp(2) = 'NOGEOM'
                    kcmp(1) = chgeom(1:8)
                    kcmp(2) = chgeo2(1:19)
                    call mecact('V', chlapl, 'MODELE', modele, 'FLAPLA  ',&
                                ncmp=2, lnomcmp=lcmp, vk=kcmp)
                    lchin(1) = chlapl
                    ifla = 1
                endif
                call gcnco2(newnom)
                resuel(10:16) = newnom(2:8)
                call corich('E', resuel, ichin_ = icha)
!
                call calcul('S', option, ligrmo, 3, lchin,&
                            lpain, 1, resuel, paout, 'V',&
                            'OUI')
                call reajre(vecele, resuel, 'V')
            else
                goto 20
            endif
        end do
 20     continue
!
    end do
!
 40 continue
    vecelz = vecele
!
    call jedema()
end subroutine
