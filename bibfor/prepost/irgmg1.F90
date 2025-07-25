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
subroutine irgmg1(numold, ima, nbord2, tabd, tabl, &
                  tabv, partie, jtype, nbno, icmp, &
                  ifi, iwri, iadmax)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cesexi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: numold(*), tabd(*), tabl(*), tabv(*), jtype
    integer(kind=8) :: icmp, ifi, ima, nbord2, iadmax, nbno
    aster_logical :: iwri
    character(len=*) :: partie
!
!     BUT: ECRITURE D'UNE CMP D'UN CHAMP "ELGA" OU "ELEM"
!     POUR UN TYPE D'ELEMENT AU FORMAT GMSH
!
!     ENTREE:
!     NUMOLD : I   : TABLEAU DE CORRESPONDANCE NOUV MAILLE ANC. MAILLE
!     IMA    : I   : NUMERO NOUVELLE MAILLE
!     NBORD2 : I   : NOMBRE DE NUM D'ORDRE
!     TABD   : I   : DECRIPTEURS DU CHMAP SIMPLE A IMMRIMER
!     TABL   : I   : DECRIPTEURS DU CHMAP SIMPLE A IMMRIMER
!     TABV   : I   : DECRIPTEURS DU CHMAP SIMPLE A IMMRIMER
!     PARTIE : K4  : IMPRESSION DE LA PARTIE COMPLEXE OU REELLE DU CHAMP
!     JTYPE  : I   : ADRESSE DU TYPE DU CHAMP ( REEL OU COMPLEXE )
!     ICMP   : I   : NUMERO COMPOSANTE CHAMP
!     IFI    : I   : NUMERO D'UNITE LOGIQUE DU FICHIER GMSH
!     IWRI    : L   : INDIQUE SI ON DOIT ECRIRE
!     SORTIE
!     IADMAX  : I   : MAX DES IAD SI >0 LE CHAMP EXISTE POUR LA MAILLE
!
!     ------------------------------------------------------------------
    integer(kind=8) :: imaold, ior, jcesd, jcesl, jcesv, nbpt, nbsp, ipt, isp, iad, ino
    real(kind=8) :: vale
!
!     ------------------------------------------------------------------
!
    call jemarq()
!
    imaold = numold(ima)
!
!     ON NE TRAITE QUE LES CHAMPS A 1 SOUS-POINT,
!     ET UNE SEULE VALEUR SCALAIRE
!
    isp = 1
    iadmax = 0
!
    do ior = 1, nbord2
        jcesd = tabd(ior)
        jcesl = tabl(ior)
        jcesv = tabv(ior)
        nbpt = zi(jcesd-1+5+4*(imaold-1)+1)
        nbsp = zi(jcesd-1+5+4*(imaold-1)+2)
        if (nbsp .ne. 1) then
            call utmess('F', 'PREPOST2_57')
        end if
        vale = 0.d0
        if (zk8(jtype-1+ior) .eq. 'R') then
            do ipt = 1, nbpt
                call cesexi('C', jcesd, jcesl, imaold, ipt, &
                            isp, icmp, iad)
                if (iad .gt. 0) then
                    iadmax = iad
                    vale = vale+zr(jcesv-1+iad)
                end if
            end do
        else if (zk8(jtype-1+ior) .eq. 'C') then
            if (partie .eq. 'REEL') then
                do ipt = 1, nbpt
                    call cesexi('C', jcesd, jcesl, imaold, ipt, &
                                isp, icmp, iad)
                    if (iad .gt. 0) then
                        iadmax = iad
                        vale = vale+dble(zc(jcesv-1+iad))
                    end if
                end do
            else if (partie .eq. 'IMAG') then
                do ipt = 1, nbpt
                    call cesexi('C', jcesd, jcesl, imaold, ipt, &
                                isp, icmp, iad)
                    if (iad .gt. 0) then
                        iadmax = iad
                        vale = vale+dimag(zc(jcesv-1+iad))
                    end if
                end do
            end if
        end if
        if (abs(vale) .le. 1.d-99) vale = 0.d0
        if (nbpt .ne. 0) vale = vale/nbpt
        if (iwri) then
            do ino = 1, nbno
                write (ifi, 1000) vale
            end do
        end if
    end do
!
    call jedema()
!
1000 format(1p, e15.7e3)
!
end subroutine
