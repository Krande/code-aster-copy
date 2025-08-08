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
subroutine vedimd(nomo, lischa, instan, vecele)
!
!
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/exisd.h"
#include "asterfort/gcnco2.h"
#include "asterfort/inical.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisico.h"
#include "asterfort/lislch.h"
#include "asterfort/lislco.h"
#include "asterfort/lisllc.h"
#include "asterfort/lisltc.h"
#include "asterfort/lisnbg.h"
#include "asterfort/lisnnb.h"
#include "asterfort/lisnol.h"
#include "asterfort/lisopt.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/reajre.h"
#include "asterfort/vemare.h"
    character(len=19) :: lischa, vecele
    character(len=8) :: nomo
    real(kind=8) :: instan
!
! ----------------------------------------------------------------------
!
! CALCUL DES VECTEURS ELEMENTAIRES
!
! MECANIQUE - DIRICHLET
!
! ----------------------------------------------------------------------
!
!
! IN  NOMO   : NOM DU MODELE
! IN  LISCHA : SD LISTE DES CHARGES
! IN  INSTAN : INSTANT DE CALCUL
! OUT VECELE : VECT_ELEM RESULTAT
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nbout, nbin
    parameter(nbout=1, nbin=3)
    character(len=8) :: lpaout(nbout), lpain(nbin)
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    character(len=8) :: nomch0
    character(len=8) :: newnom
    character(len=16) :: option
    character(len=19) :: ligcal
    character(len=13) :: prefob
    character(len=19) :: chgeom, chtime
    character(len=19) :: carte
    character(len=8) :: parain, paraou, typech
    integer(kind=8) :: iret
    integer(kind=8) :: ichar, nbchar, genrec
    aster_logical :: ldual
    character(len=24) :: nomlis
    integer(kind=8) :: jlisci, nbch, indxch
    integer(kind=8) :: nbdual
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    newnom = '.0000000'
    nomlis = '&&NOMLIS'
    call detrsd('VECT_ELEM', vecele)
!
! --- INITIALISATION DES CHAMPS POUR CALCUL
!
    call inical(nbin, lpain, lchin, nbout, lpaout, &
                lchout)
!
! --- NOMBRE DE CHARGES
!
    call lisnnb(lischa, nbchar)
    if (nbchar .eq. 0) goto 99
!
! --- PRESENCE DE CE GENRE DE CHARGEMENT
!
    nbdual = lisnbg(lischa, 'DIRI_DUAL')
    if (nbdual .eq. 0) goto 99
!
! --- ALLOCATION DU VECT_ELEM RESULTAT
!
    call vemare('V', vecele, nomo)
    call reajre(vecele, ' ', 'V')
!
! --- CHAMP DE GEOMETRIE
!
    call megeom(nomo, chgeom)
!
! --- CARTE DE L'INSTANT
!
    chtime = '&&VEDIMD.CH_INST_R'
    call mecact('V', chtime, 'MODELE', nomo, 'INST_R', &
                ncmp=1, nomcmp='INST', sr=instan)
!
! --- CHAMPS D'ENTREES STANDARDS
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PINSTR'
    lchin(2) = chtime
!
! --- LISTE DES INDEX DES CHARGES
!
    call lisnol(lischa, 'DIRI_DUAL', nomlis, nbch)
    ASSERT(nbch .eq. 1)
    call jeveuo(nomlis, 'L', jlisci)
    indxch = zi(jlisci-1+1)
!
! --- CALCUL
!
    do ichar = 1, nbchar
        call lislco(lischa, ichar, genrec)
        ldual = lisico('DIRI_DUAL', genrec)
        if (ldual) then
!
! ------- PREFIXE DE L'OBJET DE LA CHARGE
!
            call lisllc(lischa, ichar, prefob)
!
! ------- TYPE DE LA CHARGE
!
            call lisltc(lischa, ichar, typech)
!
! ------- NOM DE LA CHARGE
!
            call lislch(lischa, ichar, nomch0)
!
! ------- CALCUL SI CHARGE EXISTANTE
!
            call lisopt(prefob, nomo, typech, indxch, option, &
                        parain, paraou, carte, ligcal)
!
! ------- CARTE D'ENTREE
!
            lpain(3) = parain
            lchin(3) = carte
!
! ------- CARTE DE SORTIE
!
            lpaout(1) = paraou
            call gcnco2(newnom)
            lchout(1) = '&&VEDIMD.'//newnom(2:8)
            call corich('E', lchout(1), ichin_=ichar)
!
! ------- CALCUL
!
            call calcul('S', option, ligcal, nbin, lchin, &
                        lpain, nbout, lchout, lpaout, 'V', &
                        'OUI')
!
! ------- RESU_ELEM DANS LE VECT_ELEM
!
            call exisd('CHAMP_GD', lchout(1), iret)
            ASSERT(iret .gt. 0)
            call reajre(vecele, lchout(1), 'V')
        end if
    end do
!
99  continue
!
    call jedetr(nomlis)
!
    call jedema()
end subroutine
