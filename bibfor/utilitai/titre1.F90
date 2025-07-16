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
subroutine titre1(st, nomobj, base, nbtitr, titdon, &
                  lgdon, formr, nomsym, iordr)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/titreb.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=1) :: st
    character(len=*) :: nomobj, base, titdon(*), formr
    integer(kind=8) :: nbtitr, lgdon(*)
    character(len=*), optional, intent(in) :: nomsym
    integer(kind=8), optional, intent(in) :: iordr
!                           MXLIGS MAX DE LIGNES EN SORTIE
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icold, icols, ideb, ierx, ilig, iligd, vali(2)
    integer(kind=8) :: iligs, ldons, ldons1, lonmax, lsort, mxcold, mxligs
    aster_logical :: lDefault
!
!-----------------------------------------------------------------------
    parameter(mxligs=50)
    character(len=1) :: kavant, kcoura
!
!     ------------------------------------------------------------------
    call jemarq()
    call wkvect('&&TITRE1.TAMPON.SORTIE', 'V V K80', mxligs, ldons)
    ldons1 = ldons
!
    lDefault = ASTER_TRUE
    if (nbtitr .lt. 0) then
        lDefault = ASTER_FALSE
        nbtitr = -nbtitr
    end if
!
! -------------------------
!     ICOLS  = INDICE DE LA DERNIERE COLONNE REMPLIE DANS LA SORTIE
!     ILIGS  = INDICE DE LA LIGNE COURANTE (TABLEAU DE SORTIE)
!     ICOLD  = INDICE DE LA DERNIERE COLONNE LUE DANS LA DONNEE
!     ILIGD  = INDICE DE LA LIGNE COURANTE (TABLEAU DE DONNEE)
! -------------------------
!
!     --- TANT QU'IL Y A DES LIGNES FAIRE ---
    icols = 0
    iligd = 1
!
    if (nbtitr .gt. mxligs) then
        vali(1) = mxligs
        vali(2) = nbtitr
        call utmess('A', 'UTILITAI4_89', ni=2, vali=vali)
        nbtitr = mxligs
    end if
1000 continue
    if (iligd .le. nbtitr .and. iligd .le. mxligs) then
!
!        --- TANT QU'IL Y A DES COLONNES FAIRE ---
        icold = 1
        mxcold = lgdon(iligd)
1100    continue
        if (icold .le. mxcold) then
            if (titdon(iligd) (icold:icold) .eq. '&' .and. lDefault) then
!
!              --- ON A TROUVE UN "&",  ATTENTION DEMON ---
!               uniquement pour les titres par défaut
                call titreb(titdon, iligd, icold, nbtitr, zk80(1), &
                            ldons1, icols, formr, nomsym, iordr)
!
            else
                icols = icols+1
                if (icols .gt. 80) then
!
!                 --- ON EVITE DE COUPER LES MOTS ---
                    icols = 80
                    kavant = zk80(ldons1) (icols:icols)
                    kcoura = titdon(iligd) (icold:icold)
                    if (kavant .ne. ' ' .and. kcoura .ne. ' ') then
200                     continue
                        icols = icols-1
                        if (icols .gt. 0) then
                            kavant = zk80(ldons1) (icols:icols)
                            if (kavant .ne. ' ') goto 200
                            ideb = icols+1
                            icols = 0
                            do i = ideb, 80
                                icols = icols+1
                                kavant = zk80(ldons1) (i:i)
                                zk80(ldons1+1) (icols:icols) = kavant
                                zk80(ldons1) (i:i) = ' '
                            end do
                            icols = icols+1
                            ldons1 = ldons1+1
                        else
                            ldons1 = ldons1+1
                            icols = 1
                        end if
                    else
                        ldons1 = ldons1+1
                        icols = 1
                    end if
                end if
                if (ldons1 .gt. ldons-1+mxligs) then
                    ldons1 = ldons1-1
                    goto 1200
                end if
                zk80(ldons1) (icols:icols) = titdon(iligd) (icold:icold)
                icold = icold+1
            end if
            goto 1100
        end if
        iligd = iligd+1
        goto 1000
    end if
1200 continue
!
!     --- RECOPIE DANS L'OBJET FINAL ----
    iligs = ldons1-ldons+1
!
    call jeexin(nomobj, ierx)
    if (ierx .eq. 0) then
        call wkvect(nomobj, base(1:1)//' V K80', iligs, lsort)
        lonmax = 0
    else if (st .eq. 'C') then
        call jelira(nomobj, 'LONMAX', lonmax)
        call juveca(nomobj, lonmax+iligs)
        call jeveuo(nomobj, 'E', lsort)
    else
        call jedetr(nomobj)
        call wkvect(nomobj, base(1:1)//' V K80', iligs, lsort)
        lonmax = 0
    end if
    do ilig = 1, iligs
        if (ilig .gt. mxligs) then
            goto 2001
        end if
        zk80(lsort+lonmax+ilig-1) = zk80(ldons+ilig-1)
    end do
2001 continue
!C
!C ----- DEBUG
!C
!C    IFM = IUNIFI('MESSAGE')
!C      WRITE(IFM,*) ' ---------------------------------------------- '
!C      WRITE(IFM,*) ' TITRE ATTACHE AU CONCEPT (PRODUIT)  ',NOMOBJ
!C      WRITE(IFM,*) ' ---------------------------------------------- '
!C      DO 3000 ILIG = 1, LONMAX+ILIGS
!C         WRITE(IFM,*) ZK80(LSORT+ILIG-1)
!3000   CONTINUE
!C      WRITE(IFM,*) ' '
!C      WRITE(IFM,*) ' '
!
!     --- MENAGE ---
    call jedetr('&&TITRE1.TAMPON.SORTIE')
    call jedema()
end subroutine
