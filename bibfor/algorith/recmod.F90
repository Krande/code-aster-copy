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
subroutine recmod(modmec, nbmode, nbamor, bande, tymmec, &
                  grdmod)
    implicit none
#include "jeveux.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveut.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbmode, nbamor
    real(kind=8) :: bande(2)
    character(len=8) :: modmec, tymmec
    character(len=16) :: grdmod
!
!  BUT: RECUPERER LES INFORMATIONS DE TYPE MODE DYNAMIQUE POUR
!        LE CALCUL DYNAMIQUE ALEATOIRE
!
! OUT : MODMEC : NOMBRE DE MODES DYNAMIQUES UTILISES DANS LE CALCUL
! OUT : NBMODE : NOMBRE DE MODES DYNAMIQUES UTILISES DANS LE CALCUL
! OUT : NBAMOR : NOMBRE D AMORTISSEMENTS MODAUX DONNES
! OUT : BANDE  : LES DEUX BORNES DE L INTERVALLE DE FREQUENCES
! OUT : TYMMEC : TPYE R OU C POUR LES MODES MECAS
! OUT : GRDMOD : TYPE DE GRANDEUR A RECUPERER DANS LES MODES DYN ET STA
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ibid, nbtrou, nbmod1(1), lnumor, ilmode, iad, imod1, im, iret
    integer(kind=8) :: iadrmg, ilamod, ilamor, na1
    real(kind=8) :: rbid, freq1, amunif
    complex(kind=8) :: c16b
    character(len=8) :: k8b
    character(len=24) :: nomcha
!     ------------------------------------------------------------------
!
    call jemarq()
!
!-----MODES RETENUS
!
    call getvid('BASE_MODALE', 'MODE_MECA', iocc=1, scal=modmec, nbret=ibid)
!
    call rsorac(modmec, 'LONUTI', ibid, rbid, k8b, &
                c16b, 0.0d0, k8b, nbmod1, 1, &
                nbtrou)
    call wkvect('&&RECMOD.NUMERO.ORDRE', 'V V I', nbmod1(1), lnumor)
    call rsorac(modmec, 'TOUT_ORDRE', ibid, rbid, k8b, &
                c16b, 0.0d0, k8b, zi(lnumor), nbmod1(1), &
                nbtrou)
!
    call getvis('BASE_MODALE', 'NUME_ORDRE', iocc=1, nbval=0, nbret=nbmode)
    nbmode = -nbmode
    if (nbmode .eq. 0) then
        call getvr8('BASE_MODALE', 'BANDE', iocc=1, nbval=2, vect=bande, &
                    nbret=ibid)
        call wkvect('&&OP0131.LISTEMODES', 'V V I', nbmod1(1), ilmode)
        do im = 1, nbmod1(1)
            imod1 = zi(lnumor+im-1)
            call rsadpa(modmec, 'L', 1, 'FREQ', imod1, &
                        0, sjv=iad, styp=k8b)
            freq1 = zr(iad)
            if ((freq1-bande(1))*(freq1-bande(2)) .le. 0.d0) then
                nbmode = nbmode+1
                zi(ilmode-1+nbmode) = imod1
            end if
        end do
        if (nbmode .eq. 0) then
            call utmess('F', 'ALGORITH10_31')
        end if
    else
        call wkvect('&&OP0131.LISTEMODES', 'V V I', nbmode, ilmode)
        call getvis('BASE_MODALE', 'NUME_ORDRE', iocc=1, nbval=nbmode, vect=zi(ilmode), &
                    nbret=ibid)
        do im = 1, nbmode
            if (zi(ilmode-1+im) .gt. nbmod1(1)) then
                call utmess('F', 'ALGORITH10_32')
            end if
        end do
    end if
!
!----AMORTISSEMENTS MODAUX RETENUS
!
    call wkvect('&&OP0131.LISTEAMOR', 'V V R8', nbmode, ilamor)
    call getvr8('BASE_MODALE', 'AMOR_REDUIT', iocc=1, nbval=0, nbret=na1)
    nbamor = -(na1)
    if (nbamor .ne. 0) then
        if (na1 .ne. 0) then
            call getvr8('BASE_MODALE', 'AMOR_REDUIT', iocc=1, nbval=nbmode, vect=zr(ilamor), &
                        nbret=na1)
        end if
    else
        call getvr8('BASE_MODALE', 'AMOR_UNIF', iocc=1, scal=amunif, nbret=ibid)
        do im = 1, nbmode
            zr(ilamor-1+im) = amunif
        end do
    end if
!
!------CONSITUTION DE LA LISTE DES ADRESSES DES MODES DYNAMIQUES
!
    grdmod = 'DEPL'
    call wkvect('&&OP0131.LISTADRMODE', 'V V I', nbmode, ilamod)
    do im = 1, nbmode
        imod1 = zi(ilmode+im-1)
        call rsexch('F', modmec, grdmod, imod1, nomcha, &
                    iret)
        call jeveut(nomcha(1:19)//'.VALE', 'L', zi(ilamod+im-1))
    end do
    call jelira(nomcha(1:19)//'.VALE', 'TYPE', cval=tymmec)
!
!-----RECUPERATION DE LA MASSE GENERALISEE
!
    call wkvect('&&OP0131.MASSEGENE', 'V V R8', nbmode, iadrmg)
    do im = 1, nbmode
        imod1 = zi(ilmode+im-1)
        call rsadpa(modmec, 'L', 1, 'MASS_GENE', imod1, &
                    0, sjv=iad, styp=k8b)
        zr(iadrmg+im-1) = zr(iad)
    end do
!
    call jedetr('&&RECMOD.NUMERO.ORDRE')
    call jedema()
end subroutine
