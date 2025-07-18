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

subroutine op0081()
    implicit none
!
!     BUT:
!       OPERATEUR DE CALCUL DE MACRO-ELEMENT A PARTIR D'UNE BASE MODALE
!       ET DE MATRICES ASSEMBLEES
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
!
!      SORTIE :
!-------------
!
! ......................................................................
!
!
!
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterc/r8pi.h"
#include "asterfort/calamo.h"
#include "asterfort/calprc.h"
#include "asterfort/calpro.h"
#include "asterfort/comp81.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/impe81.h"
#include "asterfort/iner81.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedetr.h"
#include "asterfort/refe81.h"
#include "asterfort/remp81.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ioc, n1, nbvam, imod, nma
!
    real(kind=8) :: pi
!
    character(len=8) :: nomres, nomcon, nomope, mailla, basmod, blanc
    character(len=8) :: macrin
    character(len=19) :: raid, mass, amor, impe, typmat
    character(len=24) :: nommat
    integer(kind=8) :: nbmod, iocm, iocf, ioca, vali(2)
    integer(kind=8) :: lmass, lrigi, lamor, ldref, ldres
!
    data blanc/'        '/
!
!-----------------------------------------------------------------------
!
    call infmaj()
    call getres(nomres, nomcon, nomope)
!
    pi = r8pi()
!
! --- MACR_ELEM_DYNA OBTENU PAR MODELE NUMERIQUE OU EXPERIMENTAL
    call getfac('MODELE_MESURE', nbvam)

! --- MACR_ELEM_DYNA EXISTANT
    call getvid(' ', 'MACR_ELEM_DYNA', nbval=0, nbret=nma)
!
! --- RECUPERATION BASE MODALE, MATRICES ET CREATION .REFE
!     ET DETERMINATION OPTION DE CALCUL
!
    if (nma .eq. 0) then
        call refe81(nomres, basmod, raid, mass, amor, mailla)
    else
        basmod = blanc
        call getvid(' ', 'MACR_ELEM_DYNA', scal=macrin, nbret=nma)
        if (macrin .ne. nomres) call utmess('F', 'ASSEMBLA_21')
    end if
!
! --- CALCUL DES MATRICES PROJETEES (SI PAS DE MOT-CLE MODELE_MESURE)
!
    if (nbvam .eq. 0) then
        impe = blanc
        call getvid(' ', 'MATR_IMPE', scal=impe, nbret=n1)
        if (impe .ne. blanc) then
            call impe81(nomres, impe, basmod)
            typmat = 'MATR_ASSE_DEPL_R'
            goto 10
        end if
        call getvid(' ', 'MATR_IMPE_RIGI', scal=impe, nbret=n1)
        if (impe .ne. blanc) then
            call impe81(nomres, impe, basmod)
            typmat = 'MATR_ASSE_DEPL_R'
            goto 10
        end if
        call getvid(' ', 'MATR_IMPE_MASS', scal=impe, nbret=n1)
        if (impe .ne. blanc) then
            call impe81(nomres, impe, basmod)
            typmat = 'MATR_ASSE_DEPL_R'
            goto 10
        end if
        call getvid(' ', 'MATR_IMPE_AMOR', scal=impe, nbret=n1)
        if (impe .ne. blanc) then
            call impe81(nomres, impe, basmod)
            typmat = 'MATR_ASSE_DEPL_R'
            goto 10
        end if
!
        nommat = nomres//'.MAEL_RAID'
        call gettco(raid, typmat)
!
        if (typmat .eq. 'MATR_ASSE_DEPL_R') then
            call calpro(nommat, 'G', basmod, raid)
        else if (typmat .eq. 'MATR_ASSE_DEPL_C') then
            call calprc(nommat, 'G', basmod, raid)
        else
            call utmess('F', 'ALGORITH14_17', sk=typmat)
        end if
!
        nommat = nomres//'.MAEL_MASS'
        call calpro(nommat, 'G', basmod, mass)
!
        if (amor .ne. blanc) then
            nommat = nomres//'.MAEL_AMOR'
            call calpro(nommat, 'G', basmod, amor)
        else
            call getvr8(blanc, 'AMOR_REDUIT', iocc=1, nbval=0, nbret=ioc)
            if (ioc .lt. 0) then
                nommat = nomres//'.MAEL_AMOR'
                call calamo(nommat, 'G', basmod)
            end if
        end if
!
! ---   CALCUL DES FORCES D'INERTIES
!
        nommat = nomres//'.MAEL_INER'
        call iner81(nommat, 'G', basmod, mass)
10      continue
!
! ------------------------------------------------------------------- C
!
    else
!
! -----  CAS OU LES MATRICES SONT REMPLIES A LA MAIN
!
        call dismoi('NB_MODES_TOT', basmod, 'RESULTAT', repi=nbmod)
!
! ---    RECUPERATION DES VALEURS DE MASSES GENERALISEES ET VERIF :
! ---    LE NOMBRE DE VALEURS ENTREES = NOMBRE DE VECT DE LA BASE
!
! ---   MASSE GENERALISEE
        call wkvect('&&OP0081.MASS', 'V V R', nbmod, lmass)
        call getvr8('MODELE_MESURE', 'MASS_GENE', iocc=1, nbval=nbmod, vect=zr(lmass), &
                    nbret=iocm)
        if (iocm .ne. nbmod) then
            vali(1) = nbmod
            vali(2) = iocm
            if (iocm .lt. 0) vali(2) = -iocm
            call utmess('F', 'ALGORITH17_31', sk='MASSE_GENE', ni=2, vali=vali)
        end if
!
! ---   FREQUENCES PROPRES
        call wkvect('&&OP0081.RIGI', 'V V R', nbmod, lrigi)
        call getvr8('MODELE_MESURE ', 'FREQ', iocc=1, nbval=nbmod, vect=zr(lrigi), &
                    nbret=iocf)
!
        if (iocf .ne. nbmod) then
            vali(1) = nbmod
            vali(2) = iocf
            if (iocf .lt. 0) vali(2) = -iocm
            call utmess('F', 'ALGORITH17_31', sk='FREQ', ni=2, vali=vali)
        end if
!
! ---   AMORTISSEMENTS REDUITS
        call wkvect('&&OP0081.AMOR', 'V V R', nbmod, lamor)
        call getvr8('MODELE_MESURE', 'AMOR_REDUIT', iocc=1, nbval=nbmod, vect=zr(lamor), &
                    nbret=ioca)
        if (ioca .ne. 0 .and. ioca .ne. nbmod) then
            vali(1) = nbmod
            vali(2) = ioca
            if (iocf .lt. 0) vali(2) = -iocm
            call utmess('F', 'ALGORITH17_31', sk='FREQ', ni=2, vali=vali)
        end if
!
! ----- REMPLISSAGE
!
! ---   MASS_GENE : DIRECT
        call remp81(nomres//'.MAEL_MASS', lmass, basmod, nbmod)
!
! ---   RIGI_GENE : LES CALCULER A PARTIR DE MASSE_GENE ET FREQ
        do imod = 1, nbmod
            zr(lrigi-1+imod) = 4*pi**2*zr(lrigi-1+imod)**2*zr(lmass-1+imod)
        end do
        call remp81(nomres//'.MAEL_RAID', lrigi, basmod, nbmod)
!
! ---   AMOR_GENE : LES CALCULER APARTIR D'AMOR_REDUIT, MASSE ET RIGI
        if (ioca .ne. 0) then
            do imod = 1, nbmod
                zr(lamor-1+imod) = 2*zr(lamor-1+imod)*sqrt(zr(lrigi-1+imod)*zr(lmass-1+imod))
            end do
            call remp81(nomres//'.MAEL_AMOR', lamor, basmod, nbmod)
        end if
!
! ---   REMPLLISSAGE DES FORCES D'INERTIES
!
        call wkvect(nomres//'.MAEL_INER_REFE', 'G V K24', 2, ldref)
        zk24(ldref) = basmod
        zk24(ldref+1) = '        '
        call wkvect(nomres//'.MAEL_INER_VALE', 'G V R', 3*nbmod, ldres)
!
! ---   LE VEC DES MASSES EFFE EST REMPLI A 0, ON EMET UNE ALARME
!       (ON NE PEUT PAS LE REMPLIR, ON NE CONNAIT PAS LA MAT DE MASSE)
        call utmess('A', 'ALGORITH17_32', sk=' ')
        do imod = 1, 3*nbmod
            zr(ldres-1+imod) = 0.d0
        end do
!
!
    end if
!
! --- COMPATIBILITE AVEC SD MACR_ELEM_STAT
!
    if (nma .eq. 0) then
        call comp81(nomres, basmod, raid, mailla)
    end if
!
! --- MENAGE
    call jedetr('&&OP0081.MASS')
    call jedetr('&&OP0081.RIGI')
    call jedetr('&&OP0081.AMOR')
!
end subroutine
