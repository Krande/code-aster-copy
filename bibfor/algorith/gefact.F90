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
subroutine gefact(duree, nominf)
!
! GENERATION DE FCT ALEATOIRES :
!        - PREPARATION (DISCRETISATION-PROLONGEMENT) DE L INTERSPECTRE
!            (EN FONCTION DES CARACTERISTIQUES DU TEMPOREL A GENERER)
!        - FACTORISATION DE L INTERSPECTRE
!
! ----------------------------------------------------------------------
!
!      IN   :
!             DUREE : DUREE DU TEMPOREL A SIMULER (NEGATIVE SI ELLE
!                     N'EST PAS UNE DONNEE)
!      OUT  :
!             NOMINF : NOM DE L'INTERSPECTRE FACTORISE
!
! ----------------------------------------------------------------------
!
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/facint.h"
#include "asterfort/folocx.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/ordis.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    real(kind=8) :: duree
!
!
! 0.2. ==> COMMUNS
!
! 0.3. ==> VARIABLES LOCALES
!
!
    integer(kind=8) :: l, n1, nbpoin, nbpini, iret, ier
    integer(kind=8) :: dim, long, dim2, dim3, dim4
    integer(kind=8) :: nnn, nbmr, nbpt1, nbpt2, longh
    integer(kind=8) :: nbfc, nbval, nbval1, indice
    integer(kind=8) :: i, ii, jj, j, k, kf, lval2, lvalc, ipas
    integer(kind=8) :: iinf, isup, ls, lr, ld, lu, lv, lw, ix, iy
    integer(kind=8) :: lval, lval1, lchdes, inuor, lnuor, jnuor
    integer(kind=8) :: vali
!
    real(kind=8) :: prec, r8b, epsi
    real(kind=8) :: freqi, freqf, freq, frinit, fmax, fmin
    real(kind=8) :: pmin, dfreq, difpas, dt
    real(kind=8) :: pas, pas1, resure, resuim, x1, y1, x2, y2
    real(kind=8) :: pui2, pui2d, pui3d
    real(kind=8) :: valr
!
    character(len=1) :: coli
    character(len=3) :: interp
    character(len=8) :: k8b, intesp
    character(len=16) :: k16bid, nomcmd, prolgd, abscisse
    character(len=19) :: k19bid, nominf, nomint
    character(len=24) :: chvale, chdesc, chnuor, nomobj
    character(len=24) :: chnumi, chnumj, chfreq, chval, chrefe
    aster_logical :: lfreqf, lfreqi, lnbpn, linter, lprem, diag
    integer(kind=8) :: i1, lnumi, lnumj, lfreq, nbfreq, lrefe
!
!     ----------------------------------------------------------------
!     --- INITIALISATION  ---
!
    call jemarq()
    indice = 1
!
    call getres(k19bid, k16bid, nomcmd)
!
!===============
! 1. LECTURE DES DONNEES LIEES A L INTERSPECTRE ET VERIFS
!===============
!
    call getvid(' ', 'INTE_SPEC', scal=nomint, nbret=l)
!
    call getvtx(' ', 'INTERPOL', nbval=2, vect=interp, nbret=n1)
    linter = (interp .eq. 'NON')
!
    call getvis(' ', 'NB_POIN', scal=nbpoin, nbret=l)
    lnbpn = l .ne. 0
    nbpini = nbpoin
!
!=====
!  1.1  RECUPARATION DES DIMENSIONS, DES NUMEROS D'ORDRE, ...
!=====
    intesp = nomint(1:8)
    chrefe = intesp//'.REFE'
    call jeveuo(chrefe, 'L', lrefe)
    abscisse = zk16(lrefe+2)
    if (abscisse .ne. 'FREQ') then
        call utmess('F', 'UTILITAI8_72')
    end if
!
    chnumi = intesp//'.NUMI'
    chnumj = intesp//'.NUMJ'
    chfreq = intesp//'.DISC'
    chval = intesp//'.VALE'
    call jeveuo(chnumi, 'L', lnumi)
    call jeveuo(chnumj, 'L', lnumj)
    call jeveuo(chfreq, 'L', lfreq)
    call jelira(chnumi, 'LONMAX', nbmr)
    call jelira(chfreq, 'LONMAX', nbfreq)
!
    nomobj = '&&GEFACT.TEMP.NUOR'
    call wkvect(nomobj, 'V V I', nbmr, jnuor)
    do i1 = 1, nbmr
        zi(jnuor-1+i1) = zi(lnumi-1+i1)
    end do
    call ordis(zi(jnuor), nbmr)
    call wkvect('&&GEFACT.MODE', 'V V I', nbmr, inuor)
    nnn = 1
    zi(inuor) = zi(jnuor)
    do i = 2, nbmr
        if (zi(jnuor+i-1) .eq. zi(inuor+nnn-1)) goto 10
        nnn = nnn+1
        zi(inuor+nnn-1) = zi(jnuor+i-1)
10      continue
    end do
!
!=====
! 1.2 DONNEES ECHANTILLONNAGE FREQUENTIEL, DEDUCTION DES DONNEES
!     MANQUANTES (POUR RESPECTER ENTRE AUTRE LE TH. DE SHANNON)
!=====
    dim = nnn
    nbfc = dim*(dim+1)/2
    ASSERT(nbfc .eq. nbmr)
!
! 1.2.1 CAS OU ON UTILISE LA DISCRETISATION DE L INTERSPECTRE :
!     VERIFICATION DE LA COHERENCE DE LA DISCRETISATION DES FONCTIONS
!     DANS LE CAS OU CETTE DISCRETISATION EST CONSERVEE
    if (linter) then
        nbval = nbfreq
        pas = (zr(lfreq+nbval-1)-zr(lfreq))/(nbval-1)
        prec = 1.d-06
        do ii = 1, nbval-1
            pas1 = zr(lfreq+ii)-zr(lfreq+ii-1)
            difpas = abs(pas1-pas)
            if (difpas .gt. prec) then
                call utmess('F', 'ALGORITH3_78')
            end if
        end do
!
        if ((lnbpn) .and. (nbpini .lt. nbval)) then
            freqf = zr(lfreq+nbpini-1)
            valr = freqf
            call utmess('A', 'ALGORITH15_11', sr=valr)
        else
            freqf = zr(lfreq+nbval-1)
        end if
        dfreq = pas
        duree = 1.d0/dfreq
        freqi = zr(lfreq)
        frinit = mod(freqi, dfreq)
!
        nbpoin = 2**(int(log(freqf/dfreq)/log(2.d0))+1)
        if (lnbpn) then
            if (nbpoin .gt. nbpini) then
                vali = nbpoin
                r8b = 0.d0
                call utmess('A', 'ALGORITH15_12', si=vali)
            else
                pui2 = log(dble(nbpini))/log(2.d0)
                pui2d = abs(pui2-aint(pui2))
                pui3d = abs(1.d0-pui2d)
                if (pui2d .ge. 1.d-06 .and. pui3d .ge. 1.d-06) then
                    nbpini = 2**(int(pui2)+1)
                    call utmess('A', 'ALGORITH3_80')
                end if
                nbpoin = nbpini
            end if
        end if
!
    else
! 1.2.2 CAS OU ON PEUT INTERPOLER L INTERSPECTRE
!
        call getvr8(' ', 'FREQ_FIN', scal=freqf, nbret=l)
        lfreqf = l .ne. 0
!
        call getvr8(' ', 'FREQ_INIT', scal=freqi, nbret=l)
        lfreqi = l .ne. 0
!
!      RECHERCHE DES FREQUENCES MIN ET MAX ET DU PAS EN FREQUENCE MIN
!      DE L INTERSPECTRE
        pmin = 1.d+10
        fmax = 0.d0
        fmin = 1.d+10
        nbval = nbfreq
        do j = 1, nbval-1
            pas = abs(zr(lfreq+j)-zr(lfreq+j-1))
            if (pas .lt. pmin) pmin = pas
        end do
        freq = zr(lfreq+nbval-1)
        if (freq .gt. fmax) fmax = freq
        freq = zr(lfreq)
        if (freq .lt. fmin) fmin = freq
!
        if (.not. lfreqf) freqf = fmax
        if (.not. lfreqi) freqi = fmin
!
!     DETERMINATION DES PARAMETRES DE L ALGO.
        if (duree .gt. 0.d0) then
!     LA DUREE EST UNE DONNEE
            dfreq = 1.d0/duree
            if (lnbpn) then
                dt = duree/nbpoin/2.d0
                if (1.d0/dt .lt. 2.d0*freqf) then
                    nbpoin = 2**(int(log(2.d0*freqf*duree)/log(2.d0)))
                    vali = nbpoin
                    r8b = 0.d0
                    call utmess('A', 'ALGORITH15_13', si=vali)
                end if
            else
                nbpoin = 2**(int(log(2.d0*freqf*duree)/log(2.d0)))
                if ((dfreq .gt. 2*pmin) .and. (pmin .gt. 0.d0)) then
                    valr = 1.d0/pmin
                    call utmess('A', 'ALGORITH15_14', sr=valr)
                end if
            end if
        else
!     LA DUREE EST UNE INCONNUE
            if (lnbpn) then
                pui2 = log(dble(nbpoin))/log(2.d0)
                pui2d = abs(pui2-aint(pui2))
                pui3d = abs(1.d0-pui2d)
                if (pui2d .ge. 1.d-06 .and. pui3d .ge. 1.d-06) then
                    call utmess('A', 'ALGORITH3_80')
                    nbpoin = 2**(int(pui2)+1)
                end if
                dfreq = freqf/(nbpoin-1)
                frinit = freqf-dble(nbpoin)*dfreq
                if (frinit .lt. 0.d0) frinit = 0.d0
                if ((dfreq .gt. pmin) .and. (pmin .gt. 0.d0)) then
                    vali = nbpoin
                    valr = (freqf-freqi)/pmin+1
                    call utmess('A', 'ALGORITH15_15', si=vali, sr=valr)
                end if
            else
                if (pmin .gt. 0.d0) then
                    dfreq = pmin
                    nbpoin = 2**(int(log(2.d0*freqf/pmin)/log(2.d0)))
                    if (nbpoin .lt. 256) nbpoin = 256
                else
                    nbpoin = 256
                    dfreq = (freqf-freqi)/dble(nbpoin-1)
                end if
            end if
            duree = 1.d0/dfreq
        end if
!
        if (dble(nbpoin-1)*dfreq .gt. (freqf-freqi)) then
            frinit = freqf-dble(nbpoin)*dfreq
            if (frinit .lt. 0.d0) frinit = 0.d0
        else
            frinit = freqi
        end if
    end if
!
!
!===============
! 3. LECTURE DES VALEURS DES FONCTIONS ET/OU INTERPOLATION-PROLONGEMENT
!===============
!
    nbpt1 = nbpoin
    nbpt2 = nbpoin*2
    long = nbfc*nbpt2+nbpt1
    longh = dim*dim*nbpt2+nbpt1
!
!     --- CREATION D'UN VECTEUR TEMP.VALE POUR STOCKER LES VALEURS
!           DES FONCTIONS  ---
!
    call wkvect('&&GEFACT.TEMP.VALE', 'V V R', long, lval)
!
!
!     --- ON STOCKE LES FREQUENCES ET ON RECHERCHE LES INDICES AU
!         DE LA DES QUELLES LA MATRICE EST NULLE---
    lprem = .true.
    iinf = 0
    isup = nbpt1+1
    do k = 1, nbpt1
        freq = frinit+(k-1)*dfreq
        zr(lval+k-1) = freq
        if (freq .lt. freqi) iinf = k
        if ((freq .gt. freqf) .and. lprem) then
            isup = k
            lprem = .false.
        end if
    end do
!     ------------------------------------------------------------------
!     --- CHANGER LA FREQ INIT. A 0 HZ POUR LE CAS SANS INTERPOL
    if (linter) zr(lval) = 0.d0
!     ------------------------------------------------------------------
    lval1 = lval+nbpt1
!
!     --- POUR CHAQUE FONCTION CALCUL DE X,Y POUR CHAQUE FREQ.
!     (ON PROLONGE PAR 0 EN DEHORS DE (FREQI,FREQF)), PUIS ON STOCKE ---
    do kf = 1, nbfc
        call jeveuo(jexnum(chval, kf), 'L', lval2)
        diag = .false.
        if (zi(lnumi-1+kf) .eq. zi(lnumj-1+kf)) diag = .true.
!
        k8b = ' '
        do ipas = 1, nbpt1
            freq = frinit+(ipas-1)*dfreq
            ix = lval1+(kf-1)*nbpt2+ipas-1
            iy = lval1+(kf-1)*nbpt2+ipas-1+nbpt1
            if ((ipas .le. iinf) .or. (ipas .ge. isup)) then
                zr(ix) = 0.d0
                zr(iy) = 0.d0
            else
                if (linter) then
                    if (diag) then
                        resure = zr(lval2+ipas-iinf-1)
                        resuim = 0.d0
                    else
                        resure = zr(lval2+2*(ipas-iinf-1))
                        resuim = zr(lval2+2*(ipas-iinf-1)+1)
                    end if
                else
! ON INTERPOLLE
                    prolgd = 'CC      '
                    epsi = sqrt(r8prem())
                    call folocx(zr(lfreq), nbfreq, freq, prolgd, indice, &
                                epsi, coli, ier)
                    if (coli .eq. 'C') then
                        if (diag) then
                            resure = zr(lval2+indice-1)
                            resuim = 0.d0
                        else
                            resure = zr(lval2+2*(indice-1))
                            resuim = zr(lval2+2*(indice-1)+1)
                        end if
                    else if ((coli .eq. 'I') .or. (coli .eq. 'E')) then
                        x1 = zr(lfreq+indice-1)
                        x2 = zr(lfreq+indice)
                        if (diag) then
                            y1 = zr(lval2+indice-1)
                            y2 = zr(lval2+indice)
                            resure = y1+(freq-x1)*(y2-y1)/(x2-x1)
                            resuim = 0.d0
                        else
                            y1 = zr(lval2+2*(indice-1))
                            y2 = zr(lval2+2*indice)
                            resure = y1+(freq-x1)*(y2-y1)/(x2-x1)
                            y1 = zr(lval2+2*(indice-1)+1)
                            y2 = zr(lval2+2*indice+1)
                            resuim = y1+(freq-x1)*(y2-y1)/(x2-x1)
                        end if
                    else
                        call utmess('A', 'PREPOST3_6', sk=coli)
                    end if
                end if
                zr(ix) = resure
                zr(iy) = resuim
            end if
        end do
    end do
    nbval1 = nbpt1
!
!===============
! 4. FACTORISATION DES MATRICES INTERSPECTRALES (UNE PAR FREQ.)
!===============
!
!               1234567890123456789
    nominf = '&&INTESPECFACT     '
!
!     --- CREATION DE L'OBJET NOMINF//'.VALE'
    chvale = nominf//'.VALE'
    call wkvect(chvale, 'V V R', longh, lvalc)
!
!     --- CREATION DE L'OBJET NOMINF//'.DESC'
    chdesc = nominf//'.DESC'
    call wkvect(chdesc, 'V V I', 3, lchdes)
    zi(lchdes) = nbval1
    zi(lchdes+1) = dim
    zi(lchdes+2) = dim*dim
!
!     --- CREATION DE L'OBJET NOMINF//'.NUOR'
    chnuor = nominf//'.NUOR'
    call wkvect(chnuor, 'V V I', dim, lnuor)
    call jeveuo('&&GEFACT.MODE', 'L', inuor)
    do i = 1, dim
        zi(lnuor-1+i) = zi(inuor-1+i)
    end do
!
    dim2 = dim*dim
    dim3 = dim2+dim
    dim4 = 2*dim
    call wkvect('&&GEFACT.TEMP.VALS', 'V V C', dim2, ls)
    call wkvect('&&GEFACT.TEMP.VALR', 'V V C', dim2, lr)
    call wkvect('&&GEFACT.TEMP.VALD', 'V V R', dim, ld)
    call wkvect('&&GEFACT.TEMP.VALU', 'V V C', dim2, lu)
    call wkvect('&&GEFACT.TEMP.VALV', 'V V R', dim3, lv)
    call wkvect('&&GEFACT.TEMP.VALW', 'V V C', dim4, lw)
!
    call facint(nbval1, dim, longh, zr(lval), zr(lvalc), &
                long, zc(ls), zc(lr), zr(ld), zc(lu), &
                zr(lv), zc(lw))
!
    nbpt1 = nbval1
    do jj = 1, nbpt1
        zr(lvalc+jj-1) = zr(lval+jj-1)
    end do
!
    call titre()
!
    call jedetr('&&GEFACT.MODE')
    call jedetr(nomobj)
    call jeexin('&&GEFACT.FONCTION', iret)
    if (iret .ne. 0) call jedetr('&&GEFACT.FONCTION')
    call jedetr('&&GEFACT.TEMP.VALE')
    call jedetr('&&GEFACT.TEMP.VALD')
    call jedetr('&&GEFACT.TEMP.VALR')
    call jedetr('&&GEFACT.TEMP.VALS')
    call jedetr('&&GEFACT.TEMP.VALU')
    call jedetr('&&GEFACT.TEMP.VALV')
    call jedetr('&&GEFACT.TEMP.VALW')
    call jedetr('&&GEFACT.TEMP.FONC')
    call jedema()
end subroutine
