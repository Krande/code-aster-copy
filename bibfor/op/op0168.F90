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

subroutine op0168()
    implicit none
!     OPERATEUR  EXTR_MODE
!     ------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/iunifi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsorac.h"
#include "asterfort/rsvpar.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/vpcrea.h"
#include "asterfort/vprecu.h"
#include "asterfort/vpstor.h"
#include "asterfort/wkvect.h"
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, ifr, impr, iord, iprec, iret
    integer(kind=8) :: j, jadr, jme, jnom, jor, jordr, tmod(1)
    integer(kind=8) :: k, lmod, lvali, lvalk, lvalr, n1
    integer(kind=8) :: n10, n2, n3, n4, n5, n6, n7
    integer(kind=8) :: n8, n9, nbfilt, nbme, nbmode, nbmodt, nbmodu
    integer(kind=8) :: nbmr, nbpara, nbpari, nbpark, nbparr, ndimt, neq
    integer(kind=8) :: npari, npark, nparr, nume, nume1, nume2
    real(kind=8) :: cumulx, cumuly, cumulz, dx, dy, dz, fremax
    real(kind=8) :: fremin, freq, seuil, undf
!-----------------------------------------------------------------------
    parameter(nbpari=1, nbparr=15, nbpark=3, nbpara=19)
    integer(kind=8) :: lpar(3)
    integer(kind=8) :: vali(2)
    integer(kind=8) :: nbvect
    real(kind=8) :: r8b, prec, zero, mastot, valr(7)
    character(len=1) :: typmod
    character(len=3) :: ouinon
    character(len=8) :: k8b, modeou, modein, modele, mate, cara
    character(len=16) :: typcon, nomcmd, critfi, nompar(3), nomsy, nompav
    character(len=16) :: nompa2(3)
    character(len=19) :: numedd
    character(len=24) :: masse, amor, raide, refd, massi, amori, raidi, kmode
    character(len=24) :: valk(3)
    character(len=24) :: kvec, kvali, kvalr, kvalk, nopara(nbpara)
    complex(kind=8) :: c16b
!     ------------------------------------------------------------------
    data refd/'                   .REFD'/
    data kvec/'&&OP0168.VAL_PROPRE'/
    data kvali/'&&OP0168.GRAN_MODAL_I'/
    data kvalr/'&&OP0168.GRAN_MODAL_R'/
    data kvalk/'&&OP0168.GRAN_MODAL_K_'/
    data nompar/'MASS_EFFE_UN_DX', 'MASS_EFFE_UN_DY',&
     &              'MASS_EFFE_UN_DZ'/
    data nompa2/'MASS_EFFE_DX', 'MASS_EFFE_DY',&
     &              'MASS_EFFE_DZ'/
    data nopara/'NUME_MODE',&
     &  'NORME', 'TYPE_MODE', 'NOEUD_CMP',&
     &  'FREQ', 'OMEGA2', 'AMOR_REDUIT',&
     &  'MASS_GENE', 'RIGI_GENE', 'AMOR_GENE',&
     &  'MASS_EFFE_DX', 'MASS_EFFE_DY', 'MASS_EFFE_DZ',&
     &  'FACT_PARTICI_DX', 'FACT_PARTICI_DY', 'FACT_PARTICI_DZ',&
     &  'MASS_EFFE_UN_DX', 'MASS_EFFE_UN_DY', 'MASS_EFFE_UN_DZ'/
!     ------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
!
    zero = 0.d0
    ouinon = ' '
    undf = r8vide()
    call getres(modeou, typcon, nomcmd)
    ifr = iunifi('RESULTAT')
!
!     --- RECUPERATION DU NOMBRE DE MODES A EXTRAIRE ---
    kmode = '&&OP0168.MODE.RETENU'
    nbmr = 0
    ndimt = 0
    call getfac('FILTRE_MODE', nbfilt)
    if (nbfilt .ne. 0) then
        call wkvect('&&OP0168.NOM_MODE', 'V V K8', nbfilt, jnom)
        call jecrec(kmode, 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                    nbfilt)
        do i = 1, nbfilt
            call getvid('FILTRE_MODE', 'MODE', iocc=i, scal=modein, nbret=n1)
            call jelira(modein//'           .ORDR', 'LONUTI', iret)
            if (iret .eq. 0) goto 10
            nbmr = nbmr+1
!
            if (nbmr .eq. 1) then
!      --- MATRICES DE REFERENCE DES MODES ---
                call dismoi('REF_RIGI_PREM', modein, 'RESU_DYNA', repk=raide, arret='C')
                call dismoi('REF_MASS_PREM', modein, 'RESU_DYNA', repk=masse, arret='C')
                call dismoi('REF_AMOR_PREM', modein, 'RESU_DYNA', repk=amor, arret='C')
                call dismoi('NUME_DDL', modein, 'RESU_DYNA', repk=numedd)
                call vpcrea(0, modeou, masse, amor, raide, &
                            numedd, ibid)
            end if
!
            zk8(jnom+nbmr-1) = modein
            call dismoi('REF_RIGI_PREM', modein, 'RESU_DYNA', repk=raidi, arret='C')
            call dismoi('REF_MASS_PREM', modein, 'RESU_DYNA', repk=massi, arret='C')
            call dismoi('REF_AMOR_PREM', modein, 'RESU_DYNA', repk=amori, arret='C')
            if (massi .ne. masse .or. amori .ne. amor .or. raidi .ne. raide) then
                call utmess('F', 'ALGELINE3_9')
            end if
!
            call rsorac(modein, 'LONUTI', 0, r8b, k8b, &
                        c16b, r8b, k8b, tmod, 1, &
                        ibid)
            nbmodt = tmod(1)
            call wkvect('&&OP0168.NUME_ORDRE', 'V V I', nbmodt, jor)
            call rsorac(modein, 'TOUT_ORDRE', 0, r8b, k8b, &
                        c16b, r8b, k8b, zi(jor), nbmodt, &
                        ibid)
!
            call jecroc(jexnum(kmode, nbmr))
            call jeecra(jexnum(kmode, nbmr), 'LONMAX', nbmodt)
            call jeveuo(jexnum(kmode, nbmr), 'E', jordr)
!
            call getvtx('FILTRE_MODE', 'TOUT_ORDRE', iocc=i, scal=ouinon, nbret=n1)
            if (n1 .ne. 0 .and. ouinon .eq. 'OUI') then
                nbmode = nbmodt
                do j = 1, nbmode
                    zi(jordr+j-1) = zi(jor+j-1)
                end do
                call jeecra(jexnum(kmode, nbmr), 'LONUTI', nbmode)
                ndimt = ndimt+nbmode
                call jedetr('&&OP0168.NUME_ORDRE')
                goto 10
            end if
!
            call getvis('FILTRE_MODE', 'NUME_ORDRE', iocc=i, nbval=0, nbret=n2)
            if (n2 .ne. 0) then
                nbmodu = -n2
                call wkvect('&&OP0168.NUME_MODE', 'V V I', nbmodu, jme)
                call getvis('FILTRE_MODE', 'NUME_ORDRE', iocc=i, nbval=nbmodu, vect=zi(jme), &
                            nbret=n2)
                nbmode = 0
                do j = 1, nbmodu
                    do k = 1, nbmodt
                        if (zi(jme+j-1) .eq. zi(jor+k-1)) then
                            nbmode = nbmode+1
                            zi(jordr+nbmode-1) = zi(jme+j-1)
                            goto 20
                        end if
                    end do
                    valk(1) = modein
                    vali(1) = zi(jme+j-1)
                    call utmess('A', 'ALGELINE4_55', sk=valk(1), si=vali(1))
20                  continue
                end do
                call jeecra(jexnum(kmode, nbmr), 'LONUTI', nbmode)
                ndimt = ndimt+nbmode
                call jedetr('&&OP0168.NUME_MODE')
                call jedetr('&&OP0168.NUME_ORDRE')
                goto 10
            end if
!
            call getvis('FILTRE_MODE', 'NUME_MODE', iocc=i, nbval=0, nbret=n3)
            if (n3 .ne. 0) then
                nbmodu = -n3
                call wkvect('&&OP0168.NUME_MODE', 'V V I', nbmodu, jme)
                call getvis('FILTRE_MODE', 'NUME_MODE', iocc=i, nbval=nbmodu, vect=zi(jme), &
                            nbret=n3)
                nbmode = 0
                do j = 1, nbmodt
                    iord = zi(jor+j-1)
                    call rsadpa(modein, 'L', 1, 'NUME_MODE', iord, &
                                0, sjv=jadr, styp=k8b, istop=0)
                    nume = zi(jadr)
                    do k = 1, nbmodu
                        if (nume .eq. zi(jme+k-1)) then
                            nbmode = nbmode+1
                            zi(jordr+nbmode-1) = iord
                        end if
                    end do
                end do
                call jeecra(jexnum(kmode, nbmr), 'LONUTI', nbmode)
                ndimt = ndimt+nbmode
                call jedetr('&&OP0168.NUME_MODE')
                call jedetr('&&OP0168.NUME_ORDRE')
                goto 10
            end if
!
            call getvis('FILTRE_MODE', 'NUME_MODE_EXCLU', iocc=i, nbval=0, nbret=n4)
            if (n4 .ne. 0) then
                nbme = -n4
                call wkvect('&&OP0168.NUME_MODE', 'V V I', nbme, jme)
                call getvis('FILTRE_MODE', 'NUME_MODE_EXCLU', iocc=i, nbval=nbme, vect=zi(jme), &
                            nbret=n4)
                nbmode = 0
                do j = 1, nbmodt
                    iord = zi(jor+j-1)
                    call rsadpa(modein, 'L', 1, 'NUME_MODE', iord, &
                                0, sjv=jadr, styp=k8b, istop=0)
                    nume = zi(jadr)
                    do k = 1, nbme
                        if (nume .eq. zi(jme+k-1)) goto 40
                    end do
                    nbmode = nbmode+1
                    zi(jordr+nbmode-1) = iord
40                  continue
                end do
                call jeecra(jexnum(kmode, nbmr), 'LONUTI', nbmode)
                ndimt = ndimt+nbmode
                call jedetr('&&OP0168.NUME_ORDRE')
                call jedetr('&&OP0168.NUME_MODE')
                goto 10
            end if
!
            call getvr8('FILTRE_MODE', 'FREQ_MIN', iocc=i, nbval=0, nbret=n5)
            if (n5 .ne. 0) then
                call getvr8('FILTRE_MODE', 'FREQ_MIN', iocc=i, scal=fremin, nbret=n5)
                call getvr8('FILTRE_MODE', 'FREQ_MAX', iocc=i, scal=fremax, nbret=n5)
                call getvr8('FILTRE_MODE', 'PRECISION', iocc=i, scal=prec, nbret=n5)
                fremin = fremin-prec
                fremax = fremax+prec
                nbmode = 0
                do j = 1, nbmodt
                    iord = zi(jor+j-1)
                    call rsadpa(modein, 'L', 1, 'FREQ', iord, &
                                0, sjv=jadr, styp=k8b, istop=0)
                    freq = zr(jadr)
                    if (freq .ge. fremin .and. freq .le. fremax) then
                        nbmode = nbmode+1
                        zi(jordr+nbmode-1) = iord
                    end if
                end do
                call jeecra(jexnum(kmode, nbmr), 'LONUTI', nbmode)
                ndimt = ndimt+nbmode
                call jedetr('&&OP0168.NUME_ORDRE')
                goto 10
            end if
!
            call getvtx('FILTRE_MODE', 'CRIT_EXTR', iocc=i, nbval=0, nbret=n6)
            if (n6 .ne. 0) then
                call getvtx('FILTRE_MODE', 'CRIT_EXTR', iocc=i, scal=critfi, nbret=n6)
                call getvr8('FILTRE_MODE', 'SEUIL', iocc=i, scal=seuil, nbret=n7)
                call getvr8('FILTRE_MODE', 'SEUIL_X', iocc=i, scal=seuil, nbret=n8)
                call getvr8('FILTRE_MODE', 'SEUIL_Y', iocc=i, scal=seuil, nbret=n9)
                call getvr8('FILTRE_MODE', 'SEUIL_Z', iocc=i, scal=seuil, nbret=n10)
                nbmode = 0
                if (critfi .eq. 'MASS_EFFE_UN' .and. typcon(1:9) .eq. 'MODE_MECA') then
                    do j = 1, nbmodt
                        iord = zi(jor+j-1)
                        call rsadpa(modein, 'L', 3, nompar, iord, &
                                    0, tjv=lpar, styp=k8b, istop=0)
                        dx = zr(lpar(1))
                        dy = zr(lpar(2))
                        dz = zr(lpar(3))
                        if (dx .eq. undf .or. dy .eq. undf .or. dz .eq. undf) then
                            call utmess('F', 'ALGELINE3_10')
                        end if
                        if (n7 .ne. 0) then
                            if (dx .ge. seuil .or. dy .ge. seuil .or. dz .ge. seuil) then
                                nbmode = nbmode+1
                                zi(jordr+nbmode-1) = iord
                            end if
                        else if (n8 .ne. 0) then
                            if (dx .ge. seuil) then
                                nbmode = nbmode+1
                                zi(jordr+nbmode-1) = iord
                            end if
                        else if (n9 .ne. 0) then
                            if (dy .ge. seuil) then
                                nbmode = nbmode+1
                                zi(jordr+nbmode-1) = iord
                            end if
                        else if (n10 .ne. 0) then
                            if (dz .ge. seuil) then
                                nbmode = nbmode+1
                                zi(jordr+nbmode-1) = iord
                            end if
                        end if
!
                    end do
                end if
                if (critfi .eq. 'MASS_EFFE' .and. typcon(1:9) .eq. 'MODE_MECA') then
                    do j = 1, nbmodt
                        iord = zi(jor+j-1)
                        call rsadpa(modein, 'L', 3, nompa2, iord, &
                                    0, tjv=lpar, styp=k8b, istop=0)
                        dx = zr(lpar(1))
                        dy = zr(lpar(2))
                        dz = zr(lpar(3))
                        if (dx .eq. undf .or. dy .eq. undf .or. dz .eq. undf) then
                            call utmess('F', 'ALGELINE3_10')
                        end if
                        if (n7 .ne. 0) then
                            if (dx .ge. seuil .or. dy .ge. seuil .or. dz .ge. seuil) then
                                nbmode = nbmode+1
                                zi(jordr+nbmode-1) = iord
                            end if
                        else if (n8 .ne. 0) then
                            if (dx .ge. seuil) then
                                nbmode = nbmode+1
                                zi(jordr+nbmode-1) = iord
                            end if
                        else if (n9 .ne. 0) then
                            if (dy .ge. seuil) then
                                nbmode = nbmode+1
                                zi(jordr+nbmode-1) = iord
                            end if
                        else if (n10 .ne. 0) then
                            if (dz .ge. seuil) then
                                nbmode = nbmode+1
                                zi(jordr+nbmode-1) = iord
                            end if
                        end if
!
                    end do
                end if
                if (critfi .eq. 'MASS_GENE') then
                    mastot = zero
                    do j = 1, nbmodt
                        iord = zi(jor+j-1)
                        nompav = 'MASS_GENE'
                        call rsadpa(modein, 'L', 1, nompav, iord, &
                                    0, sjv=lpar(1), styp=k8b, istop=0)
                        mastot = mastot+zr(lpar(1))
                    end do
                    do j = 1, nbmodt
                        iord = zi(jor+j-1)
                        nompav = 'MASS_GENE'
                        call rsadpa(modein, 'L', 1, nompav, iord, &
                                    0, sjv=lpar(1), styp=k8b, istop=0)
                        dx = zr(lpar(1))/mastot
                        if (dx .ge. seuil) then
                            nbmode = nbmode+1
                            zi(jordr+nbmode-1) = iord
                        end if
                    end do
                end if
                call jeecra(jexnum(kmode, nbmr), 'LONUTI', nbmode)
                ndimt = ndimt+nbmode
                call jedetr('&&OP0168.NUME_ORDRE')
                goto 10
            end if
10          continue
        end do
    end if
!
!     --- STOCKAGE ---
!
    call dismoi('MODELE', modein, 'RESULTAT', repk=modele)
    call dismoi('CARA_ELEM', modein, 'RESULTAT', repk=cara)
    call dismoi('CHAM_MATER', modein, 'RESULTAT', repk=mate)

    if (ndimt .eq. 0) then
        call utmess('F', 'ALGELINE3_11')
    end if
    call rscrsd('G', modeou, typcon, ndimt)
    iprec = 0
    nomsy = 'DEPL'
    do i = 1, nbmr
        call jemarq()
        call jerecu('V')
        modein = zk8(jnom+i-1)
        call jelira(jexnum(kmode, i), 'LONUTI', nbmode)
        if (nbmode .eq. 0) then
            valk(1) = modein
            call utmess('A', 'ALGELINE4_56', sk=valk(1))
            goto 102
        end if
        call jeveuo(jexnum(kmode, i), 'L', jordr)
        nbvect = nbmode
        call vprecu(modein, nomsy, nbvect, zi(jordr), kvec, &
                    nbpara, nopara(1), kvali, kvalr, kvalk, &
                    neq, nbmode, typmod, npari, nparr, &
                    npark)
        ASSERT(npari .eq. nbpari)
        ASSERT(nparr .eq. nbparr)
        ASSERT(npark .eq. nbpark)
        call jeveuo(kvec, 'L', lmod)
        call jeveuo(kvali, 'L', lvali)
        call jeveuo(kvalr, 'L', lvalr)
        call jeveuo(kvalk, 'L', lvalk)
        if (typmod .eq. 'R') then
            call vpstor(-1, typmod, modeou, nbmode, neq, &
                        zr(lmod), [c16b], nbmode, nbpari, nbparr, &
                        nbpark, nopara, '    ', zi(lvali), zr(lvalr), &
                        zk24(lvalk), iprec, modele, mate, cara)
        else if (typmod .eq. 'C') then
            call vpstor(-1, typmod, modeou, nbmode, neq, &
                        [r8b], zc(lmod), nbmode, nbpari, nbparr, &
                        nbpark, nopara, '    ', zi(lvali), zr(lvalr), &
                        zk24(lvalk), iprec, modele, mate, cara)
        else
            call utmess('F', 'ALGELINE2_44', sk=typmod)
        end if
        iprec = iprec+nbmode
        call jedetr(kvec)
        call jedetr(kvali)
        call jedetr(kvalr)
        call jedetr(kvalk)
102     continue
        call jedema()
    end do
!
!     --- ON ALARME SI NUME_MODE IDENTIQUE ---
!
    call rsorac(modeou, 'LONUTI', 0, r8b, k8b, &
                c16b, r8b, k8b, tmod, 1, &
                ibid)
    nbmode = tmod(1)
    call wkvect('&&OP0168.NUME_ORDRE', 'V V I', nbmode, jordr)
    call rsorac(modeou, 'TOUT_ORDRE', 0, r8b, k8b, &
                c16b, r8b, k8b, zi(jordr), nbmode, &
                ibid)
    do j = 1, nbmode
        iord = zi(jordr+j-1)
        call rsadpa(modeou, 'L', 1, 'NUME_MODE', iord, &
                    0, sjv=jadr, styp=k8b, istop=0)
        nume1 = zi(jadr)
        do k = j+1, nbmode
            iord = zi(jordr+k-1)
            call rsadpa(modeou, 'L', 1, 'NUME_MODE', iord, &
                        0, sjv=jadr, styp=k8b, istop=0)
            nume2 = zi(jadr)
            if (nume1 .eq. nume2) then
                vali(1) = ibid
                vali(2) = iord
                call utmess('A', 'ALGELINE4_57', ni=2, vali=vali)
            end if
        end do
    end do
!
!     --- LES IMPRESSIONS ---
!
    call getfac('IMPRESSION', impr)
    if (impr .ne. 0) then
        valk(1) = modeou
        valk(2) = typcon
        valk(3) = nomcmd
        call utmess('I', 'ALGELINE6_8', nk=3, valk=valk)
        call getvtx('IMPRESSION', 'CUMUL', iocc=1, scal=ouinon, nbret=n1)
        call getvtx('IMPRESSION', 'CRIT_EXTR', iocc=1, scal=critfi, nbret=n2)
!
        call rsvpar(modeou, 1, 'MASS_EFFE_UN_DX', ibid, undf, &
                    k8b, iret)
        call rsadpa(modeou, 'L', 1, 'MASS_EFFE_UN_DX', 1, &
                    0, sjv=jadr, styp=k8b, istop=0)
!
        if (iret .ne. 100 .and. critfi .eq. 'MASS_EFFE_UN' .and. typcon(1:9) .eq. &
            'MODE_MECA') then
            if (ouinon .eq. 'OUI') then
                call utmess('I', 'ALGELINE6_50')
            else
                call utmess('I', 'ALGELINE6_52')
            end if
            cumulx = 0.d0
            cumuly = 0.d0
            cumulz = 0.d0
            do j = 1, nbmode
                iord = zi(jordr+j-1)
                call rsadpa(modeou, 'L', 1, 'NUME_MODE', iord, &
                            0, sjv=jadr, styp=k8b, istop=0)
                nume = zi(jadr)
                call rsadpa(modeou, 'L', 1, 'FREQ', iord, &
                            0, sjv=jadr, styp=k8b, istop=0)
                freq = zr(jadr)
                call rsadpa(modeou, 'L', 3, nompar, iord, &
                            0, tjv=lpar, styp=k8b, istop=0)
                dx = zr(lpar(1))
                dy = zr(lpar(2))
                dz = zr(lpar(3))
                if (dx .eq. undf .or. dy .eq. undf .or. dz .eq. undf) then
                    call utmess('F', 'ALGELINE3_10')
                end if
                cumulx = cumulx+dx
                cumuly = cumuly+dy
                cumulz = cumulz+dz
                vali(1) = iord
                vali(2) = nume
                valr(1) = freq
                valr(2) = dx
                valr(3) = dy
                valr(4) = dz
                valr(5) = cumulx
                valr(6) = cumuly
                valr(7) = cumulz
                if (ouinon .eq. 'OUI') then
                    call utmess('I', 'ALGELINE6_51', ni=2, vali=vali, nr=7, &
                                valr=valr)
                else
                    call utmess('I', 'ALGELINE6_53', ni=2, vali=vali, nr=4, &
                                valr=valr)
                end if
            end do
        end if
        call rsvpar(modeou, 1, 'MASS_EFFE_DX', ibid, undf, &
                    k8b, iret)
        call rsadpa(modeou, 'L', 1, 'MASS_EFFE_DX', 1, &
                    0, sjv=jadr, styp=k8b, istop=0)
!
        if (iret .ne. 100 .and. critfi .eq. 'MASS_EFFE' .and. typcon(1:9) .eq. &
            'MODE_MECA') then
            if (ouinon .eq. 'OUI') then
                call utmess('I', 'ALGELINE6_50')
            else
                call utmess('I', 'ALGELINE6_52')
            end if
            cumulx = 0.d0
            cumuly = 0.d0
            cumulz = 0.d0
            do j = 1, nbmode
                iord = zi(jordr+j-1)
                call rsadpa(modeou, 'L', 1, 'NUME_MODE', iord, &
                            0, sjv=jadr, styp=k8b, istop=0)
                nume = zi(jadr)
                call rsadpa(modeou, 'L', 1, 'FREQ', iord, &
                            0, sjv=jadr, styp=k8b, istop=0)
                freq = zr(jadr)
                call rsadpa(modeou, 'L', 3, nompa2, iord, &
                            0, tjv=lpar, styp=k8b, istop=0)
                dx = zr(lpar(1))
                dy = zr(lpar(2))
                dz = zr(lpar(3))
                if (dx .eq. undf .or. dy .eq. undf .or. dz .eq. undf) then
                    call utmess('F', 'ALGELINE3_10')
                end if
                cumulx = cumulx+dx
                cumuly = cumuly+dy
                cumulz = cumulz+dz
                vali(1) = iord
                vali(2) = nume
                valr(1) = freq
                valr(2) = dx
                valr(3) = dy
                valr(4) = dz
                valr(5) = cumulx
                valr(6) = cumuly
                valr(7) = cumulz
                if (ouinon .eq. 'OUI') then
                    call utmess('I', 'ALGELINE6_51', ni=2, vali=vali, nr=7, &
                                valr=valr)
                else
                    call utmess('I', 'ALGELINE6_53', ni=2, vali=vali, nr=4, &
                                valr=valr)
                end if
            end do
        end if
        if (critfi .eq. 'MASS_GENE') then
            if (ouinon .eq. 'OUI') then
                call utmess('I', 'ALGELINE6_54')
            else
                call utmess('I', 'ALGELINE6_56')
            end if
            cumulx = 0.d0
            do j = 1, nbmode
                iord = zi(jordr+j-1)
                call rsadpa(modeou, 'L', 1, 'NUME_MODE', iord, &
                            0, sjv=jadr, styp=k8b, istop=0)
                nume = zi(jadr)
                call rsadpa(modeou, 'L', 1, 'FREQ', iord, &
                            0, sjv=jadr, styp=k8b, istop=0)
                freq = zr(jadr)
                nompav = 'MASS_GENE'
                call rsadpa(modeou, 'L', 1, nompav, iord, &
                            0, sjv=lpar(1), styp=k8b, istop=0)
                dx = zr(lpar(1))
                cumulx = cumulx+dx
                vali(1) = iord
                vali(2) = nume
                valr(1) = freq
                valr(2) = dx
                valr(3) = cumulx
                if (ouinon .eq. 'OUI') then
                    call utmess('I', 'ALGELINE6_55', ni=2, vali=vali, nr=3, &
                                valr=valr)
                else
                    call utmess('I', 'ALGELINE6_57', ni=2, vali=vali, nr=2, &
                                valr=valr)
                end if
            end do
        end if
    end if
!
    call titre()
!
    call jedema()
end subroutine
