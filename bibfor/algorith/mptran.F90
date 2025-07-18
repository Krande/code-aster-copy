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

subroutine mptran(nombas, nommes, nbmesu, nbmode, basepr, &
                  vnoeud, vrange, vcham)
!
!
!     PROJ_MESU_MODAL : CALCUL DES CONTRIBUTIONS MODALES ET CONSTRUCTION
!                       DU TRAN_GENE OU HARM_GENE
!
!     IN  : NOMBAS : NOM DE LA BASE DE PROJECTION
!     IN  : NOMMES : NOM DE LA SD MESURE
!     IN  : NBMESU : NOMBRE DE DDL DE MESURE
!     IN  : NBMODE : NOMBRE DE VECTEURS DE BASE
!     IN  : BASEPR : NOM BASE PROJETEE SUIVANT DIRECTION MESURE
!     IN  : VNOEUD : NOM RANGEMENT NOEUD MESURE
!     IN  : VRANGE : NOM CORRESPONDANCE ORIENTATION SUIVANT LNOEUD
!
    implicit none
!
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/cnocns.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mdallo.h"
#include "asterfort/mdallr.h"
#include "asterfort/mpinv2.h"
#include "asterfort/mpinvc.h"
#include "asterfort/mpinvr.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/scalai.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nomres, nombas, nommes
    character(len=24) :: vrange, vnoeud, basepr, vcham
    integer(kind=8) :: nbmesu, nbmode, jpara, iexi
!
    character(len=1) :: typval
    character(len=8) :: k8bid, k8b, scal, kcmp, kreg
    character(len=8) :: modele, chmat, carael
    character(len=16) :: nomcmd, typres, k16bid, nomcha, kcham
    character(len=19) :: chs, chamno, sd2
!
    character(len=24) :: vabs, vmes, typba, raide
!
    aster_logical :: lfonct, zcmplx
!
    integer(kind=8) :: i, j, jabs, tmod(1)
    integer(kind=8) :: jdep, jvit, jacc, jpass, jordr, lord, imes, iret, gd
    integer(kind=8) :: labs, lmesu, lcoef, lred, jcnsd, jcnsc, jcnsv, n1
    integer(kind=8) :: ncoef, nfonc, lfonc, null, ibid, jcnsl, nbcmp
    integer(kind=8) :: lvale, lonmax, iocc, numord, ino, icmp, indice
    integer(kind=8) :: jcnsk, lrange, lnoeud, nbabs, jord, nbord
    integer(kind=8) :: jbasm, lcham, nbcham, ich, lch, jpames
!
    real(kind=8) :: r8bid, dt, pas, diff
!
    complex(kind=8) :: cbid
!
! ----------------------------------------------------------------------
!
    data nomcmd/'&PROJ_MESU_MODAL'/, k8b/'        '/, null/0/
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    kreg = ' '
!
! RECUPERATION DU NOM DU CONCEPT RESULTAT
    call getres(nomres, typres, k16bid)
! RECUPERATION DU CHAMP MESURE : NOMMES
    call getvtx('MODELE_MESURE', 'NOM_CHAM', iocc=1, nbval=0, nbret=nbcham)
    if (nbcham .ne. 0) then
        nbcham = -nbcham
    else
        call utmess('A', 'ALGORITH10_93')
    end if
!
    call wkvect('&&LISTE_CH', 'V V K16', nbcham, lch)
    call getvtx('MODELE_MESURE', 'NOM_CHAM', iocc=1, nbval=nbcham, vect=zk16(lch), &
                nbret=ibid)
!
    call getvtx('RESOLUTION', 'REGUL', iocc=1, scal=kreg, nbret=n1)
!
! RECUPERATION DU NOMBRE D ABSCISSES : NBABS
    call rsorac(nommes, 'LONUTI', 0, r8bid, k8bid, &
                cbid, r8bid, 'ABSOLU', tmod, 1, &
                ibid)
    nbabs = tmod(1)
!
    vabs = '&&ABSCISSES'
    call wkvect(vabs, 'V V R', nbabs, labs)
!
    vmes = '&&MESURE'
!
    call jeveuo(vrange, 'L', lrange)
    call jeveuo(vnoeud, 'L', lnoeud)
    call jeveuo(vcham, 'L', lcham)
!
! RECUPERATION ADRESSE DES NUMEROS D'ORDRE ET DU NOM SYMBOLIQUE
!
    call jeveuo(nommes//'           .ORDR', 'L', lord)
!
    chs = '&&MESURE.CHS'
!
! BOUCLE SUR LES CHAMPS
    do ich = 1, nbcham
        nomcha = zk16(lch-1+ich)
!
! BOUCLE SUR LES NUMEROS ORDRE
!
        do numord = 1, nbabs
!        -> EXISTENCE DES CHAMPS DANS LA STRUCTURE DE DONNEES MESURE
            call rsexch(' ', nommes, nomcha, zi(lord-1+numord), chamno, &
                        iret)
            if ((numord .le. 1) .and. (ich .eq. 1)) then
                call dismoi("NUM_GD", chamno, "CHAM_NO", repi=gd)
                scal = scalai(gd)
                typval = scal(1:1)
                if (typval .eq. 'C') then
                    zcmplx = .true.
                    call wkvect(vmes, 'V V C', nbmesu*nbabs, lmesu)
                else
                    zcmplx = .false.
                    call wkvect(vmes, 'V V R', nbmesu*nbabs, lmesu)
                end if
            end if
!
! RECUPERATION DE L ABSCISSE
            if ((typres(1:9) .eq. 'MODE_GENE') .or. (typres(1:9) .eq. 'HARM_GENE')) then
                call rsadpa(nommes, 'L', 1, 'FREQ', numord, &
                            0, sjv=jabs, styp=k8bid)
            else if (typres(1:9) .eq. 'TRAN_GENE') then
                call rsadpa(nommes, 'L', 1, 'INST', numord, &
                            0, sjv=jabs, styp=k8bid)
            end if
            zr(labs-1+numord) = zr(jabs)
!
! TRANSFORMATION DE CHAMNO EN CHAM_NO_S : CHS
            call detrsd('CHAM_NO_S', chs)
            call cnocns(chamno, 'V', chs)
            call jeveuo(chs//'.CNSK', 'L', jcnsk)
            call jeveuo(chs//'.CNSD', 'L', jcnsd)
            call jeveuo(chs//'.CNSC', 'L', jcnsc)
            call jeveuo(chs//'.CNSV', 'L', jcnsv)
            call jeveuo(chs//'.CNSL', 'L', jcnsl)
!
            nbcmp = zi(jcnsd-1+2)
!
            do imes = 1, nbmesu
                ino = zi(lnoeud-1+imes)
                kcmp = zk8(lrange-1+imes)
                kcham = zk16(lcham-1+imes)
                do icmp = 1, nbcmp
                    indice = (ino-1)*nbcmp+icmp
                    if ((zk8(jcnsc-1+icmp) .eq. kcmp) .and. (nomcha .eq. kcham)) then
                        if (zcmplx) then
                            zc(lmesu-1+(numord-1)*nbmesu+imes) = &
                                zc(jcnsv-1+indice)
                        else
                            zr(lmesu-1+(numord-1)*nbmesu+imes) = &
                                zr(jcnsv-1+indice)
                        end if
                    end if
                end do
            end do
!
! FIN BOUCLE SUR NUMERO ORDRE
        end do
!
! FIN BOUCLE SUR LES CHAMPS
    end do
!
! GESTION PARAMETRES DE REGULARISATION
    call getvr8('RESOLUTION', 'COEF_PONDER', iocc=1, nbval=0, nbret=ncoef)
    call getvid('RESOLUTION', 'COEF_PONDER_F', iocc=1, nbval=0, nbret=nfonc)
    iocc = abs(ncoef)+abs(nfonc)
    if ((ncoef .eq. 0) .and. (nfonc .eq. 0)) iocc = 0
!
    if ((iocc .eq. 0) .or. (kreg .eq. 'NON')) then
! CAS SANS REGULARISATION : PAR DEFAUT
        lfonct = .false.
        call wkvect(nomcmd//'.PONDER', 'V V R', nbmode, lcoef)
        do i = 1, nbmode
            zr(lcoef-1+i) = 0.d0
        end do
    else
        call getvr8('RESOLUTION', 'COEF_PONDER', iocc=1, nbval=0, nbret=ncoef)
        if (-ncoef .gt. 0) then
! CAS DE REGULARISATION SOUS FORME DE LISTE DE REELS
            lfonct = .false.
            if (-ncoef .gt. nbmode) then
                call utmess('F', 'ALGORITH6_27')
            end if
            if (-ncoef .gt. 0) then
                call wkvect(nomcmd//'.PONDER', 'V V R', nbmode, lcoef)
                call getvr8('RESOLUTION', 'COEF_PONDER', iocc=1, nbval=-ncoef, vect=zr(lcoef), &
                            nbret=ncoef)
            end if
            if (ncoef .lt. nbmode) then
                call utmess('I', 'ALGORITH6_28')
                do i = ncoef+1, nbmode
                    zr(lcoef-1+i) = zr(lcoef-1+ncoef)
                end do
            end if
        else
! CAS DE REGULARISATION SOUS FORME DE LISTE DE FONCTIONS
            lfonct = .true.
            call getvid('RESOLUTION', 'COEF_PONDER_F', iocc=1, nbval=0, nbret=nfonc)
            if (-nfonc .gt. nbmode) then
                call utmess('F', 'ALGORITH6_29')
            end if
            if (-nfonc .gt. 0) then
                call wkvect(nomcmd//'.FONC', 'V V K8', nbmode, lfonc)
                call getvid('RESOLUTION', 'COEF_PONDER_F', iocc=1, nbval=-nfonc, &
                            vect=zk8(lfonc), nbret=nfonc)
            end if
            if (nfonc .gt. 0 .and. nfonc .lt. nbmode) then
                call utmess('I', 'ALGORITH6_30')
                do i = nfonc+1, nbmode
                    zk8(lfonc-1+i) = zk8(lfonc-1+nfonc)
                end do
            end if
            call wkvect(nomcmd//'.PONDER', 'V V R', nbmode*nbabs, lcoef)
            do i = 1, nbmode
                call jelira(zk8(lfonc-1+i)//'           .VALE', 'LONMAX', lonmax)
                if (lonmax .ne. 2*nbabs) then
                    call utmess('F', 'ALGORITH6_31')
                end if
!
                call jeveuo(zk8(lfonc-1+i)//'           .VALE', 'L', lvale)
                do j = 1, nbabs
                    diff = zr(lvale-1+j)-zr(labs-1+j)
                    if (j .eq. 1) then
                        pas = zr(labs+1)-zr(labs)
                    else
                        pas = zr(labs-1+j)-zr(labs-1+j-1)
                    end if
                    if (abs(diff) .gt. pas*1.d-4) then
                        call utmess('F', 'ALGORITH6_32')
                    end if
!
                    zr(lcoef-1+(j-1)*nbmode+i) = zr(lvale-1+(lonmax/2)+j)
                end do
            end do
        end if
! FIN TEST SUR TYPE DE PONDERATION : REELS / LISTE DE FONCTIONS
    end if
! FIN GESTION PARAMETRES DE REGULARISATION
!
!
! INITIALISATION POUR ALLOCATION DU TRAN_GENE
!
    if (typres(1:9) .eq. 'TRAN_GENE') then
        dt = (zr(labs-1+nbabs)-zr(labs))/nbabs
    end if
!
! RECUPERATION DE LA MATRICE MODALE PROJETEE
!
    call jeveuo(basepr, 'L', lred)
!
! ALLOCATION DE TRAN_GENE OU HARM_GENE ET RESOLUTION DU SYSTEME
!
    if (.not. zcmplx) then
! SECOND MEMBRE REEL
        if (typres(1:9) .eq. 'HARM_GENE') then
            call utmess('F', 'ALGORITH6_33')
        end if
        if (typres(1:9) .eq. 'TRAN_GENE') then
! ALLOCATION
            call mdallo(nomres, 'TRAN', nbabs, sauve='GLOB', base=nombas, &
                        nbmodes=nbmode, jordr=jordr, jdisc=jabs, jdepl=jdep, jvite=jvit, &
                        jacce=jacc, dt=dt, jptem=jpass)
! RESOLUTION
            call mpinv2(nbmesu, nbmode, nbabs, zr(lred), zr(lmesu), &
                        zr(lcoef), zr(labs), lfonct, zr(jdep), zr(jvit), &
                        zr(jacc))
!
        else if (typres(1:9) .eq. 'MODE_GENE') then
            call wkvect(nomcmd//'.RETA', 'V V R', nbmode*nbabs, jdep)
            call mpinvr(nbmesu, nbmode, nbabs, zr(lred), zr(lmesu), &
                        zr(lcoef), zr(labs), lfonct, zr(jdep))
!
            call rscrsd('G', nomres, 'MODE_GENE', nbabs)
            call mdallr(nommes, nomres, nombas, nbmode, nbabs, &
                        zr(jdep), [cbid], zcmplx)
        end if
!
    else
! SECOND MEMBRE COMPLEXE
        if (typres(1:9) .eq. 'HARM_GENE') then
!
! ALLOCATION
!         -- DANS PROJ_MESU_MODAL ON REMPLIT TOUJOURS LES TROIS
!            CHAMPS, PEU IMPORTE LE TYPE DE MESURE FOURNI
            call mdallo(nomres, 'HARM', nbabs, sauve='GLOB', base=nombas, &
                        nbmodes=nbmode, jordr=jordr, jdisc=jabs, jdepl=jdep, jvite=jvit, &
                        jacce=jacc)
! RESOLUTION
            call mpinvc(nbmesu, nbmode, nbabs, zr(lred), zc(lmesu), &
                        zr(lcoef), zr(labs), lfonct, zc(jdep), zc(jvit), &
                        zc(jacc))
!
!
        else if (typres(1:9) .eq. 'MODE_GENE') then
! ALLOCATION
            call wkvect(nomcmd//'.RETA', 'V V C', nbmode*nbabs, jdep)
            call wkvect(nomcmd//'.RET1', 'V V C', nbmode*nbabs, jvit)
            call wkvect(nomcmd//'.RET2', 'V V C', nbmode*nbabs, jacc)
! RESOLUTION
            call mpinvc(nbmesu, nbmode, nbabs, zr(lred), zc(lmesu), &
                        zr(lcoef), zr(labs), lfonct, zc(jdep), zc(jvit), &
                        zc(jacc))
!
            call rscrsd('G', nomres, 'MODE_GENE', nbabs)
            call mdallr(nommes, nomres, nombas, nbmode, nbabs, &
                        [0.d0], zc(jdep), zcmplx)
        else
            call utmess('F', 'ALGORITH6_33')
        end if
!
    end if
!
!     -- REMPLISSAGE DE L'OBJET .ORDR :
!
    call jeveuo(nomres//'           .ORDR', 'E', jordr)
    do i = 1, nbabs
        if (typres(1:9) .eq. 'MODE_GENE') then
            zi(jordr-1+i) = i
        else if (typres(1:9) .eq. 'HARM_GENE') then
            zi(jordr-1+i) = i
        else
            zi(jordr-1+i) = i-1
        end if
    end do
!     -- REMPLISSAGE DE L'OBJET .PTEM :
    if (typres(1:9) .eq. 'TRAN_GENE') then
        call jeexin(nommes//'           .PTEM', iexi)
        if (iexi .gt. 0) then
            call jeveuo(nommes//'           .PTEM', 'E', jpames)
            do i = 1, nbabs
                zr(jpass-1+i) = zr(jpames-1+i)
            end do
        end if
    end if
!
!
!     -- REMPLISSAGE DE L'OBJET .NUMO :
    if (typres(1:9) .eq. 'MODE_GENE') then
        do i = 1, nbabs
            call rsadpa(nomres, 'E', 1, 'NUME_MODE', zi(jordr-1+i), &
                        0, sjv=jpara, styp=k8b)
            zi(jpara) = i
        end do
    end if
!     -- REMPLISSAGE DE "FREQ/DISC" :
    if (typres(1:9) .eq. 'TRAN_GENE' .or. typres(1:9) .eq. 'HARM_GENE') then
        call jeveuo(nomres//'           .DISC', 'E', jabs)
    else
        call jeexin(nomres//'           .FREQ', iexi)
        if (iexi .gt. 0) then
            call jeveuo(nomres//'           .FREQ', 'E', jabs)
        end if
    end if
    do i = 1, nbabs
        if (typres(1:9) .eq. 'TRAN_GENE' .or. typres(1:9) .eq. 'HARM_GENE') then
            zr(jabs-1+i) = zr(labs-1+i)
        else
            if (iexi .gt. 0) then
                zr(jabs-1+i) = zr(labs-1+i)
            else
                call rsadpa(nomres, 'E', 1, 'FREQ', zi(jordr-1+i), &
                            0, sjv=jpara, styp=k8b)
                zr(jpara) = zr(labs-1+i)
            end if
        end if
    end do
!
!
! --- STOCKAGE
    if (typres(1:9) .eq. 'MODE_GENE') then
        call jeveuo(nomres//'           .ORDR', 'L', jord)
        call jelira(nomres//'           .ORDR', 'LONUTI', nbord)
        call dismoi('TYPE_BASE', nombas, 'RESU_DYNA', repk=typba, arret='C', &
                    ier=iret)
        call dismoi('REF_RIGI_PREM', nombas, 'RESU_DYNA', repk=raide, arret='C')
        if (typba(1:1) .ne. ' ') then
!            if (raide(1:8) .eq. '        ') then
            call jeveuo(jexnum(nombas//'           .TACH', 1), 'L', jbasm)
            sd2 = zk24(jbasm) (1:8)
            call rsadpa(sd2, 'L', 1, 'MODELE', 1, &
                        0, sjv=jpara, styp=k8b)
            modele = zk8(jpara)
            call rsadpa(sd2, 'L', 1, 'CHAMPMAT', 1, &
                        0, sjv=jpara, styp=k8b)
            chmat = zk8(jpara)
            call rsadpa(sd2, 'L', 1, 'CARAELEM', 1, &
                        0, sjv=jpara, styp=k8b)
            carael = zk8(jpara)
            goto 44
!            end if
        end if
!
!       -- POUR LES BASES TYPE MODE_MECA SANS REFERENCE
        if (raide(1:8) .eq. '        ') then
            modele = '        '
            chmat = '        '
            carael = '        '
            goto 44
        end if
!
        call dismoi('NOM_MODELE', raide(1:8), 'MATR_ASSE', repk=modele)
!        call dismoi('CHAM_MATER', raide(1:8), 'MATR_ASSE', repk=chmat)
!        call dismoi('CARA_ELEM', raide(1:8), 'MATR_ASSE', repk=carael)
        chmat = '        '
        carael = '        '
44      continue
!
!
        do i = 1, nbord
            call rsadpa(nomres, 'E', 1, 'MODELE', zi(jordr-1+i), &
                        0, sjv=jpara, styp=k8b)
            zk8(jpara) = modele
            call rsadpa(nomres, 'E', 1, 'CHAMPMAT', zi(jordr-1+i), &
                        0, sjv=jpara, styp=k8b)
            zk8(jpara) = chmat
            call rsadpa(nomres, 'E', 1, 'CARAELEM', zi(jordr-1+i), &
                        0, sjv=jpara, styp=k8b)
            zk8(jpara) = carael
        end do
    end if
!
    call jedetr(vabs)
    call jedetr(vmes)
    call jedetr('&&LISTE_CH')
    call detrsd('CHAM_NO_S', chs)
!
    call jedema()
!
end subroutine
