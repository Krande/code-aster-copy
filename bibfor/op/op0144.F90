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

subroutine op0144()
    implicit none
!
!     OPERATEUR "CALC_FLUI_STRU"
!
!-----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cmpcha.h"
#include "asterfort/crprno.h"
#include "asterfort/dismoi.h"
#include "asterfort/flust1.h"
#include "asterfort/flust2.h"
#include "asterfort/flust3.h"
#include "asterfort/flust4.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/iunifi.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/pteequ.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!     UN COMMON AJOUTE POUR RESORBER UNE GLUTE ANTIQUE (VOIR HISTOR):
    character(len=8) :: typflu
    common/kop144/typflu
!
    integer(kind=8) :: ibid, nbconn, iconn
    integer(kind=8) :: ncmp
    real(kind=8) :: r8b
    complex(kind=8) :: c16b
    aster_logical :: tmode, calcul(2)
    character(len=8) :: nombm, mailla, k8b, gran, nomcmp(6)
    character(len=16) :: concep, cmd, nompar
    character(len=19) :: nomu, cham19, numeq
    character(len=24) :: desc, numo, vite, freq, masg, fact
    character(len=24) :: numoi, fsic, nomcha, matria, chrefe
    character(len=24) :: chvale
    character(len=32) :: nomvar
    integer(kind=8), pointer :: cata_to_field(:) => null()
    integer(kind=8), pointer :: field_to_cata(:) => null()
    character(len=8), pointer :: cmp_name(:) => null()
!
!-----------------------------------------------------------------------
    data nomcmp/'DX      ', 'DY      ', 'DZ      ', 'DRX     ', 'DRY     ', 'DRZ     '/
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacmp, iamor, iav, icmp, idec, idesc
    integer(kind=8) :: iec, ifact, ifm, ifr, ifreq, ifsic
    integer(kind=8) :: ii, imasg, inec, ino, inumo, io, ipar
    integer(kind=8) :: irefe, itypfl, iv, ivite
    integer(kind=8) :: j, jcrefe, jdesc, jj, long
    integer(kind=8) :: nbam, nbcomp, nbno, nbnoeu, nbocc, nbpar, nbpv
    integer(kind=8) :: nec, nivdef, nivpar, numgd
    real(kind=8) :: amor, umin, vmax, vmin, vmoy, vpas
    integer(kind=8), pointer :: prno(:) => null()
    ibid = 0

    call jemarq()
    c16b = (0.d0, 0.d0)
    r8b = 0.d0
    call infmaj()
!
    call getres(nomu, concep, cmd)
    call getvis('BASE_MODALE', 'NUME_ORDRE', iocc=1, nbval=0, nbret=nbno)
    call getvr8('BASE_MODALE', 'AMOR_REDUIT', iocc=1, nbval=0, nbret=nbam)
    call getvr8('BASE_MODALE', 'AMOR_REDUIT_CONN', iocc=1, nbval=0, nbret=nbconn)
    ifr = iunifi('RESULTAT')
    ifm = iunifi('MESSAGE')
!
!
! --- 0.VERIFICATIONS AVANT EXECUTION ---
!
    if (nbno .ne. 0 .and. nbam .ne. 0 .and. nbam .ne. nbno) then
        call utmess('F', 'ALGELINE2_82')
    end if
!
    if (nbno .ne. 0 .and. nbconn .ne. 0 .and. nbconn .ne. nbno) then
        call utmess('F', 'ALGELINE2_83')
    end if
!
! --- 1.RECUPERATION DES CARACTERISTIQUES MODALES AVANT COUPLAGE ---
!
    call getvid('BASE_MODALE', 'MODE_MECA', iocc=1, scal=nombm)
!
    tmode = .false.
    calcul(1) = .false.
    calcul(2) = .false.
    if (nbno .eq. 0) then
        tmode = .true.
        numoi = nombm//'           .ORDR'
        call jelira(numoi, 'LONUTI', nbno)
        ASSERT(nbam .eq. 0 .or. abs(nbam) .eq. nbno)
    else
        nbno = abs(nbno)
    end if
!
! --- 1.1.CREATION ET REMPLISSAGE DE L'OBJET .NUMO
!
    numo = nomu//'.NUMO'
    call wkvect(numo, 'G V I', nbno, inumo)
    if (tmode) then
        do i = 1, nbno
            zi(inumo+i-1) = i
        end do
    else
        call getvis('BASE_MODALE', 'NUME_ORDRE', iocc=1, nbval=nbno, vect=zi(inumo))
    end if
!
! --- 1.2.CREATION D'UN VECTEUR TEMPORAIRE POUR LES AMORTISSEMENTS
!
    call wkvect('&&OP0144.TEMP.AMOR', 'V V R', nbno, iamor)
    call wkvect('&&OP0144.CONNORS.AMOR', 'V V R', nbno, iconn)
    if (nbam .ne. 0) then
        call getvr8('BASE_MODALE', 'AMOR_REDUIT', iocc=1, nbval=nbno, vect=zr(iamor))
        calcul(1) = .true.
    else if (nbconn .eq. 0) then
        calcul(1) = .true.
        call getvr8('BASE_MODALE', 'AMOR_UNIF', iocc=1, scal=amor)
        do i = 1, nbno
            zr(iamor+i-1) = amor
        end do
    end if
!
    if (nbconn .ne. 0) then
        call getvr8('BASE_MODALE', 'AMOR_REDUIT_CONN', iocc=1, nbval=nbno, vect=zr(iconn))
        calcul(2) = .true.
    end if
!
!
! --- 2.RECUPERATION DE LA PLAGE DE VITESSES D'ECOULEMENT ---
!
! --- 2.1.CREATION ET REMPLISSAGE DE L'OBJET .VITE
!
    call getvr8('VITE_FLUI', 'VITE_MIN', iocc=1, scal=vmin)
    call getvr8('VITE_FLUI', 'VITE_MAX', iocc=1, scal=vmax)
    call getvis('VITE_FLUI', 'NB_POIN ', iocc=1, scal=nbpv)
    if (vmin .gt. vmax) then
        umin = vmin
        vmin = vmax
        vmax = umin
        call utmess('A', 'ALGELINE2_85')
    end if
!
    vite = nomu//'.VITE'
    call wkvect(vite, 'G V R', nbpv, ivite)
    if (nbpv .eq. 1) then
        vmoy = (vmin+vmax)/2.d0
        write (ifm, *) 'UNE SEULE VITESSE D''ECOULEMENT ETUDIEE :'// &
            ' VMOY = (VMIN+VMAX)/2'
        zr(ivite) = vmoy
    else
        vpas = (vmax-vmin)/(nbpv-1)
        do iv = 1, nbpv
            zr(ivite+iv-1) = vmin+(iv-1)*vpas
        end do
    end if
!
! --- 2.2.CREATION DE L'OBJET .FREQ
!
    freq = nomu//'.FREQ'
    call wkvect(freq, 'G V R', 2*nbno*nbpv, ifreq)
!
!
! --- 3.RECUPERATION DU CONCEPT TYPE_FLUI_STRU   ---
!
! --- 3.1.CREATION ET REMPLISSAGE DE L'OBJET .REMF
!
    call wkvect(nomu//'.REMF', 'G V K8', 2, irefe)
    call getvid(' ', 'TYPE_FLUI_STRU', scal=typflu)
    zk8(irefe) = typflu
    zk8(irefe+1) = nombm
!
! --- 3.2.DETERMINATION DU TYPE DE LA CONFIGURATION ETUDIEE
!
    fsic = typflu//'           .FSIC'
    call jeveuo(fsic, 'L', ifsic)
    itypfl = zi(ifsic)
!
!
! --- 4.RECUPERATION DES NIVEAUX D'IMPRESSION ---
!
    nivpar = 0
    nivdef = 0
    call getfac('IMPRESSION', nbocc)
    if (nbocc .ne. 0) then
        call getvtx('IMPRESSION', 'PARA_COUPLAGE', iocc=1, scal=k8b)
        if (k8b(1:3) .eq. 'OUI') nivpar = 1
        call getvtx('IMPRESSION', 'DEFORMEE', iocc=1, scal=k8b)
        if (k8b(1:3) .eq. 'OUI') nivdef = 1
    end if
!
!
! --- 5.CREATION ET REMPLISSAGE DE L'OBJET .DESC ---
!
    desc = nomu//'.DESC'
    call wkvect(desc, 'G V K16', 1, idesc)
    zk16(idesc) = 'DEPL'
!
!
! --- 6.CREATION ET REMPLISSAGE DE LA TABLE DES NOMS DES CHAMPS ---
! ---   DE DEPLACEMENTS MODAUX                                  ---
! ---   SIMULTANEMENT ON CREE LES OBJETS .REFE , .DESC ET .VALE ---
! ---   ASSOCIES A CHACUN DES CHAMPS                            ---
! ---   POUR LE PREMIER CHAMP ON CREE LE NUME_EQUA QUI VAUT     ---
! ---   ENSUITE POUR TOUS LES AUTRES CHAMPS                     ---
!
! --- 6.1.RECUPERATION D'INFORMATIONS NECESSAIRES
! ---     A LA CREATION DES OBJETS ASSOCIES AUX CHAMPS
! ---     A LA CREATION DU NUME_EQUA COMMUN
!
    call dismoi('REF_RIGI_PREM', nombm, 'RESU_DYNA', repk=matria)
!
    call dismoi('NOM_MAILLA', matria, 'MATR_ASSE', repk=mailla)
    call jelira(mailla//'.COORDO    .VALE', 'LONMAX', nbnoeu)
    nbnoeu = nbnoeu/3
    long = 6*nbnoeu
!
    gran = 'DEPL_R  '
    call jenonu(jexnom('&CATA.GD.NOMGD', gran), numgd)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', numgd), 'L', iacmp)
    call jeveuo(jexatr('&CATA.GD.NOMCMP', 'LONCUM'), 'L', iav)
    nbcomp = zi(iav+numgd)-zi(iav+numgd-1)
!
    call dismoi('NB_EC', gran, 'GRANDEUR', repi=nec)
    call wkvect('&&OP0144.DESC_NOEUD', 'V V I', nec*nbnoeu, jdesc)
    do ino = 1, nbnoeu
        do icmp = 1, 6
            j = indik8(zk8(iacmp), nomcmp(icmp), 1, nbcomp)
            if (j .ne. 0) then
                iec = (j-1)/30+1
                jj = j-30*(iec-1)
                zi(jdesc+(ino-1)*nec+iec-1) = ior(zi(jdesc+(ino-1)*nec+iec-1), 2**jj)
            end if
        end do
    end do
!
! --- 6.2.CREATION DE LA STRUCTURE TABLE
!
    if (itypfl .eq. 3) then
        nbpar = nbpv
    else
        nbpar = 1
    end if
    call tbcrsd(nomu, 'G')
    call tbajpa(nomu, 1, 'NOM_CHAM', 'K24')
!
! --- 6.31.CREATION DE L'OBJET .MASG
!
    masg = nomu//'.MASG'
    call wkvect(masg, 'G V R', nbno*nbpar, imasg)
!
! --- 6.32.CREATION DE L'OBJET .FACT
!
    fact = nomu//'.FACT'
    call wkvect(fact, 'G V R', 3*nbno*nbpar, ifact)
!
! --- 6.4.REMPLISSAGE DE LA TABLE (NOMS DES CHAMPS) ET CREATION
! ---     SIMULTANEE DES OBJETS ASSOCIES AUX CHAMPS
!
    nomcha(1:13) = nomu(1:8)//'.C01.'
!
    do io = 1, nbno
!
        write (nomvar, '(I3.3)') zi(inumo+io-1)
        nomcha(14:16) = nomvar(1:3)
!
        do ipar = 1, nbpar
!
            write (nompar, '(I3.3)') ipar
            nomcha(17:24) = nompar(1:3)//'     '
!
            call tbajli(nomu, 1, 'NOM_CHAM', [ibid], [r8b], &
                        [c16b], nomcha, 0)
!
! --------CREATION DU CHAMP
            cham19 = nomcha(1:19)
! .VALE
            chvale = cham19//'.VALE'
            call jecreo(chvale, 'G V R')
            call jeecra(chvale, 'LONMAX', long)
! .REFE
            chrefe = cham19//'.REFE'
            call wkvect(chrefe, 'G V K24', 4, jcrefe)
            call jeecra(chrefe, 'DOCU', cval='CHNO')
!
! --------AU PREMIER PASSAGE CREATION DU NUME_EQUA
            if (io .eq. 1 .and. ipar .eq. 1) then
                zk24(jcrefe+1) = cham19
                call crprno(cham19, 'G', mailla, gran, long)
                call jeveuo(cham19//'.PRNO', 'E', vi=prno)
                idec = 1
                ii = 0
                do ino = 1, nbnoeu
                    prno((nec+2)*(ino-1)+1) = idec
                    prno((nec+2)*(ino-1)+2) = 6
                    do inec = 1, nec
                        ii = ii+1
                        prno((nec+2)*(ino-1)+2+inec) = zi(jdesc+ii-1)
                    end do
                    idec = idec+6
                end do
                numeq = cham19
!
! ------------- Create object local components (field) => global components (catalog)
!
                call cmpcha(cham19, cmp_name, cata_to_field, field_to_cata, nb_cmpz=ncmp)
!
! ------------- Compute .DEEQ object
!
                call pteequ(numeq, 'G', long, numgd, ncmp, &
                            field_to_cata)
                AS_DEALLOCATE(vi=cata_to_field)
                AS_DEALLOCATE(vi=field_to_cata)
                AS_DEALLOCATE(vk8=cmp_name)
!
            else
                zk24(jcrefe+1) = numeq
            end if
!
        end do
!
    end do
!
!
! --- 7.LANCEMENT DU CALCUL EN FONCTION DU TYPE DE LA CONFIGURATION ---
! ---   ETUDIEE                                                     ---
!
    if (itypfl .eq. 1) then
!
        call flust1(nomu, typflu, nombm, zi(inumo), zr(iamor), &
                    zr(iconn), zr(ifreq), zr(imasg), zr(ifact), zr(ivite), &
                    nbno, calcul, nbpv, nivpar, nivdef)
!
    else if (itypfl .eq. 2) then
!
        call flust2(nomu, typflu, nombm, mailla, zi(inumo), &
                    zr(iamor), zr(ifreq), zr(imasg), zr(ifact), zr(ivite), &
                    nbno, nbpv, nivpar, nivdef)
!
    else if (itypfl .eq. 3) then
!
        call flust3(nomu, typflu, nombm, zi(inumo), zr(iamor), &
                    zr(ifreq), zr(imasg), zr(ifact), zr(ivite), nbno, &
                    nbpv, nivpar, nivdef)
!
    else
!
        call flust4(nomu, typflu, nombm, mailla, zi(inumo), &
                    zr(iamor), zr(ifreq), zr(imasg), zr(ifact), zr(ivite), &
                    nbno, nbpv, nivpar, nivdef)
!
    end if
!
!
    call jedema()
    call jedetc('G', '&&MEFCEN', 1)

end subroutine
