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

subroutine defint(mailla, nomres)
!    P. RICHARD     DATE 20/02/90
!-----------------------------------------------------------------------
!  BUT:  RECUPERATION DE LA DEFINITION UTLISATEUR DES INTERFACES
    implicit none
!       CREATION DE LA STRUCTURE DE DONNEES CORRESPONDANTE
!       CREATION DES LISTES GLOBALE DES NOEUDS D'INTERFACE LIBRES ET
!        BLOQUEES
!-----------------------------------------------------------------------
!
! MAILLA   /I/:NOM UT DU MAILLAGE EN AMONT
! NOMRES   /I/: NOM UT DU RESULTAT
!
!
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/defdda.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utlisi.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, ideb, iec, ifin, ioc
    integer(kind=8) :: iret, k, ldact, ldfreq, llnin, lltyp
    integer(kind=8) :: ltgui, ltlgr, ltlno, ltmas, maxgr, maxno, nball
    integer(kind=8) :: nbbid, nbcmp, nbec, nbecmx, nbgr, nbint
    integer(kind=8) :: nbno, nbtemp, nbvag, nbval, nbvan
    integer(kind=8) :: numgd
    real(kind=8) :: freq
!-----------------------------------------------------------------------
    parameter(nbecmx=10)
    character(len=8) :: nomres, mailla
    character(len=8) :: nomcou, type
    character(len=9) :: nom, no, grno
    character(len=14) :: int
    character(len=24) :: nomint, notint, typint, noeint, ddlact
    character(len=24) :: valk(2)
    character(len=24) :: temgui, temlno, temlgr, temmas
    character(len=80) :: kar80
    integer(kind=8) :: icodma(nbecmx), icodac(nbecmx)
    integer(kind=8) :: vali, ign, ii, nbis, ntrou, nuno
    integer(kind=8) :: iagm1, iagm2, ialii1, ialii2, ign1, ign2, ili1, ili2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!----------------RECUPERATION DES DONNEES GRANDEUR----------------------
!
    call jemarq()
    call dismoi('NB_CMP_MAX', nomres, 'INTERF_DYNA', repi=nbcmp)
    call dismoi('NB_EC', nomres, 'INTERF_DYNA', repi=nbec)
    call dismoi('NUM_GD', nomres, 'INTERF_DYNA', repi=numgd)
!
!-------------INITIALISATION DES NOMS TRES UTILISES---------------------
!
    int = 'INTERFACE'
    nom = 'NOM'
    no = 'NOEUD'
    grno = 'GROUP_NO'
!
!-------------- DETERMINATION DU NOMBRE D'INTERFACES -------------------
!
    call getfac(int, ioc)
!
    if (ioc .eq. 0) then
        call utmess('A', 'ALGORITH12_77')
    end if
!
    nbint = 0
!
    do i = 1, ioc
        call getvtx(int, nom, iocc=i, scal=kar80, nbret=nbval)
        if (nbval .gt. 0) nbint = nbint+1
    end do
!
!------------ALLOCATION REPERTOIRES DES NOMS D'INTERFACE ---------------
!
    nomint = nomres//'.IDC_NOMS'
    call jecreo(nomint, 'G N K8')
    call jeecra(nomint, 'NOMMAX', nbint)
!
    notint = '&&DEFINT.NOM.INTF'
    call jecreo(notint, 'V N K8')
    call jeecra(notint, 'NOMMAX', nbint)
!
!-------ALLOCATION VECTEURS DE TRAVAIL----------------------------------
!
!  NUMERO MOTCLE FACTEUR DEBUT ET FIN DESCRIPTION  INTERFACE
!
    temgui = '&&DEFINT.GUIDE.INTF'
    nball = 2*nbint
    call wkvect(temgui, 'V V I', nball, ltgui)
!
!------DETERMINATION DES OCCURENCES DEBUT ET FIN INTERFACE--------------
!          STOCKAGE NOMBRE DE GROUPES PAR INTERFACE
!   VERIFICATION EXCLUSIVITE DE NOEUD ET GROUPNO PAR OCCURENCE MOTFAC
!             ET STOCKAGE DES NOMS DES INTERFACES
!
    call getvtx(int, no, iocc=1, nbval=0, nbret=nbvan)
    call getvtx(int, grno, iocc=1, nbval=0, nbret=nbvag)
!
    call getvtx(int, nom, iocc=1, scal=kar80, nbret=nbval)
    nomcou = kar80
    call jecroc(jexnom(nomint, nomcou))
    call jecroc(jexnom(notint, nomcou))
!
    zi(ltgui) = 1
    nbgr = -nbvag
    nbno = -nbvan
    maxgr = 0
    maxno = 0
    nbtemp = 1
!
    if (ioc .gt. 1) then
!
        do i = 2, ioc
            call getvtx(int, nom, iocc=i, scal=kar80, nbret=nbval)
            nomcou = kar80
            call jeexin(jexnom(notint, nomcou), iret)
            if (iret .ne. 0) then
                vali = i
                valk(1) = nomcou
                call utmess('F', 'ALGORITH12_78', sk=valk(1), si=vali)
            end if
            call getvtx(int, no, iocc=i, nbval=0, nbret=nbvan)
            call getvtx(int, grno, iocc=i, nbval=0, nbret=nbvag)
!
            if (nbval .gt. 0) then
                maxgr = max(maxgr, nbgr)
                maxno = max(maxno, nbno)
                nbgr = 0
                nbno = 0
                zi(ltgui+nbint+nbtemp-1) = i-1
                call jecroc(jexnom(nomint, nomcou))
                call jecroc(jexnom(notint, nomcou))
                nbtemp = nbtemp+1
                zi(ltgui+nbtemp-1) = i
            end if
!
!
!
            nbgr = nbgr-nbvag
            nbno = nbno-nbvan
        end do
!
    end if
!
! TRAITEMENT DERNIERE INTERFACE
!
    maxgr = max(maxgr, nbgr)
    maxno = max(maxno, nbno)
    zi(ltgui+nbint+nbtemp-1) = ioc
!
!--------ALLOCATION DE LA LISTE DES GROUPES DE L'INTERFACE COURANTE-----
!
    temlgr = '&&DEFINT.LIST.GRP'
    if (maxgr .gt. 0) then
        call wkvect(temlgr, 'V V K24', maxgr, ltlgr)
    else
        ltlgr = 1
    end if
!
!--------ALLOCATION DE LA LISTE DES NOEUDS DE L'INTERFACE COURANTE------
!
    temlno = '&&DEFINT.LIST.NO'
    if (maxno .gt. 0) then
        call wkvect(temlno, 'V V K8', maxno, ltlno)
    else
        ltlno = 1
    end if
!
!---------ALLOCATION DU  VECTEUR DES TYPES D'INTERFACES-----------------
!
    typint = nomres//'.IDC_TYPE'
    call wkvect(typint, 'G V K8', nbint, lltyp)
!
!---------ALLOCATION DU  VECTEUR DES DDL INTERFACE ACTIFS---------------
!
    ddlact = nomres//'.IDC_DDAC'
    call jecrec(ddlact, 'G V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbint)
!
!---------ALLOCATION DU  VECTEUR DES MASQUES DDL INTERFACE--------------
!
    temmas = '&&DEFINT.MASQUE'
    call jecrec(temmas, 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbint)
!
!---------ALLOCATION DES VECTEURS LISTE DES NOEUDS INTERFACE------------
!
    noeint = nomres//'.IDC_LINO'
    call jecrec(noeint, 'G V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbint)
!
    do i = 1, nbint
        ideb = zi(ltgui+i-1)
        ifin = zi(ltgui+nbint+i-1)
!
!-----STOCKAGE DU TYPE DE L'INTERFACE COURANTE--------------------------
!
        call getvtx(int, nom, iocc=ideb, scal=kar80, nbret=nbbid)
        nomcou = kar80
        call getvtx(int, 'TYPE', iocc=ideb, scal=kar80, nbret=nbval)
        type = kar80
        zk8(lltyp+i-1) = type
!
!---------SOUS-BOUCLE COMPTAGE NOMBRE DE NOEUDS DE CHAQUE INTERFACE-----
!
        nbno = 0
        ASSERT(ideb .eq. ifin)
!
        call getvtx(int, grno, iocc=ideb, nbval=maxgr, vect=zk24(ltlgr), &
                    nbret=nbvag)
        nbgr = nbvag
        if (nbgr .gt. 0) then
! cas des groupes de no
            call jenonu(jexnom(mailla//'.GROUPENO', zk24(ltlgr+1-1)), ign1)
            call jelira(jexnum(mailla//'.GROUPENO', ign1), 'LONUTI', ili1)
            call jeveuo(jexnum(mailla//'.GROUPENO', ign1), 'L', iagm1)
            nbis = 2*ili1
            call wkvect('&&DEFINT.LII1', 'V V I', nbis, ialii1)
            call wkvect('&&DEFINT.LII2', 'V V I', nbis, ialii2)
!
            nbno = ili1
            do ii = 1, nbno
                zi(ialii1-1+ii) = zi(iagm1-1+ii)
            end do
!
            do ign = 2, nbgr
                call jenonu(jexnom(mailla//'.GROUPENO', zk24(ltlgr+ign-1)), ign2)
                call jelira(jexnum(mailla//'.GROUPENO', ign2), 'LONUTI', ili2)
                call jeveuo(jexnum(mailla//'.GROUPENO', ign2), 'L', iagm2)
                call utlisi('UNION', zi(ialii1), nbno, zi(iagm2), ili2, &
                            zi(ialii2), nbis, ntrou)
!
                if (ntrou .lt. 0) then
                    nbis = -2*ntrou
                    call jedetr('&&DEFINT.LII2')
                    call wkvect('&&DEFINT.LII2', 'V V I', nbis, ialii2)
                    call utlisi('UNION', zi(ialii1), nbno, zi(iagm2), ili2, &
                                zi(ialii2), nbis, ntrou)
                    call jedetr('&&DEFINT.LII1')
                    call wkvect('&&DEFINT.LII1', 'V V I', nbis, ialii1)
                end if
                nbno = ntrou
                do ii = 1, nbno
                    zi(ialii1-1+ii) = zi(ialii2-1+ii)
                end do
            end do
            call jedetr('&&DEFINT.LII2')
        else
! cas des noeuds
            call getvtx(int, no, iocc=ideb, nbval=maxno, vect=zk8(ltlno), &
                        nbret=nbvan)
            nbno = nbvan
            call wkvect('&&DEFINT.LII1', 'V V I', nbno, ialii1)
            do ii = 1, nbno
                nuno = char8_to_int(zk8(ltlno+ii-1))
                if (nuno .eq. 0) then
                    valk(1) = mailla
                    valk(2) = zk8(ltlno+ii-1)
                    call utmess('F', 'ALGORITH14_11', nk=2, valk=valk)
                end if
                zi(ialii1+ii-1) = nuno
            end do
        end if
!
!-------ALLOCATION DES VECTEURS NUMERO NOEUDS INTERFACES----------------
!                     MASQUE,DDL ACTIFS
!
        call jenonu(jexnom(nomint, nomcou), ibid)
        call jecroc(jexnum(noeint, ibid))
        call jeecra(jexnum(noeint, ibid), 'LONMAX', nbno)
        call jeveuo(jexnum(noeint, ibid), 'E', llnin)
!
        call jenonu(jexnom(nomint, nomcou), ibid)
        call jecroc(jexnum(ddlact, ibid))
        call jeecra(jexnum(ddlact, ibid), 'LONMAX', nbno*nbec)
        call jeveuo(jexnum(ddlact, ibid), 'E', ldact)
!
        call jenonu(jexnom(nomint, nomcou), ibid)
        call jeecra(jexnum(temmas, ibid), 'LONMAX', nbno*nbec)
        call jeveuo(jexnum(temmas, ibid), 'E', ltmas)
!
!------SOUS-BOUCLE COPIE DES NUMERO DE NOEUDS INTERFACE COURANTE-------
!        ET STOCKAGES DES MASQUES ET DDL ACTIFS POUR INTERFACES
!
        do ii = 1, nbno
            zi(llnin-1+ii) = zi(ialii1-1+ii)
        end do
        call jedetr('&&DEFINT.LII1')
!
!----RECUPERATION DES MASQUES ET DDL ACTIFS INTERFACE---------------
!
        call defdda(nbec, nbcmp, numgd, ideb, 'MASQUE', &
                    0, icodma)
        call defdda(nbec, nbcmp, numgd, ideb, 'DDL_ACTIF', &
                    1, icodac)
        do k = 1, nbno
            do iec = 1, nbec
                zi(ldact+(k-1)*nbec+iec-1) = icodac(iec)
                zi(ltmas+(k-1)*nbec+iec-1) = icodma(iec)
            end do
        end do
    end do
!
!------SAUVEGARDE DES OBJECT EN GLOBALE SAUF CEUX A MODIFIER------------
!
!
!----------------RECUPERATION FREQUENCE POUR CB_HARMO-------------------
!
    call getvr8('   ', 'FREQ', iocc=1, scal=freq, nbret=ioc)
    call wkvect(nomres//'.IDC_DY_FREQ', 'G V R', 1, ldfreq)
    zr(ldfreq) = freq
!
!-----------------DESTRUCTION DES OBJETS EN VOLATILE--------------------
!
    call jedetr(temgui)
    if (maxno .gt. 0) call jedetr(temlno)
    if (maxgr .gt. 0) call jedetr(temlgr)
!
    call jedema()
end subroutine
