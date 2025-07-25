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

subroutine sinoz2(modele, nume_equa, sigel, signo)
!   BUT :  CALCUL DES CONTRAINTES AUX NOEUDS PAR LA METHODE ZZ2
    implicit none
!
!   IN  MODELE   :   NOM DU MODELE
!   IN  NUME_EQUA   :   NUME_EQUA
!   IN  SIGEL    :   NOM DU CHAMP DE CONTRAINTES AUX POINTS DE GAUSS
!
!  OUT  SIGNO    :   NOM DU CHAMP DE CONTRAINTES AUX NOEUDS
!
! ----------------------- DECLARATIONS --------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assvec.h"
#include "asterfort/celfpg.h"
#include "asterfort/celver.h"
#include "asterfort/copisd.h"
#include "asterfort/crcnct.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mecanb.h"
#include "asterfort/mtcrou.h"
#include "asterfort/predia.h"
#include "asterfort/utelvf.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zzappa.h"
#include "asterfort/zzcala.h"
#include "asterfort/zzcalb.h"
#include "asterfort/zzpoly.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=8) :: modele, ma, typema, licmp(4), vecass, elrefe
    character(len=8) :: famil
    character(len=14) :: nu14
    character(len=19) :: nume_equa
    character(len=16) :: phen
    character(len=19) :: noeub, mo, vecel
    character(len=24) :: signo, sigel, lisvec, typmai, connex, coninv
    character(len=24) :: nomjv, elrfam
    real(kind=8) :: rcmp(4), eps, x(9), y(9), a(9, 9), b(9, 4), diag(9)
    real(kind=8) :: wk1(9, 9), wk2(9)
    integer(kind=8) :: nno, npg, ivf
    aster_logical :: app
!
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, iamav, ianew, ianob
    integer(kind=8) :: ianov, iatyma, ibid, ic, icmp
    integer(kind=8) :: ima, ino, inob, inoma, ipa
    integer(kind=8) :: jcon, jconin, jelfa, jpa
    integer(kind=8) :: jprno, jrefn, k
    integer(kind=8) :: nb, nbcmp, nbec, nbma, nbmav, nbn, nbno
    integer(kind=8) :: nbnob, nbnobp, nbnoma, nqua, ntri, num, numav
    integer(kind=8) :: numc, numel, numeq, numgr, numloc
    real(kind=8) :: xino, xinob, xinomi, xmax, xmin, xnorm, yino
    real(kind=8) :: yinob, yinomi, ymax, ymin
    integer(kind=8), pointer :: indic(:) => null()
    integer(kind=8), pointer :: longconinv(:) => null()
    integer(kind=8), pointer :: nbpatchmil(:) => null()
    aster_logical, pointer :: noeubord(:) => null()
    integer(kind=8), pointer :: numnb(:) => null()
    integer(kind=8), pointer :: celd(:) => null()
    character(len=24), pointer :: refe(:) => null()
    integer(kind=8), pointer :: repe(:) => null()
    real(kind=8), pointer :: coor(:) => null()
    real(kind=8), pointer :: sig(:) => null()
    real(kind=8), pointer :: val(:) => null()
    real(kind=8), pointer :: celv(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    nomjv = '&&SINOZ2.FFORMES        '
    elrfam = '&&SINOZ2.ELRE_FAMI'
!
!     -- ON VERIFIE QUE LE CHAM_ELEM N'EST PAS TROP DYNAMIQUE :
    call celver(sigel, 'NBVARI_CST', 'STOP', ibid)
    call celver(sigel, 'NBSPT_1', 'STOP', ibid)
!
    call dismoi('PHENOMENE', modele, 'MODELE', repk=phen)
    if (phen .ne. 'MECANIQUE') then
        call utmess('F', 'CALCULEL4_83')
    end if
!
    call dismoi('NOM_MAILLA', sigel(1:19), 'CHAM_ELEM', repk=ma)
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbno)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
    typmai = ma//'.TYPMAIL'
    connex = ma//'.CONNEX'
    call jeveuo(typmai, 'L', iatyma)
!
!
!   CONSTRUCTION DE LA CONNECTIVITE INVERSE (OBJET TEMPORAIRE)
!     --  OBJET CONINV    = FAMILLE CONTIGUE DE VECTEURS N*IS
    coninv = '&&SINOZ2.CONINV'
    call jecrec(coninv, 'V V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbno)
!
    AS_ALLOCATE(vi=longconinv, size=nbno)
    do ima = 1, nbma
        iad = iatyma-1+ima
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(iad)), typema)
        if (typema(1:4) .eq. 'TRIA' .or. typema(1:4) .eq. 'QUAD') then
            call jeveuo(jexnum(connex, ima), 'L', jcon)
            call jelira(jexnum(connex, ima), 'LONMAX', nbn)
            do ino = 1, nbn
                num = zi(jcon-1+ino)
                longconinv(num) = longconinv(num)+1
            end do
        end if
    end do
!
    do ino = 1, nbno
        call jeecra(jexnum(coninv, ino), 'LONMAX', longconinv(ino))
    end do
!
    AS_ALLOCATE(vi=indic, size=nbno)
!
    do ima = 1, nbma
        iad = iatyma-1+ima
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(iad)), typema)
        if (typema(1:4) .eq. 'TRIA' .or. typema(1:4) .eq. 'QUAD') then
            call jeveuo(jexnum(connex, ima), 'L', jcon)
            call jelira(jexnum(connex, ima), 'LONMAX', nbn)
            do ino = 1, nbn
                num = zi(jcon-1+ino)
                call jeveuo(jexnum(coninv, num), 'E', jconin)
                nb = indic(num)
                zi(jconin+nb) = ima
                indic(num) = indic(num)+1
            end do
        end if
    end do
    AS_DEALLOCATE(vi=indic)
!
!   CONSTRUCTION D'UN VECTEUR DE BOOLEENS SUR LES NOEUDS INDIQUANT
!   L'APPARTENANCE OU NON AU BORD
!
    vecel = '&&NOEUB            '
    call mecanb(modele, vecel)
    vecass = '&&VECASS'
    lisvec = vecel//'.RELR'
!
!     -- POUR POUVOIR APPELER ASSVEC, IL FAUT CREER UN "FAUX"
!        NUME_DDL AVEC UN NUME_EQUA :
    nu14 = '&&SINOZ2.NUDDL'
    call copisd('NUME_EQUA', 'V', nume_equa, nu14//'.NUME')
    call jeveuo(nu14//'.NUME.REFN', "E", jrefn)
    zk24(jrefn-1+1) = ma
    zk24(jrefn-1+2) = 'DEPL_R'
!
    call assvec('V', vecass, 1, lisvec, [1.d0], nu14)
    call detrsd('NUME_DDL', nu14)
!
    noeub = vecass
    call jeveuo(noeub//'.REFE', 'E', vk24=refe)
    refe(2) = nume_equa
!
    AS_ALLOCATE(vl=noeubord, size=nbno)
!
    call dismoi('NB_EC', 'DEPL_R', 'GRANDEUR', repi=nbec)
    call jeveuo(jexnum(nume_equa//'.PRNO', 1), 'L', jprno)
    call jeveuo(noeub//'.VALE', 'L', vr=val)
    eps = 1.d-06
    nbnob = 0
    do ino = 1, nbno
        numeq = zi(jprno-1+(nbec+2)*(ino-1)+1)
        nbcmp = zi(jprno-1+(nbec+2)*(ino-1)+2)
        xnorm = 0.d0
        do icmp = 1, nbcmp
            xnorm = xnorm+val(numeq-1+icmp)**2
        end do
        if (xnorm .le. eps) then
            noeubord(ino) = .false.
        else
            noeubord(ino) = .true.
            nbnob = nbnob+1
        end if
    end do
    AS_ALLOCATE(vi=nbpatchmil, size=nbno)
    AS_ALLOCATE(vi=numnb, size=nbnob)
    call wkvect('&&SINOZ2.NBPATCH', 'V V I', nbnob, jpa)
    inob = 0
    do ino = 1, nbno
        if (noeubord(ino)) then
            inob = inob+1
            numnb(inob) = ino
        end if
    end do
    if (inob .ne. nbnob) then
        call utmess('F', 'CALCULEL4_84')
    end if
!
!    VERIFICATION DES TYPES DE MAILLE
!
    ntri = 0
    nqua = 0
    do ima = 1, nbma
        iad = iatyma-1+ima
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(iad)), typema)
        if (typema(1:4) .eq. 'TRIA') then
            ntri = ntri+1
            if (typema(5:5) .ne. '3' .and. typema(5:5) .ne. '6') then
                call utmess('F', 'CALCULEL4_85')
            end if
        else if (typema(1:4) .eq. 'QUAD') then
            nqua = nqua+1
            if (typema(5:5) .ne. '4' .and. typema(5:5) .ne. '8' .and. typema(5:5) .ne. '9') then
                call utmess('F', 'CALCULEL4_86')
            end if
        end if
    end do
    if (ntri .ne. 0 .and. nqua .ne. 0) then
        call utmess('F', 'CALCULEL4_87')
    end if
    if (ntri .eq. 0 .and. nqua .eq. 0) then
        call utmess('F', 'CALCULEL4_88')
    end if
!
    mo = modele//'.MODELE    '
    call jeveuo(mo//'.REPE', 'L', vi=repe)
!
    rcmp(1) = 0.d0
    rcmp(2) = 0.d0
    rcmp(3) = 0.d0
    rcmp(4) = 0.d0
    licmp(1) = 'SIXX'
    licmp(2) = 'SIYY'
    licmp(3) = 'SIZZ'
    licmp(4) = 'SIXY'
    call crcnct('G', signo, ma, 'SIEF_R', 4, &
                licmp, rcmp)
    call jeveuo(signo(1:19)//'.VALE', 'E', vr=sig)
    call jeveuo(ma//'.COORDO    .VALE', 'L', vr=coor)
    call jeveuo(sigel(1:19)//'.CELD', 'L', vi=celd)
    call jeveuo(sigel(1:19)//'.CELV', 'L', vr=celv)
!
! --- RECUPERATION DE L'ELREFE ET DE LA FAMILLE POUR CHAQUE MAILLE
!
    call celfpg(sigel, elrfam, ibid)
    call jeveuo(elrfam, 'L', jelfa)
!
!   BOUCLE SUR LES NOEUDS
!   *********************
!
    ipa = 0
    do ino = 1, nbno
        if (.not. noeubord(ino)) then
            call jelira(jexnum(coninv, ino), 'LONMAX', nbmav)
!
!    TRAITEMENT DES SOMMETS
!
            if (nbmav .gt. 2) then
                ipa = ipa+1
                call jeveuo(jexnum(coninv, ino), 'L', iamav)
                call wkvect('&&SINOZ2.NOEBOPA', 'V V I', 10*nbmav, ianob)
!
!    INITIALISATION DE LA MATRICE A ET DU SECOND MEMBRE A ZERO
!
                a(:, :) = 0.d0
                b(:, :) = 0.d0
!
                nbnobp = 0
                ianew = ianob
!
!      PASSAGE EN COORDONNEES LOCALES AU PATCH (ENTRE -1. ET 1.)
!      POUR AMELIORER LE CONDITIONNEMENT DE LA MATRICE A
!
!      CALCUL DE XMIN,XMAX,YMIN,YMAX SUR LE PATCH
!
                xmin = 1.d+10
                xmax = -1.d+10
                ymin = 1.d+10
                ymax = -1.d+10
                do ima = 1, nbmav
                    numav = zi(iamav-1+ima)
                    call jelira(jexnum(connex, numav), 'LONMAX', nbnoma)
                    call jeveuo(jexnum(connex, numav), 'L', ianov)
                    do inoma = 1, nbnoma
                        num = zi(ianov-1+inoma)
                        x(inoma) = coor(3*(num-1)+1)
                        y(inoma) = coor(3*(num-1)+2)
                        if (x(inoma) .le. xmin) xmin = x(inoma)
                        if (y(inoma) .le. ymin) ymin = y(inoma)
                        if (x(inoma) .ge. xmax) xmax = x(inoma)
                        if (y(inoma) .ge. ymax) ymax = y(inoma)
                    end do
                end do
!
!      BOUCLE SUR LES MAILLES VOISINES
!
                do ima = 1, nbmav
                    numav = zi(iamav-1+ima)
                    numgr = repe(2*(numav-1)+1)
                    numel = repe(2*(numav-1)+2)
                    call jelira(jexnum(connex, numav), 'LONMAX', nbnoma)
                    call jeveuo(jexnum(connex, numav), 'L', ianov)
!
!        RECUPERATION DES COORDONNEES DES NOEUDS
!
                    do inoma = 1, nbnoma
                        num = zi(ianov-1+inoma)
!    SI NOEUD BORD
                        if (noeubord(num)) then
                            call zzappa(num, zi(ianew), nbnobp, app)
                            if (.not. app) then
                                nbnobp = nbnobp+1
                                zi(ianew-1+nbnobp) = num
                            end if
                        end if
                        x(inoma) = coor(3*(num-1)+1)
                        y(inoma) = coor(3*(num-1)+2)
                    end do
!
!        RECUPERATION DES FCTS DE FORME DE L'ELEMENT
                    elrefe = zk16(jelfa-1+numav) (1:8)
                    famil = zk16(jelfa-1+numav) (9:16)
                    call utelvf(elrefe, famil, nomjv, npg, nno)
                    call jeveuo(nomjv, 'L', ivf)
!
!    CALCUL DE LA MATRICE A
!
                    call zzcala(npg, nno, zr(ivf), x, y, &
                                xmin, xmax, ymin, ymax, a)
!
!    CALCUL DU SECOND MEMBRE B
!
                    call zzcalb(numgr, numel, npg, nno, zr(ivf), &
                                celd, celv, x, y, xmin, &
                                xmax, ymin, ymax, b)
!
                    call jedetr(nomjv)
!
                end do
                ianew = ianob+nbnobp
!
!   PRECONDITIONNEMENT PAR LA DIAGONALE ET RESOLUTION
!
                call predia(a, b, diag, nno)
                call mtcrou(a, b, 9, nno, 4, &
                            wk1, wk2)
!
                do ic = 1, 4
                    do i = 1, nno
                        b(i, ic) = b(i, ic)*diag(i)
                    end do
                end do
!
                xino = coor(3*(ino-1)+1)
                yino = coor(3*(ino-1)+2)
!
                xino = -1.d0+2.d0*(xino-xmin)/(xmax-xmin)
                yino = -1.d0+2.d0*(yino-ymin)/(ymax-ymin)
!
!   CALCUL DES CONTRAINTES LISSEES AU NOEUD INO
!
                call zzpoly(nno, ino, xino, yino, sig, &
                            b)
!
!    TRAITEMENT DES NOEUDS BORD DU PATCH
!
                do inob = 1, nbnobp
                    num = zi(ianob-1+inob)
                    do k = 1, nbnob
                        numc = numnb(k)
                        if (numc .eq. num) then
                            zi(jpa-1+k) = zi(jpa-1+k)+1
                        end if
                    end do
                    xinob = coor(3*(num-1)+1)
                    yinob = coor(3*(num-1)+2)
                    xinob = -1.d0+2.d0*(xinob-xmin)/(xmax-xmin)
                    yinob = -1.d0+2.d0*(yinob-ymin)/(ymax-ymin)
                    call zzpoly(nno, num, xinob, yinob, sig, &
                                b)
                end do
                call jedetr('&&SINOZ2.NOEBOPA')
!
!    TRAITEMENT DES NOEUDS MILIEUX (NON BORD)
!
                if (nno .ge. 6) then
                    do ima = 1, nbmav
                        numav = zi(iamav-1+ima)
                        call jelira(jexnum(connex, numav), 'LONMAX', nbnoma)
                        call jeveuo(jexnum(connex, numav), 'L', ianov)
!
!       RECUPERATION DU NOEUD MILIEU ASSOCIE AU NOEUD SOMMET DU PATCH
!             (1 SEUL PAR ELEMENT VOISIN)
!
                        do inoma = 1, nbnoma
                            num = zi(ianov-1+inoma)
                            if (num .eq. ino) then
                                numloc = inoma
                            end if
                        end do
!
                        if (nno .eq. 6) numloc = numloc+3
                        if (nno .eq. 8) numloc = numloc+4
                        if (nno .eq. 9) numloc = numloc+4
                        num = zi(ianov-1+numloc)
                        nbpatchmil(num) = nbpatchmil(num)+1
                        xinomi = coor(3*(num-1)+1)
                        yinomi = coor(3*(num-1)+2)
                        xinomi = -1.d0+2.d0*(xinomi-xmin)/(xmax-xmin)
                        yinomi = -1.d0+2.d0*(yinomi-ymin)/(ymax-ymin)
                        call zzpoly(nno, num, xinomi, yinomi, sig, &
                                    b)
                    end do
!
                end if
!
!    TRAITEMENT DU NOEUD BARYCENTRE (CAS DES Q9)
!
                if (nno .eq. 9) then
                    do ima = 1, nbmav
                        numav = zi(iamav-1+ima)
                        call jelira(jexnum(connex, numav), 'LONMAX', nbnoma)
                        call jeveuo(jexnum(connex, numav), 'L', ianov)
!
!     RECUPERATION DU NOEUD BARYCENTRE
!
                        num = zi(ianov-1+nbnoma)
                        nbpatchmil(num) = nbpatchmil(num)+1
                        xinomi = coor(3*(num-1)+1)
                        yinomi = coor(3*(num-1)+2)
                        xinomi = -1.d0+2.d0*(xinomi-xmin)/(xmax-xmin)
                        yinomi = -1.d0+2.d0*(yinomi-ymin)/(ymax-ymin)
                        call zzpoly(nno, num, xinomi, yinomi, sig, &
                                    b)
                    end do
                end if
!
            end if
        end if
    end do
!
!    MOYENNAGE SUR LES NOEUDS BORD
!
    do i = 1, nbnob
        num = numnb(i)
        if (zi(jpa-1+i) .eq. 0) then
            call utmess('F', 'CALCULEL4_89')
        end if
        do ic = 1, 4
            sig(4*(num-1)+ic) = sig(4*(num-1)+ic)/zi(jpa-1+i)
        end do
    end do
!
!    MOYENNAGE SUR LES NOEUDS NON BORD
!
    do ino = 1, nbno
!    SI PAS NOEUD BORD
        if (.not. noeubord(ino)) then
            if (nbpatchmil(ino) .eq. 0) then
                nbpatchmil(ino) = 1
            end if
            do ic = 1, 4
                sig(4*(ino-1)+ic) = sig(4*(ino-1)+ic)/nbpatchmil(ino)
            end do
        end if
    end do
    call detrsd('CHAMP_GD', '&&VECASS')
    AS_DEALLOCATE(vi=numnb)
    call jedetr('&&SINOZ2.NBPATCH')
    AS_DEALLOCATE(vi=nbpatchmil)
    AS_DEALLOCATE(vl=noeubord)
    call jedetr('&&SINOZ2.CONINV')
    AS_DEALLOCATE(vi=longconinv)
    call jedetr(elrfam)
!
    call jedema()
end subroutine
