! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine chrpel(champ1, repere, nbcmp, icham, type_cham, &
                  nomch, model, carele, lModelVariable)
! aslint: disable=W1501
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8dgrd.h"
#include "asterfort/angvxy.h"
#include "asterfort/assach.h"
#include "asterfort/assert.h"
#include "asterfort/calc_coor_elga.h"
#include "asterfort/calcul.h"
#include "asterfort/carelo.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesred.h"
#include "asterfort/cesvar.h"
#include "asterfort/chrgd.h"
#include "asterfort/chrpan.h"
#include "asterfort/copisd.h"
#include "asterfort/cylrep.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/matrot.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/normev.h"
#include "asterfort/reliem.h"
#include "asterfort/sepach.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/ddot.h"
    !
    integer :: nbcmp, icham
    character(len=*) :: champ1, repere, type_cham, nomch
    character(len=8) :: model, carele
    aster_logical, intent(in) :: lModelVariable
!
! --------------------------------------------------------------------------------------------------
!
!                       CHANGEMENT DE REPERE DANS LE CAS D'UN CHAM_ELEM
!
! --------------------------------------------------------------------------------------------------
!
!       champ1      : nom du champ a traiter (champ out)
!       repere      : type de repere (utilisateur ou cylindrique
!                         ou coque ou coque_util_intr ou coque_intr_util
!                         ou coque_util_cyl)
!       nbcmp       : nombre de composantes a traiter
!       icham       : numero d'occurrence
!       type_cham   : type du champ :'tens' 'vect' ou 'coque'
!       nomch       : nom de champ
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ii, jj, kk, ino, iad, ipt, isp
    integer :: jcesd, jcesv, jcesl, nbpt, ncmp
    integer :: ilcnx1, nbsp, inel, npain
    integer :: ibid, nbma, iret, inbno
    integer :: ndim, nbm, idmail, nbmail, imai
    integer :: inoeu, iret0, iret1, nbgno, igno, nncp
    integer :: ierk, ipaxe, ipaxe2
    integer :: nbno, nbpg, nuno, ipg
    integer :: type_pt
    integer :: iocc, nocc
    integer :: iexist, jcesd_gauss, jcesl_gauss, icoo
!
    integer, parameter :: type_unknown = 0, type_noeud = 1, type_gauss = 2
!   nb max de points (noeuds|gauss) par ??l??ment
    integer, parameter      :: nptmax = 30
    integer, dimension(6)   :: permvec
    real(kind=8)                    :: valr, xnormr, tmp
    real(kind=8), dimension(3)      :: xbary, angnot
    real(kind=8), dimension(3)      :: orig, axez, vectx, vecty, angrep
    real(kind=8), dimension(3,3)    :: pgl, pgcyl, pgu, pglelem
    real(kind=8), dimension(3, nptmax), target :: xno, xpg
    character(len=3)    :: tsca

    character(len=8)    :: ma, k8b, typmcl(2), nomgd, tych
    character(len=8)    :: lpain(5), paout, licmp(9), nomgdr, paoutc
    character(len=16)   :: option, motcle(2)
    character(len=19)   :: chams1, chams0, ligrel, canbsp
    character(len=19)   :: changl, carte, chr, chi, ch1, ch2
    character(len=19)   :: celgauss, cesgauss
    character(len=24)   :: mesmai, chgeom, lchin(5), chaout
    character(len=24)   :: valk(3), chcara(18)
!
    integer             :: jcesvrepso(3), jcesdrepso(3), adressev(3), elem_supp_num(2)
    integer             :: nbptii, nbspii, ncmpii, jmodemailsupp, mailtypel
    character(len=19)   :: chrel(3), chres(3), nomte
!
    integer, pointer                        :: connex(:)    => null()
    real(kind=8), pointer                   :: vale(:)      => null()
    real(kind=8), pointer                   :: coo_gauss(:) => null()
    real(kind=8), dimension(:,:), pointer   :: xpt          => null()
    character(len=8), pointer               :: nom_cmp(:)   => null()
    character(len=8), pointer               :: cesk(:)      => null()
!
    aster_logical :: exi_cmp, exi_local, okelem, effort_elno
! --------------------------------------------------------------------------------------------------
    call jemarq()
    ipaxe = 0
    motcle(1) = 'GROUP_MA'
    typmcl(1) = 'GROUP_MA'
    motcle(2) = 'MAILLE'
    typmcl(2) = 'MAILLE'
    canbsp = '&&CHRPEL.NBSP'
!
    mesmai = '&&CHRPEL.MES_MAILLES'
!
    if (nbcmp .gt. 0) then
        AS_ALLOCATE(vk8=nom_cmp, size=nbcmp)
        call getvtx('MODI_CHAM', 'NOM_CMP', iocc=icham, nbval=nbcmp, vect=nom_cmp, nbret=ibid)
    else
        call utmess('F', 'ALGORITH2_6')
    endif
!
    call dismoi('NOM_LIGREL', champ1, 'CHAM_ELEM', repk=ligrel)
!
!   DEFINITION ET CREATION DU CHAM_ELEM SIMPLE CHAMS1 A PARTIR DU CHAM_ELEM CHAMP1
    chams0='&&CHRPEL.CHAMS0'
    chams1='&&CHRPEL.CHAMS1'
    call celces(champ1, 'V', chams0)
    call cesred(chams0, 0, [0], nbcmp, nom_cmp, 'V', chams1)
    call detrsd('CHAM_ELEM_S', chams0)
    call jeveuo(chams1//'.CESK', 'L', vk8=cesk)
    call jeveuo(chams1//'.CESD', 'L', jcesd)
    ma = cesk(1)
    nomgd = cesk(2)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
!
!   ON EXCLUT LES MOT-CLES 'NOEUD' ET 'GROUP_NO'
    call jeveuo(ma//'.DIME   ', 'L', inbno)
    call wkvect('&&CHRPEL.NOEUDS', 'V V K8', zi(inbno), inoeu)
    call getvtx('AFFE', 'NOEUD', iocc=icham, nbval=0, nbret=iret0)
    call jeexin(ma//'.GROUPENO', ierk)
    if (ierk .ne. 0) then
        call jelira(ma//'.GROUPENO', 'NMAXOC', nbgno)
        call wkvect('&&CHRPEL.GROUP_NO', 'V V K24', nbgno, igno)
        call getvtx('AFFE', 'GROUP_NO', iocc=icham, nbval=0, nbret=iret1)
    else
        iret1=0
    endif
    if ( (iret0.lt.0) .or. (iret1.lt.0) ) then
        valk(1) = 'NOEUD ou GROUP_NO'
        valk(2) = nomch
        valk(3) = ' '
        call utmess('F', 'ALGORITH12_42', nk=3, valk=valk)
    endif
    call jedetr('&&CHRPEL.NOEUDS')
    call jedetr('&&CHRPEL.GROUP_NO')
!   nombre total de mailles dans le maillage
    nbma = zi(jcesd-1+1)
!   nombre de composantes du champ simple chams1
    ncmp = zi(jcesd-1+2)
!   chams1 a ??t?? cr???? ?? partir des composantes s??lectionn??es par l'utilisateur, on doit avoir:
    ASSERT( ncmp == nbcmp )
!   D??termination de la dimension ?? partir du maillage
    ndim = 3
    call dismoi('Z_CST', ma, 'MAILLAGE', repk=k8b)
    if (k8b .eq. 'OUI') ndim = 2
!
    call jeexin(ma//'.CONNEX', iret)
    ASSERT(iret.ne.0)
    call jeveuo(ma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
    call jeveuo(chams1//'.CESV', 'E', jcesv)
    call jeveuo(chams1//'.CESL', 'L', jcesl)
!
    effort_elno = (nomch.eq.'EFGE_ELNO').or.(nomch.eq.'SIEF_ELNO')
!
    if ( effort_elno.and.((type_cham.ne.'VECTR_3D').and.(type_cham.ne.'COQUE_GENE')) ) then
        call utmess('F', 'ALGORITH2_35', sk=type_cham)
    endif
!   Si le champ est exprim?? dans le rep??re local des ??l??ments
!       Construction du champ des rep??re locaux
    exi_local = .false.
    if ( effort_elno.and.(type_cham.eq.'VECTR_3D') ) then
!       Les noms des composantes : N, VX, VY, MT, MFY, MFZ DANS CETTE ORDRE
        if ( (nom_cmp(1).ne.'N') .or.(nom_cmp(2).ne.'VY') .or.(nom_cmp(3).ne.'VZ') .or. &
             (nom_cmp(4).ne.'MT').or.(nom_cmp(5).ne.'MFY').or.(nom_cmp(6).ne.'MFZ') ) then
            call utmess('F', 'ALGORITH2_33')
        endif
        chrel(1) = '&&CHRPEL.REPLO_1'; chrel(2) = '&&CHRPEL.REPLO_2'; chrel(3) = '&&CHRPEL.REPLO_3'
        chres(1) = '&&CHRPEL.REPSO_1'; chres(2) = '&&CHRPEL.REPSO_2'; chres(3) = '&&CHRPEL.REPSO_3'
        call carelo(model, carele, 'V', chrel(1), chrel(2), chrel(3))
        exi_local = .true.
!
        do ii = 1 , 3
            call celces(chrel(ii), 'V', chres(ii))
            call jeveuo(chres(ii)//'.CESV', 'L', jcesvrepso(ii))
            call jeveuo(chres(ii)//'.CESD', 'L', jcesdrepso(ii))
            call detrsd('CHAM_ELEM', chrel(ii))
        enddo
        ! Seulement pour les POU_D_E, POU_D_T
        call jenonu(jexnom('&CATA.TE.NOMTE', 'MECA_POU_D_T'), elem_supp_num(1))
        call jenonu(jexnom('&CATA.TE.NOMTE', 'MECA_POU_D_E'), elem_supp_num(2))
        ! Pointeur sur les ??l??ments supports du mod??le
        call jeveuo(model//'.MAILLE', 'L', jmodemailsupp)
    endif
!   Le mot-cl?? AFFE d??finit les caract??ristiques du nouveau rep??re
!   On peut d??finir un rep??re variable en d??finissant ces param??tres par mailles/groupes de mailles
    call getfac('AFFE', nocc)
!   Boucle sur les occurrences de AFFE
    do iocc = 1, nocc
!       Construction de la liste des num??ros de mailles
!       s??lectionn??es par les mots-cl??s GROUP_MA et MAILLE
        call reliem(' ', ma, 'NU_MAILLE', 'AFFE', iocc, 2, motcle, typmcl, mesmai, nbm)
        if (nbm .gt. 0) then
            nbmail = nbm
            call jeveuo(mesmai, 'L', idmail)
        else
            nbmail = nbma
        endif
!       Remise ?? z??ro des tableaux stockant les caract??ristiques du rep??re
        axez(:) = 0.0d0; orig(:) = 0.0d0; angnot(:) = 0.0d0
!       Changement de rep??re "UTILISATEUR"
        if (repere(1:11) .eq. 'UTILISATEUR') then
!           SI LE NOUVEAU REPERE EST DONNE VIA DES VECTEURS
            call getvr8('AFFE', 'VECT_X', iocc=iocc, nbval=3, vect=vectx, nbret=ibid)
            if (ibid .ne. 0) then
                call getvr8('AFFE', 'VECT_Y', iocc=iocc, nbval=3, vect=vecty, nbret=ibid)
                if (ndim .ne. 3) then
                    call utmess('F', 'ALGORITH2_4')
                endif
                call angvxy(vectx, vecty, angnot)
            else
                if (ndim .eq. 3) then
                    call getvr8('AFFE', 'ANGL_NAUT', iocc=iocc, nbval=3, vect=angnot, nbret=ibid)
                    if (ibid .ne. 3) then
                        call utmess('F', 'ALGORITH2_7')
                    endif
                else
                    call getvr8('AFFE', 'ANGL_NAUT', iocc=iocc, scal=angnot(1), nbret=ibid)
                    if (ibid .ne. 1) then
                        valr = angnot(1)
                        call utmess('A', 'ALGORITH12_43', sr=valr)
                    endif
                endif
                angnot(:) = angnot(:)*r8dgrd()
            endif
!
!           Matrice de passage du rep??re global vers le rep??re utilisateur
            call matrot(angnot, pgl)
!           Matrot retourne la transpos??e de la matrice de passage : on transpose pour avoir
!           la matrice de passage
            pgu=transpose(pgl)
!
!           Appliquer le changement de rep??re pour les mailles s??lectionn??es
!
            cinel: do inel = 1, nbmail
                if (nbm .ne. 0) then
                    imai = zi(idmail+inel-1)
                else
                    imai = inel
                endif
!               Si le champ est dans le rep??re local de l'??l??ment on va chercher la matrice
!               de passage du rep??re local au global
                if ( exi_local ) then
!                   Est-ce un ??l??ment autoris?? ?
                    okelem = .false.
                    mailtypel = zi(jmodemailsupp-1+imai)
                    cii: do ii = 1, 2
                        if ( mailtypel.eq.elem_supp_num(ii) ) then
                            okelem = .true.
                            exit cii
                        endif
                    enddo cii
                    if ( .not. okelem ) then
                        call jenuno(jexnum('&CATA.TE.NOMTE', mailtypel), nomte)
                        call utmess('A', 'ALGORITH2_34', sk=nomte)
                        exit cinel
                    endif
!                   V??rification du nombre : de point, de sous points, des composantes
!                   R??cup??ration de l'adresse des valeurs des composantes
                    do ii = 1 , 3
                        nbptii = zi(jcesdrepso(ii)-1+5+4*(imai-1)+1)
                        nbspii = zi(jcesdrepso(ii)-1+5+4*(imai-1)+2)
                        ncmpii = zi(jcesdrepso(ii)-1+5+4*(imai-1)+3)
                        ASSERT( (nbptii.eq.1).and.(nbspii.eq.1).and.(ncmpii.eq.3) )
                        adressev(ii) = jcesvrepso(ii)-1 + zi(jcesdrepso(ii)-1+5+4*(imai-1)+4)
                    enddo
!                   La matrice de changement de rep??re li?? ?? l'??l??ment (optenu par matrot)
                    do ii = 1, 3
                        pglelem(1,ii)  = zr(adressev(1)+ii)
                        pglelem(2,ii)  = zr(adressev(2)+ii)
                        pglelem(3,ii)  = zr(adressev(3)+ii)
                    enddo
!                   Passage dans le rep??re Global           Fglob = transpose(pglelem) . Floc
!                   Passage dans le rep??re utilisateur      Futil = pgl . Fglob
!                   Au Final                                Futil = pgl . transpose(pglelem) . Floc
                    do ii = 1,3
                        do jj = 1,3
                            tmp = 0.0
                            do kk=1,3
                                tmp = tmp + pgl(ii,kk)*pglelem(jj,kk)
                            enddo
                            pgu(ii,jj) = tmp
                        enddo
                    enddo
                endif
                nbpt = zi(jcesd-1+5+4*(imai-1)+1)
                nbsp = zi(jcesd-1+5+4*(imai-1)+2)
                cipt1: do ipt = 1, nbpt
                    do isp = 1, nbsp
                        exi_cmp = .false.
                        do ii = 1, nbcmp
                            call cesexi('S', jcesd, jcesl, imai, ipt, isp, ii, iad)
                            if (iad .gt. 0) then
                                exi_cmp = .true.
                            endif
                        enddo
                        if (exi_cmp) then
                            call chrgd(nbcmp, jcesd, jcesl, jcesv, imai, &
                                       ipt, isp, type_cham, tsca, pgu)
                        else
                            cycle cipt1
                        endif
                    enddo
                enddo cipt1
            enddo cinel
!
!       Changement de rep??re "CYLINDRIQUE"
        else if (repere(1:11).eq.'CYLINDRIQUE') then
!
            if (type_cham.eq.'VECTR_3D') then
                call utmess('F', 'ALGORITH2_31')
            endif
!
            call dismoi('TYPE_CHAMP', champ1, 'CHAMP', repk=tych, arret='C', ier=iret)
            if (ndim .eq. 3) then
                call getvr8('AFFE', 'ORIGINE', iocc=iocc, nbval=3, vect=orig,nbret=ibid)
                if (ibid .ne. 3) then
                    call utmess('F', 'ALGORITH2_8')
                endif
                call getvr8('AFFE', 'AXE_Z', iocc=iocc, nbval=3, vect=axez, nbret=ibid)
                if (ibid .eq. 0) then
                    call utmess('F', 'ALGORITH2_9')
                endif
            else
                call getvr8('AFFE', 'ORIGINE', iocc=iocc, nbval=2, vect=orig, nbret=ibid)
                if (ibid .ne. 2) then
                    call utmess('A', 'ALGORITH2_10')
                endif
                call getvr8('AFFE', 'AXE_Z', iocc=iocc, nbval=0, nbret=ibid)
                if (ibid .ne. 0) then
                    call utmess('A', 'ALGORITH2_11')
                endif
                axez(1) = 0.0d0
                axez(2) = 0.0d0
                axez(3) = 1.0d0
            endif
            xnormr = 0.0d0
            call normev(axez, xnormr)
            call jeveuo(ma//'.COORDO    .VALE', 'L', vr=vale)
!
!           Permutation des composantes en dimension 2
!           Initialisation ?? l'identit??
            permvec(:)=(/(ii,ii=1,6)/)
            if (ndim == 2) then
                select case (type_cham(1:4))
                    case('TENS')
                        permvec(4) = 5
                    case('VECT')
                        permvec(2) = 3
                        permvec(3) = 2
                end select
            endif
!
!           Localisation du champ : noeuds/pts de Gauss
            type_pt = type_unknown
            if (tych(1:4) == 'VECT') then
                type_pt = type_noeud
            endif
            if (tych(1:4) == 'ELNO') then
                type_pt = type_noeud
            else if (tych(1:4) == 'ELGA') then
                type_pt = type_gauss
            endif
            ASSERT(type_pt /= type_unknown )
!
!           Si le champ est un champ 'ELGA', on a besoin des
!           coordonn??es des points de Gauss dans chaque ??l??ment
            if ( type_pt == type_gauss) then
                if (lModelVariable) then
                    call utmess('F', 'RESULT4_90')
                endif
!               On utilise calc_coor_elga qui retourne un champ par ??l??ment
!               contenant les coordonn??es des points de Gauss
                celgauss='&&CHRPEL.CEL_GAUSS'
                call exisd('CHAMP', celgauss, iexist)
                if (iexist .eq. 0) then
                    call megeom(model, chgeom)
                    call calc_coor_elga(model, ligrel, chgeom, celgauss)
                endif
!               On transforme ce champ en champ simple
                cesgauss='&&CHRPEL.CES_GAUSS'
                call celces( celgauss, 'V', cesgauss)
                call jeveuo(cesgauss//'.CESD','L', jcesd_gauss)
                call jeveuo(cesgauss//'.CESL','L', jcesl_gauss)
                call jeveuo(cesgauss//'.CESV','L', vr=coo_gauss)
            endif

!           Boucle sur les mailles ?? transformer
            do inel = 1, nbmail
!               R??cup??ration de imai : indice de la maille courante
!                               nbno : nombre de noeuds de la maille courante
                if (nbm .ne. 0) then
                    imai = zi(idmail+inel-1)
                else
                    imai = inel
                endif
                nbno = zi(ilcnx1+imai) - zi(ilcnx1-1+imai)
!               Quelques caract??ristiques du champ simple ?? transformer sur cette maille :
!               nbpg : nombre de points de Gauss
                nbpg = zi(jcesd-1+5+4* (imai-1)+1)
!               nbsp : nombre de sous-points
                nbsp = zi(jcesd-1+5+4* (imai-1)+2)
!               nbcmp : nombre de composantes
                nbcmp = zi(jcesd-1+5+4* (imai-1)+3)
!
!               Coordonn??es des noeuds de la maille courante
                xno(:,:) = 0.d0
                do ino = 1, nbno
                    nuno = connex(zi(ilcnx1+imai-1)+ino-1)
                    xno(1,ino) = vale(1+3*(nuno-1)-1+1)
                    xno(2,ino) = vale(1+3*(nuno-1)-1+2)
                    if (ndim == 3) then
                        xno(3,ino) = vale(1+3*(nuno-1)-1+3)
                    endif
                enddo
!
                select case (type_pt)
                    case (type_noeud)
!                       Noeuds de la maille
                        nbpt=nbno
                        xpt => xno(:,:)
                    case (type_gauss)
!                       Points de Gauss de la maille, dont il faut r??cup??rer les coordonn??es
                        do ipg = 1, nbpg
                            do icoo = 1, 3
                                call cesexi('S', jcesd_gauss, jcesl_gauss, imai, ipg, 1, icoo, iad)
                                xpg(icoo,ipg) = coo_gauss(iad)
                            enddo
                        enddo
                        nbpt = nbpg
                        xpt => xpg(:,:)
                    case default
                        nbpt = 0
                        ASSERT(.false.)
                end select
!
!               Boucle sur les points (cette partie est commune aux champs ELNO et ELGA)
                cipt2: do ipt = 1, nbpt
!                   Calcul de la matrice de passage vers le rep??re cylindrique
                    call cylrep(ndim, xpt(:, ipt), axez, orig, pgcyl, ipaxe)
!                   Si le point x appartient ?? l'axe du rep??re cylindrique
                    if (ipaxe > 0) then
                        call utmess('A', 'ALGORITH2_13')
!                       Calcul de la matrice de passage au centre de gravit?? de l'??l??ment
                        xbary(:)=sum(xno(:,1:nbno), dim=2)
                        xbary(:) = xbary(:)/nbno
                        ipaxe2 = 0
                        call cylrep(ndim, xbary, axez, orig, pgcyl, ipaxe2)
!                       Si le centre de gravit?? de l'??l??ment est aussi sur l'axe, on s'arr??te
                        if (ipaxe2 > 0) then
                            call utmess('F', 'ALGORITH2_13')
                        endif
                    endif
!                   Boucle sur les sous-points
                    do isp = 1, nbsp
                        exi_cmp = .true.
                        do ii = 1, nbcmp
!                           la composante ii du champ existe-t-elle?
                            exi_cmp = .false.
                            call cesexi('S', jcesd, jcesl, imai, ipt, isp, ii, iad)
                            if (iad .gt. 0) then
                                exi_cmp = .true.
                            endif
                        enddo
                        if (exi_cmp) then
                            call chrgd(nbcmp, jcesd, jcesl, jcesv, imai, &
                                       ipt, isp, type_cham, tsca, pgcyl, permvec)
                        else
                            cycle cipt2
                        endif
                    enddo
                enddo cipt2
            enddo
!
            if (ipaxe .ne. 0) then
                call utmess('A', 'ALGORITH17_22', si=ipaxe)
            endif
        endif
!
        call jeexin(mesmai, iret)
        if (iret .ne. 0) call jedetr(mesmai)
    enddo
!
    if ( (repere(1:11).eq.'CYLINDRIQUE') .or. &
         (repere(1:11).eq.'UTILISATEUR') ) then
!       Champ simple -> Cham_elem
        call dismoi('NOM_OPTION', champ1, 'CHAM_ELEM', repk=option)
        call cescel(chams1, ligrel, option, ' ', 'OUI', nncp, 'G', champ1, 'F', ibid)
        call detrsd('CHAM_ELEM_S', chams1)
        if ( exi_local ) then
            do ii= 1, 3
                call detrsd('CHAM_ELEM_S', chres(ii))
            enddo
        endif
    endif
!
! --------------------------------------------------------------------------------------------------
!   Changement de rep??re sur une coque
    if ( (repere(1:5) .eq.'COQUE')           .or. &
         (repere(1:15).eq.'COQUE_INTR_UTIL') .or. &
         (repere(1:15).eq.'COQUE_UTIL_INTR') .or. &
         (repere(1:14).eq.'COQUE_UTIL_CYL') )  then
!       Pour l'instant on ne traite pas le cas de plusieurs occurrences du mot-cl?? AFFE
        if ( nocc /= 1  ) then
            call utmess('F', 'ALGORITH17_23', sk=repere, si=nocc)
        endif
!
        call megeom(model, chgeom)
        call mecara(carele, chcara)
!
        if ( (type_cham(1:10).eq.'COQUE_GENE') .and. &
             (repere(1:14).eq.'COQUE_UTIL_CYL') ) then
            call utmess('F', 'ELEMENTS5_55', nk=2, valk=(/'COQUE_UTIL_CYL', 'COQUE_GENE    '/))
        endif
        if (type_cham(1:10).eq.'COQUE_GENE') then
            option = 'REPE_GENE'
!           Nb de param??tres en entr??e de l'option
            npain = 4
        else if (type_cham(1:7).eq.'TENS_3D') then
            option = 'REPE_TENS'
            npain = 5
        else
            call utmess('F', 'ELEMENTS5_53', sk=type_cham)
        endif
!
!       GENERATION D UN CHAMP D'ANGLES (CARTE CONSTANTE)
        carte = '&&CHRPEL.ANGL_REP'
        angrep(:) = 0.0d0
!
        if (repere .eq. 'COQUE_INTR_UTIL') then
            angrep(3)=1.d0
        else if (repere.eq.'COQUE_UTIL_INTR') then
            angrep(3)=2.d0
        else if (repere.eq.'COQUE_UTIL_CYL') then
            angrep(3)=3.d0
        endif
        licmp(1) = 'ALPHA'
        licmp(2) = 'BETA'
        licmp(3) = 'REP'
        licmp(4) = 'AXE_X'
        licmp(5) = 'AXE_Y'
        licmp(6) = 'AXE_Z'
        licmp(7) = 'O_X'
        licmp(8) = 'O_Y'
        licmp(9) = 'O_Z'
        call mecact('V', carte, 'MODELE', model, 'CAORIE', ncmp=9, lnomcmp=licmp, vr=angrep)
!
!       CREATION D UN CHAM_ELEM D'ANGLES EN LISANT LES ANGL_REP
        changl = '&&CHRPEL.ANGL'
        call chrpan(model, carte, option, changl)
!
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PCACOQU'
        lchin(2) = chcara(7)
        lpain(3) = 'PANGREP'
        lchin(3) = changl
        lchin(4) = champ1
        lpain(5) = 'PNBSP_I'
        lchin(5) = chcara(16)
!
        call dismoi('NOM_GD', lchin(4), 'CHAMP', repk=nomgdr)
!
        if (type_cham .eq. 'COQUE_GENE') then
            if (nomch .eq. 'EFGE_ELGA') then
                lpain(4) = 'PEFGAIN'
                paout = 'PEFGAOUT'
                if (nomgdr(5:6) .eq. '_C') paoutc = 'PEFGAOUC'
            else if (nomch.eq.'EFGE_ELNO') then
                lpain(4) = 'PEFNOIN'
                paout = 'PEFNOOUT'
                if (nomgdr(5:6) .eq. '_C') paoutc = 'PEFNOOUC'
            else if (nomch.eq.'DEGE_ELGA') then
                lpain(4) = 'PDGGAIN'
                paout = 'PDGGAOUT'
                if (nomgdr(5:6) .eq. '_C') paoutc = 'PDGGAOUC'
            else if (nomch.eq.'DEGE_ELNO') then
                lpain(4) = 'PDGNOIN'
                paout = 'PDGNOOUT'
                if (nomgdr(5:6) .eq. '_C') paoutc = 'PDGNOOUC'
            else if (nomch.eq.'SIEF_ELGA') then
                lpain(4) = 'PEFGAIN'
                paout = 'PEFGAOUT'
                if (nomgdr(5:6) .eq. '_C') paoutc = 'PEFGAOUC'
            else
                call utmess('F', 'ELEMENTS5_51', sk=nomch)
            endif
        else if (type_cham.eq.'TENS_3D') then
            if (nomch .eq. 'SIGM_ELGA') then
                lpain(4) = 'PCOGAIN'
                paout = 'PCOGAOUT'
            else if (nomch.eq.'SIGM_ELNO') then
                lpain(4) = 'PCONOIN'
                paout = 'PCONOOUT'
            else if (nomch.eq.'EPSI_ELGA') then
                lpain(4) = 'PDEGAIN'
                paout = 'PDEGAOUT'
            else if (nomch.eq.'EPSI_ELNO') then
                lpain(4) = 'PDENOIN'
                paout = 'PDENOOUT'
            else
                call utmess('F', 'ELEMENTS5_52', sk=nomch)
            endif
            !
        endif
        call exisd('CHAM_ELEM_S', canbsp, iret1)
        if (iret1 .ne. 1) then
            call dismoi('MXNBSP', champ1(1:19), 'CHAM_ELEM', repi=nbsp)
!
!           SI LE CHAMP A DEJA ETE EXTRAIT IL FAUT APPELER CESVAR AVEC CE CHAMP
            if (nbsp .eq. 1) then
                call cesvar(champ1(1:19), ' ', ligrel, canbsp)
            else
                call cesvar(carele, ' ', ligrel, canbsp)
            endif
        endif
        chaout = chams1
        call copisd('CHAM_ELEM_S', 'V', canbsp, chaout)
!
        if (nomgdr(5:6) .eq. '_C') then
            chr='&&CHRPEL.CHR'
            chi='&&CHRPEL.CHI'
            ch1='&&CHRPEL.CH1'
            ch2='&&CHRPEL.CH2'
            call sepach(carele, lchin(4), 'V', chr, chi)
            lchin(4)=chr
            call calcul('S', option, ligrel, npain, lchin(1:npain), &
                        lpain(1:npain), 1, ch1, paout, 'V', 'OUI')
            lchin(4)=chi
            call calcul('S', option, ligrel, npain, lchin(1:npain),&
                        lpain(1:npain), 1, ch2, paout, 'V', 'OUI')
            call assach(ch1, ch2, 'V', chaout, parout=paoutc)
            call detrsd('CHAMP', chr)
            call detrsd('CHAMP', chi)
            call detrsd('CHAMP', ch1)
            call detrsd('CHAMP', ch2)
        else
            call calcul('S', option, ligrel, npain, lchin(1:npain),&
                        lpain(1:npain), 1, chaout, paout, 'V', 'OUI')
        endif
        call detrsd('CHAM_ELEM_S', chaout)
        call copisd('CHAMP_GD', 'G', chaout, champ1)
    endif
!
    AS_DEALLOCATE(vk8=nom_cmp)
    call exisd('CHAM_ELEM_S', canbsp, iret1)
    if (iret1 .ne. 0) call detrsd('CHAM_ELEM_S', canbsp)

    call jedema()
    !
end subroutine chrpel
