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
! person_in_charge: samuel.geniaut at edf.fr
! aslint: disable=W1504, W1501
!
subroutine xpomax(mo, malini, mailx, nbnoc, nbmac, &
                  prefno, nogrfi, maxfem, cns1, cns2, &
                  ces1, ces2, cesvi1, cesvi2, listgr, &
                  dirgrm, nivgrm, resuco, ngfon, comps1, &
                  comps2, pre1, iord)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/conare.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/elref2.h"
#include "asterfort/elrfvf.h"
#include "asterfort/exisd.h"
#include "asterfort/iselli.h"
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
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "asterfort/xismec.h"
#include "asterfort/xpoajc.h"
#include "asterfort/xpoajm.h"
#include "asterfort/xpocmp.h"
#include "asterfort/xpocox.h"
#include "asterfort/res2mat.h"
#include "asterfort/fointe.h"
#include "asterfort/rccome.h"
#include "asterfort/indk32.h"
#include "asterfort/rsadpa.h"
!
    integer(kind=8) :: nbnoc, nbmac, ngfon, iord
    character(len=2) :: prefno(4)
    character(len=8) :: mo, malini, maxfem, resuco
    character(len=19) :: cns1, cns2, ces1, ces2, cesvi1, cesvi2
    character(len=19) :: comps1, comps2
    character(len=24) :: mailx, listgr, dirgrm, nivgrm, gpptnn, nogrfi
!
!      TRAITEMENT DES MAILLES DE MAILX
!       - POUR POST_MAIL_XFEM : CREATION DES MAILLES CORRESPONDANTS
!                               AUX SOUS-ELEMENTS ET CREATION DES NOEUDS
!                               CORRESPONDANTS AUX SOMMETS DES SOUS
!                               ELEMENTS
!       - POUR POST_CHAM_XFEM : CALCUL DES DEPLACEMENTS AUX NOUVEAUX
!                               NOEUDS DE MAXFEM
!
!
!   IN
!       MO     : MODELE FISSURE
!       MALINI : MAILLAGE SAIN
!       MAILX  : LISTE DES NUMEROS DES MAILLES SOUS-DECOUPEES
!       NBNOC  : NOMBRE DE NOEUDS CLASSIQUES DU MAILLAGE FISSURE
!       NBMAC  : NOMBRE DE MAILLES CLASSIQUES DU MAILLAGE FISSURE
!       PREFNO : PREFERENCES POUR LE NOMAGE DES NOUVELLES ENTITES
!       NOGRFI : NOM DU GROUPE DE NOEUDS SUR LA FISSURE A CREER
!       MAXFEM : MAILLAGE FISSURE (SI POST_CHAMP_XFEM)
!       CNS1   : CHAMP_NO_S DU DEPLACEMENT EN ENTREE
!       CES1   : CHAMP_ELEM_S DE CONTRAINTES EN ENTREE
!       LISTGR : LISTE DES GROUPES CONTENANT CHAQUE MAILLE
!       DIRGRM : VECTEUR D'INDIRECTION ENTRE LES GROUP_MA
!       NIVGRM : VECTEUR DE REMPLISSAGE DES GROUP_MA DE MAXFEM
!       RESUCO : NOM DU CONCEPT RESULTAT DONT ON EXTRAIT LES CHAMPS
!       NGFON  : NOMBRE TOTAL DE FOND DE FISSURES
!       COMPS1 : CHAM_ELEM_S DU COMPORTEMENT EN ENTREE
!   OUT
!       MAXFEM : MAILLAGE FISSURE (SI POST_MAIL_XFEM)
!       CNS2   : CHAMP_NO_S DU DEPLACEMENT EN SORTIE
!       CES2   : CHAMP_ELEM_S DE CONTRAINTES EN SORTIE
!       NIVGRM : VECTEUR DE REMPLISSAGE DES GROUP_MA DE MAXFEM
!       COMPS2 : CHAM_ELEM_S DU COMPORTEMENT EN SORTIE
!
!
!
!
!
    integer(kind=8) :: nbch, ddlmax
    parameter(nbch=17, ddlmax=52)
!
    integer(kind=8) :: i, ier, jmax, nbmax, ich, ima, nse, ise, in, l, jbaslo, ncmp
    integer(kind=8) :: jcesd(nbch), jcesv(nbch), jcesl(nbch), iad, jconx1, jconx2
    integer(kind=8) :: ifisc, jfiss, kfiss
    integer(kind=8) :: j, ino, n, jdirno, jlsn, inn, nnn, nbnoma, nfiss, ifiss
    integer(kind=8) :: iacoo1, iacoo2, ndim, iad2, inntot, ndime, inm, nfimax
    integer(kind=8) :: jtypm2, inmtot, itypse(6), iad1, iadc, iadv, heavno(20, 3)
    parameter(nfimax=10)
    integer(kind=8) :: jcnse, iad4, iad3, itypel, nbelr, jhea, fisc(2*nfimax)
    integer(kind=8) :: igeom, nfh, ifh, nfe, ddlc, cmp(ddlmax), jlst, jheavn, fisco(2*nfimax)
    integer(kind=8) :: nbcmp, jcnsv1, jcnsv2, nbnofi, inofi, nfijon(3)
    integer(kind=8) :: jcnsl2, jcesv1, jcesd1, jcesl1, jcesv2, jcesd2, jcesl2, lacthm(16)
    integer(kind=8) :: jcviv1, jcvid1, jcvil1, jcviv2, jcvid2, jcvil2, ninter(4), vit(16)
    integer(kind=8) :: jnivgr, iagma, ngrm, jdirgr, iagno, iad10, iad11, npg, nlachm(2)
    integer(kind=8) :: nbar, ar(12, 3), dec, ino1, ino2, iar, iad16, nli, ncompi, nvit
    real(kind=8) :: lsno, ff(20), pinter(3)
    integer(kind=8) :: jstno
    integer(kind=8) :: nuflpg, nblfpg, nufgpg
    character(len=1) :: kbid
    character(len=8) :: k8b, typese(6), elrefp, lirefe(10), elrese(6)
    character(len=8) :: typma, noma, chmat
    character(len=16) :: tysd, k16b, nomcmd, notype
    character(len=19) :: chs(nbch), varcns
    character(len=24) :: dirno, geom, linofi, grpnoe, lsn, lst, hea, nogno, heavn, basloc, stano
    character(len=32) :: noflpg
    aster_logical :: opmail, lmeca, pre1
    integer(kind=8) :: iad9, irese, nnose, tabse(6), ncomp, ncompn, ncompn_tmp
    integer(kind=8) :: iviex, iret, jconq1, jconq2, jxc, ncompa, nfisc, nfisc2
    integer(kind=8) :: jresd1, jresv1, jresl1, nbcmpc, jresd2, jresv2, jresl2
    integer(kind=8), pointer :: maille(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
    integer(kind=8), pointer :: tmdim(:) => null()
    character(len=8), pointer :: cnsk(:) => null()
    integer(kind=8), pointer :: typm1(:) => null()
    integer(kind=8), pointer :: vtypma(:) => null()
    character(len=32), pointer :: pnlocfpg(:) => null()
    integer(kind=8), pointer :: tmfpg(:) => null()
    integer(kind=8), pointer :: nolocfpg(:) => null()
    real(kind=8) :: ka, mu, inst
    aster_logical :: cplan, lvarc, young, poiss
    integer(kind=8) :: jinst1
    character(len=8) :: nommat
    character(len=11) :: k11
    real(kind=8) :: varc(2), e, nu
    character(len=8), pointer :: cvrcvarc(:) => null()
    character(len=16), pointer :: valk(:) => null()
    integer(kind=8) :: jcesd_varc, jcesv_varc, jcesl_varc, nbvarc, k, nbf, ik, nbr, nbc, nbk

!
    data typese/'SEG2', 'TRIA3', 'TETRA4', 'SEG3', 'TRIA6', 'TETRA10'/
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data tabse/2, 3, 4, 3, 6, 10/
!
    call jemarq()
!
    call jeexin(mailx, ier)
    if (ier .eq. 0) goto 999
!
!     ------------------------------------------------------------------
!                   RECUPERATION DES OBJETS JEVEUX
!     ------------------------------------------------------------------
!
    call dismoi('DIM_GEOM', malini, 'MAILLAGE', repi=ndim)
!
!     NOM DE LA COMMANDE (POST_MAIL_XFEM OU POST_CHAM_XFEM)
    call getres(k8b, k16b, nomcmd)
    if (nomcmd .eq. 'POST_MAIL_XFEM') then
        opmail = .true.
    else if (nomcmd .eq. 'POST_CHAM_XFEM') then
        opmail = .false.
    end if
!
    call jeveuo(mailx, 'L', jmax)
    call jelira(mailx, 'LONMAX', nbmax)
!
    chs(1) = '&&XPOMAX.PINTTO'
    chs(2) = '&&XPOMAX.CNSETO'
    chs(3) = '&&XPOMAX.LONCHA'
    chs(4) = '&&XPOMAX.HEAV'
    chs(6) = '&&XPOMAX.CNSLN'
    chs(7) = '&&XPOMAX.CNSLT'
    chs(8) = '&&XPOMAX.CNSFI'
    chs(9) = '&&XPOMAX.PMILTO'
    chs(10) = '&&XPOMAX.PLONCH'
    chs(11) = '&&XPOMAX.PAIN'
    chs(12) = '&&XPOMAX.FHEAVN'
    chs(13) = '&&XPOMAX.HEAVN0'
    chs(14) = '&&XPOMAX.FISSC0'
    chs(15) = '&&XPOMAX.PINT'
    chs(16) = '&&XPOMAX.BASLOC'
    chs(17) = '&&XPOMAX.STANO'
!
!
    call celces(mo//'.TOPOSE.PIN', 'V', chs(1))
    call celces(mo//'.TOPOSE.CNS', 'V', chs(2))
    call celces(mo//'.TOPOSE.LON', 'V', chs(3))
    call celces(mo//'.TOPOSE.HEA', 'V', chs(4))
    call celces(mo//'.LNNO', 'V', chs(6))
    call celces(mo//'.LTNO', 'V', chs(7))
!
    call jeexin(mo//'.FISSNO    .CELD', ier)
    if (ier .ne. 0) then
        call celces(mo//'.FISSNO', 'V', chs(8))
        call jeveuo(chs(8)//'.CESD', 'L', jcesd(8))
        call jeveuo(chs(8)//'.CESV', 'E', jcesv(8))
        call jeveuo(chs(8)//'.CESL', 'L', jcesl(8))
    end if
!
    call jeexin(mo//'.TOPOSE.PMI.CELD', ier)
    if (ier .ne. 0) then
        call celces(mo//'.TOPOSE.PMI', 'V', chs(9))
        call jeveuo(chs(9)//'.CESD', 'L', jcesd(9))
        call jeveuo(chs(9)//'.CESV', 'E', jcesv(9))
        call jeveuo(chs(9)//'.CESL', 'L', jcesl(9))
    end if
!
!      CALL JEEXIN(MO//'.TOPOFAC.LO.CELD',IER)
!      IF (IER.NE.0) THEN
    call celces(mo//'.TOPOFAC.LO', 'V', chs(10))
    call jeveuo(chs(10)//'.CESD', 'L', jcesd(10))
    call jeveuo(chs(10)//'.CESV', 'E', jcesv(10))
    call jeveuo(chs(10)//'.CESL', 'L', jcesl(10))
!      ENDIF
!
!      CALL JEEXIN(MO//'.TOPOFAC.AI.CELD',IER)
!      IF (IER.NE.0) THEN
    call celces(mo//'.TOPOFAC.AI', 'V', chs(11))
    call jeveuo(chs(11)//'.CESD', 'L', jcesd(11))
    call jeveuo(chs(11)//'.CESV', 'E', jcesv(11))
    call jeveuo(chs(11)//'.CESL', 'L', jcesl(11))
!      ENDIF
!
    call celces(mo//'.TOPONO.HNO', 'V', chs(12))
    call jeveuo(chs(12)//'.CESD', 'L', jcesd(12))
    call jeveuo(chs(12)//'.CESV', 'E', jcesv(12))
    call jeveuo(chs(12)//'.CESL', 'L', jcesl(12))
!
    call jeexin(mo//'.HEAVNO    .CELD', ier)
    if (ier .ne. 0) then
        call celces(mo//'.HEAVNO', 'V', chs(13))
        call jeveuo(chs(13)//'.CESD', 'L', jcesd(13))
        call jeveuo(chs(13)//'.CESV', 'E', jcesv(13))
        call jeveuo(chs(13)//'.CESL', 'L', jcesl(13))
    end if
!
    call jeexin(mo//'.FISSCO    .CELD', ier)
    if (ier .ne. 0) then
        call celces(mo//'.FISSCO', 'V', chs(14))
        call jeveuo(chs(14)//'.CESD', 'L', jcesd(14))
        call jeveuo(chs(14)//'.CESV', 'E', jcesv(14))
        call jeveuo(chs(14)//'.CESL', 'L', jcesl(14))
    end if
!
    call celces(mo//'.TOPOFAC.PI', 'V', chs(15))
    call jeveuo(chs(15)//'.CESD', 'L', jcesd(15))
    call jeveuo(chs(15)//'.CESV', 'E', jcesv(15))
    call jeveuo(chs(15)//'.CESL', 'L', jcesl(15))
!
!    call imprsd('CHAMP',mo//'.BASLOC',8,'VERIF :: BASLOC')
    call celces(mo//'.BASLOC', 'V', chs(16))
    call jeveuo(chs(16)//'.CESD', 'L', jcesd(16))
    call jeveuo(chs(16)//'.CESV', 'L', jcesv(16))
    call jeveuo(chs(16)//'.CESL', 'L', jcesl(16))
!
!    call imprsd('CHAMP',mo//'.STNO',8,'VERIF :: STNO')
    call celces(mo//'.STNO', 'V', chs(17))
    call jeveuo(chs(17)//'.CESD', 'L', jcesd(17))
    call jeveuo(chs(17)//'.CESV', 'L', jcesv(17))
    call jeveuo(chs(17)//'.CESL', 'L', jcesl(17))
!
    do ich = 1, 4
        call jeveuo(chs(ich)//'.CESD', 'L', jcesd(ich))
        call jeveuo(chs(ich)//'.CESV', 'E', jcesv(ich))
        call jeveuo(chs(ich)//'.CESL', 'L', jcesl(ich))
    end do
    do ich = 6, 7
        call jeveuo(chs(ich)//'.CESD', 'L', jcesd(ich))
        call jeveuo(chs(ich)//'.CESV', 'E', jcesv(ich))
        call jeveuo(chs(ich)//'.CESL', 'L', jcesl(ich))
    end do
!
    call jeveuo(malini//'.CONNEX', 'L', jconx1)
    call jeveuo(jexatr(malini//'.CONNEX', 'LONCUM'), 'L', jconx2)
!
    call jeveuo(malini//'.COORDO    .VALE', 'L', iacoo1)
    call jeveuo(maxfem//'.COORDO    .VALE', 'E', iacoo2)
!
    call jeveuo('&CATA.TM.TMDIM', 'L', vi=tmdim)
    call jeveuo(malini//'.TYPMAIL', 'L', vi=typm1)
    call jeveuo(maxfem//'.TYPMAIL', 'E', jtypm2)
    call jeveuo(mo//'.MAILLE', 'L', vi=maille)
!
    if (.not. opmail) then
        call jeveuo(cns1//'.CNSK', 'L', vk8=cnsk)
        call jeveuo(cns1//'.CNSD', 'L', vi=cnsd)
        call jeveuo(cns1//'.CNSV', 'L', jcnsv1)
        call jeveuo(cns2//'.CNSV', 'E', jcnsv2)
        call jeveuo(cns2//'.CNSL', 'E', jcnsl2)
!
!  -----SI ON N'A PAS UNE MODE_MECA
!
        call gettco(resuco, tysd)
        if (tysd(1:9) .ne. 'MODE_MECA' .and. tysd(1:9) .ne. 'EVOL_THER') then
            call jeveuo(ces1//'.CESV', 'L', jcesv1)
            call jeveuo(ces1//'.CESD', 'L', jcesd1)
            call jeveuo(ces1//'.CESL', 'L', jcesl1)
            call jeveuo(ces2//'.CESV', 'E', jcesv2)
            call jeveuo(ces2//'.CESD', 'L', jcesd2)
            call jeveuo(ces2//'.CESL', 'E', jcesl2)
!
            call jeexin(cesvi1//'.CESV', iret)
            if (iret .ne. 0) then
                call jeveuo(cesvi1//'.CESV', 'L', jcviv1)
                call jeveuo(cesvi1//'.CESD', 'L', jcvid1)
                call jeveuo(cesvi1//'.CESL', 'L', jcvil1)
            else
                jcviv1 = 0
                jcvid1 = 0
                jcvil1 = 0
            end if
            iviex = iret
!
            call jeexin(cesvi2//'.CESV', iret)
            if (iret .ne. 0) then
                call jeveuo(cesvi2//'.CESV', 'E', jcviv2)
                call jeveuo(cesvi2//'.CESD', 'L', jcvid2)
                call jeveuo(cesvi2//'.CESL', 'E', jcvil2)
            else
                jcviv2 = 0
                jcvid2 = 0
                jcvil2 = 0
            end if
            iviex = iviex*iret
!
        end if
!
!       COMPORTEMENT
!       RECUPERATION DU CHAM_ELEM_S DU COMPORTEMENT EN ENTREE
        call exisd('CHAM_ELEM_S', comps1, iret)
        if (iret .ne. 0) then
!
!         RECUP DES INFOS SUR LE CHAM_ELEM_S DU COMPORTEMENT EN ENTREE
            call jeveuo(comps1//'.CESD', 'L', jresd1)
            call jeveuo(comps1//'.CESV', 'L', jresv1)
            call jeveuo(comps1//'.CESL', 'L', jresl1)
!
!         NB CMP DU COMPORTEMENT
            nbcmpc = zi(jresd1-1+2)
!
!         VERIF QUE LE CHAMP DE SORTIE A BIEN ETE CREE
            call exisd('CHAM_ELEM_S', comps2, iret)
            ASSERT(iret .ne. 0)
!
!         RECUP DES INFOS SUR LE CHAM_ELEM_S DU COMPORTEMENT EN SORTIE
            call jeveuo(comps2//'.CESD', 'L', jresd2)
            call jeveuo(comps2//'.CESV', 'E', jresv2)
            call jeveuo(comps2//'.CESL', 'E', jresl2)
        else
!         ON MET A ZERO NB CMP
            nbcmpc = 0
        end if
!
        if (tysd(1:4) .eq. 'EVOL') then
            call rsadpa(resuco, 'L', 1, 'INST', iord, &
                        0, sjv=jinst1, styp=kbid)
            inst = zr(jinst1)
        else
            inst = 0.
        end if
        varcns = '&&XPOMAX.VARC'
        call res2mat(resuco, inst, chmat, mu=mu, ka=ka, &
                     nommat=nommat, lvarc=lvarc, varcns=varcns, &
                     cplan=cplan)
        if (lvarc) then
            call jeveuo(varcns//'.CESD', 'L', jcesd_varc)
            call jeveuo(varcns//'.CESV', 'L', jcesv_varc)
            call jeveuo(varcns//'.CESL', 'L', jcesl_varc)
            call rccome(nommat, 'ELAS', iret, k11_ind_nomrc=k11)
            call jeexin(nommat//k11//'.VALR', iret)
            call jelira(nommat//k11//'.VALR', 'LONUTI', nbr)
            call jelira(nommat//k11//'.VALC', 'LONUTI', nbc)
            call jeveuo(nommat//k11//'.VALK', 'L', vk16=valk)
            call jelira(nommat//k11//'.VALK', 'LONUTI', nbk)
            nbf = (nbk-nbr-nbc)/2
            call jelira(chmat//'.CVRCVARC', 'LONMAX', nbvarc)
            ASSERT(nbvarc .le. 2)
            call jeveuo(chmat//'.CVRCVARC', 'L', vk8=cvrcvarc)
        end if
!
    end if
!
!     RECUP DES POINTS DE GAUSS
    call jeveuo('&CATA.TE.PNLOCFPG', 'L', vk32=pnlocfpg)
    call jelira('&CATA.TE.NOLOCFPG', 'LONMAX', nblfpg)
    call jeveuo('&CATA.TM.TMFPG', 'L', vi=tmfpg)
    call jeveuo('&CATA.TE.NOLOCFPG', 'L', vi=nolocfpg)
!
!     RECUP DES NUMEROS DES TYPE DE MAILLES DES SOUS-ELEMENTS
    call jenonu(jexnom('&CATA.TM.NOMTM', typese(1)), itypse(1))
    call jenonu(jexnom('&CATA.TM.NOMTM', typese(2)), itypse(2))
    call jenonu(jexnom('&CATA.TM.NOMTM', typese(3)), itypse(3))
    call jenonu(jexnom('&CATA.TM.NOMTM', typese(4)), itypse(4))
    call jenonu(jexnom('&CATA.TM.NOMTM', typese(5)), itypse(5))
    call jenonu(jexnom('&CATA.TM.NOMTM', typese(6)), itypse(6))
!
!     COMPTEURS DU NOMBRE DE NOUVEAUX NOEUDS ET MAILLES TOTAL
    inntot = 0
    inmtot = 0
!
!     COMPTEUR ET LISTE DE NOEUDS PORTANT DES DDLS DE CONTACT
    linofi = '&&XPOMAX.LINOFI'
!      LINOLA='&&XPOMAX.LINOLA'
    inofi = -1
    if (opmail) then
        call dismoi('NB_NO_MAILLA', maxfem, 'MAILLAGE', repi=nbnofi)
        call wkvect(linofi, 'V V I', nbnofi, inofi)
!       CALL WKVECT(LINOLA,'V V I',NBNOFI,INOLA)
    end if
    nbnofi = 0
!      NBNOLA=0
!
!     RECUPERATION DU VECTEUR DE REMPLISSAGE DES GROUP_MA
    if (opmail) then
        call jeexin(nivgrm, iret)
        if (iret .ne. 0) call jeveuo(nivgrm, 'E', jnivgr)
    end if
!
!     RECUPERATION DU VECTEUR D'INDIRECTION ENTRE LES GROUP_MA
    if (opmail) then
        call jeexin(dirgrm, iret)
        if (iret .ne. 0) call jeveuo(dirgrm, 'L', jdirgr)
    end if
!
!     RECUPERATION DES INFOS DU MAILLAGE SAIN
    if (.not. opmail) then
        noma = cnsk(1)
        call jeveuo(noma//'.TYPMAIL        ', 'L', vi=vtypma)
    end if
!
    call jeveuo(mo//'.XFEM_CONT', 'L', jxc)
!
!     LE RESULTAT EN ENTREE DE POST_CHAM_XFEM EST-IL MECANIQUE
    if (.not. opmail) lmeca = xismec()
!
!     ------------------------------------------------------------------
!                   BOUCLE SUR LES MAILLES DE MAILX
!     ------------------------------------------------------------------
!
    do i = 1, nbmax
!
        ima = zi(jmax-1+i)
!
!       COMPTEUR DES NOUVEAUX NOEUDS ET MAILLES AJOUTES (LOCAL)
        inn = 0
        inm = 0
!
!       RECUPERATION DU NOMBRE DE NOEUDS (2 METHODES)
        call jelira(jexnum(malini//'.CONNEX', ima), 'LONMAX', n)
        nbnoma = zi(jconx2+ima)-zi(jconx2+ima-1)
        ASSERT(n .eq. nbnoma)
!
!       RECUPERATION DU NOMBRE TOTAL DE SOUS ELEMENTS
!       CORRESPOD AU NOMBRE DE MAILLES À CREER
        call cesexi('C', jcesd(3), jcesl(3), ima, 1, &
                    1, 1, iad3)
        ASSERT(iad3 .gt. 0)
        nse = zi(jcesv(3)-1+iad3)
!
!       RECUPERATION DU NOMBRE DE POINTS D'INTERSECTION
!        NPI=ZI(JCESV(3)-1+IAD3+NIT+1)
!
!       RECUPERATION DU NOMBRE DE POINTS MILIEUX
!        NMI=ZI(JCESV(3)-1+IAD3+NIT+3)
!
!       RECUPERATION DU NOMBRE DE NOUVEAUX NOEUDS A CREER
        nnn = zi(jcesv(3)-1+iad3+2)
        if (nnn .eq. 0) goto 100
!        ASSERT(NNN.NE.0)
!
!       NOMBRE DE FISSURES "VUES" DANS LA MAILLE
        nfiss = zi(jcesd(6)-1+5+4*(ima-1)+2)
!
!       CORRESPONDANCE ENTRE FISSURE ET DDL HEAVISIDE
        call jeexin(mo//'.HEAVNO    .CELD', ier)
        if (pre1) then
            do ifiss = 1, nfiss
                do j = 1, n
                    if (ier .ne. 0) then
                        call cesexi('C', jcesd(13), jcesl(13), ima, j, &
                                    ifiss, 1, iad3)
                        if (iad3 .gt. 0) then
                            heavno(j, ifiss) = zi(jcesv(13)-1+iad3)
                        else
                            heavno(j, ifiss) = 1
                        end if
                    else
                        heavno(j, ifiss) = 1
                    end if
                end do
            end do
        else
            do j = 1, 20
                heavno(j, 1:3) = 0
            end do
        end if
!
!       CONNECTIVITE DES JONCTIONS DE FISSURES
        nfijon(1:3) = 0
        fisco(1:2*nfimax) = 0
        fisc(1:2*nfimax) = 0
        call jeexin(mo//'.FISSCO    .CELD', ier)
        if (ier .ne. 0) then
            do ifiss = 1, nfiss
                call cesexi('C', jcesd(14), jcesl(14), ima, 1, &
                            ifiss, 1, iad)
                if (iad .gt. 0) then
                    fisco(2*(ifiss-1)+1) = zi(jcesv(14)-1+iad)
                end if
                call cesexi('C', jcesd(14), jcesl(14), ima, 1, &
                            ifiss, 2, iad)
                if (iad .gt. 0) then
                    fisco(2*ifiss) = zi(jcesv(14)-1+iad)
                end if
            end do
!
            do ifiss = 1, nfiss
!
                ifisc = ifiss
                nfisc = 0
80              continue
                if (fisco(2*ifisc-1) .gt. 0) then
!       STOCKAGE DES FISSURES SUR LESQUELLES IFISS SE BRANCHE
                    nfisc = nfisc+1
                    fisc(2*(nfisc-1)+2) = fisco(2*ifisc)
                    ifisc = fisco(2*ifisc-1)
                    fisc(2*(nfisc-1)+1) = ifisc
                    goto 80
                end if
!
                nfisc2 = 0
                do jfiss = ifiss+1, nfiss
!       STOCKAGE DES FISSURES QUI SE BRANCHENT SUR IFISS
                    kfiss = fisco(2*jfiss-1)
                    do k = nfisc+1, nfisc+nfisc2
                        if (fisc(2*(k-1)+1) .eq. kfiss) then
                            nfisc2 = nfisc2+1
                            fisc(2*(nfisc+nfisc2-1)+1) = jfiss
                            fisc(2*(nfisc+nfisc2)) = fisco(2*jfiss)
                        end if
                    end do
                    if (kfiss .eq. ifiss) then
                        nfisc2 = nfisc2+1
                        fisc(2*(nfisc+nfisc2-1)+1) = jfiss
                        fisc(2*(nfisc+nfisc2)) = fisco(2*jfiss)
                    end if
                end do
!       ON NOTE LE NUMERO DE LA PREMIERE FISSURE BRANCHEE SUR IFISS EN CAS DE
!       JONCTION
                ASSERT(nfisc2 .le. 2)
                if (nfisc2 .gt. 0) then
                    nfijon(1) = ifiss
                    do kfiss = 1, nfisc2
                        nfijon(1+kfiss) = fisc(2*(nfisc+kfiss-1)+1)
                    end do
                end if
                if (nfisc2 .eq. 1) nfijon(3) = nfijon(2)
            end do

        end if
!
!       NOMBRE DE COMPOSANTE DE LA TOPOSE.HEA
        ncomp = zi(jcesd(4)-1+5+4*(ima-1)+3)
!
!       CREATION DU VECTEUR D'INDIRECTION DES NOEUDS (LOCAL A IMA)
        dirno = '&&XPOMAX.DIRNO'
        call wkvect(dirno, 'V V I', nnn*(2+nfiss), jdirno)
!        CALL WKVECT(DIRNO,'V V I',(N+NPI+NMI)*3,JDIRNO)
!
!       DIMENSION TOPOLOGIQUE DE LA MAILLE
        ndime = tmdim(typm1(ima))
!
!       1ER ELEMENT DE REFERENCE ASSOCIE A LA MAILLE
        itypel = maille(ima)
        call jenuno(jexnum('&CATA.TE.NOMTE', itypel), notype)
        call elref2(notype, 10, lirefe, nbelr)
        elrefp = lirefe(1)
!       NOMBRE DE POINT DE GAUSS
        if (.not. iselli(elrefp)) then
            irese = 3
        else
            irese = 0
        end if
        if (ndime .eq. ndim) then
            noflpg = notype//elrese(ndime+irese)//'XINT'
        else
            noflpg = notype//elrese(ndime+irese)//'RIGI'
        end if
        nuflpg = indk32(pnlocfpg, noflpg, 1, nblfpg)
        ASSERT(nuflpg .ne. 0)
        nufgpg = nolocfpg(nuflpg)
        npg = tmfpg(nufgpg)
!
!       CREATION DE VECTEUR DES COORDONNÉES DE LA MAILLE IMA
!       AVEC DES VALEURS CONTIGUES
        geom = '&&XPOAJD.GEOM'
        call wkvect(geom, 'V V R', ndim*n, igeom)
        do in = 1, n
            ino = zi(jconx1-1+zi(jconx2+ima-1)+in-1)
            do j = 1, ndim
                zr(igeom-1+ndim*(in-1)+j) = zr(iacoo1-1+3*(ino-1)+j)
            end do
        end do
!
        if (.not. opmail) then
!         TYPE DE LA MAILLE PARENT
            call jenuno(jexnum('&CATA.TM.NOMTM', vtypma(ima)), typma)
!         CONNECTIVITES DU MAILLAGE QUADRATIQUE
!         POUR RECUPERER LES LAGRANGES
            call jeveuo(noma//'.CONNEX', 'L', jconq1)
            call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconq2)
        end if
!
        nfe = 0
        nfh = 0
!
!       COMPOSANTES DU CHAMP DE DEPLACEMENT 1 POUR LA MAILLE IMA
!       ET RECUPERATION DES CONTRAINTES 1
        if (.not. opmail) then
            nbcmp = cnsd(2)
            ASSERT(nbcmp .le. ddlmax)
            call xpocmp(elrefp, cns1, ima, n, jconx1, &
                        jconx2, ndim, nfh, nfe, ddlc, &
                        nbcmp, cmp, lmeca, pre1)
!
            if (tysd(1:9) .ne. 'MODE_MECA' .and. tysd(1:9) .ne. 'EVOL_THER') then
                call cesexi('C', jcesd1, jcesl1, ima, 1, &
                            1, 1, iadc)
                if (iviex .ne. 0) call cesexi('C', jcvid1, jcvil1, ima, 1, &
                                              1, 1, iadv)
            end if
        end if
!
!       RECUPERATION DES COORDONNEES DES POINTS D'INTERSECTION
        call cesexi('C', jcesd(1), jcesl(1), ima, 1, &
                    1, 1, iad1)
        ASSERT(iad1 .gt. 0)
!
!       RECUPERATION DES COORDONNEES DES POINTS MILIEUX
        if (.not. iselli(elrefp)) then
            call cesexi('C', jcesd(9), jcesl(9), ima, 1, &
                        1, 1, iad9)
            ASSERT(iad9 .gt. 0)
        else
            iad9 = 0
        end if
!
!       RECUPERATION DU NOMBRE DE POINTS D'INTERSECTION
!       PAR FISSURE SUR LES MAILLES PORTEUSES DE DDL DE CONTACT
!
        do ifiss = 1, nfiss
            call cesexi('C', jcesd(10), jcesl(10), ima, 1, &
                        ifiss, 1, iad10)
            ninter(ifiss) = zi(jcesv(10)-1+iad10)
        end do
        call cesexi('C', jcesd(11), jcesl(11), ima, 1, &
                    1, 1, iad11)
        ncompa = zi(jcesd(11)-1+5+4*(ima-1)+3)
!
!       RECUPERATION DE LA CONNECTIVITE DES SOUS-ELEMENTS
        call cesexi('C', jcesd(2), jcesl(2), ima, 1, &
                    1, 1, iad2)
        ASSERT(iad2 .gt. 0)
!
!       RECUPERATION DE LA FONCTION HEAVISIDE
        call cesexi('C', jcesd(4), jcesl(4), ima, 1, &
                    1, 1, iad4)
        ASSERT(iad4 .gt. 0)
!
!       RECUPERATION DES INFOS CONCERNANT LES GROUP_MA CONTENANT IMA
        if (opmail) then
            call jelira(jexnum(listgr, ima), 'LONMAX', ngrm)
            if (ngrm .gt. 0) call jeveuo(jexnum(listgr, ima), 'L', iagma)
        end if
!
!       RECUPERATION DES LEVELS SET AUX NOEUDS
        lst = '&&XPOAJD.LST'
        lsn = '&&XPOAJD.LSN'
        hea = '&&XPOAJD.HEA'
        call wkvect(lst, 'V V R', n*nfiss, jlst)
        call wkvect(lsn, 'V V R', n*nfiss, jlsn)
        call wkvect(hea, 'V V I', nfiss, jhea)
        do ifiss = 1, nfiss
            do j = 1, n
                call cesexi('C', jcesd(6), jcesl(6), ima, j, &
                            ifiss, 1, iad)
                ASSERT(iad .gt. 0)
                zr(jlsn-1+(j-1)*nfiss+ifiss) = zr(jcesv(6)-1+iad)
                call cesexi('C', jcesd(7), jcesl(7), ima, j, &
                            ifiss, 1, iad)
                ASSERT(iad .gt. 0)
                zr(jlst-1+(j-1)*nfiss+ifiss) = zr(jcesv(7)-1+iad)
            end do
        end do
!
!       POUR LES JONCTIONS HM-XFEM, ON CONSTRUIT LACT
        nlachm(1:2) = 0
        lacthm(1:16) = 0
        vit(1:16) = 0
        if (pre1 .and. nfijon(1) .ne. 0 .and. ndim .eq. ndime .and. .not. opmail) then
            call cesexi('C', jcesd(15), jcesl(15), ima, 1, &
                        1, 1, iad16)
            ncompi = zi(jcesd(15)-1+5+4*(ima-1)+3)
            call conare(typma, ar, nbar)
            do nli = 1, ninter(nfijon(1))
                pinter(:) = 0.d0
                do j = 1, ndim
                    pinter(j) = zr(jcesv(15)-1-1+iad16+ncompi*(nfijon(1)-1)+ndim*(nli-1)+j)
                end do
                ff(:) = 0.d0
                call elrfvf(elrefp, pinter, ff, n)
                lsno = 0.d0
                do j = 1, n
                    lsno = lsno+ff(j)*zr(jlsn-1+(j-1)*nfiss+nfijon(2))
                end do
                dec = 0
                if (lsno .gt. -1.d-8) dec = 8
                iar = int(zr(jcesv(11)-1-1+iad11+ncompa*(nfijon(1)-1)+5*(nli-1)+1))
                ino = int(zr(jcesv(11)-1-1+iad11+ncompa*(nfijon(1)-1)+5*(nli-1)+2))
                nvit = int(zr(jcesv(11)-1-1+iad11+ncompa*(nfijon(1)-1)+5*(nli-1)+5))
                if (ino .gt. 0) then
                    lacthm(ino+dec) = nli
                else if (iar .gt. 0) then
                    ino1 = ar(iar, 1)
                    ino2 = ar(iar, 2)
                    if (nvit .eq. 1) then
                        lacthm(ino1+dec) = nli
                        vit(ino1+dec) = 1
                        lacthm(ino2+dec) = nli
                        vit(ino2+dec) = 1
                    else
                        if (vit(ino1+dec) .eq. 0) lacthm(ino1+dec) = nli
                        if (vit(ino2+dec) .eq. 0) lacthm(ino2+dec) = nli
                    end if
                end if
            end do
            do ino = 1, 8
                if (lacthm(ino) .ne. 0) nlachm(1) = nlachm(1)+1
                if (lacthm(ino+8) .ne. 0) nlachm(2) = nlachm(2)+1
            end do
        end if
!
!       RECUPERATION DE LA CONNECTIVITÉ DES FISSURES
        if (.not. opmail .and. nfh .gt. 0) then
!         CORRECTION DE NFH SI ON SE TROMPE DANS XPOCMP
            if (nfh .gt. nfiss) nfh = nfiss
            if (nfiss .gt. 1) nfh = zi(jcesd(8)-1+5+4*(ima-1)+2)
        end if
!
!       RECUPERATION DE LA DEFINITION DES FONCTION HEAVISIDES
        if (.not. opmail .and. nfh .gt. 0) then
!       NOMBRE DE COMPOSANTE DE LA SD TOPONO.HNO
            ncompn_tmp = zi(jcesd(12)-1+5+4*(ima-1)+3)
            heavn = '&&XPOAJD.HEAVN'
!       EN PRINCIPE, TOPONO N EST DEFINI QUE DANS LES MAILLES XH
!          MAIS IL PEUT ARRIVER QUE LA MAILLE VOIT AU MOINS
!          1 NOEUD  AVEC UNE COMPOSANTE HEAVISIDE.
!          LE DEPLACEMEMENT ASSOCIE DOIT VALOIR ZERO SUR LA MAILLE COURANTE
!          DANS CE CAS, ON NE PREND PAS DE RISQUE, ON CALCULE AUSSI UN COEF
!          HEAVISIDE EGAL ZERO PAR LE HEAVN.
            if (ncompn_tmp .eq. 0) then
                ncompn = 5
                call wkvect(heavn, 'V V I', n*ncompn, jheavn)
                do ifh = 1, ncompn
                    do j = 1, n
                        zi(jheavn-1+(j-1)*ncompn+ifh) = -99
                    end do
                end do
            else
                ncompn = ncompn_tmp
                call wkvect(heavn, 'V V I', n*ncompn, jheavn)
                do ifh = 1, ncompn
                    do j = 1, n
                        call cesexi('C', jcesd(12), jcesl(12), ima, j, &
                                    1, ifh, iad)
                        ASSERT(iad .gt. 0)
                        zi(jheavn-1+(j-1)*ncompn+ifh) = zi(jcesv(12)-1+iad)
                    end do
                end do
            end if
            ASSERT(ncompn .eq. 5)
        end if
!
!       RECUPERATION DE LA BASE LOCALE EN FOND DE FISSURE ET DU STATUT DES NOEUDS
        if (.not. opmail .and. (nfh+nfe) .gt. 0) then
            ncmp = zi(jcesd(17)-1+5+4*(ima-1)+3)
            stano = '&&XPOAJD.STANO'
            call wkvect(stano, 'V V I', n, jstno)
            if (ncmp .gt. 0) then
                do j = 1, n
                    call cesexi('C', jcesd(17), jcesl(17), ima, j, &
                                1, 1, iad)
                    ASSERT(iad .gt. 0)
                    zi(jstno-1+j) = zi(jcesv(17)-1+iad)
                end do
            else
                zi(jstno:(jstno-1+n)) = -99
            end if
        end if
        if (.not. opmail .and. nfe .gt. 0) then
            ncmp = zi(jcesd(16)-1+5+4*(ima-1)+3)
            basloc = '&&XPOAJD.BASLOC'
            call wkvect(basloc, 'V V R', 3*ndim*n, jbaslo)
            if (ncmp .gt. 0) then
                do j = 1, n
                    do l = 1, 3*ndim
                        call cesexi('C', jcesd(16), jcesl(16), ima, j, &
                                    1, l, iad)
                        ASSERT(iad .gt. 0)
                        zr(jbaslo-1+(j-1)*3*ndim+l) = zr(jcesv(16)-1+iad)
                    end do
                end do
            else
                zr(jbaslo:(jbaslo-1+n*3*ndim)) = 0.d0
            end if
        end if
!       CALCUL DES PARAMETRES MATERIAUX
!         INTERPOLATION DES VARC
        if (.not. opmail .and. lvarc) then
            ncmp = zi(jcesd_varc-1+5+4*(ima-1)+3)
            if (ncmp .eq. 0) then
                mu = 1.
                ka = 3.
                goto 15
            end if
            varc(:) = 0.
            do j = 1, n
                do k = 1, nbvarc
                    call cesexi('C', jcesd_varc, jcesl_varc, ima, j, &
                                1, k, iad)
                    ASSERT(iad .gt. 0)
                    varc(k) = varc(k)+zr(jcesv_varc-1+iad)/n
                end do
            end do
            poiss = .false.
            young = .false.
            do ik = 1, nbf
                if (valk(nbr+nbc+ik) .eq. 'NU') then
                call fointe('C', valk(nbr+nbc+nbf+ik), nbvarc, cvrcvarc(1:nbvarc), varc(1:nbvarc), &
                                nu, ier)
                    if (ier .eq. 0) poiss = .true.
                end if
                if (valk(nbr+nbc+ik) .eq. 'E') then
                call fointe('C', valk(nbr+nbc+nbf+ik), nbvarc, cvrcvarc(1:nbvarc), varc(1:nbvarc), &
                                e, ier)
                    if (ier .eq. 0) young = .true.
                end if
            end do
            if (poiss) then
                ka = 3.d0-4.d0*nu
                if (cplan) ka = (3.d0-nu)/(1.d0+nu)
                if (young) mu = e/(2.d0*(1.d0+nu))
            end if
15          continue
        end if
!
! ----- ON AJOUTE LES NOUVELLES MAILLES ET LES NOUVEAUX NOEUDS
!
        nnose = tabse(ndime+irese)
!
!         BOUCLE D'INTEGRATION SUR LES NSE SOUS-ELEMENTS
        do ise = 1, nse
            do ifiss = 1, nfiss
                zi(jhea-1+ifiss) = zi(jcesv(4)-1+iad4-1+ncomp*(ifiss-1)+ise)
            end do
            jcnse = jcesv(2)-1+iad2
            call xpoajm(maxfem, jtypm2, itypse(ndime+irese), jcnse, ise, &
                        n, nnose, prefno, jdirno, nse, &
                        inm, inmtot, nbmac, zi(jhea), jnivgr, &
                        iagma, ngrm, jdirgr, opmail, nfiss, &
                        ndim, ndime, jconx1, jconx2, jconq1, &
                        jconq2, ima, iad1+jcesv(1)-1, nnn, inn, &
                        inntot, nbnoc, nbnofi, inofi, iacoo1, &
                        iacoo2, iad9+jcesv(9)-1, ninter, jcesv(11)+iad11-1, &
                        ncompa, elrefp, jlsn, jlst, typma, igeom, jheavn, ncompn, &
                        zi(jxc), cmp, nbcmp, nfh, nfe, &
                        ddlc, jcnsv1, jcnsv2, jcnsl2, lmeca, &
                        pre1, heavno, fisco, nlachm, lacthm, &
                        jbaslo, jstno, ka, mu)
            if (.not. opmail) then
                if (tysd(1:9) .ne. 'MODE_MECA' .and. tysd(1:9) .ne. 'EVOL_THER') then
!
!             ON AJOUTE DES CONTRAINTES
                    call xpoajc(nse, inm, inmtot, nbmac, ise, &
                                npg, jcesd1, jcesd2, jcvid1, jcvid2, &
                                ima, ndim, ndime, iadc, iadv, &
                                jcesv1, jcesl2, jcesv2, jcviv1, jcvil2, &
                                jcviv2)
!
                    call xpocox(nbmac, ima, inmtot, nbcmpc, jresd1, &
                                jresv1, jresl1, jresd2, jresv2, jresl2)
!
!
                end if
            end if
        end do
!
        if (opmail) then
            ASSERT(inn .eq. nnn)
        end if
!
        call jedetr(geom)
        call jedetr(dirno)
        call jedetr(lsn)
        call jedetr(lst)
        call jedetr(hea)
        if (.not. opmail .and. nfh .gt. 0) call jedetr(heavn)
        if (.not. opmail .and. nfe .gt. 0) call jedetr(basloc)
        if (.not. opmail .and. (nfh+nfe) .gt. 0) call jedetr(stano)
!
100     continue
    end do
!
!     CREATION DU GROUPE DES NOEUDS SITUES SUR LA FISSURE
!     PORTANT DES DDLS DE CONTACT
    grpnoe = maxfem//'.GROUPENO'
    if (opmail .and. nbnofi .gt. 0) then
        nogno = nogrfi
!       ON SAIT QUE LE .GROUPENO N'EXISTE PAS
        gpptnn = maxfem//'.PTRNOMNOE'
        call jecreo(gpptnn, 'G N K24')
        call jeecra(gpptnn, 'NOMMAX', 1+ngfon)
        call jecrec(grpnoe, 'G V I', 'NO '//gpptnn, 'DISPERSE', 'VARIABLE', &
                    1+ngfon)
        call jecroc(jexnom(grpnoe, nogno))
        call jeecra(jexnom(grpnoe, nogno), 'LONMAX', max(1, nbnofi))
        call jeecra(jexnom(grpnoe, nogno), 'LONUTI', nbnofi)
        call jeveuo(jexnom(grpnoe, nogno), 'E', iagno)
        do j = 1, nbnofi
            zi(iagno-1+j) = zi(inofi-1+j)
        end do
    end if
!      IF (OPMAIL.AND.NBNOLA.GT.0) THEN
!        NOGNO=NOGRLA
!        CALL JECROC(JEXNOM(GRPNOE,NOGNO))
!        CALL JEECRA(JEXNOM(GRPNOE,NOGNO),'LONMAX',NBNOLA,K8B)
!        CALL JEVEUO(JEXNOM(GRPNOE,NOGNO),'E',IAGNOL)
!        DO 211 J = 1,NBNOLA
!          ZI(IAGNOL-1+J) = ZI(INOLA-1+J)
! 211    CONTINUE
!      ENDIF
    do ich = 1, 4
        call detrsd('CHAM_ELEM_S', chs(ich))
    end do
    do ich = 6, 17
        call detrsd('CHAM_ELEM_S', chs(ich))
    end do
!
    if (.not. opmail) call detrsd('CHAM_NO_S', varcns)
!
    if (opmail) call jedetr(mailx)
    if (opmail) call jedetr(linofi)
!
999 continue
!
!      IF (.NOT.OPMAIL) CALL IMPRSD('CHAMP',COMPS2,6,'COMPS2')
!
    call jedema()
end subroutine
