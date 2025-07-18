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

subroutine xpodim(malini, mailc, modvis, licham, nsetot, &
                  nnntot, ncotot, listno, cns1, cns2, &
                  ces1, ces2, cel2, cesvi1, cesvi2, &
                  ior, resuco, nbnoc, nbmac, logrma, &
                  dirgrm, maxfem, ngfon, comps1, comps2, &
                  pre1, mo)
!
! person_in_charge: samuel.geniaut at edf.fr
!
! aslint: disable=W1504
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/cescre.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnscre.h"
#include "asterfort/dismoi.h"
#include "asterfort/gmgnre.h"
#include "asterfort/infniv.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedup1.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "asterfort/wkvect.h"
#include "asterfort/xismec.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: nsetot, nnntot, ncotot, nbnoc, ior, ngfon
    character(len=8) :: maxfem, malini, resuco, modvis, mo
    character(len=19) :: cns1, cns2, ces1, ces2, cel2, cesvi1, cesvi2
    character(len=19) :: comps1, comps2
    character(len=24) :: mailc, listno, logrma, dirgrm, licham
!
!
!   DIMENTIONNEMENT DES OBJETS DU NOUVEAU MAILLAGE MAXFEM
!                ET DES NOUVEAUX CHAMPS RESULTAT
!
!   IN
!       MALINI : MAILLAGE SAIN
!       MAILC  : LISTE DES NUMEROS DES MAILLES NON SOUS-DECOUPEES
!       MODVIS : MODELE DE VISU (X-FEM)
!       LICHAM : LISTE DES CHAMPS A POST-TRAITER
!       NSETOT : NOMBRE TOTAL DE SOUS-ELEMENT
!       NNNTOT : NOMBRE TOTAL DE NOUVEAUX NOEUDS
!       NCOTOT : LONGUEUR DE LA CONNECTIVITE DES NOUVEAUX NOEUDS
!       CNS1   : CHAMP_NO_S DU DEPLACEMENT EN ENTREE
!       CES1   : CHAMP_ELEM_S DE CONTRAINTES EN ENTREE
!       MAXFEM : MAILLAGE FISSURE (SI POST_CHAMP_XFEM)
!       IOR    : POSITION DU NUMERO D'ORDRE (SI POST_CHAMP_XFEM)
!       RESUCO : NOM DU CONCEPT RESULTAT DONT ON EXTRAIT LES CHAMPS
!       LOGRMA : LONGUEUR DES NOUVEAUX GROUP_MA
!       DIRGRM : VECTEUR D'INDIRECTION ENTRE LES GROUP_MA
!       NGFON  : NOMBRE TOTAL DE FOND DE FISSURES
!       MO     : NOM DU MODELE
!
!   OUT
!       MAXFEM : MAILLAGE FISSURE (SI POST_MAIL_XFEM)
!       CNS2   : CHAMP_NO_S DU DEPLACEMENT EN SORTIE
!       CES2   : CHAMP_ELEM_S DE CONTRAINTES EN SORTIE
!       CEL2   : CHAMP_ELEM DE CONTRAINTES EN SORTIE
!       NBNOC  : NOMBRE DE NOEUDS CLASSIQUES DU MAILLAGE FISSURE
!       NBMAC  : NOMBRE DE MAILLES CLASSIQUES DU MAILLAGE FISSURE
!       LISTNO : LISTE DES NUMEROS DES NOEUDS CLASSIQUES
!
!
!
    integer(kind=8) :: ier, nbmac, nbma2, nbno, nbno2, iret, igeomr, nbid
    integer(kind=8) :: iadesc, iacoo2, jtypm2, jno, jmac
    integer(kind=8) :: ndim, jord, iord, i, ifm, niv, nmaxsp, nmaxcm, nbcham
    integer(kind=8) :: jdirgr, jxc, contac
    integer(kind=8) :: igma1, nbgma, n, jlogma, nbgma1, nbgma2, cptgr2, jlicha
    integer(kind=8) :: jresc1, nbcmp
    aster_logical :: pre1
    character(len=3) :: tsca
    character(len=8) :: k8b, ldep3(12), ldep2(8), ltemp(1), ldep1(12), ldep4(16)
    character(len=16) :: k16b, nomcmd
    character(len=19) :: coord2, ligrel, chn1, chsig1
    character(len=19) :: nomgd
    character(len=24) :: ordr, gpptnm, nogma
    character(len=24) :: comp1
    integer(kind=8), pointer :: litrav(:) => null()
    integer(kind=8), pointer :: nbpt(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    integer(kind=8), pointer :: cesd2(:) => null()
    integer(kind=8), pointer :: cvid1(:) => null()
    integer(kind=8), pointer :: resd1(:) => null()
    character(len=8), pointer :: cnsc(:) => null()
    character(len=8), pointer :: cnsk1(:) => null()
    character(len=8), pointer :: cnsk2(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
    data ldep3/'DX', 'DY', 'DZ', 'LAGS_C', 'LAGS_F1', 'LAGS_F2', &
        'LAG2_C', 'LAG2_F1', 'LAG2_F2', 'LAG3_C', 'LAG3_F1', 'LAG3_F2'/
    data ldep2/'DX', 'DY', 'LAGS_C', 'LAGS_F1', 'LAG2_C', 'LAG2_F1', &
        'LAG3_C', 'LAG3_F1'/
    data ldep1/'DX', 'DY', 'PRE1', 'PRE_FLU', 'LAG_FLI', 'LAG_FLS', &
        'LAGS_C', 'LAGS_F1', 'LAG2_C', 'LAG2_F1', 'LAG3_C', 'LAG3_F1'/
    data ldep4/'DX', 'DY', 'DZ', 'PRE1', 'PRE_FLU', 'LAG_FLI', 'LAG_FLS', &
        'LAGS_C', 'LAGS_F1', 'LAGS_F2', 'LAG2_C', 'LAG2_F1', &
        'LAG2_F2', 'LAG3_C', 'LAG3_F1', 'LAG3_F2'/
    data ltemp/'TEMP'/
!
    call jemarq()
    call infniv(ifm, niv)
!
!     NOM DE LA COMMANDE (POST_MAIL_XFEM OU POST_CHAM_XFEM)
    call getres(k8b, k16b, nomcmd)
!
!     NOMBRE DE MAILLES CLASSIQUES
    call jeexin(mailc, ier)
    if (ier .eq. 0) then
        nbmac = 0
        nbnoc = 0
    else
        call jeveuo(mailc, 'L', jmac)
        call jelira(mailc, 'LONMAX', nbmac)
!       RECHERCHE DE LA LISTE DE NOEUDS SOUS-JACENTE
        call dismoi('NB_NO_MAILLA', malini, 'MAILLAGE', repi=nbno)
        AS_ALLOCATE(vi=litrav, size=nbno)
        call wkvect(listno, 'V V I', nbno, jno)
        call gmgnre(malini, nbno, litrav, zi(jmac), nbmac, &
                    zi(jno), nbnoc, 'TOUS')
        AS_DEALLOCATE(vi=litrav)
    end if
!
!     NOMBRE DE MAILLES DU NOUVEAU MAILLAGE : NBMA2
    nbma2 = nbmac+nsetot
!
!     NOMBRE DE NOEUDS DU NOUVEAU MAILLAGE : NBNO2
    nbno2 = nbnoc+nnntot
    if (nomcmd .eq. 'POST_MAIL_XFEM') then
!
!       ---------------------------------------------------------------
!                   TRAITEMENT POUR POST_MAIL_XFEM
!       ---------------------------------------------------------------
!
!       CHANGEMENT DU .DIME
!       .DIME(1) : NOMBRE DE NOEUDS PHYSIQUES DU MAILLAGE
!       .DIME(3) : NOMBRE DE MAILLES DU MAILLAGE
        call jedup1(malini//'.DIME', 'G', maxfem//'.DIME')
        call jeveuo(maxfem//'.DIME', 'E', vi=dime)
        dime(1) = nbno2
        dime(3) = nbma2
        ASSERT(dime(2) .eq. 0)
        ASSERT(dime(4) .eq. 0)
        ASSERT(dime(5) .eq. 0)
!
!       CREATION DES .NOMMAI, .NOMNOE, .TYPMAIL, .COORDO  .CONNEX
!       DU MAILLAGE 2
        call jecreo(maxfem//'.NOMMAI', 'G N K8')
        call jeecra(maxfem//'.NOMMAI', 'NOMMAX', nbma2)
!
        call jecreo(maxfem//'.NOMNOE', 'G N K8')
        call jeecra(maxfem//'.NOMNOE', 'NOMMAX', nbno2)
!
        coord2 = maxfem//'.COORDO'
        call jenonu(jexnom('&CATA.GD.NOMGD', 'GEOM_R'), igeomr)
        call wkvect(coord2//'.DESC', 'G V I', 3, iadesc)
        call jeecra(coord2//'.DESC', 'DOCU', cval='CHGO')
        zi(iadesc-1+1) = igeomr
!       -- TOUJOURS 3 COMPOSANTES X, Y ET Z
        zi(iadesc-1+2) = -3
!       -- 14 = 2**1 + 2**2 + 2**3
        zi(iadesc-1+3) = 14
        call wkvect(coord2//'.VALE', 'G V R', 3*nbno2, iacoo2)
!
        call wkvect(maxfem//'.TYPMAIL', 'G V I', nbma2, jtypm2)
!
        call jecrec(maxfem//'.CONNEX', 'G V I', 'NU', 'CONTIG', 'VARIABLE', &
                    nbma2)
        call jeecra(maxfem//'.CONNEX', 'LONT', ncotot)
        call jeecra(maxfem//'.CONNEX', 'NUTIOC', nbma2)

!
!       .GROUPEMA
!       NOMBRE DE GROUP_MA A CREER DANS MAXFEM : NBGMA2
!       (PREND EN COMPTE LES FONDS DE FISSURE ET LES GROUPES DE MAILLES
!       DU MAILLAGE INITIAL)
        nbgma2 = ngfon
        call jeexin(malini//'.GROUPEMA', igma1)
        if (igma1 .ne. 0) then
            call jeveuo(logrma, 'L', jlogma)
            call jelira(logrma, 'LONMAX', nbgma)
            call jelira(malini//'.GROUPEMA', 'NUTIOC', nbgma1)
            ASSERT(nbgma .eq. nbgma1)
            do i = 1, nbgma
                if (zi(jlogma-1+i) .ne. 0) nbgma2 = nbgma2+1
            end do
        end if
        if (nbgma2 .ne. 0) then
            gpptnm = maxfem//'.PTRNOMMAI'
            call jecreo(gpptnm, 'G N K24')
            call jeecra(gpptnm, 'NOMMAX', nbgma2)
            call jecrec(maxfem//'.GROUPEMA', 'G V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE', &
                        nbgma2)
!         CREATION DU VECTEUR D'INDIRECTION DES GROUP_MA
            if (igma1 .ne. 0) then
                call wkvect(dirgrm, 'V V I', nbgma, jdirgr)
                cptgr2 = 0
                do i = 1, nbgma
                    if (zi(jlogma-1+i) .ne. 0) then
                        call jenuno(jexnum(malini//'.GROUPEMA', i), nogma)
                        call jecroc(jexnom(maxfem//'.GROUPEMA', nogma))
                        n = zi(jlogma-1+i)
                        cptgr2 = cptgr2+1
                        zi(jdirgr-1+i) = cptgr2
                        call jeecra(jexnum(maxfem//'.GROUPEMA', cptgr2), 'LONMAX', n)
                        call jeecra(jexnum(maxfem//'.GROUPEMA', cptgr2), 'LONUTI', n)
                        write (ifm, 808) nogma, n
                    end if
                end do
            end if
        end if
        write (ifm, *) ' '
!
!
        if (niv .gt. 1) then
            write (ifm, *) 'CREATION .NOMAI DE LONGUEUR ', nbma2
            write (ifm, *) 'CREATION .TYPMAIL DE LONGUEUR ', nbma2
            write (ifm, *) 'CREATION .NOMNOE DE LONGUEUR ', nbno2
            write (ifm, *) 'CREATION .COORDO.VALE DE LONGUEUR ', 3*nbno2
            write (ifm, *) 'CREATION .CONNEX DE LONGUEUR ', ncotot
            write (ifm, *) 'CREATION .GROUPEMA ', nbgma2
        end if
!
    else if (nomcmd .eq. 'POST_CHAM_XFEM') then
!
!       ---------------------------------------------------------------
!                   TRAITEMENT POUR POST_CHAM_XFEM
!       ---------------------------------------------------------------
!
!       CHAMP A POST-TRAITER
        call jelira(licham, 'LONMAX', nbcham)
        call jeveuo(licham, 'L', jlicha)
!
        call dismoi('DIM_GEOM', malini, 'MAILLAGE', repi=ndim)
        call dismoi('NB_NO_MAILLA', maxfem, 'MAILLAGE', repi=nbid)
        ASSERT(nbno2 .eq. nbid)
!
        ASSERT(nbcham .gt. 0)
!
        chn1 = '&&XPODIM.CHN1'
!
!       EXTRACTION DES DEPLACEMENTS : CNS1
        ordr = resuco//'           .ORDR'
        call jeveuo(ordr, 'L', jord)
        iord = zi(jord-1+ior)
        call rsexch('F', resuco, zk16(jlicha-1+1), iord, chn1, &
                    ier)
        call cnocns(chn1, 'V', cns1)
        call jeveuo(cns1//'.CNSK', 'L', vk8=cnsk1)
!
!       TYPE DE CONTACT
        call jeveuo(mo//'.XFEM_CONT', 'L', jxc)
        contac = zi(jxc)
!
!       CREATION D'UN CHAMP SIMPLE : CNS2
        if (xismec()) then
            if (pre1) then
!         CAS DE L'HYDRO-MECANIQUE
                if (ndim .eq. 2) then
                    call cnscre(maxfem, 'DEPL_R', 4*ndim+4, ldep1, 'V', &
                                cns2)
                else
                    ASSERT(ndim .eq. 3)
                    call cnscre(maxfem, 'DEPL_R', 4*ndim+4, ldep4, 'V', &
                                cns2)
                end if
            else
!         CAS DE LA MECANIQUE
                if (ndim .eq. 2) then
                    call cnscre(maxfem, 'DEPL_R', 4*ndim, ldep2, 'V', &
                                cns2)
                elseif (ndim .eq. 3) then
                    ASSERT(ndim .eq. 3)
                    call cnscre(maxfem, 'DEPL_R', 4*ndim, ldep3, 'V', &
                                cns2)
                end if
            end if
        else
!         CAS DE LA THERMIQUE
            call cnscre(maxfem, 'TEMP_R', 1, ltemp, 'V', &
                        cns2)
        end if
!
        call jeveuo(cns2//'.CNSK', 'E', vk8=cnsk2)
        call jeveuo(cns2//'.CNSD', 'E', vi=cnsd)
        call jeveuo(cns2//'.CNSC', 'E', vk8=cnsc)
!
!       REMPLISSAGE DES DIFFERENTS OBJETS DU CNS2
        cnsk2(1) = maxfem
        cnsk2(2) = cnsk1(2)
!
        cnsd(1) = nbno2
        if (xismec()) then
            if (pre1) then
!           CAS DE L'HYDRO-MECANIQUE
                cnsd(2) = 4*ndim+4
                do i = 1, 4*ndim+4
                    if (ndim .eq. 2) cnsc(i) = ldep1(i)
                    if (ndim .eq. 3) cnsc(i) = ldep4(i)
                end do
            else
!           CAS DE LA MECANIQUE
                cnsd(2) = 4*ndim
                do i = 1, 4*ndim
                    if (ndim .eq. 2) cnsc(i) = ldep2(i)
                    if (ndim .eq. 3) cnsc(i) = ldep3(i)
                end do
            end if
        else
!         CAS DE LA THERMIQUE
            cnsd(2) = 1
            cnsc(1) = ltemp(1)
        end if
!
!       CONTRAINTES
        if (nbcham .gt. 1) then
            chsig1 = '&&XPODIM.CHE1'
            call rsexch('F', resuco, zk16(jlicha-1+2), iord, chsig1, &
                        ier)
            call celces(chsig1, 'V', ces1)
!
!         CREATION D'UN CHAMP SIMPLE
            ligrel = modvis//'.MODELE'
            call alchml(ligrel, 'FULL_MECA', 'PCONTMR', 'V', cel2, &
                        iret, ' ')
            call celces(cel2, 'V', ces2)
!
            if (nbcham .gt. 2) then
!           VARIABLES INTERNES
                chsig1 = '&&XPODIM.CHE1'
                call rsexch(' ', resuco, zk16(jlicha-1+3), iord, chsig1, &
                            ier)
            else
                ier = 1
            end if
!
!         SI LE CHAMP EXISTE
            if (ier .eq. 0) then
                call celces(chsig1, 'V', cesvi1)
!
!           TROUVER LES VALEURS MAJORANTES DU NOMBRE DE POINTS DE GAUSS,
!           NOMBRE DES SOUS-POINTS ET NOMBRE DE COMPOSANTES
!
                call jeveuo(cesvi1//'.CESD', 'L', vi=cvid1)
                call jeveuo(ces2//'.CESD', 'L', vi=cesd2)
!
                nmaxsp = cvid1(4)
                nmaxcm = cvid1(5)
                AS_ALLOCATE(vi=nbpt, size=nbma2)
!
                do i = 1, nbma2
                    nbpt(i) = cesd2(5+4*(i-1)+1)
!
                end do
!
                call cescre('V', cesvi2, 'ELGA', maxfem, 'VARI_R', &
                            -nmaxcm, ' ', nbpt, [-nmaxsp], [-nmaxcm])
            end if
!
        end if
!
!       COMPORTEMENT
!       RECUPERATION DE LA CARTE DE COMPORTEMENT UTILISEE DANS LE CALCUL
        call rsexch(' ', resuco, 'COMPORTEMENT', iord, comp1, &
                    iret)
!
!       SI LA CARTE DE COMPORTEMENT EXISTE
        if (iret .eq. 0) then
!
!         NOM ET TYPE DE LA GRANDEUR (NORMALEMENT TYPE='K16')
            call dismoi('NOM_GD', comp1, 'CARTE', repk=nomgd)
            call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
            ASSERT(tsca .eq. 'K16')
!
!         CREATION CHAM_ELEM_S 1 A PARTIR DE LA CARTE 1
            call carces(comp1, 'ELEM', ' ', 'V', comps1, &
                        'A', iret)
            call jeveuo(comps1//'.CESD', 'L', vi=resd1)
            call jeveuo(comps1//'.CESC', 'L', jresc1)
!
!         NB CMP
            nbcmp = resd1(2)
!
!         CREATION CHAM_ELEM_S 2 : COMPS2
            call cescre('V', comps2, 'ELEM', maxfem, nomgd, &
                        nbcmp, zk8(jresc1), [0], [-1], [-nbcmp])
!
        end if
!
    end if
!
    AS_DEALLOCATE(vi=nbpt)
    call jedema()
!
808 format(30x, a8, 2x, i12)
!
!
end subroutine
