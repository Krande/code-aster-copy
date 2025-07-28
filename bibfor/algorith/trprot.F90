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

subroutine trprot(model, bamo, tgeom, imodg, iadx, &
                  iady, iadz, isst, iadrp, norm1, &
                  norm2, ndble, num, nu, ma, &
                  mateco, moint, ilires, k, icor)
    implicit none
!     AUTEUR : G.ROUSSEAU
!     ROUTINE REALISANT ,
!     A PARTIR D UN MODELE GENERALISE , ROTATION, TRANSLATION
!     DES MODES PROPRES ET CONTRAINTS CALCULES SUR UNE STRUCTURE
!     MAILLEE ET TRANSPORTANT LES CHAMPS AUX NOEUDS MODAUX SUR
!     L EMPLACEMENT
!     DES AUTRES SOUS STRUCTURES NON MAILLEES AFIN DE DETERMINER
!     LES CHAMPS DE PRESSION ENGENDRES PAR LEUR MOUVEMENT DANS LE
!     FLUIDE
!     IN: K2 : MODEL : CHARACTER TRADUISANT LA DIMENSION DU FLUIDE
!     IN: K8 : BAMO : BASE MODALE ASSOCIEE AU MACROELEMENT D UNE
!                     SOUS STRUCTURE
!     IN : R8 : TGEOM :LISTE DE REELS DEFINISSANT 3 CMP DE TRANSLATION
!     DANS LE REPERE GLOBAL PUIS 3 ANGLES NAUTIQUES PRIS PAR RAPPORT
!     A O ORIGINE DU REPERE GLOBAL
!     IN : K14 : NU :NUMEROTATION ASSOCIEE AU MODELE FLUIDE
!     IN : K14 : NUM :NUMEROTATION ASSOCIEE AU MODELE INTERFACE
!     IN : K8 : MA : MATRICE DE RAIDEUR DU FLUIDE
!     IN : K8 : MOINT: MODELE INTERFACE
!     IN : I : IMODG : INDICE DE MODE D UNE SOUS STRUCTURE
!     IN : I : ISST  : INDICE  D UNE SOUS STRUCTURE
!     IN : I : IADX,IADY,IADZ : ADRESSES DES VECTEURS DE NOMS
!              DES CHAMNOS ASSOCIES PAR CMP DE DEPLACEMENT ET
!              PAR MODE D UNE SOUS STRUCTURE DONNEE D INDICE ISST
!     IN : I : IADRP : ADRESSE DU TABLEAU D ADRESSES DES VECTEURS
!              CONTENANT LES NOMS DES CHAMPS DE PRESSION
!     IN : I : ICOR(2) : TABLEAU CONTENANT LES ADRESSES
!                        JEVEUX DE TABLEAUX D'ENTIER
!              INDIQUANT LA CORRESPONDANCE ENTRE NOEUDS DE STRUCTURE
!              ET DE FLUIDE
!     IN : R8 : NORM1,NORM2 : NORMES DES VECTEURS TRANSLATION
!               ET ROTATION DECIDANT S IL Y A TRANSPORT OU NON
!               D UN CHAMNO DANS LE MAILLAGE
!     IN : I : NDBLE : INDICATEUR DE RECHERCHE DE NOEUD DOUBLE
!
!---------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/alimrs.h"
#include "asterfort/calflu.h"
#include "asterfort/chnucn.h"
#include "asterfort/chtpcn.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/prstoc.h"
#include "asterfort/resoud.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcmbl.h"
#include "asterfort/wkvect.h"
#include "asterfort/nueq_chck.h"
!
    aster_logical :: test1, test2, test3
    integer(kind=8) :: nbvale, nbrefe, isst, iadrp, i_ligr_mesh
    integer(kind=8) :: i, iad(2), iad3d(3), icor(2), ndble
    real(kind=8) :: tgeom(6), const(2)
    real(kind=8) :: tailmi, norm1, norm2, ca(3), sa(3)
    real(kind=8) :: val(2), val3d(3), tol
    character(len=1) :: typech(2), typcst(2), base
    character(len=2) :: model
    character(len=8) :: tcorx(2), tcory(2), tcorz(2)
    character(len=8) :: moint, ma, k8bid, maflui
    character(len=8) :: bamo, mailla, gd
    character(len=14) :: nu, num
    character(len=19) :: chtmpx, chtmpy, chtmpz, chcomb, vestoc
    character(len=19) :: vesolx, vesoly, vepr, vesolz, tampon, chcmb2
    character(len=19) :: chflu, chamnx, chamny, chamnz, newcha, pchno, nomch(2)
    character(len=19), parameter :: maprec = '&&OP0152.MAPREC', chsol = '&&OP0152.SOLUTION'
    character(len=19), parameter :: solveu = '&&OP0152.SOLVEUR'
    character(len=24), parameter :: criter = '&&RESGRA_GCPC'
    character(len=24) :: nomcha
    character(len=*) :: mateco
    complex(kind=8) :: cbid
    integer(kind=8) :: iadg, iadx, iady, iadz, iaut, ichar
    integer(kind=8) :: ichav, ilires, imodg, inoe
    integer(kind=8) :: iprn, iref, iret, ival
    integer(kind=8) :: ivaleu, k, nbchar, nbchav, nbnoe, ncmp
    integer(kind=8) :: nec
    character(len=24), pointer :: refe(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    real(kind=8), parameter :: tmin = 1.d-15, epsi = 1.d-2
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
! VERIFICATION SUPPLEMENTAIRE
!
!
    call jemarq()
    if (model .eq. 'AX') then
        call utmess('F', 'ALGORITH10_99')
    end if
!
    ca(1) = cos(tgeom(4))
    sa(1) = sin(tgeom(4))
    ca(2) = cos(tgeom(5))
    sa(2) = sin(tgeom(5))
    ca(3) = cos(tgeom(6))
    sa(3) = sin(tgeom(6))
!
    call rsexch('F', bamo, 'DEPL', imodg, nomcha, &
                iret)
!
! POUR CHAQUE MODE, ON FAIT SUBIR AU CHAMNO ASSOCIE
! LA ROTATION DEFINIE POUR LA SOUS-STRUCTURE EN QUESTION
!
    call jeveuo(nomcha(1:19)//'.VALE', 'L', vr=vale)
    call jelira(nomcha(1:19)//'.VALE', 'LONMAX', nbchav)
    call jeveuo(nomcha(1:19)//'.REFE', 'L', vk24=refe)
    call jelira(nomcha(1:19)//'.REFE', 'LONMAX', nbchar)
!
!
! CHANGEMENT DE VALEUR POUR DX ET DY (OU DZ)
!
    call dismoi('NUME_EQUA', nomcha, 'CHAM_NO', repk=pchno)
    call dismoi('NOM_MAILLA', nomcha, 'CHAM_NO', repk=mailla)
    call dismoi('NOM_GD', nomcha, 'CHAM_NO', repk=gd)
    call dismoi('NB_EC', gd, 'GRANDEUR', repi=nec)
!
    call dismoi('NB_NO_MAILLA', mailla, 'MAILLAGE', repi=nbnoe)
!
!
! - Protection: no matrix shrinking
!
    call nueq_chck(pchno)
    call jenonu(jexnom(pchno//'.LILI', '&MAILLA'), i_ligr_mesh)
    call jeveuo(jexnum(pchno//'.PRNO', i_ligr_mesh), 'L', iprn)
    call jeveuo(pchno//'.NUEQ', 'L', vi=nueq)
!
    newcha = '&&TRPROT.NCHNO'
!
    call wkvect(newcha//'.VALE', 'V V R', nbchav, ival)
!
    call wkvect(newcha//'.REFE', 'V V K24', nbchar, iref)
!
!
!
    do ichav = 1, nbchav
        zr(ival+ichav-1) = 0.0d0
    end do
!
! ROTATION DU CHAMNO PROPREMENT DITE
!
    do inoe = 1, nbnoe
!
        ivaleu = zi(iprn-1+(inoe-1)*(nec+2)+1)
        ncmp = zi(iprn-1+(inoe-1)*(nec+2)+2)
        iadg = iprn-1+(inoe-1)*(nec+2)+3
!
        if (ncmp .eq. 0) goto 22
!
! SI LES CMP DX ET DY SONT PRESENTES DANS LE NOEUD :
!
        test1 = exisdg(zi(iadg), 1)
        test2 = exisdg(zi(iadg), 2)
!
        if (model .eq. '2D') then
!
            if (test1 .and. test2) then
!
                do i = 1, ncmp
!
                    if (i .le. 2) then
!
                        iad(i) = nueq(ivaleu+i-1)
                        val(i) = vale(iad(i))
!
                    else
!
! DECALAGE DE 2 INDICES CAR ON DOIT METTRE A ZERO
! LES COMPOSANTES DE ROTATION APRES LES DEUX COMPOSANTES
! DE TRANSLATION DX ET DY
!
                        iaut = nueq(ivaleu+i-1)
                        zr(ival-1+iaut+i-1-2) = 0.0d0
!
                    end if
!
                end do
!
! C EST LA ROTATION LIMITEE AU 2D
!
                zr(ival-1+iad(1)) = ca(3)*ca(1)*val(1)+val(2)*(sa(3)* &
                                                               sa(2)*ca(1)-ca(3)*sa(1))
!
!
!
                zr(ival-1+iad(2)) = sa(1)*ca(2)*val(1)+val(2)*(ca(3)* &
                                                               ca(1)+sa(2)*sa(1)*sa(3))
!
! AUCUNE OU UNE SEULE EST PRESENTE : ON MET TOUT A ZERO
!
            else
!
                do i = 1, ncmp
                    iaut = nueq(ivaleu+i-1)
                    zr(ival-1+iaut+i-1) = 0.0d0
                end do
!
            end if
!
        else if (model .eq. '3D') then
!
!      ROTATION DU CHAMNO DE DEPL MODAL EN 3D
!
            test3 = exisdg(zi(iadg), 3)
!
            if (test1 .and. test2 .and. test3) then
                do i = 1, ncmp
!
                    if (i .le. 3) then
!
                        iad3d(i) = nueq(ivaleu+i-1)
                        val3d(i) = vale(iad3d(i))
!
                    else
!
                        iaut = nueq(ivaleu+i-1)
! DEBUG
!               IF ((IAUT+I-1-1-3).GT.36576) THEN
!
!
!               ENDIF
! FIN DEBUG
!
! DECALAGE DE 3 INDICES CAR ON DOIT METTRE A ZERO
! LES COMPOSANTES DE ROTATION APRES LES TROIS COMPOSANTES
! DE TRANSLATION DX , DY ET DZ
!
                        zr(ival-1+iaut+i-1-3) = 0.0d0
!
                    end if
!
                end do
!
!
                zr(ival-1+iad3d(1)) = ca(3)*ca(1)*val3d(1)+val3d(2)*( &
                                      sa(3)*sa(2)*ca(1)-ca(3)*sa(1)) &
                                      +val3d(3)*(ca(3)*sa(2) &
                                                 *ca(1)+sa(3)*sa(1))
!
!
                zr(ival-1+iad3d(2)) = sa(1)*ca(2)*val3d(1)+val3d(2)* &
                                      (ca(3)*ca(1)+sa(2)*sa(1)*sa(3)) &
                                      +val3d(3)*(ca(3)*sa( &
                                                 1)*sa(2)-sa(3)*ca(1))
!
                zr(ival-1+iad3d(3)) = -val3d(1)*sa(2)+val3d(2)*sa(3)*ca( &
                                      2)+val3d(3)*ca(3)*ca(2)
!
!             DO 5000 IBOUCL=1,3
!               IF ((IAD3D(IBOUCL)-1).GT.36576) THEN
!
!
!               ENDIF
!5000         CONTINUE
!
! AUCUNE OU UNE SEULE EST PRESENTE : ON MET TOUT A ZERO
!
            else
!
                do i = 1, ncmp
                    iaut = nueq(ivaleu+i-1)
                    zr(ival-1+iaut+i-1) = 0.0d0
                end do
!
            end if
        end if
!
22      continue
    end do
!
!
    do ichar = 1, nbchar
        zk24(iref+ichar-1) = refe(ichar)
    end do
!
!
!
! TEST POUR SAVOIR SI LE FLUIDE ET LA SOUS-STRUCTURE REPOSENT
! SUR LE MEME MAILLAGE OU NON
!
    call dismoi('NOM_MAILLA', moint, 'MODELE', repk=maflui)
!
!
    chamnx = '&&TRPROT.CHAMNX'
    chamny = '&&TRPROT.CHAMNY'
    chamnz = '&&TRPROT.CHAMNZ'
!
    if (maflui .ne. mailla) then
        base = 'V'
        call alimrs(mateco, mailla, maflui, moint, ndble, &
                    num, newcha, chamnx, 'DX', icor)
!
        call alimrs(mateco, mailla, maflui, moint, ndble, &
                    num, newcha, chamny, 'DY', icor)
        if (model .eq. '3D') then
!
            call alimrs(mateco, mailla, maflui, moint, ndble, &
                        num, newcha, chamnz, 'DZ', icor)
        end if
        call detrsd('CHAM_NO', newcha)
!
    else
!
!
! PLONGEMENT DU MODE PROPREMENT DIT SUR L INTERFACE
!
        base = 'V'
!
!----- PLONGEMENT DE LA COMPOSANTE DX QUI DEVIENT TEMPERATURE
!
        tcorx(1) = 'DX'
        tcorx(2) = 'TEMP'
!
        call chnucn(newcha, num, 2, tcorx, 'V', &
                    chamnx)
!
!
!
!----- PLONGEMENT DE LA COMPOSANTE DY QUI DEVIENT TEMPERATURE
!
        tcory(1) = 'DY'
        tcory(2) = 'TEMP'
!
        call chnucn(newcha, num, 2, tcory, 'V', &
                    chamny)
!
!
        if (model .eq. '3D') then
!
!----- PLONGEMENT DE LA COMPOSANTE DZ QUI DEVIENT TEMPERATURE
!
            tcorz(1) = 'DZ'
            tcorz(2) = 'TEMP'
!
            call chnucn(newcha, num, 2, tcorz, 'V', &
                        chamnz)
!
!
        end if
!
        call detrsd('CHAM_NO', newcha)
!
    end if
!
!---------ON TRANSPORTE CE MODE TRANSFORME EN TEMPERATURE
!-----SUR LES CONTOURS DES AUTRES SOUS - STRUCTURES NON MAILLEES
!
    typcst(1) = 'R'
    typcst(2) = 'R'
!
    const(1) = 1.0d0
    const(2) = 1.0d0
!
    typech(1) = 'R'
    typech(2) = 'R'
!
    vesolx = '&&VESOLX'
    vesoly = '&&VESOLY'
    vesolz = '&&VESOLZ'
!
    nomch(1) = vesolx(1:8)
    nomch(2) = vesoly(1:8)
!
!
    tol = r8prem()
!
    chtmpx = '&&TRPROT.CHTMPX'
    chtmpy = '&&TRPROT.CHTMPY'
    chtmpz = '&&TRPROT.CHTMPZ'
!
    if ((norm1 .gt. tol) .or. (norm2 .ne. 0.0d0)) then
!
!        CAS OU T DIFF 0 OU R DIFF 0
!
        call chtpcn(chamnx, tgeom, tailmi, tmin, epsi, &
                    base, chtmpx)
!
        call calflu(chtmpx, moint, mateco, num, vesolx, &
                    nbrefe, nbvale, 'X')
!
        call chtpcn(chamny, tgeom, tailmi, tmin, epsi, &
                    base, chtmpy)
!
        call calflu(chtmpy, moint, mateco, num, vesoly, &
                    nbrefe, nbvale, 'Y')
!
        ilires = ilires+1
        k = k+1
!
        vestoc = '&&TRPROT.TPXSTO'
        call prstoc(chtmpx, vestoc, ilires, k, iadx, &
                    nbvale, nbrefe)
!
        vestoc = '&&TRPROT.TPYSTO'
        call prstoc(chtmpy, vestoc, ilires, k, iady, &
                    nbvale, nbrefe)
!
        if (model .eq. '3D') then
!
            call chtpcn(chamnz, tgeom, tailmi, tmin, epsi, &
                        base, chtmpz)
!
            call calflu(chtmpz, moint, mateco, num, vesolz, &
                        nbrefe, nbvale, 'Z')
            vestoc = '&&TRPROT.TPZSTO'
            call prstoc(chtmpz, vestoc, ilires, k, iadz, &
                        nbvale, nbrefe)
!
        end if
    else
!
!        'CAS OU T=0 ET R=0'
!
        call calflu(chamnx, moint, mateco, num, vesolx, &
                    nbrefe, nbvale, 'X')
        call calflu(chamny, moint, mateco, num, vesoly, &
                    nbrefe, nbvale, 'Y')
!
        ilires = ilires+1
        k = k+1
!
        vestoc = '&&TRPROT.TPXSTO'
        call prstoc(chamnx, vestoc, ilires, k, iadx, &
                    nbvale, nbrefe)
!
!           IF (VESTOC(1:12).EQ.'TPXSTO000012') THEN
!
!
!           DO 12 I=1,NBVEP
!
!12         CONTINUE
!           ENDIF
!
        vestoc = '&&TRPROT.TPYSTO'
        call prstoc(chamny, vestoc, ilires, k, iady, &
                    nbvale, nbrefe)
!
        if (model .eq. '3D') then
!
            call calflu(chamnz, moint, mateco, num, vesolz, &
                        nbrefe, nbvale, 'Z')
            vestoc = '&&TRPROT.TPZSTO'
            call prstoc(chamnz, vestoc, ilires, k, iadz, &
                        nbvale, nbrefe)
!
        end if
    end if
!
!---ON RECOMBINE LES DEUX (TROIS)CHAMPS AUX NOEUDS DE TEMP ET ON CALCULE
!-----LE FLUX FLUIDE TOTAL.....
!
!
    chcomb = '&&CHCOMB'
!
    call vtcmbl(2, typcst, const, typech, nomch, &
                'R', chcomb)
!
    chcmb2 = '&&CHCMB2'
!
    if (model .eq. '3D') then
!
        typcst(1) = 'R'
        typcst(2) = 'R'
!
        const(1) = 1.0d0
        const(2) = 1.0d0
!
        typech(1) = 'R'
        typech(2) = 'R'
!
        nomch(1) = chcomb
        nomch(2) = vesolz(1:8)
!
        call vtcmbl(2, typcst, const, typech, nomch, &
                    'R', chcmb2)
!
    end if
!
!
    chflu = '&&TRPROT.CHFLU'
!
    if (model .eq. '3D') then
        tampon = chcmb2
    else
        tampon = chcomb
    end if
!
    call chnucn(tampon, nu, 0, k8bid, 'V', &
                chflu)
!
!----ON RESOUT L EQUATION DE LAPLACE
!
    call resoud(ma, maprec, solveu, ' ', 0, &
                chflu, chsol, 'V', [0.d0], [cbid], &
                criter, .true._1, 0, iret)
    call jedupc('V', chsol(1:19), 1, 'V', chflu(1:19), &
                .false._1)
    call detrsd('CHAMP_GD', chsol)
!
!--------ON REPLONGE LA PRESSION SUR L INTERFACE
!-----------------QU 'ON STOCKE
!
    vepr = '&&TRPROT.VEPR'
    call chnucn(chflu, num, 0, k8bid, 'V', &
                vepr)
!
    vestoc = '&&TRPROT.VESTOC'
    call prstoc(vepr, vestoc, ilires, k, zi(iadrp+isst-1), &
                nbvale, nbrefe)
!
    call detrsd('CHAM_NO', vepr)
    call detrsd('CHAMP_GD', chcomb)
    call detrsd('CHAMP_GD', chcmb2)
    call detrsd('CHAM_NO', chflu)
    call detrsd('CHAM_NO', vesolx)
    call detrsd('CHAM_NO', vesoly)
    call detrsd('CHAM_NO', vesolz)
    call detrsd('CHAM_NO', chtmpx)
    call detrsd('CHAM_NO', chtmpy)
    call detrsd('CHAM_NO', chtmpz)
    call detrsd('CHAM_NO', chamnx)
    call detrsd('CHAM_NO', chamny)
    call detrsd('CHAM_NO', chamnz)
!
    call jeexin(criter(1:19)//'.CRTI', iret)
    if (iret .ne. 0) then
        call jedetr(criter(1:19)//'.CRTI')
        call jedetr(criter(1:19)//'.CRTR')
        call jedetr(criter(1:19)//'.CRDE')
    end if
!
    call jedema()
!
end subroutine
