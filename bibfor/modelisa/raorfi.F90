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

subroutine raorfi(noma, ligrel, noepou, cara, coorig, &
                  eg1, eg2, eg3, typrac, rayon)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/carcou.h"
#include "asterfort/dismoi.h"
#include "asterfort/etenca.h"
#include "asterfort/infniv.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mecact.h"
#include "asterfort/normev.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvlg.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: info, ifm
    character(len=8) :: noepou, noma, cara
    character(len=19) :: ligrel
    real(kind=8) :: coorig(3), eg1(3), eg2(3), eg3(3), r
!     POUR LE RACCORD (COQUE OU 3D)_TUYAU
!
!
    integer(kind=8) :: nma, ima, inopou, iconex, nbno, ier, i, j, ntseg3, nutyma
    integer(kind=8) :: idesc, ncmpmx, ivale, igd, idebgd
    integer(kind=8) :: inopo1, inopo2, icoude, nno, ntseg4
    character(len=8) :: nomgd, nocmp(3), noepo1, noepo2, typrac
    character(len=19) :: chcara
    real(kind=8) :: coorif(3), gpl(3), gpg(3), pgl(3, 3)
    real(kind=8) :: el1(3), el2(3), el3(3), coono1(3), coono2(3), e1(3), nore1
    real(kind=8) :: rayon
    real(kind=8) :: pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), l, omega, theta
    real(kind=8) :: pgl4(3, 3)
    integer(kind=8), pointer :: typmail(:) => null()
    real(kind=8), pointer :: coor(:) => null()
    integer(kind=8), pointer :: ptma(:) => null()
!
    call jemarq()
    call infniv(ifm, info)
!
!     RECHERCHE DE LA MAILLE IMA  CONTENANT LE NOEUD NOEPOU
!
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nma)
    ima = 0
    inopou = char8_to_int(noepou)
    do i = 1, nma
        call jeveuo(jexnum(noma//'.CONNEX', i), 'L', iconex)
        call jelira(jexnum(noma//'.CONNEX', i), 'LONMAX', nbno)
        do j = 1, nbno
            if (zi(iconex+j-1) .eq. inopou) then
                if (ima .eq. 0) then
                    ima = i
                    inopo1 = zi(iconex)
                    inopo2 = zi(iconex+1)
                else
                    call utmess('F', 'MODELISA6_36', sk=noepou)
                end if
            end if
        end do
    end do
!
!     VERIFICATION QUE LA MAILLE IMA EST UN SEG3
!
    call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG3'), ntseg3)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG4'), ntseg4)
    call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
    nutyma = typmail(ima)
    if (nutyma .eq. ntseg3) then
        nno = 3
    else if (nutyma .eq. ntseg4) then
        nno = 4
    else
        call utmess('F', 'MODELISA6_37', sk=noepou)
    end if
!
!     RECUPERATION DES ANGLES NAUTIQUES DANS LA CARTE ORIENTATION
!
    chcara = cara(1:8)//'.CARORIEN'
    call etenca(chcara, ligrel, ier)
    ASSERT(ier .eq. 0)
    nomgd = 'CAORIE_R'
    call jeveuo(chcara//'.DESC', 'L', idesc)
!
! --- NOMBRE DE COMPOSANTES ASSOCIEES A LA GRANDEUR CF NUROTA
!
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgd), 'LONMAX', ncmpmx)
!
! --- TABLEAU DE VALEURS DE LA CARTE CHCARA
!
    call jeveuo(chcara//'.VALE', 'L', ivale)
!
! --- RECUPERATION DU VECTEUR D'ADRESSAGE DANS LA CARTE CREE PAR ETENCA
!
    call jeveuo(chcara//'.PTMA', 'L', vi=ptma)
!
!     RECUPERATION DES ANGLES NAUTIQUES
!
    if (ptma(ima) .ne. 0) then
        igd = ptma(ima)
        idebgd = (igd-1)*ncmpmx
!        RECUPERATION DE L'ORIENTATION
        call carcou(zr(ivale+idebgd), l, pgl, r, theta, &
                    pgl1, pgl2, pgl3, pgl4, nno, &
                    omega, icoude)
    else
        call utmess('F', 'MODELISA6_38')
    end if
!
!     CALCUL DU VECTEUR E1 ORIENTANT LA MAILLE TUYAU
!
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=coor)
    coono1(1) = coor(3*(inopo1-1)+1)
    coono1(2) = coor(3*(inopo1-1)+2)
    coono1(3) = coor(3*(inopo1-1)+3)
    coono2(1) = coor(3*(inopo2-1)+1)
    coono2(2) = coor(3*(inopo2-1)+2)
    coono2(3) = coor(3*(inopo2-1)+3)
    e1 = coono2-coono1
    call normev(e1, nore1)
    nocmp(1) = 'X'
    nocmp(2) = 'Y'
    nocmp(3) = 'Z'
!
    call mecact('V', typrac//'.CAXE_TUY', 'LIGREL', ligrel, 'GEOM_R', &
                ncmp=3, lnomcmp=nocmp, vr=e1)
!
!     CALCUL DU VECTEUR GPL, AVEC P ORIGINE DE PHI
!
    gpl(1) = 0.d0
    gpl(2) = 0.d0
    gpl(3) = -rayon
!
!     PASSAGE DE GPL DANS LE REPERE GLOBAL ET COORDONNEES DE P
!
!      CALL MATROT (ORIEN,PGL)
    call utpvlg(1, 3, pgl, gpl, gpg)
    coorif(1) = gpg(1)+coorig(1)
    coorif(2) = gpg(2)+coorig(2)
    coorif(3) = gpg(3)+coorig(3)
!
! --- NOTATION DANS LA CARTE DE NOM '&&RAPOCO.CAORIFI' DES
! --- COORDONNEES DU POINT ORGINE DE PHI SUR LA SECTION DE RACCORD
!
    nocmp(1) = 'X'
    nocmp(2) = 'Y'
    nocmp(3) = 'Z'
!
    call mecact('V', typrac//'.CAORIFI', 'LIGREL', ligrel, 'GEOM_R', &
                ncmp=3, lnomcmp=nocmp, vr=coorif)
!
!     COORDONNEES DES VECTEURS UNITAIRES DANS LE REPERE GLOBAL
!
    el1(1) = 1.d0
    el1(2) = 0.d0
    el1(3) = 0.d0
!
    el2(1) = 0.d0
    el2(2) = 1.d0
    el2(3) = 0.d0
!
!        A CAUSE DE LA DEFINITION DU REPERE LOCAL, OU Z EST OPPOSE A
!        CELUI OBTENU PAR ROTATION DE ALPHA, BETA, GAMMA, IL FAUT
!        MODIFIER LE SIGNE DE Z (VERIF FAITE SUR LA FLEXION HORS PLAN)
!
    el3(1) = 0.d0
    el3(2) = 0.d0
!      EL3(3)=-1.D0 NON DIRECT  MAIS FTN ACTUEL
!      EL3(3)=1.D0   SERAIT LOGIQUE
    el3(3) = 1.d0
!
    call utpvlg(1, 3, pgl, el1, eg1)
    call utpvlg(1, 3, pgl, el2, eg2)
    call utpvlg(1, 3, pgl, el3, eg3)
!
    if (info .eq. 2) then
        noepo1 = int_to_char8(inopo1)
        noepo2 = int_to_char8(inopo2)
        ifm = iunifi('MESSAGE')
        write (ifm, *) 'RAYON DE LA SECTION COQUE OU 3D ', rayon
        write (ifm, *) 'BARYCENTRE DE LA SECTION COQUE OU 3D ', coorig
        write (ifm, *) 'POINT ORIGINE DE LA GENERATRICE ', coorif
        write (ifm, *) 'VECTEUR AXE DU TUYAU : E1 ', e1
        write (ifm, *) 'NOEUDS AXE DU TUYAU :  ', noepo1, noepo2
        write (ifm, *) 'VECTEURS UNITAIRES DU TUYAU : E1 ', eg1
        write (ifm, *) 'VECTEURS UNITAIRES DU TUYAU : E2 ', eg2
        write (ifm, *) 'VECTEURS UNITAIRES DU TUYAU : E3 ', eg3
    end if
!
    call jedetr(chcara//'.PTMA')
    call jedema()
end subroutine
