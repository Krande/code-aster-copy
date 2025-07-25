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

subroutine defcur(vecr1, veck1, nb, vecr2, nv, &
                  nommai, nm, prolgd, interp)
    implicit none
!     LECTURE DE LA DEFINITION D'UNE FONCTION (ABSCISSE CURVILIGNE)
!     STOCKAGE DANS UN OBJET DE TYPE FONCTION
! ----------------------------------------------------------------------
!     IN  : VECR1  : VECTEUR DE LONG. NB, CONTIENT LES VALEURS DE LA
!                    FONCTION DEFINIE AUX NOEUDS.
!     IN  : VECK1  : VECTEUR DE LONG. NB, CONTIENT LES NOMS DES NOEUDS.
!     OUT : VECR2  : VECTEUR DE LONG. NV, CONTIENT LES VALEURS DE LA
!                    FONCTION.
!     IN  : MONMAI : NOM DU MAILLAGE.
!     IN  : NM     : NOMBRE DE MAILLES.
#include "jeveux.h"
#include "asterfort/i2extf.h"
#include "asterfort/i2sens.h"
#include "asterfort/i2tgrm.h"
#include "asterfort/i2vois.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/prfcur.h"
#include "asterfort/utmess.h"
#include "asterfort/vefcur.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: ptch, pnoe, nb, nv
    real(kind=8) :: vecr1(nb), vecr2(nv)
    character(len=2) :: prolgd
    character(len=8) :: nommai, interp, veck1(nb)
    character(len=24) :: cooabs, nomnoe, connex, typmai
    character(len=8) :: typm
    character(len=24) :: conseg, typseg
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iach, iacnex, iagm, iav1, iav2, iexi
    integer(kind=8) :: ii, im, ima1, ima2, ind, ing, ino
    integer(kind=8) :: iseg2, isens, itym, itypm, jgcnx, ji, jj
    integer(kind=8) :: jp, kk, kseg, labs, lnoe, lvali, mi
    integer(kind=8) :: nbchm, nbnoma, nbpoi1, nbrma, nbrma1, nbrma2
    integer(kind=8) :: nbrse1, nbrse2, nbrseg, nbseg2, nm, numno
!
!-----------------------------------------------------------------------
    call jemarq()
    nbrma = nm
!
!     --- CONSTRUCTION DES OBJETS DU CONCEPT MAILLAGE ---
!
    nomnoe = nommai//'.NOMNOE'
    cooabs = nommai//'.ABSC_CURV .VALE'
    connex = nommai//'.CONNEX'
    typmai = nommai//'.TYPMAIL'
!
!     --- VERIFICATION DE L EXISTENCE DE L ABSCISSE CURVILIGNE ---
!
    call jeexin(cooabs, iexi)
    if (iexi .eq. 0) then
        call utmess('F', 'UTILITAI_46')
    end if
    call jeveuo(cooabs, 'L', labs)
!     --- CREATION D OBJETS TEMPORAIRES ---
!
    call wkvect('&&DEFOCU.TEMP      ', 'V V I', nbrma, iagm)
    do ii = 1, nbrma
        zi(iagm+ii-1) = ii
    end do
    nbrma2 = 2*nbrma
    nbrma1 = nbrma+1
    call wkvect('&&DEFOCU.TEMP.VOIS1', 'V V I', nbrma, iav1)
    call wkvect('&&DEFOCU.TEMP.VOIS2', 'V V I', nbrma, iav2)
    call wkvect('&&DEFOCU.TEMP.CHM  ', 'V V I', nbrma1, ptch)
    call wkvect('&&DEFOCU.TEMP.IACHM', 'V V I', nbrma2, iach)
    call wkvect('&&DEFOCU.TEMP.LNOE', 'V V I', nbrma1, lnoe)
    call wkvect('&&DEFOCU.TEMP.PNOE', 'V V I', nv, pnoe)
    call wkvect('&&DEFOCU.TEMP.IPOI1', 'V V I', nbrma, ima1)
    call wkvect('&&DEFOCU.TEMP.ISEG2', 'V V I', nbrma, ima2)
!
!     TRI DES MAILLES POI1 ET SEG2
    nbseg2 = 0
    nbpoi1 = 0
    kseg = 0
    do im = 1, nbrma
        call jeveuo(typmai, 'L', itypm)
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(itypm+im-1)), typm)
        if (typm .eq. 'SEG2') then
            kseg = zi(itypm+im-1)
            nbseg2 = nbseg2+1
            zi(ima2+nbseg2-1) = im
        else if (typm .eq. 'POI1') then
            nbpoi1 = nbpoi1+1
            zi(ima1+nbpoi1-1) = im
        else
            call utmess('F', 'MODELISA_2')
        end if
    end do
    conseg = '&&DEFOCU.CONNEX'
    typseg = '&&DEFOCU.TYPMAI'
    call wkvect(typseg, 'V V I', nbrma, itym)
    do im = 1, nbrma
        zi(itym-1+im) = kseg
    end do
!     IL FAUT CREER UNE TABLE DE CONNECTIVITE POUR LES SEG2
!
    nbnoma = 2*nbseg2
    nbrseg = nbseg2
    nbrse1 = nbseg2+1
    nbrse2 = nbseg2*2
    call jecrec(conseg, 'V V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbseg2)
    call jeecra(conseg, 'LONT', nbnoma)
    do iseg2 = 1, nbseg2
        im = zi(ima2+iseg2-1)
        call jelira(jexnum(connex, im), 'LONMAX', nbnoma)
        call jeveuo(jexnum(connex, im), 'L', iacnex)
        call jeecra(jexnum(conseg, iseg2), 'LONMAX', nbnoma)
        call jeveuo(jexnum(conseg, iseg2), 'E', jgcnx)
        do ino = 1, nbnoma
            numno = zi(iacnex-1+ino)
            zi(jgcnx+ino-1) = numno
        end do
    end do
!
    call i2vois(conseg, typseg, zi(iagm), nbrseg, zi(iav1), &
                zi(iav2))
    call i2tgrm(zi(iav1), zi(iav2), nbrseg, zi(iach), zi(ptch), &
                nbchm)
    call i2sens(zi(iach), nbrse2, zi(iagm), nbrseg, conseg, &
                typseg, zr(labs))
!
!     --- CREATION D UNE LISTE ORDONNEE DE NOEUDS ---
    do i = 1, nbrseg
        isens = 1
        mi = zi(iach+i-1)
        if (mi .lt. 0) then
            mi = -mi
            isens = -1
        end if
        call i2extf(mi, 1, conseg, typseg, ing, &
                    ind)
        if (isens .eq. 1) then
            zi(lnoe+i-1) = ing
            zi(lnoe+i) = ind
        else
            zi(lnoe+i) = ing
            zi(lnoe+i-1) = ind
        end if
    end do
!
!
!     --- VERIFICATION DE LA DEFINITION DE LA FONCTION ---
!
    call vefcur(zi(lnoe), nbrse1, veck1, zi(pnoe), nb, &
                nomnoe)
!
!
    call wkvect('&&DEFOCU.TEMP.VALE', 'V V R8', nv, lvali)
!
    do i = 1, nbrseg
        zr(lvali+2*(i-1)) = min(zr(labs+4*(i-1)), zr(labs+4*(i-1)+1))
    end do
!
    zr(lvali+2*nbrseg) = max(zr(labs+4*(nbrseg-1)), zr(labs+4*(nbrseg-1)+1))
!
    do i = 1, nb
        kk = 2*(zi(pnoe+i-1)-1)+1
        zr(lvali+kk) = vecr1(i)
    end do
!
    do i = 1, nb
        jp = zi(pnoe+i-1)
        ji = i
        do jj = i, nb
            if (zi(pnoe+jj-1) .lt. jp) then
                ji = jj
                jp = zi(pnoe+jj-1)
            end if
        end do
        zi(pnoe+ji-1) = zi(pnoe+i-1)
        zi(pnoe+i-1) = jp
    end do
!
!     ------------- INTERPOLATION DE LA FONCTION -------------
!     --- PROLONGEMENT DE LA FONCTION A GAUCHE ET A DROITE ---
!
    call prfcur(zi(pnoe), nb, zr(lvali), nv, interp, &
                prolgd)
!
!     --- REMPLISSAGE DE L OBJET .VALE ---
!
    do i = 1, nbrse1
        vecr2(i) = zr(lvali+2*(i-1))
        vecr2(nbrse1+i) = zr(lvali+2*(i-1)+1)
    end do
!
!     --- MENAGE ---
    call jedetr('&&DEFOCU.TEMP      ')
    call jedetr('&&DEFOCU.TEMP.VOIS1')
    call jedetr('&&DEFOCU.TEMP.VOIS2')
    call jedetr('&&DEFOCU.TEMP.CHM  ')
    call jedetr('&&DEFOCU.TEMP.IACHM')
    call jedetr('&&DEFOCU.TEMP.LNOE')
    call jedetr('&&DEFOCU.TEMP.PNOE')
    call jedetr('&&DEFOCU.TEMP.IPOI1')
    call jedetr('&&DEFOCU.TEMP.ISEG2')
    call jedetr('&&DEFOCU.TEMP.VALE')
    call jedetr('&&DEFOCU.CONNEX')
    call jedetr('&&DEFOCU.TYPMAI')
!
    call jedema()
end subroutine
