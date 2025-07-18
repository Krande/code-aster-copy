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

subroutine op0142()
    implicit none
!      DEFINITION D'UN PROFIL DE VITESSE LE LONG D'UNE STRUCTURE
!      EN FONCTION D'UNE ABSCISSE CURVILIGNE.
!     STOCKAGE DANS UN OBJET DE TYPE FONCTION
! ----------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/foimpr.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/i2extf.h"
#include "asterfort/i2sens.h"
#include "asterfort/i2tgrm.h"
#include "asterfort/i2vois.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/ordonn.h"
#include "asterfort/prvite.h"
#include "asterfort/reliem.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: pnoe, ptch
    character(len=2) :: prolgd
    character(len=4) :: interp(2)
    character(len=8) :: nommai
    character(len=16) :: nomcmd, tprof, typfon
    character(len=16) :: motcleini(2), motclefin(2), typmcl(2)
    character(len=19) :: nomfon
    character(len=24) :: cooabs, typmai, connex
    character(len=24) :: conseg, typseg, lisnoini, lisnofin
    character(len=8) :: typm
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iabs, iach, iacnex, iagm, iav1, iav2
    integer(kind=8) :: ibid, iexi, ifm, ij, im, ima1, ima2
    integer(kind=8) :: ind, ing, ino, iplac1, iplac2, iseg2, isens
    integer(kind=8) :: itp, itym, itypm, jgcnx, jnof, jnoi, kseg, l, labs
    integer(kind=8) :: lnoe, lpro, lval, mi, n3, nbbav, nbchm, nbno
    integer(kind=8) :: nbnoma, nbpoi1, nbrm21, nbrma, nbrma1, nbrma2, nbrse1
    integer(kind=8) :: nbrse2, nbrseg, nbseg2, niv, nnoe, num1, num2
    integer(kind=8) :: numno
    real(kind=8) :: rvale
!-----------------------------------------------------------------------
!
    data motcleini/'NOEUD_INIT', 'GROUP_NO_INIT'/
    data motclefin/'NOEUD_FIN', 'GROUP_NO_FIN'/
    data typmcl/'NOEUD', 'GROUP_NO'/
!
!=======================================================================
!
    call jemarq()
!
! --- RECUPERATION DU NIVEAU D'IMPRESSION
    call infmaj()
    call infniv(ifm, niv)
!
    call getres(nomfon, typfon, nomcmd)
!
    call getvid(' ', 'MAILLAGE', scal=nommai, nbret=l)
    lisnoini = '&&OP0142.LISTE_NO_INI'
    lisnofin = '&&OP0142.LISTE_NO_FIN'
    call reliem(' ', nommai, 'NU_NOEUD', ' ', 1, &
                2, motcleini, typmcl, lisnoini, nbno)
    ASSERT(nbno .eq. 1)
    call reliem(' ', nommai, 'NU_NOEUD', ' ', 1, &
                2, motclefin, typmcl, lisnofin, nbno)
    ASSERT(nbno .eq. 1)
!
!     --- CONSTRUCTION DES OBJETS DU CONCEPT MAILLAGE ---
!
    cooabs = nommai//'.ABSC_CURV .VALE'
    connex = nommai//'.CONNEX'
    typmai = nommai//'.TYPMAIL'
!
    call jeveuo(lisnoini, 'L', jnoi)
    num1 = zi(jnoi-1+1)
    call jeveuo(lisnofin, 'L', jnof)
    num2 = zi(jnof-1+1)
!
    call jeexin(cooabs, iexi)
    if (iexi .eq. 0) then
        call utmess('F', 'UTILITAI2_84')
    end if
    call jeveuo(cooabs, 'L', labs)
!
    call getvtx(' ', 'INTERPOL', nbval=2, vect=interp, nbret=n3)
    if (n3 .eq. 1) interp(2) = interp(1)
    call getvtx(' ', 'PROL_GAUCHE', scal=prolgd(1:1), nbret=n3)
    call getvtx(' ', 'PROL_DROITE', scal=prolgd(2:2), nbret=n3)
!
!     --- CREATION ET REMPLISSAGE DE L'OBJET NOMFON//'.PROL'
!
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', 'G V K24', 6, lpro)
!
    zk24(lpro) = 'FONCTION'
    zk24(lpro+1) = interp(1)//interp(2)
    zk24(lpro+2) = 'ABSC '
    zk24(lpro+3) = 'VITE'
    zk24(lpro+4) = prolgd
    zk24(lpro+5) = nomfon
!
!     --- LECTURE DES CARACTERISTIQUES DU GROUPE DE MAILLES : ADRESSE
!                   ET NOMBRE DE MAILLES
!
    call jelira(typmai, 'LONMAX', nbrma)
    call wkvect('&&OP0142.MAILL.TEMP', 'V V I', nbrma, iagm)
    do ij = 1, nbrma
        zi(iagm+ij-1) = ij
    end do
    nbrma2 = 2*nbrma
    nbrma1 = nbrma+1
!     --- CREATION D OBJETS TEMPORAIRES ---
!
    call wkvect('&&OP0142.TEMP.VOIS1', 'V V I', nbrma, iav1)
    call wkvect('&&OP0142.TEMP.VOIS2', 'V V I', nbrma, iav2)
    call wkvect('&&OP0142.TEMP.CHM  ', 'V V I', nbrma1, ptch)
    call wkvect('&&OP0142.TEMP.LNOE ', 'V V I', nbrma1, lnoe)
    call wkvect('&&OP0142.TEMP.NNOE ', 'V V K8', nbrma1, nnoe)
    call wkvect('&&OP0142.TEMP.PNOE ', 'V V I', nbrma1, pnoe)
    call wkvect('&&OP0142.TEMP.IABS ', 'V V R8', nbrma1, iabs)
    call wkvect('&&OP0142.TEMP.IACHM', 'V V I', nbrma2, iach)
    call wkvect('&&OP0142.TEMP.IPOI1', 'V V I', nbrma, ima1)
    call wkvect('&&OP0142.TEMP.ISEG2', 'V V I', nbrma, ima2)
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
    conseg = '&&OP0142.CONNEX'
    typseg = '&&OP0142.TYPMAI'
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
!    do i =1, nbrseg
!        write(*,*) i, ':', zi(iach+i-1)
!    end do
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
    do i = 1, nbrse1
        if (zi(lnoe+i-1) .eq. num1) then
            iplac1 = i
        end if
        if (zi(lnoe+i-1) .eq. num2) then
            iplac2 = i
        end if
    end do
    if (iplac1 .ge. iplac2) then
        call utmess('F', 'UTILITAI2_85')
    end if
!
!     --- CREATION DE L OBJET .VALE SUR LA GLOBALE ---
!
    call jeveuo(cooabs, 'L', labs)
    nbrm21 = nbrse1*2
    call wkvect(nomfon//'.VALE', 'G V R8', nbrm21, lval)
!
    do i = 1, nbrseg
        zr(lval+(i-1)) = min(zr(labs+4*(i-1)), zr(labs+4*(i-1)+1))
    end do
!
    zr(lval+nbrseg) = max(zr(labs+4*(nbrseg-1)), zr(labs+4*(nbrseg-1)+1))
!
    call getvtx('VITE ', 'PROFIL', iocc=1, scal=tprof, nbret=ibid)
    if (tprof .eq. 'UNIFORME') then
        call getvr8('VITE ', 'VALE', iocc=1, scal=rvale, nbret=ibid)
        do i = 1, nbrse1
            if (i .ge. iplac1 .and. i .le. iplac2) then
                zr(lval+nbrse1+i-1) = rvale
            else
                zr(lval+nbrse1+i-1) = 0.d0
            end if
        end do
    else
        call getvis('VITE ', 'NB_BAV', iocc=1, scal=nbbav, nbret=ibid)
        if (nbbav .eq. 0) then
            itp = 1
        else if (nbbav .eq. 2) then
            itp = 2
        else if (nbbav .eq. 3) then
            itp = 3
        end if
!
        call prvite(zr(lval), nbrm21, iplac1, iplac2, itp)
!
    end if
!
!
!     --- VERIFICATION QU'ON A BIEN CREER UNE FONCTION ---
!         ET REMISE DES ABSCISSES EN ORDRE CROISSANT
    call ordonn(nomfon, 0)
!
!     --- CREATION D'UN TITRE ---
    call titre()
!
!     --- IMPRESSIONS ---
    if (niv .gt. 1) call foimpr(nomfon, niv, ifm, 0, ' ')
    call jedema()
end subroutine
