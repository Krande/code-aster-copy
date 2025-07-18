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

subroutine xpocrf(modele, maxfem, mftot, nftot)
    implicit none
!
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/ltnotb.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: modele, maxfem
    integer(kind=8) :: mftot, nftot
!
!
!     BUT : GENERER UN MAILLAGE DESTINE UNIQUEMENT AU POST-TRAITEMENT
!           DU FOND DE FISSURE
!   IN
!       MODELE : MODELE FISSURE
!       MFTOT  : NOMBRE TOTAL DE MAILLES DE FONDS DE FISSURES
!       NFTOT  : NOMBRE TOTAL DE NOEUDS DE FONDS DE FISSURES
!       NFCOMF : NOMBRE TOTAL DE CONNEXIONS DANS LES MAILLES
!   IN/OUT
!       MAXFEM : MAILLAGE FISSURE
!     =================================================================
!     ------------------------------------------------------------------
    integer(kind=8) :: ifiss, ifon, ifon1, ifon2, ino, j, ima
    integer(kind=8) :: nfiss, nfond, nfon, ntseg2, ntpoi1
    integer(kind=8) :: icompt, ncompt, nufon, iagma, iagno
    integer(kind=8) :: ntail, ndim, nbmax
    integer(kind=8) :: ibid, nnntot, iret, jconx, igr
    integer(kind=8) :: jva00, jva0, jva1, jva2, jva3
    integer(kind=8) ::  jnom
    character(len=2) :: chn1, chn2
    character(len=6) :: chn
    character(len=8) :: fiss, mo
    character(len=19) :: nomtab, coord2
    character(len=24) :: nom, nogno, nogma
    integer(kind=8), pointer :: fondmult(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
!
!
!     INITIALISATION
!
!     RECUPERATION DES CARACTERISTIQUES DES FONDS DE FISSURES
    mo = modele
    nom = mo//'.FISS                   '
    call jeveuo(nom, 'L', jnom)
    call jelira(nom, 'LONUTI', nfiss)
!
    call dismoi('NB_NO_MAILLA', maxfem, 'MAILLAGE', repi=nnntot)
    call dismoi('NB_MA_MAILLA', maxfem, 'MAILLAGE', repi=nbmax)
    call dismoi('DIM_GEOM', mo, 'MODELE', repi=ndim)
    if (.not. (ndim .eq. 2 .or. ndim .eq. 3)) then
        call utmess('F', 'MODELISA2_6')
    end if
!
    call jeveuo(maxfem//'.TYPMAIL', 'E', vi=typmail)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), ntpoi1)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG2'), ntseg2)
!
    if ((nftot .gt. 0) .and. (mftot .gt. 0)) then
!
!       ATTRIBUTION DU NOM DES NOEUDS DU FOND DE FISSURE
        do ino = 1, nftot
            call codent(ino, 'G', chn)
            call jecroc(jexnom(maxfem//'.NOMNOE', 'NF'//chn))
        end do
!       ATTRIBUTION DU NOM DES MAILLES DU FOND DE FISSURE
        do ima = 1, mftot
            call codent(ima, 'G', chn)
            call jecroc(jexnom(maxfem//'.NOMMAI', 'MF'//chn))
        end do
!
        ncompt = 0
        icompt = 0
        coord2 = maxfem//'.COORDO'
        call jeveuo(coord2//'.VALE', 'E', vr=vale)
!
        do ifiss = 1, nfiss
!
            fiss = zk8(jnom)
            call jeexin(fiss//'.FONDFISS', iret)
            if (iret .ne. 0) then
                call jeveuo(fiss//'.FONDMULT', 'L', vi=fondmult)
                call ltnotb(fiss, 'FOND_FISS', nomtab)
                call jeveuo(nomtab//'.0001', 'L', jva0)
                call jelira(nomtab//'.0001', 'LONUTI', nfon)
                if (ndim .eq. 3) then
                    call jeveuo(nomtab//'.0002', 'L', jva00)
                    call jeveuo(nomtab//'.0004', 'L', jva1)
                    call jeveuo(nomtab//'.0005', 'L', jva2)
                    call jeveuo(nomtab//'.0006', 'L', jva3)
                else
                    call jeveuo(nomtab//'.0002', 'L', jva1)
                    call jeveuo(nomtab//'.0003', 'L', jva2)
                end if
                nfond = zi(jva0-1+nfon)
                call codent(ifiss, 'D0', chn1)
                do ifon = 1, nfond
!
                    call codent(ifon, 'D0', chn2)
                    nogma = 'MF_'//chn1//'_'//chn2
                    nogno = 'NF_'//chn1//'_'//chn2
!
!           VERIFICATION DE L'ABSENCE DE GROUPE AYANT LE MEME NOM
                    call jenonu(jexnom(maxfem//'.GROUPEMA', nogma), ibid)
                    if (ibid .gt. 0) then
                        call utmess('F', 'ALGELINE3_7', sk=nogma)
                    end if
                    call jenonu(jexnom(maxfem//'.GROUPENO', nogno), ibid)
                    if (ibid .gt. 0) then
                        call utmess('F', 'SOUSTRUC_37', sk=nogno)
                    end if
!
                    ntail = fondmult(2*ifon)-fondmult(2*ifon-1)+1
!
!           CONSTRUCTION DES GROUPES DE MAILLES DU FOND DE FISSURE
                    call jecroc(jexnom(maxfem//'.GROUPEMA', nogma))
                    call jenonu(jexnom(maxfem//'.GROUPEMA', nogma), igr)
                    call jeecra(jexnum(maxfem//'.GROUPEMA', igr), 'LONMAX', nfon)
                    call jeecra(jexnum(maxfem//'.GROUPEMA', igr), 'LONUTI', max((ntail-1), 1))
                    call jeveuo(jexnum(maxfem//'.GROUPEMA', igr), 'E', iagma)
!
!           CONSTRUCTION DES GROUPES DE NOEUDS DU FOND DE FISSURE
                    call jecroc(jexnom(maxfem//'.GROUPENO', nogno))
                    call jenonu(jexnom(maxfem//'.GROUPENO', nogno), igr)
                    call jeecra(jexnum(maxfem//'.GROUPENO', igr), 'LONMAX', nfon)
                    call jeecra(jexnum(maxfem//'.GROUPENO', igr), 'LONUTI', ntail)
                    call jeveuo(jexnum(maxfem//'.GROUPENO', igr), 'E', iagno)
!
!           COORDONNEES DES NOEUDS
                    if (ndim .eq. 3) then
                        do ifon2 = 1, ntail
                            ifon1 = ifon2+fondmult(2*ifon-1)-1
                            ino = nnntot-nftot+ifon1+ncompt
                            vale(1+3*(ino-1)-1+1) = zr(jva1-1+ifon1)
                            vale(1+3*(ino-1)-1+2) = zr(jva2-1+ifon1)
                            vale(1+3*(ino-1)-1+3) = zr(jva3-1+ifon1)
                        end do
                    else
                        ifon1 = fondmult(2*ifon-1)
                        ino = nnntot-nftot+ifon1+ncompt
                        vale(1+3*(ino-1)-1+1) = zr(jva1-1+ifon1)
                        vale(1+3*(ino-1)-1+2) = zr(jva2-1+ifon1)
                    end if
!
!           CONNEXITE DES NOEUDS
                    do ifon2 = 1, ntail
                        ifon1 = ifon2+fondmult(2*ifon-1)-1
                        ino = nnntot-nftot+ifon1+ncompt
                        if (ndim .eq. 3) then
                            nufon = zi(jva00-1+ifon1)
                            if (nufon .gt. 1) then
                                ima = nbmax-mftot+icompt+1
                                typmail(ima) = ntseg2
                                call jeecra(jexnum(maxfem//'.CONNEX', ima), 'LONMAX', 2)
                                call jeveuo(jexnum(maxfem//'.CONNEX', ima), 'E', jconx)
                                do j = 1, 2
                                    zi(jconx-1+j) = ino+j-2
                                end do
                                icompt = icompt+1
                                zi(iagma-1+ifon2-1) = ima
                            end if
                        else if (ntail .eq. 1) then
                            ima = nbmax-mftot+icompt+1
                            typmail(ima) = ntpoi1
                            call jeecra(jexnum(maxfem//'.CONNEX', ima), 'LONMAX', 1)
                            call jeveuo(jexnum(maxfem//'.CONNEX', ima), 'E', jconx)
                            zi(jconx) = ino
                            icompt = icompt+1
                            zi(iagma-1+ifon2) = ima
                        end if
                        zi(iagno-1+ifon2) = ino
                    end do
!
                end do
                ncompt = ncompt+nfon
            end if
!
        end do
!
!
    end if
!
!
end subroutine
