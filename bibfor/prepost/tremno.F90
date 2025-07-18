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

subroutine tremno(ncmp, nssche, nomsd)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/i2trgi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/numek8.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: ncmp
    character(len=19) :: nssche, nomsd
!
!   OPERATION REALISEE
!   ------------------
!
!     REORGANISATION PAR NOEUDS D' UN SOUS_CHAMELEM
!
!   ARGUMENTS
!   ---------
!
!     NOMSDL (IN) : NOM DE LA SD IMPLEMENTANT LA REORGANISATION
!
!     NCMP   (IN) : NOM DE LA COMPOSANTE TRAITEE
!
!     NSSCHE (IN) : NOM DE LA SD DE TYPE SOUS_CHAMPELEM A TRAITER
!
!     DESCRIPTION DE LA SD PRODUITES
!     ------------------------------
!
!         .VACP : XD V R8, UN OC CORRESPOND A LA LISTE DES VALEURS
!                 DE LA CMP SUR UN NOEUDS VU DES MAILLES LE CONTENANT
!
!         .NUMA : XD V I , UN OC CORRESPOND A LA LISTE DES MAILLES
!                 CONTENANT LE NOEUD DE L' OC DE VACP CORRESPONDANT
!
!         .NUND : S V I, VECTEUR DES NUMEROS DE NOEUDS CONCERNES
!
!         .NOCP : S E K8, NOM DE LA CMP
!
!         .NUCP : S E I, NUMERO DE LA CMP
!
!     L' OC NUMERO I CONTIENT LES VALEURS DE LA CMP SUR LE NOEUD
!     NUMERO NUND(I) (DANS LA NUMEROTATION DU MAILLAGE)
!
!*********************************************************************
!
!   FONCTIONS EXTERNES
!   ------------------
!
!
!   -------------------------
!
!
!   NOMS ET ADRESSES DES OJB ASSOCIES AU SOUS CHAMELEM
!   --------------------------------------------------
!
    character(len=24) :: npnbn, npadr, npcmp, nvale, nnoma, nnugd
    integer(kind=8) :: apnbn, apadr, apcmp, avale, anoma, anugd
!
!   NOMS ET ADRESSES DES OJB ASSOCIES A LA SD
!   -----------------------------------------
!
    character(len=24) :: nvacp, nnund, nnuma, nnucp, nnocp
    integer(kind=8) :: avacp, anund, anuma, anucp, anocp
    integer(kind=8) :: ptm, ptv, tco, tsp, nbco, nbsp, lngm, lngv, ico, isp
!
!   ADRESSE DE NUMERO DE CMP CONCERNEES PAR L' EXTRACTION
!   -----------------------------------------------------
!
    integer(kind=8) :: libre, nbnm, nbtcmp, numcp, gd, aconec, acmpgd
!
!   DIVERS
!   ------
!
    integer(kind=8) :: i, m, n, in, im, nbtnd, nbtmai, adrm, nbn, nbm, ndloc, nbcpac
    integer(kind=8) :: iconec, lconec, icncin, lcncin
    integer(kind=8) :: vali
    character(len=24) :: nconec, ncncin
    character(len=8) :: tk8(1), nmaila
    aster_logical :: trouve
    integer(kind=8), pointer :: entier(:) => null()
    integer(kind=8), pointer :: pnco(:) => null()
    integer(kind=8), pointer :: pnsp(:) => null()
!
!================= FIN DES DECLARATIONS ============================
!
!   RECUPERATION DES NOMS ET DES OJB DU SOUS_CHAMELEM
!   -------------------------------------------------
!
    call jemarq()
    npnbn = nssche//'.PNBN'
    npadr = nssche//'.PADR'
    npcmp = nssche//'.PCMP'
    nvale = nssche//'.VALE'
    nnoma = nssche//'.NOMA'
    nnugd = nssche//'.NUGD'
!
    call jeveuo(npnbn, 'L', apnbn)
    call jeveuo(npadr, 'L', apadr)
    call jeveuo(npcmp, 'L', apcmp)
    call jeveuo(nvale, 'L', avale)
    call jeveuo(nnoma, 'L', anoma)
    call jeveuo(nnugd, 'L', anugd)
    call jeveuo(nssche//'.PNCO', 'L', vi=pnco)
    call jeveuo(nssche//'.PNSP', 'L', vi=pnsp)
!
    nmaila = zk8(anoma)
!
    gd = zi(anugd)
!
    nconec = nmaila//'.CONNEX'
    ncncin = '&&OP0051.CONNECINVERSE'
!
!   CONSTRUCTION DES NOM DES OJB DE LA SD PRODUITE
!   ----------------------------------------------
!
    nvacp = nomsd//'.VACP'
    nnund = nomsd//'.NUND'
    nnocp = nomsd//'.NOCP'
    nnucp = nomsd//'.NUCP'
    nnuma = nomsd//'.NUMA'
!
!   CREATION REMPLISSAGE DE .NOCP
!   -----------------------------
!
    call wkvect(nnocp, 'V V K8', 1, anocp)
!
    zk8(anocp) = ncmp
!
!   CREATION REMPLISSAGE DE .NUCP
!   -----------------------------
!
    tk8(1) = ncmp
!
    call wkvect(nnucp, 'V V I', 1, anucp)
!
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', acmpgd)
    call jelira(npcmp, 'LONMAX', nbtcmp)
    call numek8(zk8(acmpgd), tk8, nbtcmp, 1, zi(anucp))
!
    numcp = zi(anucp)
!
!   RECUPERATION DU NBR TOTAL DE MAILLES
!   ------------------------------------
!
    call jelira(npnbn, 'LONMAX', nbtmai)
!
!   CALCUL DU NBR DE CMP ACTIVES
!   ----------------------------
!
    nbcpac = 0
!
    do i = 1, nbtcmp, 1
!
        nbcpac = nbcpac+min(zi(apcmp+i-1), 1)
!
    end do
!
!   CONSTRUCTION DU .NUND
!   ---------------------
!
    call dismoi('NB_NO_MAILLA', nmaila, 'MAILLAGE', repi=nbtnd)
!
    call jecreo('&&TREMNO.LISTE.ENTIER', 'V V I')
    call jeecra('&&TREMNO.LISTE.ENTIER', 'LONMAX', nbtnd)
    call jeveuo('&&TREMNO.LISTE.ENTIER', 'E', vi=entier)
    call jeveuo(nconec, 'L', iconec)
    call jeveuo(jexatr(nconec, 'LONCUM'), 'L', lconec)
    call jeveuo(ncncin, 'L', icncin)
    call jeveuo(jexatr(ncncin, 'LONCUM'), 'L', lcncin)
!
!
    libre = 1
!
    do im = 1, nbtmai, 1
!
        if (zi(apadr+im-1) .ne. 0) then
!
            adrm = iconec+zi(lconec-1+im)-1
!
            nbn = zi(apnbn+im-1)
!
            call i2trgi(entier, zi(adrm), nbn, libre)
!
        end if
!
    end do
!
    nbtnd = libre-1
!
    call jecreo(nnund, 'V V I')
    call jeecra(nnund, 'LONMAX', nbtnd)
    call jeveuo(nnund, 'E', anund)
!
    do in = 1, nbtnd, 1
!
        zi(anund+in-1) = entier(in)
!
    end do
!
    call jedetr('&&TREMNO.LISTE.ENTIER')
!
!   CONSTRUCTION DU .VACP
!   ---------------------
!
    call jecrec(nvacp, 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nbtnd)
    call jecrec(nnuma, 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbtnd)
!
    do in = 1, nbtnd, 1
!
        n = zi(anund+in-1)
!
!
        adrm = icncin+zi(lcncin-1+n)-1
        nbm = zi(lcncin+n)-zi(lcncin-1+n)
!
        lngm = 0
        lngv = 0
!
        do im = 1, nbm, 1
!
            m = zi(adrm+im-1)
!
            if (m .ne. 0) then
!
                lngm = lngm+min(zi(apadr+m-1), 1)
                lngv = lngv+min(zi(apadr+m-1), 1)*pnco(m)*pnsp(m)
!
            end if
!
        end do
!
        if (lngm .eq. 0) then
!
            vali = n
            call utmess('F', 'PREPOST5_76', si=vali)
!
        end if
!
        call jecroc(jexnum(nvacp, in))
        call jeecra(jexnum(nvacp, in), 'LONMAX', lngv)
        call jeveuo(jexnum(nvacp, in), 'E', avacp)
!
        call jecroc(jexnum(nnuma, in))
        call jeecra(jexnum(nnuma, in), 'LONMAX', lngm)
        call jeveuo(jexnum(nnuma, in), 'E', anuma)
!
        ptm = 1
        ptv = 1
!
        do im = 1, nbm, 1
!
            m = zi(adrm+im-1)
!
            if (m .ne. 0) then
!
                if (zi(apadr+m-1) .ne. 0) then
!
                    nbco = pnco(m)
                    nbsp = pnsp(m)
                    ndloc = 1
                    trouve = .false.
!
                    aconec = iconec+zi(lconec-1+m)-1
                    nbnm = zi(lconec+m)-zi(lconec-1+m)
!
!
220                 continue
                    if ((.not. trouve) .and. (ndloc .le. nbnm)) then
!
                        if (zi(aconec+ndloc-1) .eq. n) then
!
                            trouve = .true.
!
                        else
!
                            ndloc = ndloc+1
!
                        end if
!
                        goto 220
!
                    end if
!
                    tsp = nbcpac*nbsp
!
                    do ico = 1, nbco, 1
!
                        tco = nbsp*nbcpac*nbnm*(ico-1)
!
                        do isp = 1, nbsp, 1
!
                            zr(avacp+ptv-1) = zr( &
                                                avale+zi(apadr+m-1)+tco+(ndloc-1)*tsp+(i&
                                                &sp-1)*nbcpac+zi(apcmp+numcp-1)-2 &
                                                )
!
                            ptv = ptv+1
!
                        end do
!
                    end do
!
!
                    zi(anuma+ptm-1) = m
!
                    ptm = ptm+1
!
                end if
!
            end if
!
        end do
!
    end do
!
    call jedema()
end subroutine
