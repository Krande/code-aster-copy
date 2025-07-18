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

subroutine rvpost(mcf, iocc, dim, i1, i2, &
                  ncheff, xnomcp, nresu, nch19, nlsmac, &
                  nlsnac, nomtab, xnovar)
!     PILOTAGE DU POST-TRAITEMENT
!     ------------------------------------------------------------------
! IN  IOCC   : I : INDICE DE L' OCCURENCE
! IN  DIM    : K : '2D' OU '3D'
! IN  I1, I2 : I : REPERAGE DU CHAMP DANS UNE SD RESULTAT_COMPOSE
! IN  XNOMCP : K : NOM DE LA COLLECTION DES NOMS DE CMP
! IN  NCH19  : K : NOM DU CHAMP A TRAITER
! IN  NLSMAC : K : NOM DU VECTEUR DES MAILLES ACTIVES
! IN  NLSNAC : K : NOM DU VECTEUR DES NOEUDS ACTIFS
!     ------------------------------------------------------------------
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/extche.h"
#include "asterfort/extchn.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rvaffe.h"
#include "asterfort/rvaffm.h"
#include "asterfort/rvaffs.h"
#include "asterfort/rvcalq.h"
#include "asterfort/rvchgr.h"
#include "asterfort/rvcpnc.h"
#include "asterfort/rvinfo.h"
#include "asterfort/rvlieu.h"
#include "asterfort/rvpste.h"
#include "asterfort/rvpstm.h"
#include "asterfort/rvpsts.h"
#include "asterfort/tuesch.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: iocc, i1, i2
    character(len=2) :: dim
    character(len=8) :: nresu
    character(len=16) :: ncheff
    character(len=19) :: nch19, nomtab
    character(len=24) :: xnomcp, nlsmac, nlsnac, xnovar
    character(len=*) :: mcf
!
    integer(kind=8) :: gd, i, idir, niv, iret, isd, jcmpcd, jcmpnc, jdir, jlsmac
    integer(kind=8) :: jlsnac, jnomcp, jsdev, jsdli, jvec1, jvec2, n, n0, n1, n2, n3
    integer(kind=8) :: nbcac, nbcpn, nbmac, nbnac, nboper, nbsd, nr, ifm, ibid
    integer(kind=8) :: ny
    real(kind=8) :: vecty(3)
    aster_logical :: tridim
    character(len=24) :: lscpnc, quant, sdlieu, sdeval, lscpcd
    character(len=24) :: sdev, sdli, sdmoye, sdmail
    character(len=19) :: sdpost, eval, lieu, sdnewr, ssch19
    character(len=16) :: option, oper, operat(2)
    character(len=8) :: k8b, typco, courbe, mailla, repere
    character(len=4) :: docu
    character(len=1) :: ca
    aster_logical :: chok
!
!==================== CORPS DE LA ROUTINE =============================
!
    call jemarq()
    chok = .true.
    call infniv(ifm, niv)
!
    call getvtx(mcf, 'OPERATION', iocc=iocc, nbval=0, nbret=nboper)
    nboper = -nboper
    call getvtx(mcf, 'OPERATION', iocc=iocc, nbval=nboper, vect=operat, &
                nbret=n0)
!
    if (nch19(1:1) .eq. '&') then
        if (niv .gt. 1) call rvinfo(ifm, iocc, i1, i2, 'E', &
                                    ncheff)
    else
        lscpcd = '&&RVPOST.NOM.CMP.CAND.OC'
        lscpnc = '&&RVPOST.NOM.CMP.NCSR.OC'
        call dismoi("DOCU", nch19, "CHAMP", repk=docu)
        call dismoi("NUM_GD", nch19, "CHAMP", repi=gd)
        call jeveuo(jexnum(xnomcp, iocc), 'L', jnomcp)
        call jelira(jexnum(xnomcp, iocc), 'LONMAX', nbcac)
        call wkvect(lscpcd, 'V V K8', nbcac, jcmpcd)
        call wkvect('&&RVPOST.VAL.DIR', 'V V R', 3, jdir)
        do i = 1, nbcac, 1
            zk8(jcmpcd+i-1) = zk8(jnomcp+i-1)
        end do
        if (niv .gt. 1) call rvinfo(ifm, iocc, i1, i2, 'B', &
                                    ncheff)
        if (nresu(1:1) .eq. ' ') nresu = nch19(1:8)
        call rvcpnc(mcf, iocc, nch19, gd, docu, &
                    nbcac, lscpcd, lscpnc, repere, option, &
                    quant, idir, zr(jdir), iret)
!        /* POSSIBILITE AGRANDISSEMENT DE LSCPCD => ON REFAIT JEVEUO */
        call jeveuo(lscpcd, 'L', jcmpcd)
!
        if (iret .ne. 0) then
            ssch19 = '&&RVPOST.SOUS.CH.GD'
            call jelira(lscpnc, 'LONMAX', nbcpn)
            call jeveuo(lscpnc, 'L', jcmpnc)
!
            if (docu .eq. 'CHNO') then
                call jelira(nlsnac, 'LONMAX', nbnac)
                call jeveuo(nlsnac, 'L', jlsnac)
                call extchn(nch19, k8b, zi(jlsnac), zk8(jcmpnc), nbnac, &
                            nbcpn, 'NUMERO', ssch19, mcf, iocc)
!
            else
                call jeexin(nlsnac, ibid)
                if (ibid .gt. 0) then
                    call jelira(nlsnac, 'LONMAX', nbnac)
                    call jeveuo(nlsnac, 'L', jlsnac)
                else
                    jlsnac = 1
                    nbnac = 0
                end if
                call jelira(nlsmac, 'LONMAX', nbmac)
                call jeveuo(nlsmac, 'L', jlsmac)
                call extche(nch19, k8b, zi(jlsmac), zk8(jcmpnc), nbmac, &
                            nbcpn, 'NUMERO', ssch19, mcf, iocc, &
                            nbnac, zi(jlsnac))
            end if
!
            call getvr8('ACTION', 'VECT_Y', iocc=iocc, nbval=3, vect=vecty, &
                        nbret=ny)
            tridim = ny .ne. 0
!
            if (chok) then
                call dismoi('NOM_MAILLA', nch19, 'CHAMP', repk=mailla)
                typco = 'AUTRE'
                courbe = '&&YAPAS'
                sdlieu = '&&RVPOST.NOM.VECT.LIEU'
                sdeval = '&&RVPOST.NOM.VECT.EVAL'
!
                call getvtx(mcf, 'MOYE_NOEUD', iocc=iocc, scal=k8b, nbret=n)
                if (k8b(1:1) .eq. 'O') then
                    ca = 'N'
                else
                    ca = 'E'
                end if
!
                call rvlieu(mailla, typco, nlsnac, sdlieu)
                call rvpste(sdlieu, ssch19, sdeval, ca)
                call jelira(sdlieu, 'LONMAX', nbsd)
                call jeveuo(sdlieu, 'L', jsdli)
                call jeveuo(sdeval, 'L', jsdev)
                call getvtx(mcf, 'RESULTANTE', iocc=iocc, nbval=0, nbret=nr)
                sdnewr = '&&RVPOST.NEW.REPERE'
                if (repere(1:1) .ne. 'G' .and. .not. tridim) then
                    call rvchgr(mailla, nlsnac, repere, sdnewr, &
                                iret)
                else
                    iret = 1
                end if
!
                if (iret .ne. 0) then
                    sdpost = '&&RVPOST.FINAL.POST'
                    do isd = 1, nbsd, 1
                        if (repere(1:1) .ne. 'G' .and. .not. tridim) then
                            call jeveuo(jexnum(sdnewr//'.VEC1', isd), 'L', jvec1)
                            call jeveuo(jexnum(sdnewr//'.VEC2', isd), 'L', jvec2)
                        else
                            jvec1 = 0
                            jvec2 = 0
                        end if
!
                        sdev = zk24(jsdev+isd-1)
                        sdli = zk24(jsdli+isd-1)
!
                        call rvcalq(iocc, sdev, zr(jvec1), zr(jvec2), repere, &
                                    zk8(jcmpcd), nbcpn, nbcac, option, quant, &
                                    sdli, idir, zr(jdir), sdpost, courbe)
!
                        if (nr .eq. 0) then
                            if (nboper .eq. 2) then
                                sdmail = sdev(1:19)//'.MAIL'
                                sdmoye = '&&RVPOST.MOYENNE'
                                call rvpstm(sdli, sdpost, sdmoye)
                                call rvaffe(mcf, iocc, sdli, sdpost, sdmail, &
                                            ca, quant, option, repere, nomtab, &
                                            xnovar, ncheff, i1, isd)
                                oper = 'MOYENNE'
                                call rvaffm(mcf, iocc, sdli, sdpost, sdmoye, &
                                            oper, quant, option, repere, nomtab, &
                                            xnovar, ncheff, i1, isd)
                                call jedetr(sdmoye)
!
                            else
                                oper = operat(1)
                                if (oper .eq. 'EXTRACTION') then
                                    sdmail = sdev(1:19)//'.MAIL'
                                    call rvaffe(mcf, iocc, sdli, sdpost, sdmail, &
                                                ca, quant, option, repere, nomtab, &
                                                xnovar, ncheff, i1, isd)
!
                                else
                                    sdmoye = '&&RVPOST.MOYENNE'
                                    call rvpstm(sdli, sdpost, sdmoye)
                                    call rvaffm(mcf, iocc, sdli, sdpost, sdmoye, &
                                                oper, quant, option, repere, nomtab, &
                                                xnovar, ncheff, i1, isd)
                                    call jedetr(sdmoye)
                                end if
                            end if
                        else
                            sdmoye = '&&RVPOST.SOMME'
                            call rvpsts(iocc, sdli, sdpost, sdmoye)
                            call rvaffs(mcf, iocc, sdli, sdpost, sdmoye, &
                                        quant, option, repere, nomtab, ncheff, &
                                        i1, isd)
                            sdmoye(20:24) = '.VALE'
                            call jedetr(sdmoye)
                            sdmoye(20:24) = '.NOCP'
                            call jedetr(sdmoye)
                        end if
!
                        call jedetr(sdpost//'.VALE')
                        call jedetr(sdpost//'.PADR')
                        call jedetr(sdpost//'.PNBN')
                        call jedetr(sdpost//'.NOCP')
                        call jedetr(sdpost//'.PNCO')
                        call jedetr(sdpost//'.PNSP')
                    end do
!
                end if
                call jeexin(sdnewr//'.VEC1', n1)
                if (n1 .ne. 0) then
                    call jedetr(sdnewr//'.VEC1')
                    call jedetr(sdnewr//'.VEC2')
                end if
                call jelira(sdlieu, 'LONMAX', n)
                call jeveuo(sdlieu, 'L', n1)
                call jeveuo(sdeval, 'L', n2)
                do i = 1, n, 1
                    lieu = zk24(n1+i-1) (1:19)
                    eval = zk24(n2+i-1) (1:19)
                    call jedetr(lieu//'.ABSC')
                    call jedetr(lieu//'.REFE')
                    call jedetr(lieu//'.DESC')
                    call jedetr(lieu//'.NUME')
                    call jedetr(lieu//'.COOR')
                    call tuesch(eval)
                    call jeexin(eval//'.MAIL', n3)
                    if (n3 .ne. 0) then
                        call jedetr(eval//'.MAIL')
                    end if
                end do
                call jedetr(sdlieu)
                call jedetr(sdeval)
            end if
            call jedetr(lscpnc)
            call tuesch(ssch19)
        end if
        call jedetr(lscpcd)
        call jedetr('&&RVPOST.VAL.DIR')
    end if
!
    call jedema()
end subroutine
