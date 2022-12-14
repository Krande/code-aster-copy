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
!
subroutine asecon(nomsy, neq, mome, resu)
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/r8vide.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/codent.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsexis.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/vtdefs.h"
#include "asterfort/wkvect.h"
!
    integer :: neq
    character(len=16) :: nomsy
    character(len=*) :: mome, resu
!     COMMANDE : COMB_SISM_MODAL
!        CALCUL DES TERMES D'ENTRAINEMENT
!     ------------------------------------------------------------------
! IN  : NOMSY  : OPTION DE CALCUL
! IN  : NEQ    : NOMBRE D'EQUATIONS
! IN  : NOME   : MODES MECANIQUES
! IN  : RESU   : NOM UTILISATEUR DE LA COMMANDE
!     ------------------------------------------------------------------
    integer :: iad, ibid, icas, idep, idir, ier, ii, in, ino, ioc, iocc, iordr
    integer :: iorst, iret, jcas, jdir, jno, jord
    integer :: jvale, jval1, lnod, nbmode, nbno, nboc, nbtrou
    integer :: ncas, ndep, nucas, nume, tordr(1)
    real(kind=8) :: r8b, epsmac, xxx, xx1, xx2, xx3
    complex(kind=8) :: cbid
    character(len=8) :: k8b, noeu, cmp, nomcmp(3), knum, kdir, stat
    character(len=8) :: meca, occur
    character(len=16) :: monacc, concep, nomcmd, def
    character(len=19) :: chextr, champ, moncha
    character(len=24) :: vale, noms2, valk(3)
    real(kind=8), pointer :: vabs(:) => null()
    real(kind=8), pointer :: aux(:) => null()
    real(kind=8), pointer :: cumul(:) => null()
    real(kind=8), pointer :: line(:) => null()
    real(kind=8), pointer :: quad(:) => null()
    real(kind=8), pointer :: rep(:) => null()
    integer, pointer :: type(:) => null()
    character(len=8), pointer :: vstat(:) => null()
!     ------------------------------------------------------------------
    data  nomcmp / 'DX' , 'DY' , 'DZ' /
    data  vale / '                   .VALE' /
!     ------------------------------------------------------------------
!
    call jemarq()
    call getfac('COMB_DEPL_APPUI', nboc)
    call getfac('DEPL_MULT_APPUI', ndep)
    call getvid(' ', 'MODE_MECA', scal=meca, nbret=ibid)
!
    AS_ALLOCATE(vr=cumul, size=neq)
    AS_ALLOCATE(vr=aux, size=neq*nboc)
!
    epsmac = r8vide()
    nbmode = 10
    ii = 0
!
! -- PREPARATION DU STOCKAGE DE LA REPONSE SECONDAIRE
!
    call getres(k8b, concep, nomcmd)
!     --- CREATION DE LA STRUCTURE D'ACCUEIL ---
    call rsexis(resu, ier)
    if (ier .eq. 0) call rscrsd('G', resu, concep, nbmode)
    noms2 = nomsy
    if (nomsy(1:4) .eq. 'VITE') noms2 = 'DEPL'
!
    call rsorac(mome, 'TOUT_ORDRE', ibid, r8b, k8b,&
                cbid, r8b, k8b, tordr, 1,&
                nbtrou)
    iordr=tordr(1)
!
    call rsexch('F', meca, noms2, iordr, moncha,&
                ier)
    def = 'SECONDAIRE'
    iordr = 200
!           --- CHAMP PAR OCCURENCE DE COMB_DPL_APPUI ---
!
    call jeveuo('&&ASENAP.TYPE', 'L', vi=type)
    AS_ALLOCATE(vr=rep, size=neq)
    AS_ALLOCATE(vr=quad, size=neq)
    AS_ALLOCATE(vr=line, size=neq)
    AS_ALLOCATE(vr=vabs, size=neq)
    call jeexin('&&ASECON.NORD', iret)
    if (iret .eq. 0) then
        call wkvect('&&ASECON.NORD', 'V V I', nboc+1, jord)
    else
        call jeveuo('&&ASECON.NORD', 'E', jord)
    endif
!
    do iocc = 1, nboc
!
! POUR CHAQUE OCCURENCE ON STOQUE LE CHAMP
!
        call rsexch(' ', resu, nomsy, iordr, champ,&
                    ier)
        if (ier .eq. 100) then
            call vtdefs(champ, moncha, 'G', 'R')
        else
            valk (1) = nomsy
            valk (2) = champ
            call utmess('F', 'SEISME_25', nk=2, valk=valk, si=iocc)
        endif
        vale(1:19) = champ
        call jeexin(vale(1:19)//'.VALE', ibid)
        if (ibid .gt. 0) then
            vale(20:24)='.VALE'
        else
            vale(20:24)='.CELV'
        endif
        call jeveuo(vale, 'E', jvale)
!
        do in = 1, neq
            quad(in)= 0.0d0
            line(in)= 0.0d0
            vabs(in)= 0.0d0
        end do
        call jelira(jexnum('&&ASENAP.LISTCAS', iocc), 'LONMAX', ncas)
        call jeveuo(jexnum('&&ASENAP.LISTCAS', iocc), 'L', jcas)
        do icas = 1, ncas
            nucas = zi(jcas+icas-1)
            do idep = 1, ndep
                call getvis('DEPL_MULT_APPUI', 'NUME_CAS', iocc=idep, scal=nume, nbret=ibid)
                if (nume .eq. nucas) then
                    knum = 'N       '
                    call codent(nucas, 'D0', knum(2:8))
                    kdir = 'D       '
                    call codent(nucas, 'D0', kdir(2:8))
                    call jelira(jexnom('&&ASENAP.LINOEU', knum), 'LONMAX', nbno)
                    call jeveuo(jexnom('&&ASENAP.LINOEU', knum), 'L', jno)
                    lnod = 3*nbno
                    call jelira(jexnom('&&ASENAP.LIDIR', kdir), 'LONMAX', lnod)
                    call jeveuo(jexnom('&&ASENAP.LIDIR', kdir), 'L', jdir)
                    call jeveuo('&&ASENAP.STAT', 'L', vk8=vstat)
                    stat = vstat(icas)
                    do ino = 1, nbno
                        noeu =zk8(jno+ino-1)
                        do idir = 1, 3
                            if (zr(jdir+3*(ino-1)+idir-1) .ne. epsmac) then
                                cmp = nomcmp(idir)
                                monacc = noeu//cmp
                                xx1 = zr(jdir+3*(ino-1)+idir-1)
                                call rsorac(stat, 'NOEUD_CMP', ibid, r8b, monacc,&
                                            cbid, r8b, k8b, tordr, 1,&
                                            nbtrou)
                                iorst=tordr(1)
                                call rsexch('F', stat, nomsy, iorst, chextr,&
                                            iret)
                                call jeexin(chextr//'.VALE', ibid)
                                if (ibid .gt. 0) then
                                    call jeveuo(chextr//'.VALE', 'L', jval1)
                                else
                                    call jeveuo(chextr//'.CELV', 'L', jval1)
                                endif
                                do in = 1, neq
                                    rep(in) = zr(jval1+in-1) * xx1
                                end do
                                if (type(iocc) .eq. 1) then
!                 --- COMBINAISON QUADRATIQUE ---
                                    do in = 1, neq
                                        xxx = rep(in)
                                        quad(in)= quad(in)+&
                                        xxx*xxx
                                    end do
                                else if (type(iocc).eq.2) then
!               --- COMBINAISON LINEAIRE ---
                                    do in = 1, neq
                                        line(in)= line(in)+&
                                        rep(in)
                                    end do
                                else
!              --- COMBINAISON VALEUR ABSOLUE ---
                                    do in = 1, neq
                                        xx1 = abs(rep(in))
                                        vabs(in)= vabs(in)+&
                                        xx1
                                    end do
                                endif
                            endif
                        end do
                    end do
                endif
            end do
        end do
        do in = 1, neq
            xx1 = line(in)
            xx2 = vabs(in)
            xx3 = sqrt(quad(in))
            zr(jvale+in-1) = xx1 + xx2 + xx3
            ii = ii + 1
            aux(ii) = zr(jvale+in-1)
        end do
!
        call rsnoch(resu, nomsy, iordr)
        call rsadpa(resu, 'E', 1, 'NOEUD_CMP', iordr,&
                    0, sjv=iad, styp=k8b)
        call codent(iocc, 'D', occur)
        zk16(iad) = 'COMBI'// occur
        call rsadpa(resu, 'E', 1, 'TYPE_DEFO', iordr,&
                    0, sjv=iad, styp=k8b)
        zk16(iad) = def
        call jelibe(vale)
!
        zi(jord+iocc-1) = iordr
        iordr = iordr + 1
    end do
    zi(jord+nboc) = iordr
!
    call rsexch(' ', resu, nomsy, iordr, champ,&
                ier)
    if (ier .eq. 100) then
        call vtdefs(champ, moncha, 'G', 'R')
    else
        valk(1) = nomsy
        valk(2) = champ
        call utmess('F', 'SEISME_25', nk=2, valk=valk, si=iordr)
    endif
    vale(1:19) = champ
    call jeexin(vale(1:19)//'.VALE', ibid)
    if (ibid .gt. 0) then
        vale(20:24)='.VALE'
    else
        vale(20:24)='.CELV'
    endif
    call jeveuo(vale, 'E', jvale)
!
    do ioc = 1, nboc
        do in = 1, neq
            xx1 = aux(1+(ioc-1)*neq+in-1)
            cumul(in) = cumul(in)+xx1*xx1
        end do
    end do
! STOCKAGE DU CUMUL QUADRATIQUE
    do in = 1, neq
        zr(jvale+in-1) = sqrt( abs ( cumul(in) ) )
    end do
    call jelibe(vale)
    call rsnoch(resu, nomsy, iordr)
!
!        --- PARAMETRE ---
    call rsadpa(resu, 'E', 1, 'NOEUD_CMP', iordr,&
                0, sjv=iad, styp=k8b)
    zk16(iad) = 'CUMUL'//' QUAD'
    call rsadpa(resu, 'E', 1, 'TYPE_DEFO', iordr,&
                0, sjv=iad, styp=k8b)
    zk16(iad) = def
!
    AS_DEALLOCATE(vr=cumul)
    AS_DEALLOCATE(vr=aux)
    AS_DEALLOCATE(vr=rep)
    AS_DEALLOCATE(vr=quad)
    AS_DEALLOCATE(vr=line)
    AS_DEALLOCATE(vr=vabs)
    call jedema()
end subroutine
