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
!
subroutine op0196()
    implicit none
! person_in_charge: samuel.geniaut at edf.fr
!     =================================================================
!                      OPERATEUR POST_CHAM_XFEM
!                      ------------------------
!     BUT : DETERMNATION DES CHAMPS DE DEPLACEMENTS, DE CONTRAINTES
!           ET DE VARIABLES INTERNES SUR LE MAILLAGE FISSURE X-FEM
!     =================================================================
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cescar.h"
#include "asterfort/cescel.h"
#include "asterfort/cnscno.h"
#include "asterfort/detrsd.h"
#include "asterfort/gettco.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/xpodim.h"
#include "asterfort/xpofon.h"
#include "asterfort/xpoini.h"
#include "asterfort/xpomac.h"
#include "asterfort/xpomax.h"
#include "asterfort/xposep.h"
#include "asterfort/xxishm.h"
    integer(kind=8) :: jlicha, nbordr, ior, jord, iord, jinst1, jinst2, nbcham
    integer(kind=8) :: ibid, iret, nsetot, nnntot, ncotot, nbnoc, nbmac, ifm, niv, ic
    integer(kind=8) :: jmod, mftot, nftot, nfcomf, ngfon
    character(len=1) :: kbid
    character(len=2) :: k2b(4)
    character(len=8) :: maxfem, mo, malini, resuco, resux, modvis
    character(len=16) :: tysd, nomcha
    character(len=19) :: cns1, cns2, ces1, ces2, cel2, ch, cesvi1, cesvi2, k19
    character(len=19) :: comps1, comps2
    character(len=24) :: mailx, mailc, licham, ordr, listno, logrma, k24, listgr, k24b
    aster_logical :: pre1
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
!
!     ------------------------------------------------------------------
!     1. RECUPERATION DES CONCEPTS UTILISATEURS
!     ------------------------------------------------------------------
!
    if (niv .gt. 1) write (ifm, *) ' '
    if (niv .gt. 1) write (ifm, *) '1. XPOINI'
    licham = '&&OP0196.LICHAM'
    call xpoini(maxfem, mo, malini, modvis, licham, &
                resuco, resux, k2b, k24b)
    call xpofon(mo, mftot, nftot, nfcomf, ngfon)
!
!     ------------------------------------------------------------------
!     2. SEPARATION DES MAILLES DE MALINI EN 2 GROUPES
!              - MAILC : MAILLES NON AFFECTEES D'UN MODELE
!                        OU NON SOUS-DECOUPEES (CLASSIQUE)
!              - MAILX : MAILLES SOUS-DECOUPEES (X-FEM)
!     ------------------------------------------------------------------
!
    if (niv .gt. 1) write (ifm, *) ' '
    if (niv .gt. 1) write (ifm, *) '2. XPOSEP'
    mailc = '&&OP0196.MAILC'
    mailx = '&&OP0196.MAILX'
    logrma = '&&OP0196.LOGRMA'
    listgr = '&&OP0196.LISTGR'
    call xposep(mo, malini, mailc, mailx, nsetot, &
                nnntot, ncotot, logrma, listgr)
!
!     PRE1=.FALSE. :
!     ON NE POST TRAITE PAS SUR DES ELEMENTS HM-XFEM
!     PRE1=.TRUE. :
!     ON POST TRAITE SUR DES ELEMENTS HM-XFEM --> ATTENTION AUX DDLs DE
!     PRESSION
!
    pre1 = xxishm(mailc, mailx, mo)
!
!     CREATION DE LA NOUVELLE SD RESULTAT
    ordr = resuco//'           .ORDR'
    call jeveuo(ordr, 'L', jord)
    call jelira(ordr, 'LONUTI', nbordr)
    call gettco(resuco, tysd)
    call rscrsd('G', resux, tysd, nbordr)
!
!     BOUCLE SUR LES NBORDR NUMEROS D'ORDRE
    do ior = 1, nbordr
!
        iord = zi(jord-1+ior)
!       ----------------------------------------------------------------
!       3. DIMENSIONNEMENT DES OBJETS DU RESU X-FEM
!       ----------------------------------------------------------------
!
        if (niv .gt. 1) write (ifm, *) ' '
        if (niv .gt. 1) write (ifm, *) '3. XPODIM'
        cns1 = '&&OP0196.CNS1'
        ces1 = '&&OP0196.CES1'
        cesvi1 = '&&OP0196.CESVI1'
        cns2 = '&&OP0196.CNS2'
        ces2 = '&&OP0196.CES2'
        cesvi2 = '&&OP0196.CESVI2'
        cel2 = '&&OP0196.CEL2'
        listno = '&&OP0196.LISTNO'
        comps1 = '&&OP0196.COMPOR1'
        comps2 = '&&OP0196.COMPOR2'
        call xpodim(malini, mailc, modvis, licham, nsetot+mftot, &
                    nnntot+nftot, ncotot+nfcomf, listno, cns1, cns2, &
                    ces1, ces2, cel2, cesvi1, cesvi2, &
                    ior, resuco, nbnoc, nbmac, logrma, &
                    k24, maxfem, ibid, comps1, comps2, &
                    pre1, mo)
!
!       ----------------------------------------------------------------
!       4. TRAITEMENT DES MAILLES DE MAILC
!       ----------------------------------------------------------------
!
        if (niv .gt. 1) write (ifm, *) ' '
        if (niv .gt. 1) write (ifm, *) '4. XPOMAC'
        call xpomac(malini, mailc, listno, nbnoc, nbmac, &
                    maxfem, k24, cns1, cns2, ces1, &
                    ces2, cesvi1, cesvi2, resuco, comps1, &
                    comps2, pre1)
!
!       ----------------------------------------------------------------
!       5. TRAITEMENT DES MAILLES DE MAILX
!       ----------------------------------------------------------------
!
        if (niv .gt. 1) write (ifm, *) ' '
        if (niv .gt. 1) write (ifm, *) '5. XPOMAX'
        call xpomax(mo, malini, mailx, nbnoc, nbmac, &
                    k2b, k24b, maxfem, cns1, cns2, &
                    ces1, ces2, cesvi1, cesvi2, listgr, &
                    k24, k24, resuco, ibid, comps1, &
                    comps2, pre1, iord)
!
!       ----------------------------------------------------------------
!       6. ENREGISTREMENT DES CHAMPS DE SORTIES
!       ----------------------------------------------------------------
!
        call jelira(licham, 'LONMAX', nbcham)
        call jeveuo(licham, 'L', jlicha)
        do ic = 1, nbcham
            nomcha = zk16(jlicha-1+ic)
!
!         SI LE CHAMP N'EXISTE PAS DANS LE RESU EN ENTREE, ON PASSE
            call rsexch(' ', resuco, nomcha, iord, k19, &
                        iret)
            if (iret .gt. 0) goto 20
!
            if (niv .gt. 1) write (ifm, *) '6. ENREGISTREMENT DE ', nomcha
!
!         RECUPERATION DU NOM DU CHAMP A ECRIRE : CH
            call rsexch(' ', resux, nomcha, iord, ch, &
                        iret)
!
            if (nomcha .eq. 'DEPL' .or. nomcha .eq. 'TEMP') then
!
                call cnscno(cns2, ' ', 'NON', 'G', ch, &
                            'F', ibid)
!
            else if (nomcha .eq. 'SIEF_ELGA') then
!
                call cescel(ces2, modvis//'.MODELE', 'FULL_MECA', 'PCONTMR', 'NON', &
                            ibid, 'G', ch, 'F', ibid)
!
            else if (nomcha .eq. 'VARI_ELGA') then
!
                call cescel(cesvi2, modvis//'.MODELE', 'FULL_MECA', 'PVARIMR', 'NON', &
                            ibid, 'G', ch, 'F', ibid)
!
            end if
            call rsnoch(resux, nomcha, iord)
20          continue
        end do
!
!       CARTE DU COMPORTEMENT
        call rsexch(' ', resuco, 'COMPORTEMENT', iord, k19, &
                    iret)
        if (iret .eq. 0) then
!         RECUPERATION DU NOM DU CHAMP A ECRIRE : CH
            call rsexch(' ', resux, 'COMPORTEMENT', iord, ch, &
                        iret)
            call cescar(comps2, ch, 'G')
            call rsnoch(resux, 'COMPORTEMENT', iord)
        end if
!
        if (tysd(1:4) .eq. 'EVOL') then
            call rsadpa(resuco, 'L', 1, 'INST', iord, &
                        0, sjv=jinst1, styp=kbid)
            call rsadpa(resux, 'E', 1, 'INST', iord, &
                        0, sjv=jinst2, styp=kbid)
            zr(jinst2) = zr(jinst1)
            call rsadpa(resux, 'E', 1, 'MODELE', iord, &
                        0, sjv=jmod, styp=kbid)
            zk8(jmod) = modvis
        else if (tysd(1:9) .eq. 'MODE_MECA') then
            call rsadpa(resuco, 'L', 1, 'FREQ', iord, &
                        0, sjv=jinst1, styp=kbid)
            call rsadpa(resux, 'E', 1, 'FREQ', iord, &
                        0, sjv=jinst2, styp=kbid)
            zr(jinst2) = zr(jinst1)
!
            call rsadpa(resux, 'E', 1, 'MODELE', iord, &
                        0, sjv=jmod, styp=kbid)
            zk8(jmod) = mo
!
            call refdcp(resuco, resux)
        end if
!
        call detrsd('CHAM_NO_S', cns1)
        call detrsd('CHAM_NO_S', cns2)
        call detrsd('CHAM_ELEM_S', ces1)
        call detrsd('CHAM_ELEM_S', ces2)
        call detrsd('CHAM_ELEM', cel2)
        call detrsd('CHAM_ELEM_S', cesvi1)
        call detrsd('CHAM_ELEM_S', cesvi2)
        call detrsd('CHAM_ELEM_S', comps1)
        call detrsd('CHAM_ELEM_S', comps2)
!
    end do
!
    call jeexin(mailc, iret)
    if (iret .ne. 0) call jedetr(mailc)
    call jeexin(mailx, iret)
    if (iret .ne. 0) call jedetr(mailx)
    call jedetr(licham)
    call jeexin(listgr, iret)
    if (iret .ne. 0) call jedetr(listgr)
!
    if (niv .gt. 1) write (ifm, *) 'FIN DE POST_CHAM_XFEM'
!
!
    call jedema()
!
end subroutine
