! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine rc3600_momeq()

! ======================================================================
!     COMMANDE :  POST_RCCM / B3600 / OPTION MOMENT_EQUIVALENT
! ----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/chmima.h"
#include "asterfort/dismoi.h"
#include "asterfort/detrsd.h"
#include "asterfort/exlima.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rc3600_chtotab.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rslesd.h"
#include "asterfort/rsnopa.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/utmess.h"
!#include "asterfort/imprsd.h"
    character(len=8) :: resu, champ

    integer :: ifm, niv
    integer :: n0, n1, iret
    integer :: nbordr, i, iadin, iadou, j, jordr, nuordr
    
    real(kind=8) :: prec
    character(len=3) :: type
    character(len=16) :: crit, typesd, nomsym
    character(len=19) ::resuTmp19
    
    integer ::ibid, nbac, nbpa, nbpara, jnompa, nbma
    character(len=8) :: modele, resuTmp, noma, nomtab, conceptin
    character(len=8) :: modeav, lpain(1), lpaout(1), typmax, typresu
    character(len=4) :: tsca
    character(len=24) :: nompar
    character(len=16) :: nomcmd, concep, nopara
    character(len=19) :: chin, chextr, ligrel, resu19, lchin(1), lchout(1)
    character(len=19) :: noch19, tychlu, mcf
    integer :: vali(2), iexi
    aster_logical :: lreel, lnoeu, ldetli
    

    
!     ------------------------------------------------------------------
!
    call jemarq()
!
!
    call infmaj()
    call infniv(ifm, niv)
    
    call getres(nomtab, concep, nomcmd)
    
    mcf = 'RESU_MECA'
    call getvid(mcf, 'RESULTAT', iocc=1, scal=resu, nbret=n0)
    call getvid(mcf, 'CHAM_GD', iocc=1,  scal=champ, nbret=n1)
    
    if (n0.ge.1) then
        conceptin = resu
        resu19=resu
        call dismoi('TYPE_RESU', resu, 'RESULTAT', repk=typesd)
        call getvid(mcf, 'NOM_CHAM', iocc=1, scal=nomsym, nbret=ibid)
        lnoeu=nomsym.eq.'EFGE_ELNO'
!
!        -- SELECTION DES NUMERO D'ORDRE :
!        ---------------------------------
        prec=-1.d0
        crit=' '
        call getvr8(mcf, 'PRECISION', iocc=1, scal=prec, nbret=iret)
        call getvtx(mcf, 'CRITERE', iocc=1, scal=crit, nbret=iret)
        call rsutnu(resu19, mcf, 1, '&&RC3600.NUME_ORDRE', nbordr,&
                    prec, crit, iret)
        ASSERT(iret.eq.0)
        ASSERT(nbordr.gt.0)
        call jeveuo('&&RC3600.NUME_ORDRE', 'L', jordr)
!
!
!        -- 1. ON CREE LA SD_RESULTAT TEMPORAIRE :
!        ---------------------------------------------
        resuTmp = '&&RC3600'
        call rscrsd('V', resuTmp, typesd, nbordr)
!
!        -- 2. : BOUCLE SUR LES NUMERO D ORDRE
!        --------------------------------------------------
        
        modeav=' '
        ldetli=.false.
        do  i = 1, nbordr
            nuordr=zi(jordr-1+i)
            call rsexch(' ', resu19, nomsym, nuordr, chin,&
                        iret)
            if (iret .eq. 0) then
!
!         -- 3.1 : MODELE, LIGREL :
                call rslesd(resu, nuordr, model_ = modele)
                if (modele .ne. modeav) then
                    if (ldetli) call detrsd('LIGREL', ligrel)
                    call exlima('ZONE_ANALYSE', 1, 'V', modele, ligrel)
                    modeav=modele
!             -- SI ON CREE UN LIGREL, IL FAUT VERIFIER QUE L'ON S'EN
!                SERT VRAIMENT. SINON, IL FAUT LE DETRUIRE:
                    ldetli=.false.
                    if (ligrel(1:8) .ne. modele) ldetli=.true.
                endif
!
                call rsexch(' ', resuTmp, nomsym, nuordr, chextr,&
                            iret)
                ASSERT(iret.eq.100)
!
                call jelira(chin//'.CELV', 'TYPE', cval=tsca)
                ASSERT(tsca .eq. 'R')
!
                if (lnoeu) then
                    lpain(1)='PEFFONR'
                    lpaout(1)='PEFFOENR'
                else
                    ASSERT(.false.)
                endif
!
                lchin(1)=chin
                lchout(1)=chextr
!
                call calcul('C', 'EFGE_EQUIV', ligrel, 1, lchin,&
                            lpain, 1, lchout, lpaout, 'V',&
                            'OUI')
!
                call jeexin(lchout(1)//'.CELV', iexi)
                if (iexi .eq. 0) then
                    call utmess('A', 'CALCULEL2_18', si=nuordr)
                else
!                    call imprsd('CHAMP',chextr, 6, 'test')
                    ldetli=.false.
                    call rsnoch(resuTmp, nomsym, nuordr)
                endif
            endif
        end do
    

!
!       -- 3. RECOPIE DES PARAMETRES DE RESU VERS resuTmp :
!       --------------------------------------------------
        nompar='&&RC3600'//'.NOMS_PARA'
        call rsnopa(resu, 2, nompar, nbac, nbpa)
        nbpara=nbac+nbpa
        call jeveuo(nompar, 'L', jnompa)
        resuTmp19 = resuTmp
        call jeveuo(resuTmp19//'.ORDR', 'L', jordr)
        call jelira(resuTmp19//'.ORDR', 'LONUTI', nbordr)
!
        do i = 1, nbordr
            nuordr=zi(jordr-1+i)
            do j = 1, nbpara
                nopara=zk16(jnompa-1+j)
                call rsadpa(resu, 'L', 1, nopara, nuordr,&
                            1, sjv=iadin, styp=type, istop=0)
                call rsadpa(resuTmp, 'E', 1, nopara, nuordr,&
                            1, sjv=iadou, styp=type)
                if (type(1:1) .eq. 'I') then
                    zi(iadou)=zi(iadin)
                else if (type(1:1).eq.'R') then
                    zr(iadou)=zr(iadin)
                else if (type(1:1).eq.'C') then
                    zc(iadou)=zc(iadin)
                else if (type(1:3).eq.'K80') then
                    zk80(iadou)=zk80(iadin)
                else if (type(1:3).eq.'K32') then
                    zk32(iadou)=zk32(iadin)
                else if (type(1:3).eq.'K24') then
                    zk24(iadou)=zk24(iadin)
!                    if (nopara(1:5) .eq. 'EXCIT' .and. zk24(iadin)(1:2) .ne. '  ') then
!                        zk24(iadou)=resuTmp//zk24(iadin)(9:)
!                        call copisd(' ', 'V', zk24(iadin)(1:19), zk24(iadou)( 1:19))
!                    endif
                else if (type(1:3).eq.'K16') then
                    zk16(iadou)=zk16(iadin)
                else if (type(1:2).eq.'K8') then
                    zk8(iadou)=zk8(iadin)
                endif
            enddo
        enddo
        call jedetr(nompar)
!
!
!       -- 6. RECOPIE DE L'OBJET .REFD (SI NECESSAIRE) :
!       --------------------------------------------------
!        call refdcp(resu19, resuTmp19)
!
        call jedetr('&&RC3600.NUME_ORDRE')

!       sélection des maximas
!       ---------------------

        typmax = 'MAXI_ABS'
        tychlu = 'ELNO_SIEF_R'
        noch19 = '&&RC3600_CHMAXI'
        typresu = 'VALE_ABS'
        
        call chmima(resuTmp, nomsym, tychlu, typmax, noch19,&
                    typresu=typresu, mcfz = mcf)
        
        if (ldetli) call detrsd('LIGREL', ligrel)
        
    else
        nomsym = ' '
        conceptin = champ
        call dismoi('NOM_MODELE', champ, 'CHAMP', repk=modele)
    endif


!    creation de la table

    call rc3600_chtotab(nomtab, conceptin, nomsym, modele, noch19)
!
    call jedema()
end subroutine
