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

subroutine modirepresu(resuou, resuin)
!
    implicit none
    character(len=19) :: resuou, resuin
!
! ----------------------------------------------------------------------
!
!     COMMANDE : MODI_REPERE / RESULTAT
!
! ----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterc/getfac.h"
#include "asterfort/gettco.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesfus.h"
#include "asterfort/chrpel.h"
#include "asterfort/chrpno.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlima.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/medome_once.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsinfo.h"
#include "asterfort/rslesd.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsnopa.h"
#include "asterfort/rsorac.h"
#include "asterfort/rsutnu.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: n0, n1, nbordr, iret, nocc, i, j, np, iordr, ndim
    integer(kind=8) :: iord, ioc, ibid, nc
    integer(kind=8) :: jordr, nbnosy, jpa, iadin, iadou
    integer(kind=8) :: nbpara, nbac, nbpa, ifm, niv, nncp
    real(kind=8) :: prec
    real(kind=8) :: lcoer(2)
    real(kind=8) :: r8b
    complex(kind=8) :: c16b
    complex(kind=8) :: lcoec(2)
    character(len=8) :: k8b
    character(len=8) :: crit, tych, nomma, model, modelRefe
    character(len=8) :: carele, exipla, exicoq
    character(len=16) :: option, tysd, type, type_cham, repere, option1, optinit
    character(len=16) :: cham_resu
    character(len=19) :: knum
    character(len=19) :: chams1, chams0, chafus, chs(2), ligrel, ligrel1
    character(len=24) :: nompar, champ0, champ01, champ1, champ2
    character(len=24) :: valk(2)
    integer(kind=8), pointer :: nume_ordre(:) => null()
!
    aster_logical :: lreuse, lcumu(2), lcoc(2), lModelVariable, check
!
    data lcumu/.false., .false./
    data lcoc/.false., .false./
    data lcoer/1.d0, 1.d0/
! ---------------------------------------------------------------------
    call jemarq()
!
    call infmaj()
    call infniv(ifm, niv)
!
!   LE CONCEPT EST REENTRANT SI REPERE = 'COQUE_INTR_UTIL' OU 'COQUE_UTIL_INTR' ou 'COQUE_UTIL_CYL'
!   DANS CE CAS ON CREE UNE SD RESULTAT TEMPORAIRE POUR LES CALCULS ET ENSUITE ON SURCHARGE
!   RESUIN PAR LES CHAMPS MODIFIES STOCKES DANS RESUOU
    lreuse = .false.
    if (resuin .eq. resuou) then
        lreuse = .true.
        resuou = 'MODIREPE'
    end if
!
    call jelira(resuin//'.DESC', 'NOMMAX', nbnosy)
    if (nbnosy .eq. 0) goto 999
!
    call gettco(resuin, tysd)
!   RECUPERATION DU NOMBRE DE CHAMPS SPECIFIE
    call getfac('MODI_CHAM', nocc)
!
!   DEFINITION DU REPERE UTILISE
    call getvtx(' ', 'REPERE', scal=repere, nbret=i)

    if (lreuse) then
        if (i .eq. 0) then
            call utmess('F', 'MODELISA3_14')
        end if
        if ((repere .ne. 'COQUE_INTR_UTIL') .and. &
            (repere .ne. 'COQUE_UTIL_INTR') .and. &
            (repere .ne. 'UTILISATEUR')) then
            call utmess('F', 'MODELISA3_15', nk=1, valk=repere)
        end if
    end if
!
!   RECUPERATION DES NUMEROS D'ORDRE DE LA STRUCTURE DE DONNEES DE TYPE RESULTAT RESU A PARTIR
!   DES VARIABLES D'ACCES UTILISATEUR 'NUME_ORDRE','FREQ','INST','NOEUD_CMP'
!   (VARIABLE D'ACCES 'TOUT_ORDRE' PAR DEFAUT)
    knum = '&&OP0191.NUME_ORDRE'
    call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
    call rsutnu(resuin, ' ', 1, knum, nbordr, prec, crit, iret)
    if (iret .eq. 10) then
        call utmess('F', 'CALCULEL4_8', sk=resuin)
    end if
    if (iret .ne. 0) then
        call utmess('F', 'ALGORITH3_41')
    end if
    call jeveuo(knum, 'L', jordr)
    call rscrsd('G', resuou, tysd, nbordr)
!
    ! QUELQUES INITIALISATIONS
    lModelVariable = ASTER_FALSE
    modelRefe = " "
    ligrel1 = " "
    option1 = 'EGRU_ELNO       '
    optinit = ' '
    n1 = 0

    ! POUR REPERE = COQUE_* VERIFIER QU'ON N'EST PAS EN MULTI MODELES
    if (repere(1:5) .eq. 'COQUE') then
        AS_ALLOCATE(vi=nume_ordre, size=nbordr)
        call rsorac(resuin, 'TOUT_ORDRE', ibid, r8b, k8b, &
                    c16b, r8b, k8b, nume_ordre, nbordr, &
                    ndim)
        call medome_once(resuin, nume_ordre, nbordr)
        AS_DEALLOCATE(vi=nume_ordre)
    end if

    do ioc = 1, nocc
        call getvtx('MODI_CHAM', 'NOM_CHAM', iocc=ioc, scal=option, nbret=n0)
        call getvtx('MODI_CHAM', 'TYPE_CHAM', iocc=ioc, scal=type_cham, nbret=n0)
        call getvtx('MODI_CHAM', 'NOM_CHAM_RESU', iocc=ioc, scal=cham_resu, nbret=n1)
        optinit = option
        check = .false.
        do iord = 1, nbordr
            ! SI OPTION CHANGE DANS LE CASE DE EGRU_ELNO
            option = optinit
            !
            call jemarq()
            call jerecu('V')
            iordr = zi(jordr-1+iord)
            call rsexch('F', resuin, option, iordr, champ0, iret)
            call dismoi('NOM_MAILLA', champ0(1:19), 'CHAMP', repk=nomma)
            call dismoi('TYPE_CHAMP', champ0, 'CHAMP', repk=tych, arret='C', ier=iret)

            !
            if ((option(1:9) .eq. option1(1:9)) .and. (repere(1:15) .eq. 'COQUE_UTIL_INTR')) then
                valk(1) = option(1:9)
                valk(2) = 'COQUE_UTIL_INTR'
                call utmess('F', 'ALGORITH5_87', nk=2, valk=valk)
            end if

            ! TRAITEMENT SPECIFIQUE POUR EFGE_ELNO
            if ((option(1:9) .eq. 'EFGE_ELNO') .and. (n1 .gt. 0)) then
                ASSERT(cham_resu(1:9) .eq. 'EGRU_ELNO')
                ! CREER LE CHAMP EGRU
                call rsexch(' ', resuou, option, iordr, champ2, iret)
                call copisd('CHAMP_GD', 'G', champ0, champ2)
                call rsnoch(resuou, option, iordr)
                if (lreuse) then
                    call rsexch(' ', resuin, option1, iordr, champ01, iret)
                    check = .true.
                end if
                option = option1
            end if

            ! CHAMP1 SERA ENSUITE RECREE SUR LA BASE GLOBALE
            call rsexch(' ', resuou, option, iordr, champ1, iret)
            call copisd('CHAMP_GD', 'G', champ0, champ1)

!           RECUPERATION DU MODELE ASSOCIE AU CHAMP
            model = ''; carele = ''
            call rslesd(resuin(1:8), iordr, model_=model, cara_elem_=carele)
            if (iord .eq. 1) then
                modelRefe = model
            else
                if (modelRefe .ne. model) then
                    lModelVariable = ASTER_TRUE
                end if
            end if
            if (model .ne. '') then
                call dismoi('EXI_PLAQUE', model, 'MODELE', repk=exipla)
                call dismoi('EXI_COQUE', model, 'MODELE', repk=exicoq)
                if (((exipla(1:3) .eq. 'OUI') .or. (exicoq(1:3) .eq. 'OUI')) .and. &
                    ((type_cham .eq. 'TENS_2D') .or. (type_cham .eq. 'TENS_3D')) .and. &
                    (repere .eq. 'UTILISATEUR')) then
                    call utmess('F', 'ALGORITH3_7')
                end if
            end if
!               Obligatoire : modèle , cara_elem
!                             repere = UTILISATEUR

            if (type_cham .eq. '1D_GENE') then
                if ((model .eq. '') .or. (carele .eq. '') .or. (repere .ne. 'UTILISATEUR')) then
                    call utmess('F', 'ALGORITH2_32')
                end if
            end if
!
!           RECUPERATION DE LA NATURE DES CHAMPS (CHAM_NO OU CHAM_ELEM)
            if (tych(1:4) .eq. 'NOEU') then
                call chrpno(champ1, repere, option, type_cham)
            else if (tych(1:2) .eq. 'EL') then
                if (iord .eq. 1 .or. modelRefe .ne. model) then
                    call exlima('MODI_CHAM', ioc, 'G', model, ligrel1)
                end if
                call chrpel(champ1, repere, option, ioc, type_cham, &
                            option, model, carele, ligrel1, lModelVariable)
            else
                valk(1) = tych
                valk(2) = champ1
                call utmess('A', 'ALGORITH9_69', nk=2, valk=valk)
            end if
            call rsnoch(resuou, option, iordr)
            ! COPIER CHAMP1 SI EGRU_ELNO ET REUSE
            if (check) then
                call copisd('CHAMP_GD', 'G', champ1, champ01)
                call rsnoch(resuin, option, iordr)
            end if
            call jedema()
        end do
    end do

!
    nompar = '&&OP0191.NOMS_PARA'
    call rsnopa(resuin, 2, nompar, nbac, nbpa)
    nbpara = nbac+nbpa
    call jeveuo(nompar, 'L', jpa)
    do iord = 1, nbordr
        iordr = zi(jordr-1+iord)
        do j = 1, nbpara
            call rsadpa(resuin, 'L', 1, zk16(jpa+j-1), iordr, 1, sjv=iadin, styp=type, istop=0)
            call rsadpa(resuou, 'E', 1, zk16(jpa+j-1), iordr, 1, sjv=iadou, styp=type)
            if (type(1:1) .eq. 'I') then
                zi(iadou) = zi(iadin)
            else if (type(1:1) .eq. 'R') then
                zr(iadou) = zr(iadin)
            else if (type(1:1) .eq. 'C') then
                zc(iadou) = zc(iadin)
            else if (type(1:3) .eq. 'K80') then
                zk80(iadou) = zk80(iadin)
            else if (type(1:3) .eq. 'K32') then
                zk32(iadou) = zk32(iadin)
            else if (type(1:3) .eq. 'K24') then
                zk24(iadou) = zk24(iadin)
            else if (type(1:3) .eq. 'K16') then
                zk16(iadou) = zk16(iadin)
            else if (type(1:2) .eq. 'K8') then
                zk8(iadou) = zk8(iadin)
            end if
        end do
    end do
!
    call titre()
    if (niv .eq. 2) call rsinfo(resuou, ifm)
!
999 continue
!
!   CREATION DE L'OBJET .REFD SI NECESSAIRE :
    call refdcp(resuin, resuou)
!
!   Traitement du cas ou il y a reentrance
!       REPERE = 'COQUE_INTR_UTIL' ou 'COQUE_UTIL_INTR' ou 'UTILISATEUR'
    if (lreuse) then
        do ioc = 1, nocc
            call getvtx('MODI_CHAM', 'NOM_CHAM', iocc=ioc, scal=option, nbret=n0)
            do iord = 1, nbordr
                call jemarq()
                call jerecu('V')
                iordr = zi(jordr-1+iord)
                call rsexch('F', resuin, option, iordr, champ0, iret)
                call rsexch(' ', resuou, option, iordr, champ1, iret)
                chams0 = '&&CHRPEL.CHAMS0'
                chams1 = '&&CHRPEL.CHAMS1'
                chafus = '&&CHRPEL.CHAFUS'
                chs(1) = chams0
                chs(2) = chams1
                call celces(champ0, 'V', chams0)
                call celces(champ1, 'V', chams1)
                call cesfus(2, chs, lcumu, lcoer, lcoec, lcoc(1), 'V', chafus)
                call dismoi('NOM_LIGREL', champ0, 'CHAM_ELEM', repk=ligrel)
                call cescel(chafus, ligrel, option, ' ', 'NAN', nncp, 'G', champ0, 'F', ibid)
                call detrsd('CHAMP', champ1)
                call jedema()
            end do
        end do
        call detrsd('CHAMP', chams0)
        call detrsd('CHAMP', chams1)
        call detrsd('CHAMP', chafus)
        call detrsd('RESULTAT', resuou)
    end if
!
    call jedema()
end subroutine
