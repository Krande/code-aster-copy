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

subroutine thcalr(newcal, tysd, knum, lload_name, resuco, &
                  resuc1, nbordr, modele, mate, cara, &
                  nb_load)
    implicit none
!
! ------------------------------------------------------------------
! IN  NEWCAL : TRUE POUR UN NOUVEAU CONCEPT RESULTAT, FALSE SINON
! IN  TYSD   : TYPE DU CONCEPT ATTACHE A RESUCO
! IN  KNUM   : NOM D'OBJET DES NUMERO D'ORDRE
! IN  KCHA   : NOM JEVEUX OU SONT STOCKEES LES CHARGES
! IN  PHENO  : PHENOMENE (MECA,THER,ACOU)
! IN  RESUCO : NOM DE CONCEPT RESULTAT
! IN  RESUC1 : NOM DE CONCEPT DE LA COMMANDE CALC_ERREUR
! IN  CONCEP : TYPE DU CONCEPT ATTACHE A RESUC1
! IN  NBORDR : NOMBRE DE NUMERO D'ORDRE
! IN  MODELE : NOM DU MODELE
! IN  MATE   : NOM DU CHAMP MATERIAU
! IN  CARA   : NOM DU CHAMP DES CARACTERISTIQUES ELEMENTAIRES
! IN  NCHAR  : NOMBRE DE CHARGES
! ----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/calcop.h"
#include "asterfort/dismoi.h"
#include "asterfort/erglth.h"
#include "asterfort/exlima.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mecara.h"
#include "asterfort/mecham.h"
#include "asterfort/medom1.h"
#include "asterfort/modopt.h"
#include "asterfort/reslgn.h"
#include "asterfort/resth2.h"
#include "asterfort/resthe.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexc1.h"
#include "asterfort/rsexc2.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsnopa.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbordr, nb_load
    integer(kind=8) :: vali
    character(len=8) :: resuco, resuc1, modele, cara
    character(len=16) :: tysd
    character(len=19) :: knum, lload_name
    character(len=24) :: mate
    aster_logical :: newcal
!
!     --- VARIABLES LOCALES ---
    integer(kind=8) :: iaux, jordr, iordr, jcha, iret1, iret, bufin1, iad
    integer(kind=8) :: ifm, niv, linst, niveau, n2
    integer(kind=8) :: nuord, nh, nbac, nbpa, jpa, nbpara
    integer(kind=8) :: iadin, iadou, iopt, nbopt
    integer(kind=8) :: jopt, j
    real(kind=8) :: valthe, insold, inst
    character(len=4) :: type
    character(len=6) :: nompro
    parameter(nompro='THCALR')
    character(len=8) :: ma, k8b
    character(len=8) :: psourc
    character(len=16) :: option, nomcmd, k16b
    character(len=19) :: cartef, nomgdf, carteh, nomgdh, cartet, nomgdt, cartes
    character(len=19) :: nomgds, leres1
    character(len=24) :: chcara(18), chelem, chtemm, chtemp
    character(len=24) :: chflum, chsour, chflup, cherre, cherrn
    character(len=24) :: chgeom, chharm, nompar, mateco
    character(len=24) :: lesopt, blan24
    character(len=24) :: ligrel, ligrmo
    aster_logical :: evol
    real(kind=8) :: zero
    parameter(zero=0.d0)
!
    call jemarq()
    call getres(k8b, k16b, nomcmd)
    call jerecu('V')
!
!            '123456789012345678901234'
    blan24 = '                        '
    k8b = '        '
    nh = 0
    chgeom = blan24
    chtemp = blan24
    chharm = blan24
    chelem = blan24
    lesopt = '&&'//nompro//'.LES_OPTION     '
!
!
    call infmaj()
    call infniv(ifm, niv)
!
    call getvtx(' ', 'OPTION', nbval=0, nbret=n2)
    nbopt = -n2
    call wkvect(lesopt, 'V V K16', nbopt, jopt)
    call getvtx(' ', 'OPTION', nbval=nbopt, vect=zk16(jopt), nbret=n2)
    call modopt(resuco, modele, lesopt, nbopt)
    call jeveuo(lesopt, 'L', jopt)
!
    call jeveuo(lload_name//'.LCHA', 'L', jcha)

    call jeveuo(knum, 'L', jordr)
!
    if (newcal) then
        call rscrsd('G', resuc1, tysd, nbordr)
        call titre()
    end if
!
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrmo)
!
    call jenonu(jexnom(resuco//'           .NOVA', 'INST'), iret)
!
    call exlima(' ', 0, 'V', modele, ligrel)
!
!
    leres1 = resuc1
!
!
    if (newcal) then
        call rscrsd('G', leres1, tysd, nbordr)
        call titre()
    end if
!
!
!============ DEBUT DE LA BOUCLE SUR LES OPTIONS A CALCULER ============
    do iopt = 1, nbopt
        option = zk16(jopt+iopt-1)
!
        call jeveuo(knum, 'L', jordr)
!
        call calcop(option, lesopt, resuco, resuc1, knum, &
                    nbordr, tysd, iret)
        if (iret .eq. 0) goto 120
!
!
        nuord = zi(jordr)
        call medom1(modele, mate, mateco, cara, lload_name, nb_load, &
                    resuco, nuord)
        call jeveuo(lload_name//'.LCHA', 'L', jcha)
!
        call mecham(option, modele, cara, nh, chgeom, &
                    chcara, chharm, iret)
        if (iret .ne. 0) goto 190
!
!    ------------------------------------------------------------------
!    -- OPTION "ERTH_ELEM"
!    ------------------------------------------------------------------
!
        if (option .eq. 'ERTH_ELEM') then
!
!
! PAR DECRET EDA DU 22/08/01 ON SUPPRIME LE PARAMETRE NIVEAU ET ON LE
! FIXE A 2 (15 VALEURS DE PARAMETRES).
            niveau = 2
!
! RECUPERATION NIVEAU AFFICHAGE
            call infniv(ifm, niv)
!
! BOUCLE SUR LES INSTANTS CHOISIS PAR LE USER
            insold = zero
            chtemm = ' '
            chtemp = ' '
            chflum = ' '
            chflup = ' '
!
! PREPARATION DES CALCULS D'INDICATEUR (CONNECTIVITE INVERSE, CHARGE)
            call jeveuo(lload_name//'.LCHA', 'L', jcha)
            call resth2(modele, ligrmo, zk8(jcha), nb_load, ma, &
                        cartef, nomgdf, carteh, nomgdh, cartet, &
                        nomgdt, cartes, nomgds, chgeom, chsour, &
                        psourc)
!
            if (niv .ge. 1) then
                write (ifm, *)
                write (ifm, *)&
     &       '*********************************************'
                write (ifm, *) '  CALCUL DE CARTES D''ERREURS EN RESIDU'
                write (ifm, *) '       POUR LE PROBLEME THERMIQUE'
                write (ifm, *)
                write (ifm, *) '  OPTION DE CALCUL   ERTH_ELEM'
                write (ifm, *) '  MODELE                ', modele
                write (ifm, *) '  SD EVOL_THER DONNEE   ', resuco
                write (ifm, *) '             RESULTAT   ', resuc1
                write (ifm, *)
                write (ifm, *)&
     &        '* CONTRAIREMENT AUX CALCULS THERMIQUES, POUR *'
                write (ifm, *)&
     &        '* UN TYPE DE CHARGEMENT DONNE, ON NE RETIENT *'
                write (ifm, *)&
     &        '* QUE LA DERNIERE OCCURENCE DE AFFE_CHAR_THER*'
                write (ifm, *) '  LISTE DES CHARGEMENTS :'
                do bufin1 = 1, nb_load
                    write (ifm, *) '                        ', zk8(jcha+bufin1-1)
                end do
                write (ifm, *) '  CL DE FLUX RETENUE      ', nomgdf
                write (ifm, *) '  CL D''ECHANGE RETENUE    ', nomgdh
                write (ifm, *) '  SOURCE RETENUE          ', nomgds
                write (ifm, *) '  MATERIAU PRIS EN COMPTE ', mate(1:8)
                write (ifm, *) '  NOMBRE DE NUMERO D''ORDRE ', nbordr
            end if
!
! BOUCLE SUR LES PAS DE TEMPS
            do iaux = 1, nbordr
                call jemarq()
                call jerecu('V')
                iordr = zi(jordr+iaux-1)
                call medom1(modele, mate, mateco, cara, lload_name, nb_load, &
                            resuco, iordr)
                call mecara(cara, chcara)
! RECUPERATION DU PARM_THETA CORRESPONDANT A IORDR
                call jenonu(jexnom(resuco//'           .NOVA', 'PARM_THETA'), iad)
                if (iad .eq. 0) then
                    valthe = 0.57d0
                    call utmess('A', 'CALCULEL4_98', sk=resuco)
                else
                    call rsadpa(resuco, 'L', 1, 'PARM_THETA', iordr, &
                                0, sjv=iad, styp=k8b)
                    valthe = zr(iad)
                    if ((valthe .gt. 1.d0) .or. (valthe .lt. 0.d0)) then
                        call utmess('F', 'INDICATEUR_5', sk=resuco)
                    end if
                end if
                if (niv .ge. 1) then
                    write (ifm, *) '   PARAM-THETA/IORDR ', valthe, iordr
                    if (iaux .eq. nbordr) then
                        write (ifm, *) &
                            '*************************************'// &
                            '*********'
                        write (ifm, *)
                    end if
                end if
!
! CALCUL DU CRITERE D'EVOLUTION LEVOL (TRUE=TRANSITOIRE)
! CAS PARTICULIER DE L'INSTANT INITIAL D'UN CALCUL TRANSITOIRE
! ON ESTIME SON ERREUR COMME EN STATIONNAIRE
                if (iaux .eq. 1) then
                    evol = .false.
                else
                    evol = .true.
                end if
!
! RECUPERATION DU NOM DES CHAMP_GD = RESUCO('FLUX_ELNO',I)
! ET RESUCO('TEMP',I) POUR I=IORDR. POUR IORDR-1 ILS SONT STOCKES
! DANS CHFLUM/CHTEMM DEPUIS LA DERNIERE ITERATION.
! RESUCO = NOM USER DE LA SD DESIGNEE PAR LE MOT-CLE RESULTAT
                call rsexc2(1, 1, resuco, 'TEMP', iordr, &
                            chtemp, option, iret)
                if (iret .gt. 0) then
                    vali = iordr
                    call utmess('F', 'CALCULEL6_46', si=vali)
                end if
                call rsexc2(1, 1, resuco, 'FLUX_ELNO', iordr, &
                            chflup, option, iret)
                if (iret .gt. 0) then
                    vali = iordr
                    call utmess('F', 'CALCULEL6_47', si=vali)
                end if
!
! RECUPERATION DE L'INSTANT CORRESPONDANT A IORDR
                call rsadpa(resuco, 'L', 1, 'INST', iordr, &
                            0, sjv=linst, styp=k8b)
                inst = zr(linst)
!
! IMPRESSIONS NIVEAU 2 POUR DIAGNOSTIC...
                if (niv .eq. 2) then
                    write (ifm, *) nompro, ' **********'
                    write (ifm, *) 'EVOL/I/IORDR', evol, iaux, iordr
                    write (ifm, *) 'INST/INSOLD', inst, insold
                    write (ifm, *) 'CHTEMM/CHTEMP', chtemm, ' / ', chtemp
                    write (ifm, *) 'CHFLUM/CHFLUP', chflum, ' / ', chflup
                end if
!
! RECUPERATION DU NOM DU CHAMP_GD = RESUC1('ERTH_ELEM',IORDR)
! RESUC1 = NOM USER DE LA SD CORRESPONDANT AU RESULTAT DE CALC_ERREUR
                call rsexc1(leres1, option, iordr, chelem)
! PREPARATION DES DONNEES/LANCEMENT DU CALCUL DES INDICATEURS
                call resthe(ligrmo, evol, chtemm, chtemp, chflum, &
                            chflup, mateco, valthe, insold, inst, &
                            chelem, niveau, ifm, niv, ma, &
                            cartef, nomgdf, carteh, nomgdh, cartet, &
                            nomgdt, cartes, nomgds, chgeom, chsour, &
                            psourc, iaux)
! CALCUL DE L'ESTIMATEUR GLOBAL
                call erglth(chelem, inst, niveau, iordr, resuco)
! NOTATION DE LA SD RESULTAT LERES1
                call rsnoch(leres1, option, iordr)
!
! INIT. POUR LE NUMERO D'ORDRE SUIVANT
                if (nbordr .ne. 1 .and. iaux .ne. nbordr) then
                    chtemm = chtemp
                    chflum = chflup
                    insold = inst
                end if
                call jedema()
            end do
! DESTRUCTION DES OBJETS JEVEUX VOLATILES
            call jedetr(cartef//'.PTMA')
            call jedetr(carteh//'.PTMA')
            call jedetr(cartet//'.PTMA')
            call jedetr(cartef//'.PTMS')
            call jedetr(carteh//'.PTMS')
            call jedetr(cartet//'.PTMS')
!
!    ------------------------------------------------------------------
!    -- OPTION "ERTH_ELNO"
!    ------------------------------------------------------------------
!
        else if (option .eq. 'ERTH_ELNO') then
!
            do iaux = 1, nbordr
                call jemarq()
                call jerecu('V')
                iordr = zi(jordr+iaux-1)
! RECUPERATION DU NOM DU CHAMP_GD = RESUCO('ERTH_ELEM',IORDR)
                call rsexc2(1, 1, resuco, 'ERTH_ELEM', iordr, &
                            cherre, option, iret1)
                if (iret1 .gt. 0) goto 40
! RECUPERATION DU NOM DU CHAMP_GD = RESUC1('ERTH_ELNO',IORDR)
! RESUC1 = NOM USER DE LA SD CORRESPONDANT AU RESULTAT DE CALC_ERREUR
                call rsexc1(leres1, option, iordr, cherrn)
                call reslgn(ligrmo, option, cherre, cherrn)
! NOTATION DE LA SD RESULTAT LERES1
                call rsnoch(leres1, option, iordr)
40              continue
                call jedema()
            end do
!
!    ------------------------------------------------------------------
        else
            call utmess('A', 'CALCULEL3_22', sk=option)
        end if
!
120     continue
    end do
!       ====== FIN DE LA BOUCLE SUR LES OPTIONS A CALCULER =======
!
    if (newcal) then
        nompar = '&&'//nompro//'.NOMS_PARA '
        call rsnopa(resuco, 2, nompar, nbac, nbpa)
        nbpara = nbac+nbpa
        call jeveuo(nompar, 'L', jpa)
        do iaux = 1, nbordr
            iordr = zi(jordr+iaux-1)
            do j = 1, nbpara
                call rsadpa(resuco, 'L', 1, zk16(jpa+j-1), iordr, &
                            1, sjv=iadin, styp=type)
                call rsadpa(leres1, 'E', 1, zk16(jpa+j-1), iordr, &
                            1, sjv=iadou, styp=type)
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
    end if
!
!
190 continue
!
    call jedema()
end subroutine
