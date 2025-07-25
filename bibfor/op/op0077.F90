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

subroutine op0077()
    implicit none
!
!
!     OPERATEUR REST_SOUS_STRUC
!
! ----------------------------------------------------------------------
!
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/dcapno.h"
#include "asterfort/excygl.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/harm75.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/recyec.h"
#include "asterfort/recygl.h"
#include "asterfort/refdaj.h"
#include "asterfort/regeec.h"
#include "asterfort/rege2c.h"
#include "asterfort/regegl.h"
#include "asterfort/regegc.h"
#include "asterfort/regene.h"
#include "asterfort/regres.h"
#include "asterfort/rehaec.h"
#include "asterfort/rehagl.h"
#include "asterfort/retrec.h"
#include "asterfort/retrgl.h"
#include "asterfort/rsadpa.h"
#include "asterfort/tran77.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!
!
    character(len=1) :: base
    character(len=8) :: k8b, nomres, resin, nomsst, mailsk, mode
    character(len=8) :: k8bid, result, blanc, param(3)
    character(len=16) :: concep, nomcmd, typres, typrep, champ(4)
    character(len=19) :: nume_equa
    character(len=24) :: matgen, numgen, vbl24(1)
    integer(kind=8) :: ioc1, nbord, i, iord, lpaout(3)
!
!     -----------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid, ir, ir1, iret, isk, j
    integer(kind=8) :: j2refe, j3refe, jrefn, jrefnb, lmacr, lmodge
    integer(kind=8) ::  n1, n2, nbcham, numsec, neq, neq2
    integer(kind=8), pointer :: ordr(:) => null()
    character(len=3) :: typesca
    character(len=8), pointer :: refm(:) => null()
    character(len=16) :: depl
    character(len=24) :: chamol
!-----------------------------------------------------------------------
    data depl/'DEPL            '/
    data k8b/'        '/
    data param/'MODELE', 'CHAMPMAT', 'CARAELEM'/
    data vbl24/'                        '/
!
    call jemarq()
    call infmaj()
    k8bid = '        '
    blanc = '        '
!
!     -----------------------------------------------------------------
!
!
    call getres(nomres, typres, nomcmd)
!
! --- PHASE DE TEST SUR LES CHAMPS A RESTITUER
    call getvtx(' ', 'NOM_CHAM', nbval=4, vect=champ, nbret=nbcham)
    if (nbcham .lt. 0) then
        call utmess('E', 'ALGORITH9_44')
    else
        do i = 1, nbcham
            do j = i+1, nbcham
                if (champ(i) .eq. champ(j)) then
                    call utmess('E', 'ALGORITH9_30')
                end if
            end do
        end do
    end if
!
!
! --- CREATION DU PROFIL :
!     ---------------------------
    call getvid(' ', 'SQUELETTE', scal=k8b, nbret=ir)
!     --- SI RESTITUTION SUR UNE SQUELETTE, ALORS ATTACHER UN NUME_EQUA
!         AU RESULTAT
    if (ir .eq. 0) then
        nume_equa = '&&OP0077.PROFC.NUME'
        base = 'V'
    else
        base = 'G'
        nume_equa = nomres//'.PROFC.NUME'
    end if
! --- CREATION D'UN OBJET REFN DU PROFIL SUR BASE VOLATILE
    call wkvect(nume_equa//'.REFN', base//' V K24', 5, jrefn)
    zk24(jrefn+1) = 'DEPL_R'
!
!
!
! --- LE RESULTAT EST-IL GENERALISE OU PAS :
!     ---------------------------
    call getvid(' ', 'RESULTAT', nbval=0, nbret=ir)
    if (ir .eq. 0) then
        call getvid(' ', 'RESU_GENE', scal=resin, nbret=ir1)
        call gettco(resin, concep)
    else
!      --- PROJECTION RESULTAT SUR UN SQUELETTE ENRICHI ---
        call getvid(' ', 'SQUELETTE', scal=mailsk, nbret=ibid)
        call getvid(' ', 'RESULTAT', scal=result, nbret=ibid)
        zk24(jrefn) = mailsk
        call getfac('CYCLIQUE', ioc1)
        if (ioc1 .gt. 0) then
            call excygl(nomres, typres, result, mailsk, nume_equa)
            call jeveuo(nume_equa//'.REFN', 'E', jrefnb)
            zk24(jrefnb) = mailsk
            zk24(jrefnb+1) = 'DEPL_R'
            concep(1:9) = '         '
            resin = result
            goto 30
!
        else
            if (typres .eq. 'MODE_MECA') then
                call regres(nomres, mailsk, result, nume_equa)
                call jeveuo(nume_equa//'.REFN', 'E', jrefnb)
                zk24(jrefnb) = mailsk
                zk24(jrefnb+1) = 'DEPL_R'
                concep(1:9) = '         '
                resin = result
                goto 30
!
            else
                call utmess('E', 'ALGORITH9_46')
            end if
        end if
    end if
!
!
! INDICATEUR CALCUL SANS MATRICE GENERALISEE (PROJ_MESU_MODAL)
!      PROMES=.FALSE.
    if ((concep(1:9) .eq. 'TRAN_GENE') .or. (concep(1:9) .eq. 'MODE_GENE') .or. &
        (concep(1:9) .eq. 'HARM_GENE')) then
        call dismoi('REF_RIGI_PREM', resin, 'RESU_DYNA', repk=matgen, arret='C')
        call dismoi('NUME_DDL', resin, 'RESU_DYNA', repk=numgen)
! LE RESU_GENE VIENT DE PROJ_MESU_MODAL
        if ((matgen(1:8) .eq. blanc) .and. (numgen(1:8) .eq. blanc)) then
!          PROMES=.TRUE.
            typrep = blanc
        else
            if (numgen(1:8) .eq. blanc) then
                call jeveuo(matgen(1:8)//'           .REFA', 'L', j2refe)
                numgen = zk24(j2refe+1) (1:14)
            end if
            call jeveuo(numgen(1:14)//'.NUME.REFN', 'L', j3refe)
            call gettco(zk24(j3refe), typrep)
        end if
    end if
!
!     --- DYNAMIQUE TRANSITOIRE ---
!
    if (concep(1:9) .eq. 'TRAN_GENE') then
!
        if (typrep(1:11) .eq. 'MODELE_GENE') then
            call getvid(' ', 'SQUELETTE', nbval=0, nbret=isk)
            if (isk .eq. 0) then
                call getvtx(' ', 'SOUS_STRUC', scal=nomsst, nbret=ibid)
                call retrec(nomres, resin, nomsst)
            else
                call getvid(' ', 'SQUELETTE', scal=mailsk, nbret=ibid)
                call retrgl(nomres, resin, mailsk, nume_equa)
                call jeveuo(nume_equa//'.REFN', 'E', jrefnb)
                zk24(jrefnb) = mailsk
                zk24(jrefnb+1) = 'DEPL_R'
            end if
!
!
!
        else if (typrep(1:9) .eq. 'MODE_GENE') then
            call getvtx(' ', 'SOUS_STRUC', scal=nomsst, nbret=n1)
            call getvid(' ', 'SQUELETTE', scal=mailsk, nbret=n2)
            if ((n1 .ne. 0 .and. n2 .ne. 0)) then
                call utmess('F', 'ALGORITH9_47')
            end if
            call getvid(' ', 'MODE_MECA', scal=mode, nbret=ibid)
            if (ibid .eq. 0) then
                call utmess('F', 'ALGORITH9_48')
            end if
            call tran77(nomres, typres, resin, mode)
        end if
!
!
!
!     --- CALCUL MODAL PAR SOUS-STRUCTURATION CLASSIQUE ---
!                  OU SANS SOUS-STRUCTURATION
!
    else if (concep(1:9) .eq. 'MODE_GENE') then
!
! --- CAS DE LA SOUS-STRUCTURATION MODALE
        if (typrep(1:11) .eq. 'MODELE_GENE') then
!
            call getvid(' ', 'SQUELETTE', nbval=0, nbret=isk)
            call jeveuo(resin(1:8)//'           .ORDR', 'L', ibid)
            call dcapno(resin, depl, zi(ibid), chamol)
            call dismoi('TYPE_SCA', chamol(1:19), 'CHAMP', repk=typesca)

            if (isk .eq. 0) then
                call getvtx(' ', 'SOUS_STRUC', scal=nomsst, nbret=ibid)
                !-- les routines REGEEC et REGE2C font la meme chose, une en reel,
                !-- l'autre en complexe. En cas de modification d'une des routines,
                !-- ne pas oublier de reporter le changement dans l'autre.
                if (typesca .eq. "R") then
                    call regeec(nomres, resin, nomsst)
                elseif (typesca .eq. "C") then
                    call rege2c(nomres, resin, nomsst)
                else
                    ASSERT(.false.)
                end if

            else
                call getvid(' ', 'SQUELETTE', scal=mailsk, nbret=ibid)

                !-- les routines REGEGL et REGEGC font la meme chose, une en reel,
                !-- l'autre en complexe. En cas de modification d'une des routines,
                !-- ne pas oublier de reporter le changement dans l'autre.
                if (typesca .eq. "R") then
                    call regegl(nomres, resin, mailsk, nume_equa)
                elseif (typesca .eq. "C") then
                    call regegc(nomres, resin, mailsk, nume_equa)
                else
                    ASSERT(.false.)
                end if
                call jeveuo(nume_equa//'.REFN', 'E', jrefnb)
                zk24(jrefnb) = mailsk
                zk24(jrefnb+1) = 'DEPL_R'
            end if
        else
!
!     --- CALCUL MODAL SANS SOUS-STRUCTURATION ---
            call regene(nomres, resin, nume_equa)
        end if
!
!     --- CALCUL MODAL PAR SOUS-STYRUCTURATION CYCLIQUE ---
!
    else if (concep(1:9) .eq. 'MODE_CYCL') then
        call getvid(' ', 'SQUELETTE', nbval=0, nbret=isk)
        if (isk .eq. 0) then
            call getvis(' ', 'SECTEUR', scal=numsec, nbret=ibid)
            call recyec(nomres, resin, numsec, 'MODE_MECA')
        else
            call getvid(' ', 'SQUELETTE', scal=mailsk, nbret=ibid)
            call recygl(nomres, 'MODE_MECA', resin, mailsk, nume_equa)
            call jeveuo(nume_equa//'.REFN', 'E', jrefnb)
            zk24(jrefnb) = mailsk
            zk24(jrefnb+1) = 'DEPL_R'
        end if
!
!     --- CALCUL HARMONIQUE PAR SOUS-STRUCTURATION CLASSIQUE ---
!
    else if (concep(1:9) .eq. 'HARM_GENE') then
!
! --- CAS DE LA SOUS-STRUCTURATION HARMONIQUE
        if (typrep(1:11) .eq. 'MODELE_GENE') then
            call getvid(' ', 'SQUELETTE', nbval=0, nbret=isk)
            if (isk .eq. 0) then
                call getvtx(' ', 'SOUS_STRUC', scal=nomsst, nbret=ibid)
                call rehaec(nomres, resin, nomsst)
            else
                call getvid(' ', 'SQUELETTE', scal=mailsk, nbret=ibid)
                call rehagl(nomres, resin, mailsk, nume_equa)
                call jeveuo(nume_equa//'.REFN', 'E', jrefnb)
                zk24(jrefnb) = mailsk
                zk24(jrefnb+1) = 'DEPL_R'
            end if
!
        else if (typrep(1:9) .eq. 'MODE_GENE') then
            call getvtx(' ', 'SOUS_STRUC', scal=nomsst, nbret=n1)
            call getvid(' ', 'SQUELETTE', scal=mailsk, nbret=n2)
            if ((n1 .ne. 0 .and. n2 .ne. 0)) then
                call utmess('F', 'ALGORITH9_47')
            end if
            call getvid(' ', 'MODE_MECA', scal=mode, nbret=ibid)
            if (ibid .eq. 0) then
                call utmess('F', 'ALGORITH9_48')
            end if
            call harm75(nomres, typres, resin, mode)
        end if
!
    end if
!
30  continue
!
! --- STOCKAGE
    call gettco(resin, concep)
    if ((concep(1:9) .ne. 'TRAN_GENE') .and. (concep(1:9) .ne. 'MODE_CYCL') .and. &
        (concep(1:9) .ne. 'HARM_GENE')) then
        call jeveuo(nomres//'           .ORDR', 'L', vi=ordr)
        call jelira(nomres//'           .ORDR', 'LONUTI', nbord)
        call jelira(nomres//'           .ORDR', 'LONUTI', nbord)
!
        call getvid(' ', 'SQUELETTE', nbval=0, nbret=isk)
        if (isk .eq. 0) then
            call getvtx(' ', 'SOUS_STRUC', scal=nomsst, nbret=ibid)
!
            call dismoi('NUME_DDL', resin, 'RESU_DYNA', repk=numgen)
            call jeveuo(numgen(1:14)//'.NUME.REFN', 'L', lmodge)
            call jenonu(jexnom(zk24(lmodge) (1:8)//'      .MODG.SSNO', nomsst), iret)
            call jeveuo(jexnum(zk24(lmodge) (1:8)//'      .MODG.SSME', iret), 'L', lmacr)
!-- RECUPERATION DES INFOS CARA_ELEM / MATER / MODELE POUR LES SST
!-- DANS LE .REFM DANS LE MACRO ELEMENT CORRESPONDANT
            call jeveuo(zk8(lmacr)//'.REFM', 'L', vk8=refm)
            do iord = 1, nbord
                call rsadpa(nomres, 'E', 3, param, ordr(iord), &
                            0, tjv=lpaout, styp=k8b)
                zk8(lpaout(1)) = refm(1)
                zk8(lpaout(2)) = refm(3)
                zk8(lpaout(3)) = refm(4)
            end do
        end if
!
    end if
!
!     -- CREATION DE L'OBJET .REFD SI NECESSAIRE:
!     -------------------------------------------
    call jeexin(nomres//'           .REFD', iret)
    if (iret .eq. 0) call refdaj(' ', nomres, -1, nume_equa, 'INIT', &
                                 ' ', iret)

    if (base == 'G') then
        call jeexin(nume_equa//'.NUEQ', iret)
        if (iret > 0) then
            call jelira(nume_equa//'.NUEQ', 'LONMAX', neq)
        else
            neq = 1
        end if
        call jeexin(nume_equa//'.NEQU', iret)
        if (iret > 0) then
            call jeveuo(nume_equa//'.NEQU', 'E', jrefn)
            ASSERT(zi(jrefn) == neq)
        else
            call wkvect(nume_equa//'.NEQU', 'G V I', 2, jrefn)
            zi(jrefn) = neq
            zi(jrefn+1) = neq
        end if

        call jeexin(nume_equa//'.DELG', iret)
        if (iret > 0) then
            call jeveuo(nume_equa//'.DELG', 'E', jrefn)
            call jelira(nume_equa//'.DELG', 'LONMAX', neq2)
        else
            call wkvect(nume_equa//'.DELG', 'G V I', neq, jrefn)
            neq2 = neq
        end if
        do iord = 1, neq2
            zi(jrefn-1+iord) = 0
        end do
    end if
!
    call jedema()
end subroutine
