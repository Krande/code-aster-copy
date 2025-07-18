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
subroutine op0153()
    implicit none
!
!     OPERATEUR  "POST_USURE"
!
! ----------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/r8prem.h"
#include "asterfort/detrsd.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/motubn.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/tbexp2.h"
#include "asterfort/tbextb.h"
#include "asterfort/tbexv1.h"
#include "asterfort/tbexve.h"
#include "asterfort/tbliva.h"
#include "asterfort/titre.h"
#include "asterfort/usupru.h"
#include "asterfort/usupus.h"
#include "asterfort/usuvu2.h"
#include "asterfort/usuvus.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, idangt, idvcob, idvctu, ifires
    integer(kind=8) :: indic, info, iobst, ipoupr, ipourp, iprfuo, iprfut
    integer(kind=8) :: ipus, ire1, ire2, iret, itube, ivuso
    integer(kind=8) :: ivusob, ivust, ivustu, jfn, jins2, jinst
    integer(kind=8) :: jprut, jsect, jusuo, jusut, jvg, k
    integer(kind=8) :: n0, n1, n5, na, nbinst, nbpair, nbpar
    integer(kind=8) :: nbpar2, nbpmr, nbpt, nbsec2, nbsect, nbtota, nbv
    integer(kind=8) :: nbvpu, ni1, nis, npu, ntn
    real(kind=8) :: coinst, haut, puusur, rayoo, rayot
!
!-----------------------------------------------------------------------
    parameter(nbpar=16, nbpar2=12, nbpmr=5)
    real(kind=8) :: pmoye, insdeb, epsil, dinst, dinst0
    character(len=8) :: k8b
    character(len=8) :: typar(nbpar), typpmr(nbpmr)
    real(kind=8) :: valer(nbpar)
    character(len=16) :: nopar(nbpar), nompmr(nbpmr), nopar2(nbpar2)
    character(len=16) :: concep, nomcmd, valek(2)
    character(len=19) :: resu, linst, kforn, kvgli
    character(len=19) :: tabpus, nomta, newtab
    character(len=24) :: type, valk(2)
    complex(kind=8) :: c16b
    real(kind=8), pointer :: ins3(:) => null()
    real(kind=8), pointer :: ins5(:) => null()
    data nopar/'PUIS_USUR_GLOBAL',&
     &              'INST', 'DUREE', 'ORIG_INST',&
     &              'V_USUR_TUBE', 'V_USUR_OBST', 'P_USUR_TUBE',&
     &              'SECTEUR', 'ANGLE_DEBUT', 'ANGLE_FIN',&
     &              'V_USUR_TUBE_SECT', 'V_USUR_OBST_SECT',&
     &              'P_USUR_TUBE_SECT', 'P_USUR_OBST_SECT',&
     &              'V_USUR_TUBE_CUMU', 'V_USUR_OBST_CUMU'/
    data typar/'R', 'R', 'R', 'R', 'R', 'R', 'R', 'I', 'R', 'R', 'R', 'R', 'R',&
     &             'R', 'R', 'R'/
    data nopar2/'INST', 'DUREE', 'ORIG_INST',&
     &              'SECTEUR', 'ANGLE_DEBUT', 'ANGLE_FIN',&
     &              'V_USUR_TUBE_SECT', 'V_USUR_OBST_SECT',&
     &              'P_USUR_TUBE_SECT', 'P_USUR_OBST_SECT',&
     &              'V_USUR_TUBE_CUMU', 'V_USUR_OBST_CUMU'/
    data nompmr/'PUIS_USUR_GLOBAL',&
     &             'INST', 'V_USUR_TUBE', 'V_USUR_OBST', 'P_USUR_TUBE'/
    data typpmr/'R', 'R', 'R', 'R', 'R'/
!     ------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
!
    indic = 0
    ifires = iunifi('RESULTAT')
    kforn = '&&OP0153.FORC_N'
    kvgli = '&&OP0153.VITE_G'
!
    call getres(resu, concep, nomcmd)
!
!     ------------------------------------------------------------------
!            REMPLACEMENT DU TUBE PERCE PAR UN TUBE NEUF
!     ------------------------------------------------------------------
!
    dinst = 0.d0
    call getvtx(' ', 'TUBE_NEUF', scal=k8b, nbret=ntn)
    if (ntn .ne. 0) then
        call exisd('TABLE', resu, iret)
        if (iret .eq. 0) then
            call utmess('F', 'PREPOST4_7')
        end if
        call getvid(' ', 'TABL_USURE', scal=k8b, nbret=n1)
        if (k8b .ne. resu(1:8)) then
            call utmess('F', 'PREPOST4_7')
        end if
        call tbexv1(resu, 'INST', '&&OP0153.INST', 'V', nbv, &
                    k8b)
        call jeveuo('&&OP0153.INST', 'L', jinst)
        call getvr8(' ', 'INST', scal=dinst, nbret=nis)
        if (nis .eq. 0) then
            dinst = zr(jinst+nbv-1)
        end if
        call tbexv1(resu, 'SECTEUR', '&&OP0153.SECT', 'V', nbv, &
                    k8b)
        call jeveuo('&&OP0153.SECT', 'L', jsect)
        nbsect = zi(jsect+nbv-1)
        call jedetr('&&OP0153.SECT')
        call jelira('&&OP0153.INST', 'LONUTI', nbinst)
        do i = 0, nbinst-1
            dinst0 = zr(jinst+i)
            if (dinst0 .ge. dinst) then
                call motubn(resu, dinst0, nbsect)
            end if
        end do
        goto 888
    end if
!
!     ------------------------------------------------------------------
    call getvis(' ', 'INFO', scal=info, nbret=n0)
    if (info .gt. 1) then
        write (ifires, 100)
        write (ifires, *)
        write (ifires, *) resu
    end if
!
!     --- CALCUL DE LA PUISSANCE D'USURE ---
    call usupus(puusur, kforn, kvgli, nbpt)
    call jeexin(kforn, iret)
    jfn = 1
    jvg = 1
    if (iret .gt. 0) then
        call jeveuo(kforn, 'E', jfn)
        call jeveuo(kvgli, 'E', jvg)
    end if
!
!     --- RECUPERATION DES INSTANTS DE CALCUL ---
    call getvr8(' ', 'INST', nbval=0, nbret=ni1)
    if (ni1 .ne. 0) then
        nbinst = -ni1
        call wkvect('&&OP0153.INSTANT', 'V V R', nbinst, jinst)
        call getvr8(' ', 'INST', nbval=nbinst, vect=zr(jinst), nbret=n1)
    else
        call getvid(' ', 'LIST_INST', scal=linst, nbret=n1)
        call jelira(linst//'.VALE', 'LONUTI', nbinst)
        call jeveuo(linst//'.VALE', 'L', jinst)
    end if
    call wkvect('&&OP0153.INSTAN2', 'V V R', nbinst, jins2)
    do i = 0, nbinst-1
        zr(jins2+i) = zr(jinst+i)
    end do
    call getvr8(' ', 'COEF_INST', scal=coinst, nbret=n1)
    if (n1 .ne. 0) then
        do i = 0, nbinst-1
            zr(jins2+i) = zr(jins2+i)*coinst
        end do
    end if
!
    call wkvect('&&OP0153.USURE_TUBE', 'V V R', nbinst, jusut)
    call wkvect('&&OP0153.USURE_OBST', 'V V R', nbinst, jusuo)
    call wkvect('&&OP0153.PRONF_TUBE', 'V V R', nbinst, jprut)
!
    call getfac('SECTEUR', nbsect)
    if (nbsect .ne. 0) then
        nbpair = nbsect+1
        nbtota = nbsect*nbinst
        call wkvect('&&OP0153.ANGT', 'V V R', nbpair, idangt)
        call wkvect('&&OP0153.VUSTUB', 'V V R', nbtota, ivustu)
        call wkvect('&&OP0153.VUSOB', 'V V R', nbtota, ivusob)
        call wkvect('&&OP0153.PRFUST', 'V V R', nbtota, iprfut)
        call wkvect('&&OP0153.PRFUSO', 'V V R', nbtota, iprfuo)
        call wkvect('&&OP0153.VUST', 'V V R', nbsect, ivust)
        call wkvect('&&OP0153.VUSO', 'V V R', nbsect, ivuso)
        call wkvect('&&OP0153.PUS', 'V V R', nbsect, ipus)
        call wkvect('&&OP0153.POUPRE', 'V V R', nbsect, ipoupr)
        call wkvect('&&OP0153.POURPU', 'V V R', nbsect, ipourp)
        call wkvect('&&OP0153.VCTU', 'V V R', nbsect, idvctu)
        call wkvect('&&OP0153.VCOB', 'V V R', nbsect, idvcob)
        epsil = 1.d-4
        do i = 1, nbsect
            if (i .eq. 1) then
                call getvr8('SECTEUR', 'ANGL_INIT', iocc=1, scal=zr(idangt), nbret=na)
!
!              LES ANGLES SONT CROISSANTS ENTRE -180. ET +180. :
!              -----------------------------------------------
                if ((zr(idangt) .lt. (-180.d0-epsil)) .or. (zr(idangt) .gt. (-180.d0+epsil))) then
                    call utmess('F', 'PREPOST4_8')
                end if
            end if
            call getvr8('SECTEUR', 'ANGL_FIN', iocc=i, scal=zr(idangt+i), nbret=na)
            if (zr(idangt+i) .lt. zr(idangt+i-1)) then
                call utmess('F', 'PREPOST4_9')
            end if
            if (i .eq. nbsect) then
                if ((zr(idangt+i) .lt. (180.d0-epsil)) .or. (zr(idangt+i) .gt. (180.d0+epsil))) then
                    call utmess('F', 'PREPOST4_10')
                end if
            end if
            call getvr8('SECTEUR', 'COEF_USUR_MOBILE', iocc=i, scal=zr(idvctu+i-1), nbret=n5)
            call getvr8('SECTEUR', 'COEF_USUR_OBST', iocc=i, scal=zr(idvcob+i-1), nbret=n5)
        end do
    else
        indic = 1
    end if
!
    if (puusur .le. r8prem()) goto 777
!
!     --- CALCUL DU VOLUME D'USURE TUBE ---
    itube = 1
!
    if (indic .eq. 0) then
        call usuvu2(puusur, zr(jusut), nbinst, zr(jins2), itube, &
                    nbpt, nbsect, zr(idvctu), zr(idangt), zr(jfn), &
                    zr(jvg), iret, zr(ivustu), zr(ivusob), zr(ipus), &
                    pmoye, zr(ipourp), zr(ipoupr))
    else if (indic .eq. 1) then
        call usuvus(puusur, zr(jusut), nbinst, zr(jins2), itube, &
                    nbpt, zr(jfn), zr(jvg), iret)
    end if
    if (iret .ne. 0) goto 999
!
!     --- CALCUL DU VOLUME D'USURE OBSTABLE ---
    iobst = 2
!
    if (indic .eq. 0) then
        call usuvu2(puusur, zr(jusuo), nbinst, zr(jins2), iobst, &
                    nbpt, nbsect, zr(idvcob), zr(idangt), zr(jfn), &
                    zr(jvg), iret, zr(ivustu), zr(ivusob), zr(ipus), &
                    pmoye, zr(ipourp), zr(ipoupr))
    else if (indic .eq. 1) then
        call usuvus(puusur, zr(jusuo), nbinst, zr(jins2), iobst, &
                    nbpt, zr(jfn), zr(jvg), iret)
    end if
    if (iret .ne. 0) goto 999
!
    if (indic .eq. 0) then
!
        call getvr8(' ', 'ANGL_ISTHME', scal=haut, nbret=n1)
        if (n1 .ne. 0) call utmess('A', 'PREPOST4_3', sk='ANGL_ISTHME')
        call getvr8(' ', 'ANGL_INCLI', scal=haut, nbret=n1)
        if (n1 .ne. 0) call utmess('A', 'PREPOST4_3', sk='ANGL_INCLI')
        call getvr8(' ', 'ANGL_IMPACT', scal=haut, nbret=n1)
        if (n1 .ne. 0) call utmess('A', 'PREPOST4_3', sk='ANGL_IMPACT')
!
        call getvr8(' ', 'LARGEUR_OBST', scal=haut, nbret=n1)
        if (n1 .le. 0) then
            haut = 0.011d0
        end if
        if (info .gt. 1) write (ifires, 200)
        call getvr8(' ', 'RAYON_MOBILE', scal=rayot, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('F', 'PREPOST4_11')
        end if
        do i = 1, nbsect
            do k = 1, nbinst
                if (rayot*rayot-2.d0*zr( &
                    ivustu+(k-1)*nbsect+i-1)/(haut*(zr(idangt+i)-zr(idangt+i-1))) &
                    .lt. 0.d0) then
                    call utmess('F', 'PREPOST4_4')
                end if
!
                zr(iprfut+(k-1)*nbsect+i-1) = rayot-sqrt( &
                                              rayot*rayot-2.d0*zr( &
                                              ivustu+(k-1)*nbsect+i-1)/(haut*(zr(idangt+i)-zr(i&
                                              &dangt+i-1)) &
                                              ) &
                                              )
            end do
        end do
        call getvr8(' ', 'RAYON_OBST', scal=rayoo, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('F', 'PREPOST4_12')
        end if
        do i = 1, nbsect
            do k = 1, nbinst
                zr(iprfuo+(k-1)*nbsect+i-1) = rayoo-sqrt( &
                                              rayoo*rayoo-2.d0*zr( &
                                              ivusob+(k-1)*nbsect+i-1)/(haut*(zr(idangt+i)-zr(i&
                                              &dangt+i-1)) &
                                              ) &
                                              )
            end do
        end do
    end if
!
!      --- IMPRESSIONS DES RESULTATS ---
!
!
    if (indic .ne. 0) goto 666
    do i = 1, nbsect
        if (info .gt. 1) then
            write (ifires, *)
            write (ifires, *)
            write (ifires, 190) 'SECTEUR : ', zr(idangt+i-1), ' / ', &
                zr(idangt+i)
            write (ifires, *)
            write (ifires, 140) 'COEF USURE TUBE     ', ':',&
     &    zr(idvctu+i-1)
            write (ifires, 140) 'COEF USURE OBSTACLE ', ':', zr(idvcob+ &
                                                                i-1)
            write (ifires, 130) 'PRESENCE DU CRAYON  ', ':',&
     &                       zr(ipoupr+i-1)*100.d0, '%'
            write (ifires, 160) 'PUISSANCE D USURE   ', ':',&
     &    zr(ipus+i-1), 'W'
            write (ifires, 130) '% PU DANS CE SECTEUR', ':', zr(ipourp+ &
                                                                i-1), '%'
            write (ifires, *)
            write (ifires, 110) 'ANNEES', 'V_USUR_TUBE', 'V_USUR_OBST', &
                'P_USUR_TUBE', 'P_USUR_OBST'
        end if
        do k = 1, nbinst
            if (info .gt. 1) write (ifires, 180) (zr(jins2+k-1)/coinst), &
                zr(ivustu+(k-1)*nbsect+i-1), zr(ivusob+(k-1)*nbsect+i-1), &
                zr(iprfut+(k-1)*nbsect+i-1), zr(iprfuo+(k-1)*nbsect+i-1)
        end do
    end do
666 continue
!     --- CALCUL DE PROFONDEUR D'USURE ---
    if (indic .eq. 1) call usupru(zr(jusut), zr(jusuo), nbinst, zr(jprut))
777 continue
!
    if (indic .eq. 1) then
!        --- CREATION DE LA TABLE ---
        call tbcrsd(resu, 'G')
        call tbajpa(resu, nbpmr, nompmr, typpmr)
        call tbajli(resu, 1, 'PUIS_USUR_GLOBAL', [ibid], [puusur], &
                    [c16b], k8b, 0)
        do k = 1, nbinst
            valer(1) = zr(jins2+k-1)/coinst
            valer(2) = zr(jusut+k-1)
            valer(3) = zr(jusuo+k-1)
            valer(4) = zr(jprut+k-1)
            call tbajli(resu, 4, nompmr(2), [ibid], valer, &
                        [c16b], k8b, 0)
        end do
        goto 888
    end if
!
!     REPRISE EVENTUELLE ET STOCKAGE DE LA TABLE POST_USURE :
!     -----------------------------------------------------
!
    call getvid('ETAT_INIT', 'TABL_USURE', iocc=1, scal=tabpus, nbret=npu)
    if (npu .eq. 0) then
        dinst = 0.d0
        call tbcrsd(resu, 'G')
        call tbajpa(resu, nbpar, nopar, typar)
    else
        if (tabpus .ne. resu) then
            call utmess('F', 'PREPOST4_13')
        end if
!   ON REPREND UNE TABLE EXISTANTE
        nomta = tabpus
        call tbexp2(nomta, 'INST')
        call tbexp2(nomta, 'SECTEUR')
        call tbexp2(nomta, 'V_USUR_OBST_CUMU')
        call tbexp2(nomta, 'V_USUR_TUBE_CUMU')
        call tbexve(nomta, 'INST', '&&OP0153.INS3', 'V', nbvpu, &
                    type)
        call jeveuo('&&OP0153.INS3', 'L', vr=ins3)
        insdeb = ins3(nbvpu)
        call tbexve(nomta, 'SECTEUR', '&&OP0153.SECT', 'V', nbvpu, &
                    type)
        call jeveuo('&&OP0153.SECT', 'L', jsect)
        nbsec2 = zi(jsect+nbvpu-1)
        if (nbsec2 .ne. nbsect) then
            call utmess('F', 'PREPOST4_14')
        end if
        call getvr8('ETAT_INIT', 'INST_INIT', iocc=1, scal=dinst, nbret=nis)
        if (nis .eq. 0) then
            dinst = insdeb
        else if (dinst .gt. insdeb) then
            dinst = insdeb
        else
            newtab = '&&OP0153.NEWTAB'
            call tbextb(nomta, 'V', newtab, 1, 'INST', &
                        'LE', [ibid], [dinst], [c16b], k8b, &
                        [1.d-03], 'RELA', iret)
            if (iret .eq. 10) then
                valk(1) = 'INST'
                valk(2) = nomta
                call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
            else if (iret .eq. 20) then
                valk(1) = nomta
                valk(2) = 'INST'
                call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
            end if
            call detrsd('TABLE', nomta)
            call tbextb(newtab, 'G', nomta, 1, 'INST', &
                        'LE', [ibid], [dinst], [c16b], k8b, &
                        [1.d-03], 'RELA', iret)
            if (iret .eq. 10) then
                valk(1) = 'INST'
                valk(2) = newtab
                call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
            else if (iret .eq. 20) then
                valk(1) = newtab
                valk(2) = 'INST'
                call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
            end if
            call tbexve(nomta, 'INST', '&&OP0153.INS5', 'V', nbvpu, &
                        type)
            call jeveuo('&&OP0153.INS5', 'L', vr=ins5)
            dinst = ins5(nbvpu)
        end if
!
!        DETERMINATION PAR SECTEUR DES VOLUS PAR TUBE ET OBST A DINST
!        ------------------------------------------------------------
!
        valek(1) = 'INST'
        valek(2) = 'SECTEUR'
        do i = 1, nbsect
            call tbliva(nomta, 2, valek, [i], [dinst], &
                        [c16b], k8b, 'RELA', [1.d-03], 'V_USUR_TUBE_CUMU', &
                        k8b, ibid, zr(ivust+i-1), c16b, k8b, &
                        ire1)
            call tbliva(nomta, 2, valek, [i], [dinst], &
                        [c16b], k8b, 'RELA', [1.d-03], 'V_USUR_OBST_CUMU', &
                        k8b, ibid, zr(ivuso+i-1), c16b, k8b, &
                        ire2)
            if ((ire1+ire2) .gt. 0) then
                call utmess('F', 'PREPOST4_15', sk=nomta)
            end if
        end do
    end if
!
    call tbajli(resu, 1, 'PUIS_USUR_GLOBAL', [ibid], [puusur], &
                [c16b], k8b, 0)
!
    do k = 1, nbinst
!        -INST-
        valer(1) = zr(jins2+k-1)/coinst+dinst
!        -DUREE-
        valer(2) = zr(jins2+k-1)/coinst
!        -ORIG_INST-
        valer(3) = dinst
!        -V_USUR_TUBE-
        valer(4) = zr(jusut+k-1)
!        -V_USUR_OBST-
        valer(5) = zr(jusuo+k-1)
!        -P_USUR_TUBE-
        valer(6) = zr(jprut+k-1)
        call tbajli(resu, 6, nopar(2), [ibid], valer, &
                    [c16b], k8b, 0)
        do i = 1, nbsect
!           -ANGLE_DEBUT-
            valer(4) = zr(idangt+i-1)
!           -ANGLE_FIN-
            valer(5) = zr(idangt+i)
!           -V_USUR_TUBE_SECT-
            valer(6) = zr(ivustu+(k-1)*nbsect+i-1)
!           -V_USUR_OBST_SECT-
            valer(7) = zr(ivusob+(k-1)*nbsect+i-1)
!           -P_USUR_TUBE_SECT-
            valer(8) = zr(iprfut+(k-1)*nbsect+i-1)
!           -P_USUR_OBST_SECT-
            valer(9) = zr(iprfuo+(k-1)*nbsect+i-1)
!           -V_USUR_TUBE_CUMU-
            valer(10) = zr(ivustu+(k-1)*nbsect+i-1)+zr(ivust+i-1)
!           -V_USUR_OBST_CUMU-
            valer(11) = zr(ivusob+(k-1)*nbsect+i-1)+zr(ivuso+i-1)
            call tbajli(resu, nbpar2, nopar2, [i], valer, &
                        [c16b], k8b, 0)
        end do
    end do
    if (nbsect .ne. 0 .and. info .gt. 1) then
        write (ifires, *)
        write (ifires, *) 'PUISSANCE D USURE MOYENNE'
        write (ifires, 120) pmoye, 'W'
    end if
!
888 continue
!
    call titre()
!
100 format(/, 80('-'))
110 format(a11, 2x, a15, 2x, a15, 2x, a15, 2x, a15)
120 format(1pe12.5, 1x, a1)
130 format(a20, 1x, a1, 1x, f6.2, 1x, a1)
140 format(a20, 1x, a1, 1x, 1pe11.4)
160 format(a20, 1x, a1, 1x, 1pe12.5, 1x, a1)
180 format(1pe12.5, 2x, 1pe16.9, 2x, 1pe16.9, 2x, 1pe16.9, 2x, 1pe16.9)
190 format(a10, 1x, f7.2, a3, f7.2)
200 format(&
     &'LES PROFONDEURS USEES PAR SECTEUR SONT DES APPROXIMATIONS')
!
999 continue
    call jedema()
end subroutine
