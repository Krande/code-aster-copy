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

subroutine mdrecf(nexci, nexcir, idescf, nomfon, coefm, &
                  iadvec, inumor, fondep, fonvit, fonacc, &
                  neq, typbas, basemo, nbmode, riggen, &
                  nommot, nomres)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codlet.h"
#include "asterfort/copy_field_with_numbering.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeveut.h"
#include "asterfort/rsexch.h"
#include "asterfort/trmult.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
    integer(kind=8) :: nexci, nexcir, neq, nbmode
    integer(kind=8) :: idescf(*), inumor(*), iadvec(*)
    real(kind=8) :: coefm(*), riggen(nbmode)
    character(len=8) :: nomfon(2*nexci)
    character(len=8) :: fondep(2*nexci), fonvit(2*nexci)
    character(len=8) :: fonacc(2*nexci)
    character(len=8) :: basemo, nommot, nomres
    character(len=16) :: typbas
!
!     CALCULE LES FORCES EXTERIEURES A CHAQUE PAS DE TEMPS
!     ------------------------------------------------------------------
! IN  : NEXCI  : NOMBRE D'EXCITATIONS MODALES
! IN  : NEXCIR : NOMBRE D'EXC MODALES DONNEES SOUS FORME DE TRAN_GENE
! IN  : IDESCF  :
! IN  : NOMFON  : NOM DE LA FONCTION EXCITATION
! IN  : COEFM   :
! IN  : IADVEC  :
! IN  : INUMOR  :
! OUT : FONDEP  : NOM DE LA FONCTION DEPLACEMENT
! OUT : FONVIT  : NOM DE LA FONCTION VITESSE
! OUT : FONACC  : NOM DE LA FONCTION ACCELERATION
! OUT : PSDEL  : VECTEUR DES PSI*DELTA OU CORRECTIONS MODALES
! IN  : NEQ    : NOMBRE D'EQUATIONS
! IN  : BASEMO : NOM DE LA BASE MODALE DE PROJECTION
! OUT : NOMMOT :OUI SI MULTI APPUIS OU CORRECTION MODALE
! ----------------------------------------------------------------------
!
!
!
!
    integer(kind=8) :: i, ier, ninst, jprol, nexcit
    real(kind=8) :: alpha
    real(kind=8) :: coef
    character(len=2) :: ires
    character(len=5) :: imo
    character(len=8) :: modsta, modcor
    character(len=8) :: matass, mailla, monmot(2)
    character(len=14) :: numddl
    character(len=19) :: veasge, fonct, facce, tmpcha
    character(len=19) :: chamn2, chamno, nofk19, resu
    character(len=24) :: deeq, typeba, mesh
    integer(kind=8) :: jpsdel, npsdel, iipsdl
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iddeeq, ieq, ii, iinst
    integer(kind=8) ::  inum, iret, jdepl, jmod, jdesc
    integer(kind=8) :: jvale, l1, lprol, m1, n1, n2, n3
    integer(kind=8) :: n4, n5, na, nbv, nf, nm
    real(kind=8), pointer :: modco(:) => null()
    integer(kind=8), pointer :: ordr(:) => null()
    real(kind=8), pointer :: disc(:) => null()
    real(kind=8), pointer :: base(:) => null()
    real(kind=8), pointer :: mod(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    ier = 0
! ---    CALCUL TRANSITOIRE CLASSIQUE
    call dismoi('TYPE_BASE', basemo, 'RESU_DYNA', repk=typeba, arret='C', &
                ier=ier)
!
!
    call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=numddl)
    call dismoi('NOM_MAILLA', numddl, 'NUME_DDL', repk=mailla)
    deeq = numddl//'.NUME.DEEQ'
    call jeveuo(deeq, 'L', iddeeq)
    nommot = 'NON'
!
    nexcit = nexci+nexcir*nbmode
    npsdel = nexcit*neq
!
! --- EXCITATIONS SOUS LE MOT-CLE EXCIT
!
    do i = 1, nexci
!
        call getvis('EXCIT', 'NUME_ORDRE', iocc=i, scal=inum, nbret=nf)
        call getvid('EXCIT', 'VECT_ASSE_GENE', iocc=i, scal=veasge, nbret=l1)
        call getvid('EXCIT', 'FONC_MULT', iocc=i, scal=fonct, nbret=n1)
        call getvr8('EXCIT', 'COEF_MULT', iocc=i, scal=alpha, nbret=m1)
        call getvid('EXCIT', 'ACCE', iocc=i, scal=facce, nbret=na)
        call getvtx('EXCIT', 'MULT_APPUI', iocc=i, scal=monmot(1), nbret=n2)
        call getvtx('EXCIT', 'CORR_STAT', iocc=i, scal=monmot(2), nbret=n3)
!
        if (n1 .ne. 0) then
!         CAS D'UNE FONC_MULT
            nomfon(i) = fonct(1:8)
            call jeveuo(fonct//'.PROL', 'L', lprol)
            nomfon(i+nexcit) = zk24(lprol) (1:8)
            if (l1 .ne. 0) then
!           CAS D'UN VECT_ASSE_GENE
                call jeveut(veasge//'.VALE', 'L', jvale)
                iadvec(i) = jvale
                idescf(i) = 1
            else
!           CAS D'UN NUME_ORDRE
!           VERIF : LE NUME_ORDRE EST INFERIEUR AU NUME_ORDRE MAX
                if (inum .gt. neq) then
                    call utmess('F', 'ALGORITH5_76')
                end if
                inumor(i) = inum
                idescf(i) = 2
            end if
        else if (m1 .ne. 0) then
!         CAS D'UN COEF MULT
            coefm(i) = alpha
            if (l1 .ne. 0) then
!           CAS D'UN VECT_ASSE_GENE
                call jeveut(veasge//'.VALE', 'L', jvale)
                iadvec(i) = jvale
                idescf(i) = 3
            else
!           CAS D'UN NUME_ORDRE
                if (inum .gt. neq) then
                    call utmess('F', 'ALGORITH5_76')
                end if
                inumor(i) = inum
                idescf(i) = 4
            end if
        else if (na .ne. 0) then
!         CAS D'UN ACCELEROGRAMME
            nomfon(i) = facce(1:8)
            fonacc(i) = facce(1:8)
            call jeveuo(facce//'.PROL', 'L', lprol)
            nomfon(i+nexcit) = zk24(lprol) (1:8)
            fonacc(i+nexcit) = zk24(lprol) (1:8)
            if (l1 .ne. 0) then
!           CAS D'UN VECT_ASSE_GENE
                call jeveut(veasge//'.VALE', 'L', jvale)
                iadvec(i) = jvale
                idescf(i) = 1
            else
!           CAS D'UN NUME_ORDRE
                if (inum .gt. neq) then
                    call utmess('F', 'ALGORITH5_76')
                end if
                inumor(i) = inum
                idescf(i) = 2
            end if
        end if
        if (n2 .ne. 0) then
            if (monmot(1) .eq. 'OUI') then
                nommot = 'OUI'
                call jeexin(nomres//'           .IPSD', iipsdl)
                if (iipsdl .eq. 0) then
                    call wkvect(nomres//'           .IPSD', 'G V R8', npsdel, jpsdel)
                else
                    call jeveuo(nomres//'           .IPSD', 'E', jpsdel)
                end if
!               Add the identifier .DESC[6] to 1 , i.e. MULT_APPUI exists
                call jeexin(nomres//'           .DESC', iret)
                if (iret .eq. 0) then
                    call wkvect(nomres//'           .DESC', 'G V I', 6, jdesc)
                else
                    call jeveuo(nomres//'           .DESC', 'E', jdesc)
                end if
                zi(jdesc+6-1) = 1
!
                call getvid(' ', 'MODE_STAT', scal=modsta, nbret=nbv)
                if (nbv .eq. 0) then
                    ier = ier+1
                    call utmess('F', 'ALGORITH13_46')
                    goto 10
                end if

                call trmult(modsta, i, mailla, neq, iddeeq, &
                            zr(jpsdel+(i-1)*neq), numddl)
                call getvid('EXCIT', 'VITE', iocc=i, scal=fonvit(i), nbret=n4)
                fonct = fonvit(i)
                call jeveuo(fonct//'.PROL', 'L', lprol)
                fonvit(i+nexcit) = zk24(lprol) (1:8)
                call getvid('EXCIT', 'DEPL', iocc=i, scal=fondep(i), nbret=n5)
                fonct = fondep(i) (1:8)
                call jeveuo(fonct//'.PROL', 'L', lprol)
                fondep(i+nexcit) = zk24(lprol) (1:8)
            else
                ASSERT(.false.)
            end if
        end if
        if (n3 .ne. 0) then
            if (monmot(2) .eq. 'OUI') then
                nommot = 'OUI'
                call jeexin(nomres//'           .IPSD', iipsdl)
                if (iipsdl .eq. 0) then
                    call wkvect(nomres//'           .IPSD', 'G V R8', npsdel, jpsdel)
                else
                    call jeveuo(nomres//'           .IPSD', 'E', jpsdel)
                end if
                call getvid(' ', 'MODE_CORR', scal=modcor, nbret=nbv)
                if (nbv .eq. 0) then
                    ier = ier+1
                    call utmess('F', 'ALGORITH13_47')
                    goto 10
                end if
!

!               Add the identifier .DESC[5] to 1 , i.e. CORR_STAT exists
                call jeexin(nomres//'           .DESC', iret)
                if (iret .eq. 0) then
                    call wkvect(nomres//'           .DESC', 'G V I', 6, jdesc)
                else
                    call jeveuo(nomres//'           .DESC', 'E', jdesc)
                end if
                zi(jdesc+5-1) = 1
!
                call getvid('EXCIT', 'D_FONC_DT', iocc=i, scal=fonvit(i), nbret=n4)
                fonct = fonvit(i) (1:8)
                call jeveuo(fonct//'.PROL', 'L', lprol)
                fonvit(i+nexcit) = zk24(lprol) (1:8)
                call getvid('EXCIT', 'D_FONC_DT2', iocc=i, scal=fonacc(i), nbret=n5)
                fonct = fonacc(i) (1:8)
                call jeveuo(fonct//'.PROL', 'L', lprol)
                fonacc(i+nexcit) = zk24(lprol) (1:8)
                fondep(i) = nomfon(i) (1:8)
                fondep(i+nexcit) = nomfon(i+nexcit) (1:8)
!
                call rsexch('F', modcor, 'DEPL', i, chamno, iret)
                tmpcha = '&&COPMOD.CHAMP'
                call dismoi('NOM_MAILLA', numddl, 'NUME_DDL', repk=mesh)
                call copy_field_with_numbering(chamno, tmpcha, mesh, numddl//'.NUME', 'V')
                call jeveuo(tmpcha//'.VALE', 'L', vr=modco)
                do ieq = 1, neq
                    zr(jpsdel+ieq-1+(i-1)*neq) = modco(ieq)
                end do
                do nm = 1, nbmode
                    coef = zr(iadvec(i)+nm-1)/riggen(nm)
                    call rsexch('F', basemo, 'DEPL', nm, chamn2, iret)
                    call jeveuo(chamn2//'.VALE', 'L', vr=mod)
                    do ieq = 1, neq
                        zr(jpsdel+ieq-1+(i-1)*neq) = zr( &
                                                     jpsdel+ieq-1+(i-1)*neq &
                                                     )-coef*mod(1+ieq-1 &
                                                                )
                    end do
                    call jelibe(chamn2//'.VALE')
                end do

                call jelibe(chamno//'.VALE')
                call detrsd('CHAMP', tmpcha)
!
!           --- MISE A ZERO DES DDL DE LAGRANGE
                call zerlag(neq, zi(iddeeq), vectr=zr(jpsdel+(i-1)*neq))
!
            else
                ASSERT(.false.)
            end if
        end if
10      continue
    end do
!
!
! --- EXCITATIONS SOUS LE MOT-CLE EXCIT_RESU
!
! ON TRANSFORME CES EXCITATIONS SOUS LA FORME D'EXCTATIONS MODALES
! (MC NUME_ORDRE), ASSOCIEES A UNE FONCTION MULTIPLICATRICE
!
    call jeveuo(basemo//'           .ORDR', 'L', vi=ordr)

    !if (nbmode .gt. 999) then
    !   call utmess('F', 'ALGORITH5_77')
    !endif

    do i = 1, nexcir
!
        call getvid('EXCIT_RESU', 'RESULTAT', iocc=i, scal=resu, nbret=l1)
        call getvr8('EXCIT_RESU', 'COEF_MULT', iocc=i, scal=alpha, nbret=m1)
! ----- NOMBRE DE PAS DE TEMPS DU RESULTAT
        call jelira(resu//'.DISC', 'LONMAX', ninst)
        call jeveuo(resu//'.DISC', 'L', vr=disc)
! ----- EXCITATION STOCKEE DANS LE CHAMP DEPL (IDEM QUE DYNA_VIBRA//TRAN/PHYS)
        call jeveuo(resu//'.DEPL', 'L', jdepl)
!
        ii = nexci+nbmode*(i-1)
        call codlet(i, 'D0', ires)
!
        do jmod = 1, nbmode
            call codlet(jmod, 'D0', imo)
            nomfon(ii+jmod) = '_'//ires//imo
            nofk19 = nomfon(ii+jmod)
            call wkvect(nofk19//'.VALE', 'V V R8', 2*ninst, jvale)
            do iinst = 1, ninst
                zr(jvale-1+iinst) = disc(iinst)
                zr(jvale-1+ninst+iinst) = alpha*zr(jdepl-1+nbmode*(iinst-1)+jmod)
            end do
!
            call wkvect(nofk19//'.PROL', 'V V K24', 6, jprol)
            zk24(jprol-1+1) = 'FONCTION'
            zk24(jprol-1+2) = 'LIN LIN'
            zk24(jprol-1+3) = 'INST'
            zk24(jprol-1+4) = 'TOUTRESU'
            zk24(jprol-1+5) = 'EE'
            zk24(jprol-1+6) = nomfon(ii)
!
            nomfon(nexcit+ii+jmod) = zk24(jprol) (1:8)
!
            inumor(ii+jmod) = ordr(jmod)
            idescf(ii+jmod) = 2
!
        end do
!
    end do
!
    call jedema()
end subroutine
