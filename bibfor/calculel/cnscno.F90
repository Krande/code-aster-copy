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

subroutine cnscno(cnsz, nume_equaz, prol0, basez, cnoz, &
                  kstop, iret, nbz, vchamz, lprofconst, prolong)
!
! aslint: disable=
    use proj_champ_module
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/cmpcha.h"
#include "asterfort/codent.h"
#include "asterfort/nume_equa_crsd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/copisd.h"
#include "asterfort/gcncon.h"
#include "asterfort/gnomsd.h"
#include "asterfort/idensd.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/pteequ.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcreb.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*) :: cnsz, cnoz, basez, nume_equaz, prol0
    character(len=1) :: kstop
    integer(kind=8), optional :: nbz
    character(len=24), optional :: vchamz
    aster_logical, optional :: lprofconst
    type(prolongation), optional :: prolong
!
! ------------------------------------------------------------------
! BUT : TRANSFORMER UN CHAM_NO_S (CNSZ) EN CHAM_NO (CNOZ)
! ------------------------------------------------------------------
!     ARGUMENTS:
! CNSZ    IN/JXIN  K19 : SD CHAM_NO_S A TRANSFORMER
! NUME_EQUAZ  IN/JXVAR K19 : SD NUME_EQUA  (OU ' ')
!          SI NUME_EQUAZ EXISTE ON CREE CNOZ CONFORMEMENT A NUME_EQUAZ :
!             => SI CNSZ CONTIENT DES VALEURS QUE L'ON NE SAIT PAS
!                STOCKER DANS NUME_EQUAZ, ON LES "OUBLIE"
!             => SI NUME_EQUAZ EXIGE DES VALEURS QUE L'ON NE TROUVE PAS
!                DANS CNSZ :
!                  - SI PROL0='OUI' : ON PRENDS LA VALEUR "ZERO"
!                  - SI PROL0='NON' : ERREUR <F>
!
!          SI NUME_EQUAZ N'EXISTE PAS ON CREE CNOZ EN FONCTION
!             DU CONTENU DE CNSZ
!             SI NUME_EQUAZ  = ' ' ON CREE UN NUME_EQUA "SOUS-TERRAIN"
!             SI NUME_EQUAZ /= ' ' ON CREE UN NUME_EQUA DE NOM NUME_EQUAZ
! PROL0   IN   K3  :  POUR PROLONGER (OU NON) LE CHAMP PAR "ZERO"
!        /OUI /NON  ( CET ARGUMENT N'EST UTILISE QUE SI NUME_EQUAZ /= ' ')
!        "ZERO" : / 0       POUR LES CHAMPS NUMERIQUES (R/C/I)
!                 / ' '     POUR LES CHAMPS "KN"
!                 / .FALSE. POUR LES CHAMPS DE "L"
!
! BASEZ   IN       K1  : BASE DE CREATION POUR CNOZ : G/V/L
! CNOZ    IN/JXOUT K19 : SD CHAM_NO A CREER
! KSTOP   IN       K1  : COMPORTEMENT EN CAS DE PROBLEME :
!              / 'A' : ON EMET UNE ALARME ET ON REND IRET > 0
!              / 'F' : ON EMET UNE ERREUR FATALE
!              / ' ' : ON N'EMET PAS DE MESSAGE
! IRET    OUT       I  : CODE DE RETOUR :
!              / 0 : OK
!              / 1 : LE CHAM_NO N'A PAS PU ETRE CREE
!----------------------------------------------------------------------
!
    character(len=24) :: noojb
    character(len=24) :: valk(3)
!     -----------------------------------------------------------------
    integer(kind=8) :: icmp, nec, jcnsv, jcnsl, gd, iexi, ncmp, jvcham, jrefk
    integer(kind=8) :: reste, iec, code, nbno, nb
    integer(kind=8) :: ncmpmx, jrefe, ncmp1, nb_equa, jcmpgd, icmp1, k, ieq2, nbec
    integer(kind=8) :: jprn2, ino, idg2, ico, jvale, iret
    integer(kind=8) :: lshift, nuprf, nb_equa_gl, jrefn
    character(len=1) :: base
    character(len=8) :: ma, nomgd, nomno, nomcmp, prolv
    aster_logical :: l_crea_nume_equa, l_chck_nume_equa, ldist, l_pmesh
    aster_logical :: lprof_const
    character(len=3) :: tsca
    character(len=19) :: cns, cno, nume_equa, messag, prnoav, nume_equa_tmp
    integer(kind=8), pointer :: deeq(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
    character(len=8), pointer :: cnsc(:) => null()
    character(len=24) :: vcham
    character(len=24), pointer :: refn(:) => null()
    character(len=8), pointer :: cnsk(:) => null()
    integer(kind=8), pointer :: tmp_nucm1(:) => null()
    integer(kind=8), pointer :: tmp_nucmp(:) => null()
    integer(kind=8), pointer :: cata_to_field(:) => null()
    integer(kind=8), pointer :: field_to_cata(:) => null()
    character(len=8), pointer :: cmp_name(:) => null()
    integer(kind=8), pointer :: v_nequ(:) => null()
!     -----------------------------------------------------------------
    call jemarq()
!   EVENTUEL PARALLELISME EN TEMPS
    if (present(nbz) .and. present(vchamz)) then
        ldist = .True.
        ASSERT(nbz .ge. 2)
        nb = nbz
        vcham = vchamz
    else
        ldist = .False.
        nb = 0
        vcham = ' '
    end if
!
!   le profchno est constant par défaut
    lprof_const = ASTER_TRUE
    if (present(lprofconst)) then
        lprof_const = lprofconst
    end if
!
    base = basez
    ASSERT((base .eq. 'G') .or. (base .eq. 'V'))
    cns = cnsz
    cno = cnoz
!
    l_chck_nume_equa = .false.
    l_crea_nume_equa = .false.
!     CALL UTIMSD(6,2,.TRUE.,.TRUE.,CNS,1,' ')
!
    call jeveuo(cns//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(cns//'.CNSD', 'L', vi=cnsd)
    call jeveuo(cns//'.CNSC', 'L', vk8=cnsc)
    call jeveuo(cns//'.CNSV', 'L', jcnsv)
    call jeveuo(cns//'.CNSL', 'L', jcnsl)
!
    ma = cnsk(1)
    nomgd = cnsk(2)
    nbno = cnsd(1)
    ncmp1 = cnsd(2)

    l_pmesh = isParallelMesh(ma)
!
    call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nec)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    call dismoi('NB_CMP_MAX', nomgd, 'GRANDEUR', repi=ncmpmx)
    call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nbec)
    call dismoi('NUM_GD', nomgd, 'GRANDEUR', repi=gd)
!
!   Si 'prolong' est donné, 'prol0' ne sert pas
    prolv = '???'
    if (present(prolong)) then
        prol0 = '???'
        ! si prolong%prol_vale_r sur tsca : R C
        !   Si C ==> c'est seulement 0
        if (prolong%prol_vale_r .eq. 'OUI') then
            ASSERT((tsca .eq. 'R') .or. (tsca .eq. 'C'))
        end if
        prolv = prolong%prol_vale_r
    end if
!   SI CNO EXISTE DEJA, ON LE DETRUIT :
    call detrsd('CHAM_NO', cno)
!
! - NUME_EQUA name
!
    if (nume_equaz .eq. ' ') then
        if (base .eq. 'G') then
            noojb = '12345678.NUME000000.PRNO'
            call gnomsd(' ', noojb, 14, 19)
            noojb(1:8) = cno(1:8)
            nume_equa = noojb(1:19)
        else
            call gcncon('.', nume_equa)
        end if
        l_chck_nume_equa = .false.
    else
        nume_equa = nume_equaz
        l_chck_nume_equa = .true.
    end if
!
! - Create NUME_EQUA ?
!
    call jeexin(nume_equa//'.PRNO', iexi)
    l_crea_nume_equa = (iexi .eq. 0)
    if (l_crea_nume_equa) then
        l_chck_nume_equa = .false.
    end if
!
!
!   1- REMPLISSAGE DE .TMP_NUCMP ET .TMP_NUCM1 :
!   --------------------------------------------
    AS_ALLOCATE(vi=tmp_nucmp, size=ncmpmx)
    if (ncmp1 .ne. 0) then
        AS_ALLOCATE(vi=tmp_nucm1, size=ncmp1)
    end if
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', jcmpgd)
    do icmp1 = 1, ncmp1
        nomcmp = cnsc(icmp1)
        icmp = indik8(zk8(jcmpgd), nomcmp, 1, ncmpmx)
        ASSERT(icmp .gt. 0)
        tmp_nucmp(icmp) = icmp1
        tmp_nucm1(icmp1) = icmp
    end do
!
! - Check NUME_EQUA
!
    if (l_chck_nume_equa) then
!       ON PEUT VERIFIER maillage et grandeur
        call jeveuo(nume_equa(1:19)//'.REFN', 'L', vk24=refn)
        if ((refn(1) .ne. ma) .or. (refn(2) .ne. nomgd)) then
            if (nomgd(1:5) == refn(2) (1:5)) then
                if (base .eq. 'G') then
                    noojb = '12345678.NUME000000.PRNO'
                    call gnomsd(' ', noojb, 14, 19)
                    noojb(1:8) = cno(1:8)
                    nume_equa_tmp = noojb(1:19)
                else
                    call gcncon('.', nume_equa_tmp)
                end if
                call copisd("NUME_EQUA", base, nume_equa, nume_equa_tmp)
                nume_equa = nume_equa_tmp
                call jeveuo(nume_equa//".REFN", 'E', jrefn)
                zk24(jrefn+1) = nomgd
            else
                valk(1) = cno
                valk(2) = nume_equa
                call utmess('F', 'CALCULEL4_6', nk=2, valk=valk)
            end if
        end if
    end if
!
!
!     2- ON CREE (SI NECESSAIRE) LE NUME_EQUA  :
!     ------------------------------------------
    if (l_crea_nume_equa) then
!
!       2.1 ON COMPTE LES CMPS PORTEES PAR CNS :
        nb_equa = 0
        do k = 1, nbno*ncmp1
            if (zl(jcnsl-1+k)) then
                nb_equa = nb_equa+1
            else if (prol0 .eq. 'OUI') then
                zl(jcnsl-1+k) = .true.
                if (tsca .eq. 'R') then
                    zr(jcnsv-1+k) = 0.d0
                else if (tsca .eq. 'C') then
                    zc(jcnsv-1+k) = (0.d0, 0.d0)
                else if (tsca .eq. 'I') then
                    zi(jcnsv-1+k) = 0
                else if (tsca .eq. 'L') then
                    zl(jcnsv-1+k) = .false.
                else if (tsca .eq. 'K8') then
                    zk8(jcnsv-1+k) = ' '
                else
                    ASSERT(.false.)
                end if
                nb_equa = nb_equa+1
            else if (prolv .eq. 'OUI') then
                zl(jcnsl-1+k) = .true.
                if (tsca .eq. 'R') then
                    zr(jcnsv-1+k) = prolong%vale_r
                else if (tsca .eq. 'C') then
                    zc(jcnsv-1+k) = (0.d0, 0.d0)
                end if
                nb_equa = nb_equa+1
            end if
        end do
!
        nb_equa_gl = nb_equa
        if (l_pmesh) then
            call asmpi_comm_vect("MPI_SUM", "I", sci=nb_equa_gl)
        end if
!
        if (nb_equa_gl .eq. 0) then
            valk(1) = cns
            valk(2) = cno
            messag = 'CALCULEL2_12'
            goto 70
!
        end if
!
!       2.2 ALLOCATION DES OBJETS :
        call nume_equa_crsd(nume_equa, base, nb_equa, meshz=ma, &
                            gran_namez=nomgd, l_coll_constz=lprof_const)
        call jecroc(jexnum(nume_equa//'.PRNO', 1))
!
!       2.3 REMPLISSAGE DE .PRNO :
        call jeveuo(jexnum(nume_equa//'.PRNO', 1), 'E', jprn2)
        do ino = 1, nbno
            do icmp1 = 1, ncmp1
                if (zl(jcnsl-1+(ino-1)*ncmp1+icmp1)) then
                    icmp = tmp_nucm1(icmp1)
                    iec = (icmp-1)/30+1
                    reste = icmp-30*(iec-1)
                    code = lshift(1, reste)
                    idg2 = jprn2-1+((2+nec)*(ino-1))+2+iec
                    zi(idg2) = ior(zi(idg2), code)
                    zi(jprn2-1+((2+nec)*(ino-1))+2) = zi(jprn2-1+((2+nec)*(ino-1))+2)+1
                end if
            end do
        end do
!
        ico = 0
        do ino = 1, nbno
            zi(jprn2-1+((2+nec)*(ino-1))+1) = ico+1
            ico = ico+zi(jprn2-1+((2+nec)*(ino-1))+2)
        end do
        call jelibe(nume_equa//'.PRNO')
!
!       2.4 CREATION DE .DEEQ : POUR DES RAISONS DE PERFORMANCES, IL VAUT MIEUX LE FAIRE PLUS TARD.
    end if
!
!   Get number of equations
    call jeexin(nume_equa//'.NEQU', iexi)
    if (iexi .eq. 0) then
        call jelira(nume_equa//'.NUEQ', 'LONUTI', nb_equa)
    else
        call jeveuo(nume_equa//'.NEQU', 'L', vi=v_nequ)
        nb_equa = v_nequ(1)
    end if
!
! - Create node field
!
!   CREATION DES SDS CHAM_NOS SIMPLE OU SIMULTANES
    if (ldist .and. (nb .ge. 2)) then
        call vtcreb(cno, base, tsca, &
                    meshz=ma, nume_equaz=nume_equa, idx_gdz=gd, nb_equa_inz=nb_equa, &
                    nbz=nb, vchamz=vcham)
    else
        call vtcreb(cno, base, tsca, &
                    meshz=ma, nume_equaz=nume_equa, idx_gdz=gd, nb_equa_inz=nb_equa)
    end if
!
!
!   5-BIS ON CREE SI NECESSAIRE LE .DEEQ DU NUME_EQUA
!   -------------------------------------------------
    if (l_crea_nume_equa) then
        !
        ! Create object local components (field) => global components (catalog)
        call cmpcha(cno, cmp_name, cata_to_field, field_to_cata, nb_cmpz=ncmp)
        ! Pour économiser la mémoire (pendant pteequ) on libère temporairement .cnsv et .cnsl
        call jelibe(cns//'.CNSV')
        call jelibe(cns//'.CNSL')
        call pteequ(nume_equa, base, nb_equa, gd, ncmp, field_to_cata)
        AS_DEALLOCATE(vi=cata_to_field)
        AS_DEALLOCATE(vi=field_to_cata)
        AS_DEALLOCATE(vk8=cmp_name)
        call jeveuo(cns//'.CNSV', 'L', jcnsv)
        call jeveuo(cns//'.CNSL', 'L', jcnsl)
    end if
!
!   6- ON REMPLIT LE .VALE
!   ----------------------
    call jeveuo(nume_equa//'.DEEQ', 'L', vi=deeq)
    call jeveuo(cno//'.VALE', 'E', jvale)
!
    do ieq2 = 1, nb_equa
        ino = deeq(2*(ieq2-1)+1)
        icmp = deeq(2*(ieq2-1)+2)
        if (ino*icmp .gt. 0) then
            nomcmp = zk8(jcmpgd-1+icmp)
            icmp1 = tmp_nucmp(icmp)
            !
            if (icmp1 .eq. 0) then
                if ((prol0 .eq. 'NON') .or. (prolv .eq. 'NON')) then
                    nomno = int_to_char8(ino)
                    valk(1) = nomcmp
                    valk(2) = nomno
                    valk(3) = cno
                    messag = 'CALCULEL2_13'
                    goto 70
                else
                    if (prol0 .eq. 'OUI') then
                        if (tsca .eq. 'R') then
                            zr(jvale-1+ieq2) = 0.d0
                        else if (tsca .eq. 'C') then
                            zc(jvale-1+ieq2) = (0.d0, 0.d0)
                        else if (tsca .eq. 'I') then
                            zi(jvale-1+ieq2) = 0
                        else if (tsca .eq. 'L') then
                            zl(jvale-1+ieq2) = .false.
                        else if (tsca .eq. 'K8') then
                            zk8(jvale-1+ieq2) = ' '
                        else
                            ASSERT(.false.)
                        end if
                    else if (prolv .eq. 'OUI') then
                        if (tsca .eq. 'R') then
                            zr(jvale-1+ieq2) = prolong%vale_r
                        else if (tsca .eq. 'C') then
                            zc(jvale-1+ieq2) = (0.d0, 0.d0)
                        end if
                    end if
                    goto 60
                end if
            end if
!
            if (zl(jcnsl-1+(ino-1)*ncmp1+icmp1)) then
!
                if (tsca .eq. 'R') then
! ----------------- Test for protect when nb_equa.ne.nb_dof
                    if (zr(jvale-1+ieq2) .ne. 0.d0) then
                        if (zr(jcnsv-1+(ino-1)*ncmp1+icmp1) .ne. zr(jvale-1+ieq2)) then
                            ASSERT(.false.)
                        end if
                    end if
                    zr(jvale-1+ieq2) = zr(jcnsv-1+(ino-1)*ncmp1+icmp1)
!
                else if (tsca .eq. 'C') then
! ----------------- Test for protect when nb_equa.ne.nb_dof
                    if (zc(jvale-1+ieq2) .ne. (0.d0, 0.d0)) then
                        if (zc(jvale-1+ieq2) .ne. zc(jcnsv-1+(ino-1)*ncmp1+icmp1)) then
                            ASSERT(.false.)
                        end if
                    end if
                    zc(jvale-1+ieq2) = zc(jcnsv-1+(ino-1)*ncmp1+icmp1)
!
                else if (tsca .eq. 'I') then
! ----------------- Test for protect when nb_equa.ne.nb_dof
                    if (zi(jvale-1+ieq2) .ne. 0) then
                        if (zi(jvale-1+ieq2) .ne. zi(jcnsv-1+(ino-1)*ncmp1+icmp1)) then
                            ASSERT(.false.)
                        end if
                    end if
                    zi(jvale-1+ieq2) = zi(jcnsv-1+(ino-1)*ncmp1+icmp1)
!
                else if (tsca .eq. 'L') then
! ----------------- Test for protect when nb_equa.ne.nb_dof
                    if (.not. zl(jvale-1+ieq2)) then
                        if (zl(jvale-1+ieq2) .neqv. zl(jcnsv-1+(ino-1)*ncmp1+icmp1)) then
                            ASSERT(.false.)
                        end if
                    end if
                    zl(jvale-1+ieq2) = zl(jcnsv-1+(ino-1)*ncmp1+icmp1)
!
                else if (tsca .eq. 'K8') then
! ----------------- Test for protect when nb_equa.ne.nb_dof
                    if (zk8(jvale-1+ieq2) .ne. ' ') then
                        if (zk8(jvale-1+ieq2) .ne. zk8(jcnsv-1+(ino-1)*ncmp1+icmp1)) then
                            ASSERT(.false.)
                        end if
                    end if
                    zk8(jvale-1+ieq2) = zk8(jcnsv-1+(ino-1)*ncmp1+icmp1)
!
                else
                    ASSERT(.false.)
                end if
!
            else
                if ((prol0 .eq. 'NON') .or. (prolv .eq. 'NON')) then
                    nomno = int_to_char8(ino)
                    valk(1) = nomcmp
                    valk(2) = nomno
                    valk(3) = cno
                    messag = 'CALCULEL2_13'
                    goto 70
                else
                    if (prol0 .eq. 'OUI') then
                        if (tsca .eq. 'R') then
                            zr(jvale-1+ieq2) = 0.d0
                        else if (tsca .eq. 'C') then
                            zc(jvale-1+ieq2) = (0.d0, 0.d0)
                        else if (tsca .eq. 'I') then
                            zi(jvale-1+ieq2) = 0
                        else if (tsca .eq. 'L') then
                            zl(jvale-1+ieq2) = .false.
                        else if (tsca .eq. 'K8') then
                            zk8(jvale-1+ieq2) = ' '
                        else
                            ASSERT(.false.)
                        end if
                    else if (prolv .eq. 'OUI') then
                        if (tsca .eq. 'R') then
                            zr(jvale-1+ieq2) = prolong%vale_r
                        else if (tsca .eq. 'C') then
                            zc(jvale-1+ieq2) = (0.d0, 0.d0)
                        end if
                    end if
                    goto 60
                end if
            end if
        end if
60      continue
    end do
!
!   7 - POUR ECONOMISER LES NUME_EQUA, ON REGARDE SI LE PRECEDENT NE CONVIENDRAIT PAS
!   ---------------------------------------------------------------------------------
    if (nume_equaz .eq. ' ' .and. base .eq. 'G') then
        read (nume_equa(15:19), '(I5)') nuprf
        if (nuprf .gt. 0) then
            prnoav = nume_equa
            call codent(nuprf-1, 'D0', prnoav(15:19))
            if (idensd('NUME_EQUA', nume_equa, prnoav)) then
                call detrsd('NUME_EQUA', nume_equa)
                call jeveuo(cno//'.REFE', 'E', jrefe)
                zk24(jrefe-1+2) = prnoav
                ! SI PARALLELISME EN TEMPS: ON PREND LA MEME DECISION POUR TOUS LES CHAM_NOS
                ! CALCULES AU MEME PAS PARALLELE (COHERENCE AVEC LE VTCREB PLUS HAUT)
                if (ldist) then
                    call jeveuo(vcham, 'L', jvcham)
                    do k = 1, nb
                        call jeveuo(zk24(jvcham+k-1) (1:19)//'.REFE', 'E', jrefk)
                        zk24(jrefk-1+2) = prnoav
                    end do
                end if
            end if
        end if
    end if
!
!
    iret = 0
    goto 80
!
!   MESSAGES D'ERREUR
70  continue
    ASSERT(kstop .eq. 'F' .or. kstop .eq. 'A' .or. kstop .eq. ' ')
    iret = 1
    call detrsd('CHAMP', cno)
    if (kstop .eq. ' ') goto 80
!
    if (messag .eq. 'CALCULEL2_12') then
        call utmess(kstop, 'CALCULEL2_12', nk=2, valk=valk)
    else if (messag .eq. 'CALCULEL2_13') then
        call utmess(kstop, 'CALCULEL2_13', nk=3, valk=valk)
    else
        ASSERT(.false.)
    end if
!
80  continue
    AS_DEALLOCATE(vi=tmp_nucmp)
    AS_DEALLOCATE(vi=tmp_nucm1)
    call jedema()
!     CALL UTIMSD(6,2,.TRUE.,.TRUE.,CNO,1,' ')
end subroutine
