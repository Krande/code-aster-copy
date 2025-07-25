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

subroutine op0014()
    implicit none
!     OPERATEUR FACTORISER
!     BUT: - FACTORISE UNE MATRICE ASSEMBLEE EN 2 MATRICES TRIANGULAIRES
!            (SOLVEUR MUMPS,MULT_FRONT,LDLT), OU
!          - DETERMINE UNE MATRICE DE PRECONDITIONNEMENT POUR L'ALGO DU
!            GRADIENT CONJUGUE PRCONDITIONNE (SOLVEUR GCPC)
!     ------------------------------------------------------------------

#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8prem.h"
#include "asterfort/apetsc.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/crsolv.h"
#include "asterfort/gcncon.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/mtcopy.h"
#include "asterfort/mtdefs.h"
#include "asterfort/mtdsc2.h"
#include "asterfort/mtdscr.h"
#include "asterfort/mtexis.h"
#include "asterfort/pcldlt.h"
#include "asterfort/pcmump.h"
#include "asterfort/titre.h"
#include "asterfort/tldlg3.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpu.h"
#include "asterfort/vrrefe.h"

    character(len=3) :: kstop
    character(len=5) :: klag2
    character(len=24) :: valk(2)
    character(len=8) :: matass, matfac, type, ktypr, ktyps, precon, mixpre, kacmum
    character(len=24) :: usersm
    character(len=16) :: concep, nomcmd, metres, renum
    character(len=19) :: mass, mfac, solveu, solvbd
    integer(kind=8) :: lch, i, lslvo
    integer(kind=8) :: nprec, iatfac, ibdeb, ibfin, ibid, ier1, ifm, ildeb, ilfin
    integer(kind=8) :: iret, isingu, istop, jadia, pcpiv, niremp
    integer(kind=8) :: ldtblo, lfnblo, ndeci, neq, niv, npvneg
    integer(kind=8) :: jslvk, jslvr, jslvi, reacpr, jslvo
    real(kind=8) :: fillin, epsmat, eps, blreps
    character(len=24), pointer :: refa(:) => null()
    character(len=2500) :: myopt
    aster_logical :: lreuse
!   ------------------------------------------------------------------
    call jemarq()

    call infmaj()
    call infniv(ifm, niv)

    call getres(matfac, concep, nomcmd)
    mfac = matfac
    call getvid('  ', 'MATR_ASSE', scal=matass, nbret=ibid)
    mass = matass
    lreuse = mfac .eq. mass

!   -- recuperation de certains mots cles :
!   ----------------------------------------
    call getvtx(' ', 'METHODE', scal=metres)
    call getvtx(' ', 'RENUM', scal=renum)

    if (metres .eq. 'MUMPS') then
        call getvtx(' ', 'ACCELERATION', iocc=1, scal=kacmum)
        call getvr8(' ', 'LOW_RANK_SEUIL', iocc=1, scal=blreps)
    else
        kacmum = 'XXXX'
        blreps = 0.d0
    end if

    if (metres .eq. 'GCPC' .or. metres .eq. 'PETSC') then
        call getvtx(' ', 'PRE_COND', scal=precon)
    else
        precon = ' '
    end if

    if (metres .eq. 'LDLT' .or. metres .eq. 'MUMPS' .or. metres .eq. 'MULT_FONT') then
        call getvis('  ', 'NPREC', scal=nprec)
        call getvtx('  ', 'STOP_SINGULIER', scal=kstop)
        if (kstop .eq. 'OUI') then
            istop = 0
        else if (kstop .eq. 'NON') then
            istop = 1
        end if
    else
        nprec = 0
        kstop = ' '
        istop = 0
    end if

!   -- la commande peut-elle ou doit-elle etre reentrante ?
!      .not.reuse <=> PCPC + LDLT_SP
!   --------------------------------------------------------
    precon = ' '
    call getvtx(' ', 'PRE_COND', scal=precon, nbret=ibid)
    if (metres .eq. 'GCPC' .and. precon .eq. 'LDLT_INC') then
        if (lreuse) then
            call utmess('F', 'ALGELINE5_56')
        end if
    else
        if (.not. lreuse) then
            call utmess('F', 'ALGELINE5_56')
        end if
    end if

!   -- on cree un solveur minimal pour retenir les infos entre FACTORISER et RESOUDRE:
!   ----------------------------------------------------------------------------------
    solveu = mfac(1:8)//'.SOLVEUR'
    call crsolv(metres, renum, kacmum, blreps, solveu, 'G')
    call jeveuo(mass//'.REFA', 'E', vk24=refa)
    refa(7) = solveu
    call jeveuo(solveu//'.SLVK', 'E', jslvk)
    call jeveuo(solveu//'.SLVR', 'E', jslvr)
    call jeveuo(solveu//'.SLVI', 'E', jslvi)
    call jeveuo(solveu//'.SLVO', 'E', jslvo)

    call uttcpu('CPU.RESO.1', 'DEBUT', ' ')
    call uttcpu('CPU.RESO.4', 'DEBUT', ' ')

!   -- CAS DU SOLVEUR  GCPC :
!   -------------------------
    if (metres .eq. 'GCPC') then
        if (concep(16:16) .eq. 'C') call utmess('F', 'ALGELINE5_57')
        zk24(jslvk-1+2) = precon

        if (precon .eq. 'LDLT_INC') then
            ASSERT(.not. lreuse)
            call copisd('MATR_ASSE', 'G', mass, mfac)
!           -- on ecrit dans la sd solveur le type de preconditionneur

            call getvis(' ', 'NIVE_REMPLISSAGE', scal=niremp, nbret=iret)
            call pcldlt(mfac, mass, niremp, 'G')
            call jeveuo(mfac//'.REFA', 'E', vk24=refa)
            refa(7) = solveu

        else if ((precon .eq. 'LDLT_SP') .or. (precon .eq. 'LDLT_DP')) then
!           -- nom du solveur pour mumps simple precision
            call gcncon('.', solvbd)
            zk24(jslvk-1+3) = solvbd

            call getvis(' ', 'REAC_PRECOND', scal=reacpr, nbret=ibid)
            call getvis(' ', 'PCENT_PIVOT', scal=pcpiv, nbret=ibid)
            call getvtx(' ', 'GESTION_MEMOIRE', scal=usersm, nbret=ibid)
            call getvr8(' ', 'LOW_RANK_SEUIL', iocc=1, scal=blreps, nbret=ibid)
            if (abs(blreps) < r8prem()) then
                kacmum = 'AUTO'
            else
                kacmum = 'LR'
            end if
            zi(jslvi-1+1) = -9999
            zi(jslvi-1+5) = 0
            zi(jslvi-1+6) = reacpr
            zi(jslvi-1+7) = pcpiv
            zk24(jslvk-1+5) = kacmum
            zk24(jslvk-1+9) = usersm
            zr(jslvr-1+4) = blreps

!           -- appel a la construction du preconditionneur
            call pcmump(mass, solveu, iret)
            if (iret .ne. 0) then
                call utmess('F', 'ALGELINE5_76', sk=precon)
            end if
        end if
        goto 999
    end if

!   -- CAS DU SOLVEUR MUMPS :
!   --------------------------
    if (metres .eq. 'MUMPS') then
        call getvtx(' ', 'TYPE_RESOL', scal=ktypr)
        call getvtx(' ', 'PRETRAITEMENTS', scal=ktyps)
        call getvtx(' ', 'ELIM_LAGR', scal=klag2)
        ASSERT(klag2 .eq. 'NON' .or. klag2 .eq. 'LAGR2')
        call getvtx(' ', 'GESTION_MEMOIRE', scal=usersm)
        mixpre = 'NON'
        epsmat = -1.d0
        eps = -1.d0
        call getvis(' ', 'PCENT_PIVOT', scal=pcpiv, nbret=ibid)

        zi(jslvi-1+1) = nprec
        zi(jslvi-1+2) = pcpiv
        zi(jslvi-1+3) = istop
        zi(jslvi-1+6) = 1
        zi(jslvi-1+7) = -9999
        zk24(jslvk-1+2) = ktyps
        zk24(jslvk-1+3) = ktypr
        zk24(jslvk-1+6) = klag2
        zk24(jslvk-1+7) = mixpre
        zk24(jslvk-1+8) = 'NON'
        zk24(jslvk-1+9) = usersm
        zk24(jslvk-1+10) = 'XXXX'
        zk24(jslvk-1+11) = 'XXXX'
        zk24(jslvk-1+12) = 'XXXX'
        zr(jslvr-1+1) = epsmat
        zr(jslvr-1+2) = eps
        ildeb = 1
        ilfin = 0
    end if

!   -- CAS DU SOLVEUR PETSC :
!   --------------------------
    if (metres .eq. 'PETSC') then
!       -- avec PETSC, on est forcement en reuse :
        mfac = mass
        zk24(jslvk-1+2) = precon

        call getvtx(' ', 'OPTION_PETSC', scal=myopt, nbret=ibid)
        ASSERT(ibid .eq. 1)
        lch = lxlgut(myopt)
        ASSERT(lch .lt. 2500)
        lslvo = int(lch/80)+1
        do i = 1, lslvo
            zk80(jslvo-1+i) = myopt(80*(i-1):80*i)
        end do
        !

        if (precon .eq. 'LDLT_INC') then
            call getvis(' ', 'NIVE_REMPLISSAGE', scal=niremp, nbret=ibid)
            call getvr8(' ', 'REMPLISSAGE', scal=fillin, nbret=ibid)
            zr(jslvr-1+3) = fillin
            zi(jslvi-1+4) = niremp

        else if ((precon .eq. 'LDLT_SP') .or. (precon .eq. 'LDLT_DP')) then
!           -- nom du solveur pour mumps simple precision/low_rank
            call gcncon('.', solvbd)
            zk24(jslvk-1+3) = solvbd
            kacmum = 'XXXX'
            blreps = 0.d0
            call getvis(' ', 'REAC_PRECOND', scal=reacpr, nbret=ibid)
            ASSERT(ibid .eq. 1)
            call getvis(' ', 'PCENT_PIVOT', scal=pcpiv, nbret=ibid)
            ASSERT(ibid .eq. 1)
            call getvtx(' ', 'GESTION_MEMOIRE', scal=usersm)
            ASSERT(ibid .eq. 1)
            call getvr8(' ', 'LOW_RANK_SEUIL', iocc=1, scal=blreps, nbret=ibid)
            ASSERT(ibid .eq. 1)
            if (abs(blreps) < r8prem()) then
                kacmum = 'AUTO'
            else
                kacmum = 'LR'
            end if
            zi(jslvi-1+1) = -9999
            zi(jslvi-1+5) = 0
            zi(jslvi-1+6) = reacpr
            zi(jslvi-1+7) = pcpiv
            zi(jslvi-1+9) = lslvo
            zk24(jslvk-1+5) = kacmum
            zk24(jslvk-1+9) = usersm
            zr(jslvr-1+4) = blreps
        end if

        call apetsc('DETR_MAT', ' ', mfac, [0.d0], ' ', &
                    0, ibid, iret)
        call apetsc('PRERES', solveu, mfac, [0.d0], ' ', &
                    0, ibid, iret)
        ASSERT(iret .eq. 0)
        goto 999
    end if

!   -- CAS DES SOLVEURS LDLT/MULT_FRONT/MUMPS :
!   -------------------------------------------

!   -- RECUPERATION DES INDICES DE DEBUT ET FIN DE LA FACTORISATION -
    if (metres .ne. 'MUMPS') then
!       -- 1) AVEC DDL_XXX
        ildeb = 1
        ilfin = 0
        call getvis('  ', 'DDL_DEBUT', scal=ildeb, nbret=ibid)
        call getvis('  ', 'DDL_FIN', scal=ilfin, nbret=ibid)
!       -- 2) AVEC BLOC_XXX
        ibdeb = 1
        ibfin = 0
        call getvis('  ', 'BLOC_DEBUT', scal=ibdeb, nbret=ldtblo)
        call getvis('  ', 'BLOC_FIN', scal=ibfin, nbret=lfnblo)

!       -- EXISTENCE / COMPATIBILITE DES MATRICES ---
        call mtexis(mfac, iret)
        if (iret .ne. 0) then
            call vrrefe(mass, mfac, ier1)
            if (ier1 .ne. 0) then
                valk(1) = matass
                valk(2) = matfac
                call utmess('F', 'ALGELINE2_18', nk=2, valk=valk)
            else if (mfac .ne. mass) then
                if (ildeb .eq. 1 .and. ibdeb .eq. 1) then
                    call mtcopy(mass, mfac, iret)
                    ASSERT(iret .eq. 0)
                end if
            end if
        else
            type = ' '
            call mtdefs(mfac, mass, 'GLOBALE', type)
            call mtcopy(mass, mfac, iret)
            ASSERT(iret .eq. 0)
        end if
    end if

!   -- CHARGEMENT DES DESCRIPTEURS DE LA MATRICE A FACTORISER ---
    call mtdscr(mfac)
    call jeveuo(mfac(1:19)//'.&INT', 'E', iatfac)
    if (iatfac .eq. 0) then
        call utmess('F', 'ALGELINE2_19', sk=matfac)
    end if
    call mtdsc2(zk24(zi(iatfac+1)), 'SXDI', 'L', jadia)

!   -- NEQ : NOMBRE D'EQUATIONS (ORDRE DE LA MATRICE) ---
    neq = zi(iatfac+2)

    if (metres .ne. 'MUMPS') then

!       -- VERIFICATION DES ARGUMENTS RELATIF A LA PARTIE A FACTORISER
!       -- 1) AVEC DDL_XXX
        if (ilfin .lt. ildeb .or. ilfin .gt. neq) ilfin = neq

!       -- 2) AVEC BLOC_XXX
        if (ldtblo .ne. 0) then
            if (ibdeb .lt. 1) then
                call utmess('A', 'ALGELINE2_1')
                ibdeb = 1
            else if (ibdeb .gt. zi(iatfac+13)) then
                call utmess('F', 'ALGELINE2_20')
            end if
            ildeb = zi(jadia+ibdeb-2)+1
        end if
        if (lfnblo .ne. 0) then
            if (ibfin .lt. 1) then
                call utmess('F', 'ALGELINE2_21')
            else if (ibdeb .gt. zi(iatfac+13)) then
                call utmess('A', 'ALGELINE2_8')
                ibfin = zi(iatfac+13)
            end if
            ilfin = zi(jadia+ibfin-1)
        end if

!       -- IMPRESSION SUR LE FICHIER MESSAGE ----------------------------
        if (niv .eq. 2) then
            write (ifm, *) ' +++ EXECUTION DE "', nomcmd, '"'
            write (ifm, *) '       NOM DE LA MATRICE ASSEMBLEE  "', matass, '"'
            write (ifm, *) '       NOM DE LA MATRICE FACTORISEE "', matfac, '"'
            if (ildeb .eq. 1 .and. ilfin .eq. neq) then
                write (ifm, *) '     FACTORISATION COMPLETE DEMANDEE'
            else
                write (ifm, *) '     FACTORISATION PARTIELLE DE LA LIGNE',&
     &        ildeb, ' A LA LIGNE ', ilfin
            end if
            write (ifm, *) '     NOMBRE TOTAL D''EQUATIONS  ', neq
            write (ifm, *) '     NB. DE CHIFFRES SIGNIF. (NPREC) ', nprec
            write (ifm, *) ' +++ -------------------------------------------'
        end if
    end if

!   ------------------ FACTORISATION EFFECTIVE -------------------
    call tldlg3(metres, renum, istop, iatfac, ildeb, &
                ilfin, nprec, ndeci, isingu, npvneg, &
                iret, ' ')
!   --------------------------------------------------------------

999 continue

    call uttcpu('CPU.RESO.1', 'FIN', ' ')
    call uttcpu('CPU.RESO.4', 'FIN', ' ')

    call titre()

    call jedema()
end subroutine
