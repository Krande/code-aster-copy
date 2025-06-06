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

subroutine amumpu(option, type, kxmps, usersm, nprec, &
                  lresol, nbfact)
!
!
    implicit none
!--------------------------------------------------------------
! BUT : UTILITAIRE POUR LES TRAITEMENTS CONNEXES LORS DU LANCEMENT
!       DES DIFFERENTES ETAPES DE MUMPS.
!
! OPTION=1 GESTION DE LA STRATEGIE MEMOIRE MUMPS (APRES ANALYSE)
!       CETTE ROUTINE DOIT ETRE APPELLEE ENTRE LA PHASE D'ANALYSE ET
!       CELLE DE FACTORISATION NUMERIQUE
! OPTION=2 DETECTION DES SINGULARITES (APRES FACTO) ET STOCKAGE DE CES
!          INFOS DS L'OBJET JEVEUX '&&AMUMP.PIVNUL' (V V I DIM=N+2)
!
! OPTION=4 RECUPERE LE DETERMINANT ET ON LE STOCKE DS L'OBJET JEVEUX
!          '&&AMUMP.DETERMINANT' (V V R DIM=3)
!
! IN  KXMPS  :   IN   : INDICE DE L'INSTANCE MUMPS DANS XMPS
!                       (INUTILE POUR OPTION=31)
! IN  TYPE   :   K1   : TYPE DU POINTEUR R OU C
!
! SI OPTION=1
! IN  USERSM :   K12  : STRATEGIE MEMOIRE DE L'UTILISATEUR
! IN  NBFACT :   IN   : NBRE DE FACTORISEES EN SIMULTANE (SI GESTION_MEMOIRE='AUTO')
!                 (INFORMATION SOUVENT ISSUE DE SD_SOLVEUR.SLVK(8))
! SI OPTION=2
! IN  NPREC  :   IN   : NBRE DE DIGITS POUR DETECTION DE SINGULARITE
! IN LRESOL  :  LOG   : .TRUE. SI ON FAIT LE SOLVE, .FALSE. SINON
!
! SI OPTION=4
! RAS
!---------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
!
#include "asterf_types.h"
#include "asterf.h"
#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_comm_jev.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jjldyn.h"
#include "asterfort/utgtme.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "mumps/cmumps.h"
#include "mumps/dmumps.h"
#include "mumps/smumps.h"
#include "mumps/zmumps.h"
    integer(kind=4) :: option
    integer(kind=8) :: kxmps, nprec, nbfact
    character(len=1) :: type
    character(len=12) :: usersm
    aster_logical :: lresol
!
#ifdef ASTER_HAVE_MUMPS
#include "asterf_mumps.h"
    type(smumps_struc), pointer :: smpsk => null()
    type(cmumps_struc), pointer :: cmpsk => null()
    type(dmumps_struc), pointer :: dmpsk => null()
    type(zmumps_struc), pointer :: zmpsk => null()
    real(kind=8) :: rval(3), rval1, rval2, rval3, rval2b, rval3b, rinf12
    real(kind=8) :: rinf13
    integer(kind=8) :: infog16, infog26, infog36, infog38, icntl35
    integer(kind=8) :: maxmem_ic, maxmem_ooc, vali(10), icoefm, icn22, icn23, rang, n, iaux1
    integer(kind=8) :: info3, nbproc, ifm, niv, ibid, ipiv, info28, info12, i
    integer(kind=8) :: tmax, tmaxb, ltot, iret, isizemu, nsizemu, nsizema, execmu
    integer(kind=8) :: info34, icnt33
    integer(kind=8) :: pid
!    integer ::  NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
    mpi_int :: mpicou
    aster_logical :: lpara, lpbmem, lpb1
    character(len=2) :: fstring
    character(len=8) :: k8tab(3)
    character(len=10) :: strpid
    character(len=24) :: kpiv, ksizemu
    character(len=80) :: nvers
!
    call jemarq()
! --- COMMUNICATEUR MPI DE TRAVAIL
    call asmpi_comm('GET', mpicou)
    call infdbg('SOLVEUR', ifm, niv)
!
!       ------------------------------------------------
! ---   INITS
!       ------------------------------------------------
    nvers(1:80) = ''
! --- OCCURENCE DE MUMPS EXISTE DEJA DS UN VECTEUR XMPS
    select case (type)
    case ('S')
        smpsk => smps(kxmps)
        lpara = (smpsk%nprocs .gt. 1)
        nbproc = smpsk%nprocs
        rang = smpsk%myid
        nvers = smpsk%version_number
        n = smpsk%n
    case ('C')
        cmpsk => cmps(kxmps)
        lpara = (cmpsk%nprocs .gt. 1)
        nbproc = cmpsk%nprocs
        rang = cmpsk%myid
        nvers = cmpsk%version_number
        n = cmpsk%n
    case ('D')
        dmpsk => dmps(kxmps)
        lpara = (dmpsk%nprocs .gt. 1)
        nbproc = dmpsk%nprocs
        rang = dmpsk%myid
        nvers = dmpsk%version_number
        n = dmpsk%n
    case ('Z')
        zmpsk => zmps(kxmps)
        lpara = (zmpsk%nprocs .gt. 1)
        nbproc = zmpsk%nprocs
        rang = zmpsk%myid
        nvers = zmpsk%version_number
        n = zmpsk%n
    case default
        ASSERT(.false.)
    end select
!
!       ------------------------------------------------
! ---   GESTION STRATEGIE MEMOIRE MUMPS (APRES ANALYSE)
!       ------------------------------------------------
    if (option .eq. 1) then
!
! ---   INITS GENERALES
        lpb1 = .false.
        lpbmem = .false.
        maxmem_ic = -9999
        maxmem_ooc = -9999
        info3 = -9999
        tmax = -9999
        tmaxb = -9999
        icn22 = -9999
        icn23 = -9999
        nsizema = -9999
        rval1 = -9999.0d0
        rval2 = -9999.0d0
        rval3 = -9999.0d0
        rval2b = -9999.0d0
        rval3b = -9999.0d0
!
! ---   INITS. PROPRE A L'OPTION
!       On récupère dans la structure de données MUMPS
!       quelques paramètres :
!       - icntl(35) => low rank/full rank
!       - estimations MUMPS de la mémoire nécessaire,
!         qui sont stockées dans des cases différentes
!         du vecteur infog selon que l'on est en full rank
!         ou bien en low rank
        select case (type)
        case ('S')
            infog16 = smpsk%infog(16)
            infog26 = smpsk%infog(26)
            infog36 = smpsk%infog(36)
            infog38 = smpsk%infog(38)
            icntl35 = smpsk%icntl(35)
            info3 = smpsk%infog(3)*4
        case ('C')
            infog16 = cmpsk%infog(16)
            infog26 = cmpsk%infog(26)
            infog36 = cmpsk%infog(36)
            infog38 = cmpsk%infog(38)
            icntl35 = cmpsk%icntl(35)
            info3 = cmpsk%infog(3)*8
        case ('D')
            infog16 = dmpsk%infog(16)
            infog26 = dmpsk%infog(26)
            infog36 = dmpsk%infog(36)
            infog38 = dmpsk%infog(38)
            icntl35 = dmpsk%icntl(35)
            info3 = dmpsk%infog(3)*8
        case ('Z')
            infog16 = zmpsk%infog(16)
            infog26 = zmpsk%infog(26)
            infog36 = zmpsk%infog(36)
            infog38 = zmpsk%infog(38)
            icntl35 = zmpsk%icntl(35)
            info3 = zmpsk%infog(3)*16
        end select
!
        select case (icntl35)
!       Full rank
        case (0, 3)
            maxmem_ic = infog16
            maxmem_ooc = infog26
!       Low rank
        case (1, 2)
            maxmem_ic = infog36
            maxmem_ooc = infog38
        case default
            ASSERT(.false.)
        end select
!
        ASSERT(nbproc > 0)
        if (info3 .lt. 0) then
            info3 = -info3/nbproc
        else
            info3 = info3/(1024*1024*nbproc)
        end if
!
! ---   nsizema: TAILLE CUMULEE EN MO DES OBJETS MUMPS A,IRN,RHS..
! ---   EXECMU:  TAILLE EN MO DE L'EXECUTABLE MUMPS
        execmu = 30
        ksizemu = '&&TAILLE_OBJ_MUMPS'
        call jeveuo(ksizemu, 'L', isizemu)
        do i = 1, nbproc
            nsizemu = zi(isizemu+i-1)
            if (nsizemu .gt. nsizema) nsizema = nsizemu
        end do
!
! ---   MARGES POUR LES ESTIMATIONS (EN %) DE MUMPS IC ET OOC
! ---     MARGE DU AU PARALLELISME
        if (lpara) then
            icoefm = 30
        else
! ---     MARGE EN SEQUENTIEL
            icoefm = 10
        end if
! ---   MARGE POUR LES TRES PETITS CAS
        if (n .lt. 100) icoefm = 50

#ifdef ASTER_PLATFORM_MINGW
        if (type .eq. 'Z' .and. (usersm(1:4) .eq. 'AUTO' .or. usersm(1:11) .eq. 'OUT_OF_CORE')) then
            usersm = 'IN_CORE'
            call utmess('I', 'FACTOR_86')
        end if
#endif
!
! ---   CONSO MUMPS  + VERIFICATION DE SA VALIDITE
        maxmem_ic = int(maxmem_ic*((icoefm+100)*1.d0/100.d0))
        maxmem_ooc = int(maxmem_ooc*((icoefm+100)*1.d0/100.d0))
        if ((maxmem_ic .lt. 0) .or. (maxmem_ooc .lt. 0) .or. (nsizema .lt. 0)) then
            lpbmem = .true.
            call utmess('A', 'FACTOR_83')
            if (usersm(1:4) .eq. 'AUTO') then
                lpb1 = .true.
                usersm = 'OUT_OF_CORE'
            end if
        end if
!
! ---   MEM ASTER DISPONIBLE + VERIFICATION DE SA VALIDITE
        k8tab(1) = 'VMPEAK  '
        k8tab(2) = 'MEM_TOTA'
        k8tab(3) = 'VMSIZE'
        call utgtme(3_8, k8tab, rval, iret)
        rval1 = rval(1)
        rval2 = rval(2)
        rval3 = rval(3)
        if (rval1 .le. 0.d0 .or. rval2 .le. 0.d0 .or. rval3 .le. 0.d0 &
            .or. rval2 .le. rval3 .or. iret .ne. 0) then
            lpbmem = .true.
            call utmess('I', 'FACTOR_82')
            if (usersm(1:4) .eq. 'AUTO') then
                lpb1 = .true.
                usersm = 'OUT_OF_CORE'
            end if
        else
            ASSERT((nbfact .ge. 1) .and. (nbfact .le. nmxins))
            tmax = max(int(0.95d0*(rval2-rval3)/nbfact), 1)
        end if
!
        if (niv .ge. 2) write (ifm, *) '<AMUMPU> RVAL1/2/3, maxmem_ic/26, NSIZEMA, TMAX ', rval1, &
            rval2, rval3, maxmem_ic, maxmem_ooc, nsizema, tmax
!

        select case (usersm)
        case ('IN_CORE')
! --------------
! ---   IN-CORE
! --------------
            icn22 = 0
            icn23 = 0
            if ((tmax .lt. maxmem_ic) .and. (.not. lpbmem)) then
                vali(1) = maxmem_ic
                vali(2) = tmax
                vali(3) = nsizema+execmu
                vali(4) = nbfact
                call utmess('A', 'FACTOR_74', ni=4_8, vali=vali)
            end if
        case ('OUT_OF_CORE')
! ------------------
! ---   OUT-OF-CORE
!-------------------
            icn22 = 1
            icn23 = 0
            if ((tmax .lt. maxmem_ooc) .and. (.not. lpbmem)) then
                vali(1) = maxmem_ooc
                vali(2) = tmax
                vali(3) = nsizema+execmu
                vali(4) = nbfact
                call utmess('A', 'FACTOR_75', ni=4_8, vali=vali)
            end if
        case ('AUTO')
! -----------------------------------------------------------------
! ----- STRATEGIE DECIDEE EN FONCTION DES CAPACITES MACHINES ET DES
! ----- CONSOMMATIONS REQUISES PAR MUMPS
! -----------------------------------------------------------------
            ASSERT((tmax .gt. 0) .and. (.not. lpbmem) .and. (maxmem_ic .ge. 0))
            ASSERT((maxmem_ooc .ge. 0) .and. (nsizema .gt. 0))
            if (tmax .ge. maxmem_ic) then
                icn22 = 0
                icn23 = max(min(3*maxmem_ic, tmax), 1_8)
            else
                call jjldyn(0_8, -1_8, ltot)
                k8tab(1) = 'MEM_TOTA'
                k8tab(2) = 'VMSIZE'
                call utgtme(2_8, k8tab, rval, iret)
                rval2b = rval(1)
                rval3b = rval(2)
                if ((rval2b .le. 0) .or. (rval3b .le. 0) .or. (rval2b .le. rval3b) .or. &
                    (iret .ne. 0)) then
                    lpbmem = .true.
                    call utmess('A', 'FACTOR_82')
                else
                    tmaxb = max(int(0.95d0*(rval2b-rval3b)/nbfact, 8), 1_8)
                end if
                if (niv .ge. 2) then
                    write (ifm, *) '<AMUMPU> RVALB2/3, TMAXB ', rval2b, rval3b, tmaxb
                    if (.not. lpbmem) then
                        vali(1) = int(rval3-rval3b)
                        call utmess('I', 'FACTOR_51', si=vali(1))
                    end if
                end if
                if ((tmaxb .gt. maxmem_ic) .and. (.not. lpbmem)) then
                    icn22 = 0
                    icn23 = max(min(3*maxmem_ic, tmaxb), 1_8)
                else if ((tmaxb .gt. maxmem_ooc) .and. (tmaxb .lt. maxmem_ic) .and. &
                         (.not. lpbmem)) then
                    icn22 = 1
                    icn23 = max(min(3*maxmem_ooc, tmaxb), 1_8)
                else
                    icn22 = 1
                    icn23 = 0
                    vali(1) = tmax
                    vali(2) = tmaxb
                    vali(3) = maxmem_ic
                    vali(4) = maxmem_ooc
                    vali(5) = nsizema+execmu
                    vali(6) = nbfact
                    if (.not. lpbmem) then
                        call utmess('A', 'FACTOR_76', ni=6_8, vali=vali)
                    else
                        call utmess('A', 'FACTOR_69', ni=6_8, vali=vali)
                    end if
                end if
            end if
        case ('EVAL')
! --------------------------------------------------
! ---   OPTION DE PRE-EVALUATION DES CONSOS MEMOIRE
! --------------------------------------------------
            icn22 = -1
            icn23 = -1
            k8tab(1) = 'CUSE_JV'
            k8tab(2) = 'RLQ_MEM'
            call utgtme(2_8, k8tab, rval, iret)
            rval1 = rval(1)
            rval2 = rval(2)
            if ((rval1 .le. 0) .or. (rval2 .le. 0) .or. (iret .ne. 0)) then
                call utmess('A', 'FACTOR_82')
            end if
            iaux1 = int(nbfact*rval1+rval2)
            vali(1) = n
            vali(2) = max(iaux1, 1_8)
            vali(3) = max((maxmem_ic+nsizema)*nbfact+execmu, 1_8)
            vali(4) = max((maxmem_ooc+nsizema)*nbfact+execmu, 1_8)
            vali(5) = max(info3*nbfact, 1_8)
            vali(6) = vali(2)+vali(3)
            vali(7) = vali(2)+vali(4)
            vali(8) = nbfact
            call utmess('I', 'FACTOR_81', ni=8_8, vali=vali)
        case default
            ASSERT(.false.)
        end select
! --- CORRECTIF POUR BENEFICIER DES BOUCLES DE RATTRAPAGE SI PB DS L'EVALUATION MEMOIRE
! --- ET GESTION_MEMOIRE='AUTO'
        if (lpb1) then
            ASSERT(usersm(1:11) .eq. 'OUT_OF_CORE')
            usersm = 'AUTO'
        end if
!
! ---  MODIFICATION DU PARAMETRAGE MUMPS POUR LA SUITE DU PROCESSUS
! ---- (FACTORISATION NUMERIQUE + SOLVE)
        select case (type)
        case ('S')
            smpsk%icntl(22) = to_mumps_int(icn22)
            smpsk%icntl(23) = to_mumps_int(icn23)
            smpsk%ooc_tmpdir = '.'
        case ('C')
            cmpsk%icntl(22) = to_mumps_int(icn22)
            cmpsk%icntl(23) = to_mumps_int(icn23)
            cmpsk%ooc_tmpdir = '.'
        case ('D')
            dmpsk%icntl(22) = to_mumps_int(icn22)
            dmpsk%icntl(23) = to_mumps_int(icn23)
            dmpsk%ooc_tmpdir = '.'
        case ('Z')
            zmpsk%icntl(22) = to_mumps_int(icn22)
            zmpsk%icntl(23) = to_mumps_int(icn23)
            zmpsk%ooc_tmpdir = '.'
        end select
!
        if (niv .ge. 2) then
! ---  NIVEAU DEVELOPPEUR
! ---  AFFICHAGE DE CONTROLE POUR DIAGNOSTIC MEMOIRE FIN
! ---  RECUPERATION DE L'AFFICHAGE DES CONSOS SYSTEMES
! ---  (VMPEAK, VMSIZE, VMDATA) + FREE DS LE FICHIER FORT.11
! ---  SI ON DECOMMENTARISE LES LIGNES 'CALL SYSTEM()' + 'GETPID'
            pid = 0
!          PID=getpid()
            if (abs(pid) < 10) then
                fstring = 'I1'
            else if (pid < 100) then
                fstring = 'I2'
            else if (pid < 1000) then
                fstring = 'I3'
            else if (pid < 10000) then
                fstring = 'I4'
            else if (pid < 100000) then
                fstring = 'I5'
            else if (pid < 1000000) then
                fstring = 'I6'
            else
                write (6, *) 'READ_VMPEAK : PB FORMAT CHOICE !'
            end if
            write (strpid, '('//fstring//')') pid
!          str=""
!          str="/proc/"//trim(adjustl(strpid))//"/status"info
!          CALL SYSTEM("cat "//str//" > fort.11")
!          CALL SYSTEM('free -m >> fort.11')
            write (ifm, *)
            write (ifm, *) '*********************************************'
            write (ifm, *) '<AMUMPU> GESTION MEMOIRE USERSM/ICN22/ICN23/NBFACT: ',&
     &      usersm, icn22, icn23, nbfact
            write (ifm, *) '<AMUMPU> CONSO MUMPS EXEC/OBJET_AIRNJCN/IC/OOC ',&
     &                 execmu, nsizema, maxmem_ic, maxmem_ooc
            write (ifm, *) '<AMUMPU> 1ERE ESTIMATION VMSIZE/MEM_TOTA/TMAX: ',&
     &                 rval3, rval2, tmax
            write (ifm, *) '<AMUMPU> 2NDE ESTIMATION VMSIZE/MEM_TOTA/TMAX: ',&
     &                 rval3b, rval2b, tmaxb
            write (ifm, *) '*********************************************'
        end if
!
!       ------------------------------------------------
! ---   DETECTION DES SINGULARITES (APRES FACTO)
!       ------------------------------------------------
    else if (option .eq. 2) then
!
! ---   INITS. PROPRE A L'OPTION
        select case (type)
        case ('S')
            info28 = smpsk%infog(28)
            info12 = smpsk%infog(12)
        case ('C')
            info28 = cmpsk%infog(28)
            info12 = cmpsk%infog(12)
        case ('D')
            info28 = dmpsk%infog(28)
            info12 = dmpsk%infog(12)
        case ('Z')
            info28 = zmpsk%infog(28)
            info12 = zmpsk%infog(12)
        end select
!
        if (nprec .ge. 0) then
            kpiv = '&&AMUMP.PIVNUL'
            call jeexin(kpiv, ibid)
            if (ibid .ne. 0) then
                ASSERT(.false.)
            else
                call wkvect(kpiv, 'V V I', n+2, ipiv)
                if (lresol) then
! ---   KPIV(1)= NOMBRE DE PIVOTS QUASI NULS (TOUS LE PROCS)
                    if (info28 .gt. n) then
                        ASSERT(.false.)
                    else
                        zi(ipiv) = info28
                    end if
! ---   KPIV(2)= NOMBRE DE PIVOTS NEGATIFS (TOUS LE PROCS)
                    if (info12 .gt. n) then
                        ASSERT(.false.)
                    else
                        zi(ipiv+1) = info12
                    end if
                    if (rang .eq. 0) then
! ---   KPIV(3..) LES PIVOTS QUASI NULS (ONLY PROC 0)
                        select case (type)
                        case ('S')
                            do i = 1, info28
                                zi(ipiv+1+i) = smpsk%pivnul_list(i)
                            end do
                        case ('C')
                            do i = 1, info28
                                zi(ipiv+1+i) = cmpsk%pivnul_list(i)
                            end do
                        case ('D')
                            do i = 1, info28
                                zi(ipiv+1+i) = dmpsk%pivnul_list(i)
                            end do
                        case ('Z')
                            do i = 1, info28
                                zi(ipiv+1+i) = zmpsk%pivnul_list(i)
                            end do
                        end select
                    end if
! ---   BCAST POUR COMMUNIQUER L'INFO AUX AUTRES PROCS
                    call asmpi_comm_jev('BCAST', kpiv)
                end if
            end if
! ---  AFFICHAGE DE CONTROLE
            if (niv .ge. 2) then
                write (ifm, *)
                write (ifm, *) &
                    '*********************************************'
                write (ifm, *) '<AMUMPU> TEST KPIV', zi(ipiv), zi(ipiv+1), &
                    zi(ipiv+2)
                write (ifm, *) &
                    '*********************************************'
            end if
!
        end if
!
    else if (option .eq. 4) then
        !
! ---   INITS. PROPRE A L'OPTION
        select case (type)
        case ('S')
            rinf12 = smpsk%rinfog(12)
            rinf13 = smpsk%rinfog(13)
            info34 = smpsk%infog(34)
            icnt33 = smpsk%icntl(33)
        case ('C')
            rinf12 = cmpsk%rinfog(12)
            rinf13 = cmpsk%rinfog(13)
            info34 = cmpsk%infog(34)
            icnt33 = cmpsk%icntl(33)
        case ('D')
            rinf12 = dmpsk%rinfog(12)
            rinf13 = dmpsk%rinfog(13)
            info34 = dmpsk%infog(34)
            icnt33 = dmpsk%icntl(33)
        case ('Z')
            rinf12 = zmpsk%rinfog(12)
            rinf13 = zmpsk%rinfog(13)
            info34 = zmpsk%infog(34)
            icnt33 = zmpsk%icntl(33)
        end select
        if (icnt33 .eq. 1) then
            kpiv = '&&AMUMP.DETERMINANT'
            call jeexin(kpiv, ibid)
            if (ibid .ne. 0) then
                call jeveuo(kpiv, 'E', ipiv)
            else
                call wkvect(kpiv, 'V V R', 3_8, ipiv)
            end if
! --- ON STOCKE LE CALCUL DU DET: MANTISSE * (2**EXP)
! --- MANTISSE=DCMPLX(RINF12,RINF13)
! --- EXP     =INFO34
            zr(ipiv) = rinf12
            zr(ipiv+1) = rinf13
            zr(ipiv+2) = info34
        end if

!       ------------------------------------------------
! ---   TESTS DIVERS SUR OPENMP (CF. AMUMPD)
!       ------------------------------------------------
!    else if (option.eq.99) then
!
!        call system('echo $OMP_NUM_THREADS')
!        call system('echo $MKL_NUM_THREADS')
!        call system('env')
!!$OMP PARALLEL PRIVATE(NTHREADS, TID)
!        TID = OMP_GET_THREAD_NUM()
!        write(ifm,*) 'Hello World from thread = ', TID
!       if (TID .EQ. 0) then
!            NTHREADS = OMP_GET_NUM_THREADS()
!            write(ifm,*) 'Number of threads = ', NTHREADS
!        endif
!!$OMP END PARALLEL
!
! --- CASE SUR LA VARIABLE OPTION
!
!
    else
        ASSERT(.false.)
    end if
!
    call jedema()
#endif
end subroutine
