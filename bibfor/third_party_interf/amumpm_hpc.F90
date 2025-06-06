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

subroutine amumpm_hpc(kxmps, kmonit, impr, ifmump, &
                      klag2, type, epsmat, ktypr, &
                      lpreco, lbloc)
!
!
    implicit none
!--------------------------------------------------------------
! BUT : REMPLIR LES OBJETS F90 DE MUMPS REPRESENTANT LA MATRICE A PARTIR
!       DE CELLE DE CODE_ASTER.
!
! IN  KXMPS  :   IN   : INDICE DE L'INSTANCE MUMPS DANS DMPS
! IN  KMONIT :  K24   : VECTEUR DE NOMS DES OBJ JEVEUX
! IN  IMPR   :  K14   : FLAG POUR IMPRESSION MATRICE
! IN  IFMUMP :   IN   : UNITE LOGIQUE POUR IMPRESSION FICHIER
! IN  KLAG2  :   K5   : PARAMETRE DE SOLVEUR/ELIM_LAGR
! IN  TYPE   :   K1   : TYPE DU POINTEUR R OU C
! IN  EPSMAT :   R8   : SEUIL DE FILTRAGE DES TERMES DE LA MATRICE
! IN  KTYPR  :   K8   : TYPE DE RESOLUTION MUMPS (SYMDEF...)
! IN  LPRECO :  LOG   : MUMPS EST-IL UTILISE COMME PRECONDITIONNEUR ?
! IN  LBLOC  :  LOG   : LOGIQUE PRECISANT SI ON EFFECTUE L ANALYSE PAR BLOCS
!---------------------------------------------------------------
! aslint: disable=W1501
! person_in_charge: olivier.boiteau at edf.fr
!
#include "asterf.h"
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r4maem.h"
#include "asterc/r4miem.h"
#include "asterc/r8maem.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/asmpi_comm_jev.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibd.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/crnustd.h"
#include "asterfort/filtering.h"
    integer(kind=8) :: kxmps, ifmump
    aster_logical :: lpreco, lbloc
    real(kind=8) :: epsmat
    character(len=1) :: type
    character(len=5) :: klag2
    character(len=8) :: ktypr
    character(len=14) :: impr
    character(len=24) :: kmonit(12)
!
#ifdef ASTER_HAVE_MUMPS
#include "asterf_mumps.h"
!
!
    type(smumps_struc), pointer :: smpsk => null()
    type(cmumps_struc), pointer :: cmpsk => null()
    type(dmumps_struc), pointer :: dmpsk => null()
    type(zmumps_struc), pointer :: zmpsk => null()
    integer(kind=8) :: nsmdi, nsmhc, neql, n1, nzl, nvale, jvale
    integer(kind=8) :: nlong, jvale2, nzloc, iterm, ifm, niv, k
    integer(kind=8) :: sym, iret, jcoll, iligl, kzero, ltot
    integer(kind=8) :: iad, vali(2), nbproc, nfilt1, nfilt2
    integer(kind=8) :: nfilt3, isizemu, nsizemu, rang, esizemu
    integer(kind=8) :: nuno1, nuno2, procol, prolig, nucmp1, nucmp2
    integer(kind=8) ::  nzdeb, nzfin, numno1, numno2, iligg, jcolg, nloc
    mumps_int :: neqg, nz2, iligg4, jcolg4
    character(len=8) :: k8bid
    character(len=14) :: nonu
    character(len=16) :: k16bid, nomcmd
    character(len=19) :: nomat
    character(len=24) :: kfiltr, kpiv, kpiv2, ksizemu
    real(kind=8) :: raux, rfiltr, epsmac, rmax, rmin, rtest, raux2
    complex(kind=8) :: caux, caux2
    aster_logical :: lmnsy, ltypr, lnn, lfiltr, eli2lg, lsimpl, lcmde
    aster_logical :: lgive, lnn2, lspd
    aster_logical, parameter :: ldebug = ASTER_FALSE
    integer(kind=8), pointer :: smdi(:) => null()
    integer(kind=8), pointer :: nequ(:) => null()
    integer(kind=8), pointer :: nuls(:) => null()
    integer(kind=8), pointer :: deeg(:) => null()
    integer(kind=8), pointer :: owner(:) => null()
    integer(kind=8), pointer :: delg(:) => null()
    integer(kind=8), pointer :: nulg(:) => null()
    integer(kind=8), pointer :: pddl(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    integer(kind=4), pointer :: smhc(:) => null()
    integer(kind=4), pointer :: iok(:) => null()
    integer(kind=4), pointer :: iok2(:) => null()
    real(kind=8), pointer :: filter(:) => null()
!
!-----------------------------------------------------------------------
    call jemarq()
    call infdbg('SOLVEUR', ifm, niv)
!
!       ------------------------------------------------
!        INITS
!       ------------------------------------------------
!
    ASSERT(.not. lbloc)
    epsmac = r8prem()
    nomat = nomats(kxmps)
    nonu = nonus(kxmps)

!   Adresses needed to get the stiffness matrix wrt nodes and dof numbers (see below)
    if (ldebug) then
        print *, "DEBUG IN AMUMPM_HPC"
        call jeexin(nonu//'.NUME.NULS', iret)
        if (iret == 0) then
            call crnustd(nonu)
        end if
        call jeveuo(nonu//'.NUME.NULS', 'L', vi=nuls)
        call jeveuo(nonu//'.NUME.DEEG', 'L', vi=deeg)
    end if

! --- REMPLISSAGE DE DIFFERENTS OBJETS SUIVANT LE TYPE DU POINTEUR
! --- DE MUMPS: DMUMPS_STRUC OU ZMUMPS_STRUC
    if (type == 'S') then
        smpsk => smps(kxmps)
        ltypr = ASTER_TRUE
        sym = smpsk%sym
        rmax = r4maem()*real(0.5d0, 4)
        rmin = r4miem()*real(2.0d0, 4)
        nbproc = smpsk%nprocs
        rang = smpsk%myid
        esizemu = 4
    else if (type == 'C') then
        cmpsk => cmps(kxmps)
        ltypr = ASTER_FALSE
        sym = cmpsk%sym
        rmax = r4maem()*real(0.5d0, 4)
        rmin = r4miem()*real(2.0d0, 4)
        nbproc = cmpsk%nprocs
        rang = cmpsk%myid
        esizemu = 8
    else if (type == 'D') then
        dmpsk => dmps(kxmps)
        ltypr = ASTER_TRUE
        sym = dmpsk%sym
        rmax = r8maem()*0.5d0
        rmin = r8miem()*2.0d0
        nbproc = dmpsk%nprocs
        rang = dmpsk%myid
        esizemu = 8
    else if (type == 'Z') then
        zmpsk => zmps(kxmps)
        ltypr = ASTER_FALSE
        sym = zmpsk%sym
        rmax = r8maem()*0.5d0
        rmin = r8miem()*2.0d0
        nbproc = zmpsk%nprocs
        rang = zmpsk%myid
        esizemu = 16
    else
        ASSERT(ASTER_FALSE)
    end if
!
!
    if (ktypr(1:6) == 'SYMDEF') then
        if ((type == 'C') .or. (type == 'Z')) then
            call utmess('F', 'FACTOR_80')
        end if
        lspd = ASTER_TRUE
    else
        lspd = ASTER_FALSE
    end if
!
    lsimpl = (type == 'S') .or. (type == 'C')
!
!     ------------------------------------------
!     DETERMINATION DE LA SYMETRIE DE LA MATRICE
!     ------------------------------------------
    call jelira(nomat//'.VALM', 'NMAXOC', nvale)
! --- LMNSY EST INDEPENDANT DE XMPSK%SYM POUR POUVOIR TRAITER
! --- DES CAS SYMETRIQUES EN MUMPS NON SYMETRIQUE
    if (nvale == 1) then
        lmnsy = ASTER_FALSE
    else if (nvale == 2) then
        lmnsy = ASTER_TRUE
    else
        ASSERT(ASTER_FALSE)
    end if
!
!
!   ------------------------------------------------
!    lecture d'adresses et de parametres preliminaires
!   ------------------------------------------------
    call jeveuo(nonu//'.SMOS.SMDI', 'L', vi=smdi)
    call jelira(nonu//'.SMOS.SMDI', 'LONMAX', nsmdi)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', vi4=smhc)
    call jelira(nonu//'.SMOS.SMHC', 'LONMAX', nsmhc)
    call jelira(nonu//'.NUME.DELG', 'LONMAX', n1)
    call jeveuo(nonu//'.NUME.DELG', 'L', vi=delg)
    call jeveuo(nonu//'.NUME.NEQU', 'L', vi=nequ)
    call jeveuo(nonu//'.NUME.PDDL', 'L', vi=pddl)
    call jeveuo(nonu//'.NUME.DEEQ', 'L', vi=deeq)
    call jeveuo(nonu//'.NUME.NULG', 'L', vi=nulg)
    neqg = to_mumps_int(nequ(2))
    nloc = nequ(1)
    ASSERT(n1 == nsmdi)
! --- CALCUL DE N
    neql = nsmdi
!
! --- GESTION ELIM_LAGR='LAGR2'
! --- on a essaye une basule automatique ELIM_LAGR='LAGR2'/'NON'
! --- en fonction de la proportion de lagranges. en fait, 'LAGR2'
! --- offre la plupart du temps le meilleur compromis:
! ---     cpu x memoire x qualite --> on ne programme pas de
! --- bascule, on laisse 'LAGR2' par defaut (pour l'instant)
    select case (klag2(1:5))
    case ('LAGR2')
        eli2lg = ASTER_TRUE
    case ('OUI', 'NON', 'XXXX')
        eli2lg = ASTER_FALSE
    case default
        ASSERT(ASTER_FALSE)
    end select
!
! --- CALCUL DE NZ2
    nzl = smdi(neql)
    ASSERT(nzl .le. nsmhc)
!
    call jeveuo(jexnum(nomat//'.VALM', 1_8), 'L', jvale)
    call jelira(jexnum(nomat//'.VALM', 1_8), 'LONMAX', nlong)
    ASSERT(nlong == nzl)
    if (lmnsy) then
        call jeveuo(jexnum(nomat//'.VALM', 2_8), 'L', jvale2)
        call jelira(jexnum(nomat//'.VALM', 2_8), 'LONMAX', nlong)
        ASSERT(nlong == nzl)
    end if
!
! ---  DETERMINATION DES TERMES DE REFERENCE POUR LE FILTRAGE
! ---  NFILT1 : TERMES RIGOUREUSEMENT NULS
! ---  NFILT2 : DEPASSEMENT DE CAPACITE PAR VALEUR MAX
! ---  NFILT3 : DEPASSEMENT DE CAPACITE PAR VALEUR MIN
    lfiltr = ASTER_FALSE
    kfiltr = '&&AMUMPM.FILTRAGE'
    nfilt1 = 0
    nfilt2 = 0
    nfilt3 = 0
    if (epsmat > 0.d0) then
        lfiltr = ASTER_TRUE
        call jedetr(kfiltr)
        call wkvect(kfiltr, 'V V R', neql, vr=filter)
        if (ltypr) then
            do k = 1, neql
                filter(k) = epsmat*abs(zr(jvale-1+smdi(k)))
            end do
        else
            do k = 1, neql
                filter(k) = epsmat*abs(zc(jvale-1+smdi(k)))
            end do
        end if
! --- SEUILLAGE DES TERMES DE FILTRAGE POUR EVITER LES VALEURS ABBERANTE
        do k = 1, neql
            rtest = filter(k)
            if (rtest < rmin) filter(k) = 0.d0
            if (rtest > rmax) filter(k) = rmax
        end do
    else
        lfiltr = ASTER_FALSE
    end if
!
! --- VECTEURS AUXILIAIRES POUR FILTRER LES TERMES MATRICIELS
! --- KPIV: PROFIL TRIANGULAIRE INF
! --- KPIV2: IDEM SUP
    kpiv = '&&AMUMPM.TERMEOK'
    kpiv2 = '&&AMUMPM.TERMEOK2'
    call jedetr(kpiv)
    call jedetr(kpiv2)
    call wkvect(kpiv, 'V V S', nzl, vi4=iok)
    iok(1:nzl) = 0
    call wkvect(kpiv2, 'V V S', nzl, vi4=iok2)
    iok2(1:nzl) = 0
!
! --- On commence par compter le nombre de terme pour le filtrage
    nzloc = 0
    jcoll = 1
    rfiltr = -1.d0
    do jcoll = 1, nloc
        if (jcoll == 1) then
            nzdeb = 1
        else
            nzdeb = smdi(jcoll-1)+1
        end if
        nzfin = smdi(jcoll)
        procol = pddl(jcoll)
        jcolg = nulg(jcoll)
        nuno2 = 0
        if (deeq((jcoll-1)*2+1) > 0) then
            nuno2 = 1
        end if
        do k = nzdeb, nzfin
            iligl = smhc(k)
            prolig = pddl(iligl)
            iligg = nulg(iligl)
            nuno1 = 0
            if (deeq((iligl-1)*2+1) > 0) then
                nuno1 = 1
            end if
!
! --- Dans le cas symétrique, on ne garde que la triangulaire supérieur
!     de la matrice globale. Attention, ces termes ne viennent pas
!     forcément de la triangulaire supérieur
!     Attention MUMPS cumule les valeurs si les indices sont répétés en local et global
!     C'est un peu l'équivalent de ADD_VALUES de Petsc
!
! --- Filtrage
!
            if (lfiltr) then
                rfiltr = filter(iligl)+filter(jcoll)
            end if
!
            if (nuno1 /= 0 .and. nuno2 /= 0) then
                if (prolig == rang) then
                    if (sym == 0 .or. jcolg >= iligg) then
                        call filtering(ltypr, k, ASTER_FALSE, jvale, jvale2, rfiltr, rmin, rmax, &
                                       nfilt1, nfilt2, nfilt3, nzloc, iok)
                    end if
                    if (procol == rang .and. iligg /= jcolg) then
                        if (sym == 0 .or. jcolg < iligg) then
                            call filtering(ltypr, k, lmnsy, jvale, jvale2, rfiltr, rmin, rmax, &
                                           nfilt1, nfilt2, nfilt3, nzloc, iok2)
                        end if
                    end if
                else if (procol == rang) then
                    if (sym == 0 .or. jcolg < iligg) then
                        call filtering(ltypr, k, lmnsy, jvale, jvale2, rfiltr, rmin, rmax, &
                                       nfilt1, nfilt2, nfilt3, nzloc, iok2)
                    end if
                end if
            else
!               Si on est sur un ddl de Lagrange et qu'on possede le ddl d'en face
!               ou que les deux ddl sont de Lagrange, on doit donner le terme
                lgive = (nuno1 == 0 .and. procol == rang) .or. &
                        (nuno2 == 0 .and. prolig == rang) .or. &
                        (nuno1 == 0 .and. nuno2 == 0)
                if (lgive) then
                    if (sym == 0 .or. jcolg >= iligg) then
                        call filtering(ltypr, k, ASTER_FALSE, jvale, jvale2, rfiltr, rmin, rmax, &
                                       nfilt1, nfilt2, nfilt3, nzloc, iok)
                    end if
                    if (iligg /= jcolg) then
                        if (sym == 0 .or. jcolg < iligg) then
                            call filtering(ltypr, k, lmnsy, jvale, jvale2, rfiltr, rmin, rmax, &
                                           nfilt1, nfilt2, nfilt3, nzloc, iok2)
                        end if
                    end if
                end if
            end if
        end do
    end do
    nz2 = to_mumps_int(nzloc+1000)
    if (niv >= 2 .or. ldebug) then
        write (ifm, *) '<AMUMPM>     NZLOC: ', nzloc
        write (ifm, *) '       TERMES NULS: ', nfilt1
        write (ifm, *) '   UNDER/OVERFLOWS: ', nfilt3, '/', nfilt2
    end if
!
!       ------------------------------------------------
!       ALLOCATION DES OBJETS MUMPS F90
!       ------------------------------------------------
!
! ---   OBJET JEVEUX STOCKANT, PAR PROC, LA TAILLE DES OBJETS ALLOUES
! ---   POUR MUMPS. UTILE A AMUMPU. ON SUPPOSE UN SEUL RHS.
! ---   EN FIN DE ROUTINE, TOUS LES PROCS CONNAISSENT CE VECTEUR
! ----  VIA UN MPI_ALLREDUCE + SUM.
    ksizemu = '&&TAILLE_OBJ_MUMPS'
    call jeexin(ksizemu, iret)
    if (iret == 0) then
        call wkvect(ksizemu, 'V V I', nbproc, isizemu)
    else
        call jeveuo(ksizemu, 'E', isizemu)
    end if
    do k = 1, nbproc
        zi(isizemu+k-1) = 0
    end do
    nsizemu = 0
    if (type == 'S') then
        smpsk%n = neqg
        smpsk%nnz_loc = nz2
        allocate (smpsk%irn_loc(nz2))
        allocate (smpsk%jcn_loc(nz2))
        allocate (smpsk%a_loc(nz2))
    else if (type == 'C') then
        cmpsk%n = to_mumps_int(neqg)
        cmpsk%nnz_loc = nz2
        allocate (cmpsk%irn_loc(nz2))
        allocate (cmpsk%jcn_loc(nz2))
        allocate (cmpsk%a_loc(nz2))
    else if (type == 'D') then
        dmpsk%n = neqg
        dmpsk%nnz_loc = nz2
        allocate (dmpsk%irn_loc(nz2))
        allocate (dmpsk%jcn_loc(nz2))
        allocate (dmpsk%a_loc(nz2))
    else if (type == 'Z') then
        zmpsk%n = to_mumps_int(neqg)
        zmpsk%nnz_loc = nz2
        allocate (zmpsk%irn_loc(nz2))
        allocate (zmpsk%jcn_loc(nz2))
        allocate (zmpsk%a_loc(nz2))
    else
        ASSERT(ASTER_FALSE)
    end if
    nsizemu = nz2*(4+4+esizemu)+esizemu*neql
    zi(isizemu+rang) = 1+(nsizemu/(1024*1024))
!
    if (ldebug) allocate (owner(nz2))
!
!       ------------------------------------------------
!       INTERPRETATION DES PBS RENCONTRES LORS DU FILTRAGE
!       REMPLISSAGE DE LA MATRICE
!       EN CAS D'OVERFLOW:
!         SI SOLVEUR DIRECT: UTMESS_F
!         SI SIMPLE PRECISION : ON SEUILLE LES TERMES IMPACTES
!                               LORS DU REMPLISSAGE EFFECTIF.
!       EN CAS D'UNDERFLOW: ON FIXE A ZERO.
!       ------------------------------------------------
!
    if (.not. lsimpl) then
        vali(1) = nfilt3
        vali(2) = nfilt2
        do k = 1, nzl
            if (iok(k) == -1) then
                call utmess('F', 'FACTOR_78', ni=2_8, vali=vali)
            end if
        end do
        do k = 1, nzl
            if (iok2(k) == -1) then
                call utmess('F', 'FACTOR_78', ni=2_8, vali=vali)
            end if
        end do
    end if
!
! --- On remplit la matrice comme avec Petsc (voir apmamh.F90)
    raux = 0.d0
    raux2 = 0.d0
    caux = dcmplx(0.d0, 0.d0)
    caux2 = dcmplx(0.d0, 0.d0)
    iterm = 0
    do jcoll = 1, nloc
        if (jcoll == 1) then
            nzdeb = 1
        else
            nzdeb = smdi(jcoll-1)+1
        end if
        nzfin = smdi(jcoll)
        procol = pddl(jcoll)
        jcolg = nulg(jcoll)
        jcolg4 = to_mumps_int(jcolg+1)
        nuno2 = 0
        if (deeq((jcoll-1)*2+1) > 0) then
            nuno2 = 1
        else
            if (deeq((jcoll-1)*2+1) < 0 .and. deeq((jcoll-1)*2+2) > 0) then
                nuno2 = 1
            end if
        end if
        do k = nzdeb, nzfin
            iligl = smhc(k)
            prolig = pddl(iligl)
            iligg = nulg(iligl)
            iligg4 = to_mumps_int(iligg+1)
            nuno1 = 0
            if (deeq((iligl-1)*2+1) > 0) then
                nuno1 = 1
            else
                if (deeq((iligl-1)*2+1) < 0 .and. deeq((iligl-1)*2+2) > 0) then
                    nuno2 = 1
                end if
            end if
!
! --- Filtrage
!
            lnn = ASTER_FALSE
            lnn2 = ASTER_FALSE
            if (ltypr) then
                raux = zr(jvale-1+k)
                if (iok(k) == 1) then
                    lnn = ASTER_TRUE
                else if (iok(k) == -1) then
                    lnn = ASTER_TRUE
                    raux = rmax*sign(1.d0, raux)
                end if
                iok(k) = 0
!
                raux2 = zr(jvale-1+k)
                if (lmnsy) raux2 = zr(jvale2-1+k)
                if (iok2(k) == 1) then
                    lnn2 = ASTER_TRUE
                else if (iok2(k) == -1) then
                    lnn2 = ASTER_TRUE
                    raux2 = rmax*sign(1.d0, raux2)
                end if
                iok2(k) = 0
            else
                caux = zc(jvale-1+k)
                if (iok(k) == 1) then
                    lnn = ASTER_TRUE
                else if (iok(k) == -1) then
                    lnn = ASTER_TRUE
                    caux = rmax*dcmplx(1.d0*sign(1.d0, dble(caux)), &
                                       1.d0*sign(1.d0, imag(caux)))
                end if
                iok(k) = 0
!
                caux2 = zc(jvale-1+k)
                if (lmnsy) caux2 = zc(jvale2-1+k)
                if (iok2(k) == 1) then
                    lnn2 = ASTER_TRUE
                else if (iok2(k) == -1) then
                    lnn2 = ASTER_TRUE
                    caux2 = rmax*dcmplx(1.d0*sign(1.d0, dble(caux2)), &
                                        1.d0*sign(1.d0, imag(caux2)))
                end if
                iok2(k) = 0
            end if
!
! --- Pour ELIM_LAGR
            kzero = 0
            if (eli2lg) then
! ------      ON ELIMINE LE DERNIER TERME DE A/A_loc SI LAG2
                if (delg(iligl) == -1) then
                    if (jcoll == iligl) kzero = 1
                    if (delg(jcoll) == -2) kzero = 1
                end if
                if (delg(iligl) == -2 .or. delg(jcoll) == -2) then
                    if (jcoll /= iligl) kzero = 1
                end if
                if (kzero == 1) then
                    lnn = ASTER_FALSE
                    lnn2 = ASTER_FALSE
                end if
            end if
!
! --- SI RESOLUTION SPD DEMANDEE ET TERME NEGATIF OU NUL ON S'ARRETE
            if ((lspd) .and. (kzero == 0)) then
                if ((iligg == jcolg) .and. (raux .lt. epsmac)) then
                    call utmess('F', 'FACTOR_84')
                end if
            end if
!
            if (nuno1 /= 0 .and. nuno2 /= 0) then
                if (prolig == rang) then
                    if (lnn) then
                        iterm = iterm+1
                        if (type == 'S') then
                            smpsk%irn_loc(iterm) = iligg4
                            smpsk%jcn_loc(iterm) = jcolg4
                            smpsk%a_loc(iterm) = real(raux, kind=4)
                        else if (type == 'C') then
                            cmpsk%irn_loc(iterm) = iligg4
                            cmpsk%jcn_loc(iterm) = jcolg4
                            cmpsk%a_loc(iterm) = cmplx(caux, kind=4)
                        else if (type == 'D') then
                            dmpsk%irn_loc(iterm) = iligg4
                            dmpsk%jcn_loc(iterm) = jcolg4
                            dmpsk%a_loc(iterm) = raux
                        else if (type == 'Z') then
                            zmpsk%irn_loc(iterm) = iligg4
                            zmpsk%jcn_loc(iterm) = jcolg4
                            zmpsk%a_loc(iterm) = caux
                        else
                            ASSERT(ASTER_FALSE)
                        end if
!                   Writings to get the stiffness matrix wrt nodes and dof numbers
                        if (ldebug .and. raux /= 0.d0) then
                            owner(iterm) = rang
                            numno1 = deeg(2*(iligl-1)+1)
                            nucmp1 = deeg(2*(iligl-1)+2)
                            numno2 = deeg(2*(jcoll-1)+1)
                            nucmp2 = deeg(2*(jcoll-1)+2)
                            write (11+rang, *) numno1, nucmp1, numno2, nucmp2, raux, &
                                iligg4, jcolg4, iligl, jcoll, prolig, procol, rang
                        end if
                    end if
                    if (procol == rang) then
                        if (iligg /= jcolg) then
                            if (lnn2) then
                                iterm = iterm+1
                                if (type == 'S') then
                                    smpsk%irn_loc(iterm) = jcolg4
                                    smpsk%jcn_loc(iterm) = iligg4
                                    smpsk%a_loc(iterm) = real(raux2, kind=4)
                                else if (type == 'C') then
                                    cmpsk%irn_loc(iterm) = jcolg4
                                    cmpsk%jcn_loc(iterm) = iligg4
                                    cmpsk%a_loc(iterm) = cmplx(caux2, kind=4)
                                else if (type == 'D') then
                                    dmpsk%irn_loc(iterm) = jcolg4
                                    dmpsk%jcn_loc(iterm) = iligg4
                                    dmpsk%a_loc(iterm) = raux2
                                else if (type == 'Z') then
                                    zmpsk%irn_loc(iterm) = jcolg4
                                    zmpsk%jcn_loc(iterm) = iligg4
                                    zmpsk%a_loc(iterm) = caux2
                                else
                                    ASSERT(ASTER_FALSE)
                                end if
!                               Writings to get the stiffness matrix wrt nodes and dof numbers
                                if (ldebug .and. raux2 /= 0.d0) then
                                    owner(iterm) = rang
                                    numno1 = deeg(2*(iligl-1)+1)
                                    nucmp1 = deeg(2*(iligl-1)+2)
                                    numno2 = deeg(2*(jcoll-1)+1)
                                    nucmp2 = deeg(2*(jcoll-1)+2)
                                    write (11+rang, *) numno2, nucmp2, numno1, nucmp1, raux2, &
                                        jcolg4, iligg4, jcoll, iligl, prolig, procol, rang
                                end if
                            end if
                        end if
                    end if
                else if (procol == rang) then
                    if (lnn2) then
                        iterm = iterm+1
                        if (type == 'S') then
                            smpsk%irn_loc(iterm) = jcolg4
                            smpsk%jcn_loc(iterm) = iligg4
                            smpsk%a_loc(iterm) = real(raux2, kind=4)
                        else if (type == 'C') then
                            cmpsk%irn_loc(iterm) = jcolg4
                            cmpsk%jcn_loc(iterm) = iligg4
                            cmpsk%a_loc(iterm) = cmplx(caux2, kind=4)
                        else if (type == 'D') then
                            dmpsk%irn_loc(iterm) = jcolg4
                            dmpsk%jcn_loc(iterm) = iligg4
                            dmpsk%a_loc(iterm) = raux2
                        else if (type == 'Z') then
                            zmpsk%irn_loc(iterm) = jcolg4
                            zmpsk%jcn_loc(iterm) = iligg4
                            zmpsk%a_loc(iterm) = caux2
                        else
                            ASSERT(ASTER_FALSE)
                        end if
!                       Writings to get the stiffness matrix wrt nodes and dof numbers
                        if (ldebug .and. raux2 /= 0.d0) then
                            owner(iterm) = rang
                            numno1 = deeg(2*(iligl-1)+1)
                            nucmp1 = deeg(2*(iligl-1)+2)
                            numno2 = deeg(2*(jcoll-1)+1)
                            nucmp2 = deeg(2*(jcoll-1)+2)
                            write (11+rang, *) numno2, nucmp2, numno1, nucmp1, raux2, &
                                jcolg4, iligg4, jcoll, iligl, prolig, procol, rang
                        end if
                    end if
                end if
            else
!               Si on est sur un ddl de Lagrange et qu'on possede le ddl d'en face
!               ou que les deux ddl sont de Lagrange, on doit donner le terme
                lgive = (nuno1 == 0 .and. procol == rang) .or. &
                        (nuno2 == 0 .and. prolig == rang) .or. &
                        (nuno1 == 0 .and. nuno2 == 0)
                if (lgive) then
                    if (lnn) then
                        iterm = iterm+1
                        if (type == 'S') then
                            smpsk%irn_loc(iterm) = iligg4
                            smpsk%jcn_loc(iterm) = jcolg4
                            smpsk%a_loc(iterm) = real(raux, kind=4)
                        else if (type == 'C') then
                            cmpsk%irn_loc(iterm) = iligg4
                            cmpsk%jcn_loc(iterm) = jcolg4
                            cmpsk%a_loc(iterm) = cmplx(caux, kind=4)
                        else if (type == 'D') then
                            dmpsk%irn_loc(iterm) = iligg4
                            dmpsk%jcn_loc(iterm) = jcolg4
                            dmpsk%a_loc(iterm) = raux
                        else if (type == 'Z') then
                            zmpsk%irn_loc(iterm) = iligg4
                            zmpsk%jcn_loc(iterm) = jcolg4
                            zmpsk%a_loc(iterm) = caux
                        else
                            ASSERT(ASTER_FALSE)
                        end if

!                       Writings to get the stiffness matrix wrt nodes and dof numbers
                        if (ldebug .and. raux /= 0.d0) then
                            owner(iterm) = rang
                            numno1 = deeg(2*(iligl-1)+1)
                            nucmp1 = deeg(2*(iligl-1)+2)
                            numno2 = deeg(2*(jcoll-1)+1)
                            nucmp2 = deeg(2*(jcoll-1)+2)
                            write (11+rang, *) numno1, nucmp1, numno2, nucmp2, raux, &
                                iligg4, jcolg4, iligl, jcoll
                        end if
                    end if
                    if (iligg /= jcolg) then
                        if (lnn2) then
                            iterm = iterm+1
                            if (type == 'S') then
                                smpsk%irn_loc(iterm) = jcolg4
                                smpsk%jcn_loc(iterm) = iligg4
                                smpsk%a_loc(iterm) = real(raux2, kind=4)
                            else if (type == 'C') then
                                cmpsk%irn_loc(iterm) = jcolg4
                                cmpsk%jcn_loc(iterm) = iligg4
                                cmpsk%a_loc(iterm) = cmplx(caux2, kind=4)
                            else if (type == 'D') then
                                dmpsk%irn_loc(iterm) = jcolg4
                                dmpsk%jcn_loc(iterm) = iligg4
                                dmpsk%a_loc(iterm) = raux2
                            else if (type == 'Z') then
                                zmpsk%irn_loc(iterm) = jcolg4
                                zmpsk%jcn_loc(iterm) = iligg4
                                zmpsk%a_loc(iterm) = caux2
                            else
                                ASSERT(ASTER_FALSE)
                            end if
!                        Writings to get the stiffness matrix wrt nodes and dof numbers
                            if (ldebug .and. raux2 /= 0.d0) then
                                owner(iterm) = rang
                                numno1 = deeg(2*(iligl-1)+1)
                                nucmp1 = deeg(2*(iligl-1)+2)
                                numno2 = deeg(2*(jcoll-1)+1)
                                nucmp2 = deeg(2*(jcoll-1)+2)
                                write (11+rang, *) numno2, nucmp2, numno1, nucmp1, raux2, &
                                    jcolg4, iligg4, jcoll, iligl
                            end if
                        end if
                    end if
                end if
            end if
        end do
    end do
!
    ASSERT(iterm .le. nz2)
    nz2 = to_mumps_int(iterm)
    if (type == 'S') then
        smpsk%nnz_loc = nz2
    else if (type == 'C') then
        cmpsk%nnz_loc = nz2
    else if (type == 'D') then
        dmpsk%nnz_loc = nz2
    else if (type == 'Z') then
        zmpsk%nnz_loc = nz2
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jelibe(nonu//'.SMOS.SMDI')
    call jelibe(nonu//'.SMOS.SMHC')
    call jelibe(jexnum(nomat//'.VALM', 1_8))
    if (lmnsy) call jelibe(jexnum(nomat//'.VALM', 2_8))
    call jelibe(nonu//'.NUME.DELG')
!
    if (niv >= 2) then
! --- TEST POUR EVITER LE MONITORING DES CMDES ECLATEES ET DE LDLT_SP
!     LES OBJETS TEMPORAIRES DE MONITORING SONT EFFACES A CHAQUE
!     FIN DE COMMANDE (NUM_DDL/FACTORISER/RESOUDRE)
        call getres(k8bid, k16bid, nomcmd)
        if ((nomcmd(1:8) == 'NUME_DDL') .or. (nomcmd(1:10) == 'FACTORISER') .or. &
            (nomcmd(1:8) == 'RESOUDRE')) then
            lcmde = ASTER_TRUE
        else
            lcmde = ASTER_FALSE
        end if
        if ((.not. lcmde) .and. (.not. lpreco)) then
            call jeexin(kmonit(1), iret)
            if (iret /= 0) then
! --- CAS CMDE STD AVEC MUMPS SOLVEUR DIRECT
                call jeveuo(kmonit(1), 'E', iad)
                zi(iad+rang) = nz2
            else
! --- L'OBJET KMONIT(1) DEVRAIT EXISTER
                ASSERT(ASTER_FALSE)
            end if
        end if
    end if
!
!
!       ------------------------------------------------
!       IMPRESSION DE LA MATRICE (SI DEMANDEE) :
!       ------------------------------------------------
    if (impr(1:3) == 'OUI') then
        write (ifmump+rang, *) sym, '   : SYM', rang, '   : RANG'
        write (ifmump+rang, *) neql, '   : N'
        write (ifmump+rang, *) nz2, '   : NZ_loc'
        if (type == 'S') then
            do k = 1, nz2
                write (ifmump, *) smpsk%irn_loc(k), smpsk%jcn_loc(k), smpsk%a_loc(k)
            end do
        else if (type == 'C') then
            do k = 1, nz2
                write (ifmump, *) cmpsk%irn_loc(k), cmpsk%jcn_loc(k), cmpsk%a_loc(k)
            end do
        else if (type == 'D') then
            do k = 1, nz2
                write (ifmump+rang, *) dmpsk%irn_loc(k), dmpsk%jcn_loc(k), dmpsk%a_loc(k)
            end do
        else if (type == 'Z') then
            do k = 1, nz2
                write (ifmump, *) zmpsk%irn_loc(k), zmpsk%jcn_loc(k), zmpsk%a_loc(k)
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
        write (ifmump+rang, *) 'MUMPS FIN A_loc'

!  -------   VIDANGE DES BUFFERS D'IMPRESSION
        flush (ifmump+rang)

    end if
!
    if (ldebug) then
        do k = 1, iterm
            write (41+rang, *) dmpsk%irn_loc(k), dmpsk%jcn_loc(k), dmpsk%a_loc(k), owner(k)
        end do
        flush (41+rang)
    end if
!
! --- VIDANGE DES BUFFERS D'IMPRESSION
    if (ldebug) flush (11+rang)
!
! --- COMMUNICATION DU VECTEUR KSIZEMU A TOUS LES PROCS
    call asmpi_comm_jev('MPI_SUM', ksizemu)
!
! --- NETTOYAGE VECTEURS TEMPORAIRES LOCAUX
    call jedetr(kfiltr)
    call jedetr(kpiv)
    call jedetr(kpiv2)
    if (ldebug) deallocate (owner)
!
! --- DECHARGEMENT CIBLE D'OBJETS JEVEUX
    call jelibd(nonu//'.SMOS.SMDI', ltot)
    call jelibd(nonu//'.SMOS.SMHC', ltot)
    call jelibd(nonu//'.NUME.DEEQ', ltot)
    call jelibd(nonu//'.NUME.NUEQ', ltot)
    call jelibd(nonu//'.NUME.LILI', ltot)
    call jelibd(nonu//'.NUME.DELG', ltot)
    call jelibd(nonu//'.NUME.NEQU', ltot)
    call jelibd(nonu//'.NUME.PDDL', ltot)
    call jelibd(nonu//'.NUME.NULG', ltot)
    call jelibd(jexnum(nomat//'.VALM', 1_8), ltot)
    if (lmnsy) call jelibd(jexnum(nomat//'.VALM', 2_8), ltot)
!
!
    call jedema()
#endif
end subroutine
