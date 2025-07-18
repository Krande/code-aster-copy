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

subroutine tldlg3(metrez, renum, istop, lmat, ildeb, &
                  ilfin, ndigit, ndeci, isingu, npvneg, &
                  iret, solvop)
    implicit none
! person_in_charge: jacques.pellet at edf.fr
!    BUT  FACTORISER UNE MATRICE ASSEMBLEE
!         DIAGNOSTIQUER LES SINGULARITES OU LE NBRE DE PIVOTS NEGATIFS
!         POUR LES SOLVEURS LINEAIRES: LDLT, MULT_FRONT, MUMPS
!
!     IN  METRES :  /'LDLT' /'MULT_FRONT'/'MUMPS'
!     IN  RENUM :  /'MD' /'MDA' /' ' (sert a MULT_FRONT)
!     IN  ISTOP :  /0 -> SI IRET>0 : ERREUR <F>
!                  /1 -> SI IRET=1 : ALARME <A>
!                        SI IRET=2 : ERREUR <F>
!                  /2 -> LE PROGRAMME NE S'ARRETE PAS
!                        SI IRET>0 : INFO <I>
!     IN  LMAT  : DESCRIPTEUR DE LA MATRICE A FACTORISER
!     IN  ILDEB : NUMERO DE LA LIGNE DE DEPART DE FACTORISATION
!     IN  ILFIN : NUMERO DE LA LIGNE DE FIN    DE FACTORISITION
!    OUT  IRET : CODE RETOUR :
!                  /0 -> OK
!                  /1 -> LE NOMBRE DE DECIMALES PERDUES SUR LE
!                        TERME DIAGONAL DE L'EQUATION ISINGU
!                        EST SUPERIEUR A NDIGIT
!                  /2 -> LA FACTORISATION N'A PAS PU SE FAIRE
!                        JUSQU'AU BOUT.(ARRET A LA LIGNE ISINGU)
!                        SI UN PIVOT DEVIENT (EN MODULE) INFERIEUR
!                        A EPS=/1./R8GAEM()
!    OUT  NPVNEG : NOMBRE DE PIVOTS NEGATIFS SUR LA MATRICE
!                  FACTORISEE.
!                  CE NOMBRE N'A DE SENS QUE SI LA MATRICE EST
!                  DE TYPE REEL ET QUE IRET<2
!     IN  NDIGIT: NOMBRE MAX DE DECIMALES A PERDRE SUR LES TERMES
!                 DIAGONAUX DE LA MATRICE
!              SI NDIGIT <0 ON NE TESTE PAS LA SINGULARITE AVEC MUMPS
!                 SI NDIGIT=0 ON PREND LA VALEUR PAR DEFAUT :8
!                 SI NDIGIT EST GRAND (99 PAR EX.) ON N'AURA
!                    JAMAIS D'ALARME.
!    OUT  NDECI : NOMBRE MAX DE DECIMALES PERDUES SUR LES TERMES
!                 DIAGONAUX DE LA MATRICE (OUTPUT ACTIVE SI
!                 NDIGIT >=0 ET SI METRES NON MUMPS)
!    OUT  ISINGU: NUMERO DE L'EQUATION OU LA PERTE DE DECIMALES
!                 EST MAXIMUM OU BIEN NUMERO DE L'EQUATION POUR
!                 LA QUELLE LA FACTORISATION S'EST ARRETEE
!    IN SOLVOP: SD_SOLVEUR DE L'OPERATEUR (PARFOIS DIFFERENT DE CELUI
!               ASSOCIE AU MATRICE). CELA SERT UNIQUEMENT A
!               MUMPS POUR LEQUEL SEUL CE JEU DE PARAMETRES FAIT FOI SI
!               IL EST DIFFERENT DE CELUI DES MATRICES.
!     ------------------------------------------------------------------
!
!
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/amumph.h"
#include "asterfort/assert.h"
#include "asterfort/diagav.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mlfc16.h"
#include "asterfort/mtmchc.h"
#include "asterfort/mulfr8.h"
#include "asterfort/rgndas.h"
#include "asterfort/tldlc8.h"
#include "asterfort/tldlr8.h"
#include "asterfort/tlduc8.h"
#include "asterfort/tldur8.h"
#include "asterfort/ualfcr.h"
#include "asterfort/utmess.h"
#include "asterfort/isParallelMatrix.h"

    character(len=1) :: codmes
    character(len=19) :: noma19, stolci
    character(len=*) :: solvop
    character(len=14) :: nu
    character(len=*) :: metrez, renum
    character(len=16) :: metres
    character(len=24) :: kpiv
    character(len=40) ::  valk(2)
    character(len=3) :: mathpc
    integer(kind=8) :: istop, lmat, ildeb, ilfin, ndigit, ndigi2, iret, npvneg, iretz
    integer(kind=8) :: ifm, niv, nom, neq, iretp, npvnez, neqg, jnequ
    integer(kind=8) :: typvar, typsym, nbbloc, ilfin1, iexi
    integer(kind=8) :: ieq3, isingu, ieq, ndeci, jdigs, npivot
    integer(kind=8) :: ndeci1, ndeci2, ieq4, nzero, vali(6), ipiv
    real(kind=8) :: eps, dmax, dmin, d1
    aster_logical :: l_parallel_matrix
    complex(kind=8) :: cbid
    integer(kind=8), pointer :: schc(:) => null()
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer :: scbl(:) => null()
    integer(kind=8), pointer :: scdi(:) => null()
    integer(kind=8), pointer :: lc2m(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
    nom = zi(lmat+1)
    neq = zi(lmat+2)
    typvar = zi(lmat+3)
    cbid = dcmplx(0.d0, 0.d0)
    typsym = zi(lmat+4)
    noma19 = zk24(nom) (1:19)
    metres = metrez
    ieq4 = 0
    iretz = 0
    npivot = 0
!
!
!
    if (metres .ne. 'LDLT' .and. metres .ne. 'MULT_FRONT' .and. metres .ne. 'MUMPS') then
        call utmess('F', 'ALGELINE4_1')
    end if
!
!   -- DDLS ELIMINES :
    call jeveuo(noma19//'.REFA', 'L', vk24=refa)
    ASSERT(refa(3) .ne. 'ELIMF')
    if (refa(3) .eq. 'ELIML') call mtmchc(noma19, 'ELIMF')
    ASSERT(refa(3) .ne. 'ELIML')
!
    call dismoi('NOM_NUME_DDL', noma19, 'MATR_ASSE', repk=nu)
    ASSERT(nu .ne. ' ')
    ASSERT(refa(2) (1:14) .eq. nu)
!
    call infdbg('FACTOR', ifm, niv)
    if (niv .eq. 2) then
        valk(1) = noma19
        valk(2) = metres
        call utmess('I', 'FACTOR3_9', nk=2, valk=valk)
        if (typsym .eq. 1) then
            call utmess('I', 'FACTOR3_2')
        else
            call utmess('I', 'FACTOR3_10')
        end if
        if (typvar .eq. 1) then
            call utmess('I', 'FACTOR3_4')
        else
            call utmess('I', 'FACTOR3_11')
        end if
    end if
!
!   -- EST-ON EN HPC :
    l_parallel_matrix = isParallelMatrix(noma19)
    neqg = -1
    if (l_parallel_matrix) then
        call jeveuo(nu//'.NUME.NEQU', 'L', jnequ)
        neqg = zi(jnequ+1)
    end if
!
!   -- VALEUR DE NDIGIT PAR DEFAUT : 8
    if (ndigit .eq. 0) then
        ndigi2 = 8
    else
        ndigi2 = ndigit
    end if

!   -- ON NE PERMET PAS LE DEBRANCHEMENT DE LA RECHERCHE DE SINGU
!        LARITE AVEC LDLT ET MULT_FRONT (POUR L'INSTANT)
    if (metres .ne. 'MUMPS') ndigi2 = abs(ndigi2)
!
    if (ilfin .lt. ildeb .or. ilfin .gt. neq) then
        ilfin1 = neq
    else
        ilfin1 = ilfin
    end if
!
!   ON ALLOUE (SI NECESSAIRE) UN VECTEUR QUI CONTIENDRA
!   LA DIAGONALE "AVANT" ET LA DIAGONALE "APRES" :
    if (metres .ne. 'MUMPS') call diagav(noma19, neq, ilfin1, typvar, eps)
!
!
    if (metres .eq. 'LDLT') then
!   ---------------------------------------
!       -- ALLOCATION DE LA MATRICE FACTORISEE (.UALF)  ET RECOPIE
!          DE .VALM DANS .UALF
        if ((noma19 .ne. '&&OP0070.RESOC.MATC') .and. (noma19 .ne. '&&OP0070.RESUC.MATC')) then
            call ualfcr(noma19, ' ')
        end if
        call jelira(noma19//'.UALF', 'NMAXOC', nbbloc)
!
        stolci = nu//'.SLCS'
        call jeveuo(stolci//'.SCDI', 'L', vi=scdi)
        call jeveuo(stolci//'.SCBL', 'L', vi=scbl)
        call jeveuo(stolci//'.SCHC', 'L', vi=schc)
        if (typvar .eq. 1) then
            if (typsym .eq. 1) then
                call tldlr8(noma19, schc, scdi, scbl, npivot, &
                            neq, nbbloc, ildeb, ilfin1, eps)
            else if (typsym .eq. 0) then
                call tldur8(noma19, schc, scdi, scbl, npivot, &
                            neq, nbbloc/2, ildeb, ilfin1, eps)
            end if
!
        else if (typvar .eq. 2) then
            if (typsym .eq. 1) then
                call tldlc8(noma19, schc, scdi, scbl, npivot, &
                            neq, nbbloc, ildeb, ilfin1, eps)
            else if (typsym .eq. 0) then
                call tlduc8(noma19, schc, scdi, scbl, npivot, &
                            neq, nbbloc/2, ildeb, ilfin1, eps)
            end if
        end if

    else if (metres .eq. 'MULT_FRONT') then
!   ---------------------------------------
        if (typvar .eq. 1) then
            call mulfr8(noma19, npivot, neq, typsym, eps, &
                        renum)
        else if (typvar .eq. 2) then
            call mlfc16(noma19, npivot, neq, typsym, eps, &
                        renum)
        end if
!
!
    else if (metres .eq. 'MUMPS') then
!   ---------------------------------------
        call amumph('DETR_OCC', solvop, noma19, [0.d0], [cbid], &
                    ' ', 0, iretz, .true._1)
        call amumph('PRERES', solvop, noma19, [0.d0], [cbid], &
                    ' ', 0, iretz, .true._1)
!
!       -- mumps ne nous dit pas le nombre de decimales reellement perdues :
        ndeci = -999
!
        nzero = -999
        iretp = 0
        kpiv = '&&AMUMP.PIVNUL'
        if (iretz .eq. 2) then
!     -- LA FACTORISATION NE S'EST PAS BIEN PASSEE. PEUT IMPORTE LA VALEUR
!        NPREC ET L'ACTIVATION OU NON DE LA RECHERCHE DE SINGULARITE.
!        MATRICE SINGULIERE NUMERIQUEMENT OU EN STRUCTURE (DETECTE EN
!        AMONT DS AMUMPH. ON NE SAIT PAS PRECISER ISINGU CONTRAIREMENT A
!        MF/LDLT. ON MET ISINGU=-999 ARBITRAIREMENT)
            isingu = -999
            npivot = -999
        else
            if (ndigi2 .gt. 0) then
!     -- ON RECUPERE LE TABLEAU DONNANT SUR LA LISTE DES PIVOTS NON-NULS
!        IL N'EXISTE QUE SI NPREC>=0 ET SI IRETZ=0
                call jeexin(kpiv, iretp)
                if (iretp .ne. 0) then
                    call jeveuo(kpiv, 'L', ipiv)
                else
                    ASSERT(.false.)
                end if
!    -- LE PREMIER ELEMENT DU TABLEAU CORRESPOND A INFOG(28)
!    -- IL INDIQUE LE NOMBRE DE PIVOTS INFERIEUR A UN CERTAIN SEUIL
!       DEFINI DANS AMUMPR
!    -- SI INFOG(28) > 0 ALORS IRETZ=1 ( LE NOMBRE DE DECIMALES PERDUES
!         SUR LE TERME DIAGONAL DE L'EQUATION ISINGU> A NDIGIT)
!    -- ATTENTION ON N'EST PAS RIGOUREUSEMENT IDENTIQUE AU CRITERE
!       HABITUEL EMPLOYE AVEC MF ET LDLT. AVEC MUMPS, LE CRITERE
!          - UTILISE LA NORME INFINIE DE LA LIGNE OU DE LA COLONNE
!            DU PIVOT ET NON PAS EXPLICITEMENT LE RAPPORT DE TERMES
!            DIAGONAUX
!          - EST GLOBAL A TOUTE LA MATRICE ET NON LOCAL PAR LIGNE
!          - ON NE DETECTE PAS LE NUMERO DE LIGNE DE PIVOT VRAIMENT NUL
!            (CAS IRETZ=2 PAS EXPLOIE ICI MAIS DIRECTEMENT DS AMUMPR/C
!             AVEC LES ERREURS MUMPS INFO(1)=-10)
!
!               -- LA FACTORISATION S'EST BIEN PASSEE. ON CHERCHE LES SINGULARITES
                if (zi(ipiv) .eq. 0) then
!                   -- PAS DE SINGULARITE
                    iretz = 0
                    isingu = -999
                    npivot = -zi(ipiv+1)
                else if (zi(ipiv) .gt. 0) then
!                   -- AU MOINS UNE SINGULARITE
                    iretz = 1
                    isingu = zi(ipiv+2)
                    if (l_parallel_matrix) then
                        ASSERT(isingu .gt. 0 .and. isingu .le. neqg)
                    else
                        ASSERT(isingu .gt. 0 .and. isingu .le. neq)
                    end if
                    npivot = -zi(ipiv+1)
                else
                    ASSERT(.false.)
                end if
            else
!               -- LA FACTO S'EST BIEN PASSEE ET ON NE CHERCHE PAS A TESTER LES
!                  EVENTUELLES SINGULARITES
                isingu = -999
                npivot = -999
            end if
        end if
        if (iretp .ne. 0) call jedetr(kpiv)
    end if
!
!
!     -- CALCUL DE NPVNEG :
!     ---------------------
    if (npivot .lt. 0) then
        npvnez = npivot
    else
        npvnez = 0
    end if
!
!
!   -- Calcul du code retour: iretz, ndeci et isingu:
!   ------------------------------------------------
    if (metres(1:5) .ne. 'MUMPS') then
        if (npivot .gt. 0) then
            iretz = 2
            ndeci = -999
            isingu = npivot
            ASSERT(isingu .gt. 0 .and. isingu .le. neq)
        else
!
!           -- On regarde ce que sont devenus les termes diagonaux :
!           -------------------------------------------------------
            call jeveuo(noma19//'.DIGS', 'L', jdigs)
            dmax = 0.d0
            dmin = r8maem()
            nzero = 0
            do ieq = ildeb, ilfin1
                if (typvar .eq. 1) then
                    d1 = abs(zr(jdigs-1+ieq)/zr(jdigs+neq-1+ieq))
                else
                    d1 = abs(zc(jdigs-1+ieq)/zc(jdigs+neq-1+ieq))
                end if
                if (d1 .gt. dmax) then
                    dmax = d1
                    ieq3 = ieq
                end if
                if (d1 .eq. 0.d0) then
                    nzero = nzero+1
                else
                    if (d1 .lt. dmin) then
                        dmin = d1
                        ieq4 = ieq
                    end if
                end if
            end do
            ASSERT(dmax .gt. 0.d0)
            ndeci1 = int(log10(dmax))
            ndeci2 = int(log10(1.d0/dmin))
            ndeci = ndeci1
            isingu = ieq3
            ASSERT(isingu .gt. 0 .and. isingu .le. neq)
            if (ndeci .ge. ndigi2) then
                iretz = 1
            else
                iretz = 0
            end if
        end if
    end if

!   -- Emission eventuelle d'un message d'erreur :
!   ----------------------------------------------
    if ((ndigi2 .lt. 0) .and. (metres .eq. 'MUMPS')) goto 30

    if (iretz .eq. 0) then
        goto 20
    else if (istop .eq. 2) then
        codmes = 'I'
    else if (istop .eq. 1) then
        if (iretz .eq. 1) then
            codmes = 'A'
        else if (iretz .eq. 2) then
            codmes = 'F'
        else
            ASSERT(.false.)
        end if
    else if (istop .eq. 0) then
        codmes = 'F'
    end if

!   -- si LDLT et si stolci.LC2M existe et si isingu > 0, il faut en tenir compte :
    if (metres .eq. 'LDLT') then
        if (isingu .gt. 0) then
            call jeexin(stolci//'.LC2M', iexi)
            if (iexi .gt. 0) then
                call jeveuo(stolci//'.LC2M', 'L', vi=lc2m)
                isingu = lc2m(isingu)
            end if
        end if
    end if

    ASSERT(isingu .eq. -999 .or. isingu .gt. 0)
    vali(1) = isingu

    ASSERT(ndeci .eq. -999 .or. ndeci .ge. 0)
    if (isingu .eq. -999) then
        ASSERT(ndeci .eq. -999)
    end if
    vali(2) = ndeci
!
! - Error
!
    if (isingu .eq. -999) then
        call utmess(codmes, 'FACTOR_12')
    end if
    if (isingu .gt. 0) then
        if (refa(20) == '') then
            if (l_parallel_matrix) then
                call utmess('I', 'FACTOR2_7')
            else
                call rgndas(nu, isingu, l_print=.true.)
            end if
        end if
        if (ndeci .eq. -999) then
            call utmess(codmes, 'FACTOR_11', si=isingu)
        else
            call utmess(codmes, 'FACTOR_10', ni=2, vali=vali)
        end if
    end if
20  continue

!   -- impressions info=2 :
!   ------------------------
    if (niv .eq. 2) then
        call utmess('I', 'FACTOR3_12', sk=noma19)
        if (nzero .gt. 0) then
            call utmess('I', 'FACTOR3_13', si=nzero)
        end if
        vali(1) = ndigi2
        vali(2) = ndeci
        vali(3) = isingu
        vali(4) = -npvnez
        vali(5) = istop
        vali(6) = iretz
        call utmess('I', 'FACTOR3_14', ni=6, vali=vali)
!
!       -- ALARME EVENTUELLE SI LE PIVOT DEVIENT TROP GRAND :
        if ((metres .ne. 'MUMPS') .and. (ndeci2 .ge. ndigi2)) then
            ASSERT(ieq4 .gt. 0 .and. ieq4 .le. neq)
            if (refa(20) == '') then
                call rgndas(nu, ieq4, l_print=.true.)
            end if
            vali(1) = ieq4
            vali(2) = ndeci2
            call utmess('I', 'FACTOR3_15', ni=2, vali=vali)
        end if
    end if
!
!
30  continue
!
!
    iret = iretz
    npvneg = npvnez
    call jedema()
!
end subroutine
