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
subroutine ircmva(numcmp, indcmp, ncmpve, ncmprf, nvalec, &
                  nbpg, nbsp, adsv, adsd, adsl, &
                  adsk, cplxFormatZ, tymast, modnum, nuanom, &
                  fieldSupport, val, profas, ideb, ifin, &
                  codret)
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterc/r8rddg.h"
#include "asterfort/cesexi.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
#include "asterfort/jeveuo.h"
!
    integer(kind=8) :: ncmpve, ncmprf, nvalec, nbpg, nbsp
    integer(kind=8) :: numcmp(ncmprf)
    integer(kind=8) :: adsv, adsd, adsl, adsk
    integer(kind=8) :: tymast, codret
    integer(kind=8) :: modnum(MT_NTYMAX), nuanom(MT_NTYMAX, *)
    integer(kind=8) :: profas(*)
    integer(kind=8) :: ideb, ifin
    real(kind=8) :: val(ncmpve, nbsp, nbpg, nvalec)
    character(len=8), intent(in) :: fieldSupport
    character(len=*), intent(in) :: cplxFormatZ
    character(len=24), intent(in) :: indcmp
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE D'UN CHAMP -  FORMAT MED - CREATION DES VALEURS
!
! --------------------------------------------------------------------------------------------------
!
!     ENTREES :
!       NUMCMP : NUMEROS DES COMPOSANTES
!       NCMPVE : NOMBRE DE COMPOSANTES VALIDES EN ECRITURE
!       NVALEC : NOMBRE DE VALEURS A ECRIRE
!       NBPG   : NOMBRE DE POINTS DE GAUSS (1 POUR DES CHAMNO)
!       NBSP   : NOMBRE DE SOUS-POINTS (1 POUR DES CHAMNO)
!       TYPECH : TYPE DE CHAMP (ELEM,ELNO,ELGA,NOEU)
!       ADSV,D,L,K : ADRESSES DES TABLEAUX DES CHAMPS SIMPLIFIES
! In  cplxFormat       : format of complex numbers (IMAG, REAL, PHASE, MODULE or ' ')
!       TYMAST : TYPE ASTER DE MAILLE QUE L'ON VEUT (0 POUR LES NOEUDS)
!       MODNUM : INDICATEUR SI LA SPECIFICATION DE NUMEROTATION DES
!                NOEUDS DES MAILLES EST DIFFERENTES ENTRE ASTER ET MED:
!                     MODNUM = 0 : NUMEROTATION IDENTIQUE
!                     MODNUM = 1 : NUMEROTATION DIFFERENTE
!       NUANOM : TABLEAU DE CORRESPONDANCE DES NOEUDS MED/ASTER.
!                NUANOM(ITYP,J): NUMERO DANS ASTER DU J IEME NOEUD DE LA
!                MAILLE DE TYPE ITYP DANS MED.
!       PROFAS : PROFIL ASTER. C'EST LA LISTE DES NUMEROS ASTER
!                DES NOEUDS OU DES ELEMENTS POUR LESQUELS LE CHAMP
!                EST DEFINI
!       IDEB   : INDICE DE DEBUT DANS PROFAS
!       IFIN   : INDICE DE FIN DANS PROFAS
!     SORTIES :
!       VAL    : VALEURS EN MODE ENTRELACE
!       CODRET : CODE RETOUR, S'IL VAUT 100, IL Y A DES COMPOSANTES
!                 MISES A ZERO
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: gd, valk(2), typcha, cplxFormat
    integer(kind=8) :: iaux, jaux, kaux, fieldScalar
    integer(kind=8) :: adsvxx, adslxx
    integer(kind=8) :: ino, ima, nrcmp, nrcmpr, nrpg, nrsp
    integer(kind=8) :: ifm, niv, jindcm, poscmp, cntcmp
    aster_logical :: logaux, lprolz
    real(kind=8) :: end, start
!
! --------------------------------------------------------------------------------------------------
!
    lprolz = .false.
    gd = zk8(adsk-1+2)
    codret = 0
    cplxFormat = cplxFormatZ
!
    call dismoi('TYPE_SCA', gd, 'GRANDEUR', repk=typcha)
    call jeveuo(indcmp, 'L', jindcm)
!
    if (cplxFormat .eq. 'REEL') then
        fieldScalar = 2
    else if (cplxFormat .eq. 'IMAG') then
        fieldScalar = 3
    else if (cplxFormat .eq. 'MODULE') then
        fieldScalar = 4
    else if (cplxFormat .eq. 'PHASE') then
        fieldScalar = 5
    else
        if (typcha .eq. 'R') then
            fieldScalar = 1
        else if (typcha .eq. 'C') then
            call utmess('F', 'RESULT3_69')
        end if
    end if

    if ((typcha .ne. 'R') .and. (typcha .ne. 'C')) then
        valk(1) = gd
        valk(2) = 'IRCMVA'
        call utmess('F', 'DVP_3', nk=2, valk=valk)
    end if
!
! 1.1. ==> RECUPERATION DU NIVEAU D'IMPRESSION
!
    call infniv(ifm, niv)
    if (niv .gt. 1) then
        call cpu_time(start)
        write (ifm, *) '<IRCMVA> DEBUT DE CREATION DES VALEURS'
    end if
!
! 1.2. ==> INFORMATION
!
    if (niv .gt. 1) then
        call utmess('I', 'MED_47')
        write (ifm, 130) nvalec, ncmpve, nbpg, nbsp, fieldSupport
    end if
130 format('  NVALEC =', i8, ', NCMPVE =', i8, ', NBPG   =', i8, &
           ', NBSP   =', i8, /, '  TYPECH =', a8)
!
!====
! 2. CREATION DU CHAMP DE VALEURS AD-HOC
!    LE TABLEAU DE VALEURS EST UTILISE AINSI :
!        TV(NCMPVE,NBSP,NBPG,NVALEC)
!    EN FORTRAN, CELA CORRESPOND AU STOCKAGE MEMOIRE SUIVANT :
!    TV(1,1,1,1), TV(2,1,1,1), ..., TV(NCMPVE,1,1,1),
!    TV(1,2,1,1), TV(2,2,1,1), ..., TV(NCMPVE,2,1,1),
!            ...     ...     ...
!    TV(1,NBSP,NBPG,NVALEC), TV(2,NBSP,NBPG,NVALEC), ... ,
!                                      TV(NCMPVE,NBSP,NBPG,NVALEC)
!    C'EST CE QUE MED APPELLE LE MODE ENTRELACE
!    ATTENTION : LE CHAMP SIMPLIFIE EST DEJA PARTIELLEMENT FILTRE ...
!    ATTENTION ENCORE : LE CHAMP SIMPLIFIE N'A PAS LA MEME STRUCTURE
!    POUR LES NOEUDS ET LES ELEMENTS. IL FAUT RESPECTER CE TRAITEMENT
!    REMARQUE : SI UNE COMPOSANTE EST ABSENTE, ON AURA UNE VALEUR NULLE
!    REMARQUE : ATTENTION A BIEN REDIRIGER SUR LE NUMERO DE
!    COMPOSANTE DE REFERENCE
!====
!
! 2.1. ==> POUR LES NOEUDS : ON PREND TOUT CE QUI FRANCHIT LE FILTRE
!
    if (tymast .eq. 0) then
        cntcmp = 0
        do nrcmp = 1, ncmpve
            poscmp = zi(jindcm+nrcmp-1)
            if (poscmp .eq. 0) cycle
            cntcmp = cntcmp+1
            adsvxx = adsv-1+numcmp(cntcmp)-ncmprf
            adslxx = adsl-1+numcmp(cntcmp)-ncmprf
            jaux = 0
            do iaux = ideb, ifin
                ino = profas(iaux)
                jaux = jaux+1
                kaux = ino*ncmprf
                if (zl(adslxx+kaux)) then
                    if (fieldScalar .eq. 1) then
                        val(poscmp, 1, 1, jaux) = zr(adsvxx+kaux)
                    else if (fieldScalar .eq. 2) then
                        if (typcha .eq. 'R') then
                            val(poscmp, 1, 1, jaux) = zr(adsvxx+kaux)
                        else if (typcha .eq. 'C') then
                            val(poscmp, 1, 1, jaux) = dble(zc(adsvxx+kaux))
                        end if
                    else if (fieldScalar .eq. 3) then
                        if (typcha .eq. 'R') then
                            val(poscmp, 1, 1, jaux) = 0.d0
                        else if (typcha .eq. 'C') then
                            val(poscmp, 1, 1, jaux) = dimag(zc(adsvxx+kaux))
                        end if
                    else if (fieldScalar .eq. 4) then
                        if (typcha .eq. 'R') then
                            val(poscmp, 1, 1, jaux) = abs(zr(adsvxx+kaux))
                        else if (typcha .eq. 'C') then
                            val(poscmp, 1, 1, jaux) = abs(zc(adsvxx+kaux))
                        end if
                    else if (fieldScalar .eq. 5) then
                        if (typcha .eq. 'R') then
                            val(poscmp, 1, 1, jaux) = atan2(0.d0, zr(adsvxx+kaux)) &
                                                      *r8rddg()
                        else if (typcha .eq. 'C') then
                            val(poscmp, 1, 1, jaux) = atan2(dimag(zc(adsvxx+kaux)), &
                                                            dble(zc(adsvxx+kaux))) &
                                                      *r8rddg()
                        end if
                    end if
                else
                    lprolz = .true.
                    val(poscmp, 1, 1, jaux) = 0.d0
                end if
            end do
        end do
        if (lprolz) then
            codret = 100
        end if
    else
!
! 2.2. ==> POUR LES MAILLES : ON PREND TOUT CE QUI FRANCHIT LE FILTRE
!          ET QUI EST DU TYPE EN COURS
!          REMARQUE : ON NE REDECODE PAS LES NOMBRES DE POINTS DE GAUSS
!          NI DE SOUS-POINT CAR ILS SONT INVARIANTS POUR UNE IMPRESSION
!          DONNE
!          REMARQUE : DANS LE CAS DE CHAMPS AUX NOEUDS PAR ELEMENTS,
!          L'ORDRE DE STOCKAGE DES VALEURS DANS UNE MAILLE DONNEE EST
!          CELUI DE LA CONNECTIVITE LOCALE DE LA MAILLE. OR POUR
!          CERTAINES MAILLES, CET ORDRE CHANGE ENTRE ASTER ET MED. IL
!          FAUT DONC RENUMEROTER.
!
!GN        PRINT *,'PREMIERE MAILLE : ',PROFAS(IDEB)
!GN        PRINT *,'DERNIERE MAILLE : ',PROFAS(IFIN)
!
! 2.2.1. ==> A-T-ON BESOIN DE RENUMEROTER ?
!            REMARQUE : LE MODE DE RANGEMENT FAIT QUE CELA NE FONCTIONNE
!            QUE POUR LES CHAMPS AVEC 1 SEUL SOUS-POINT.
!
        logaux = .false.
        if (fieldSupport(1:4) .eq. 'ELNO') then
            if (modnum(tymast) .eq. 1) then
                logaux = .true.
            end if
        end if
!
        if (logaux) then
            if (nbsp .gt. 1) then
                write (ifm, 130) nvalec, ncmpve, nbpg, nbsp
                call utmess('F', 'MED_48')
            end if
        end if
!
! 2.2.2. ==> TRANSFERT
!            ON FAIT LE TEST AVANT LA BOUCLE 211. IL EST DONC FAIT
!            AUTANT DE FOIS QUE DE COMPOSANTES A TRANSFERER. AU-DELA, CE
!            SERAIT AUTANT DE FOIS QUE DE MAILLES, DONC COUTEUX
!
        cntcmp = 0
        do nrcmp = 1, ncmpve
            poscmp = zi(jindcm+nrcmp-1)
            if (poscmp .eq. 0) cycle
            cntcmp = cntcmp+1
            nrcmpr = numcmp(cntcmp)
            jaux = 0
            if (logaux) then
                nrsp = 1
                do iaux = ideb, ifin
                    ima = profas(iaux)
                    jaux = jaux+1
                    do nrpg = 1, nbpg
                        call cesexi('C', adsd, adsl, ima, nrpg, nrsp, nrcmpr, kaux)
                        if ((kaux .gt. 0)) then
                            if (fieldScalar .eq. 1) then
                                val(poscmp, nrsp, nuanom(tymast, nrpg), jaux) = &
                                    zr(adsv-1+kaux)
                            else if (fieldScalar .eq. 2) then
                                if (typcha .eq. 'R') then
                                    val(poscmp, nrsp, nuanom(tymast, nrpg), jaux) = &
                                        zr(adsv-1+kaux)
                                else if (typcha .eq. 'C') then
                                    val(poscmp, nrsp, nuanom(tymast, nrpg), jaux) = &
                                        dble(zc(adsv-1+kaux))
                                end if
                            else if (fieldScalar .eq. 3) then
                                if (typcha .eq. 'R') then
                                    val(poscmp, nrsp, nuanom(tymast, nrpg), jaux) = &
                                        0.d0
                                else if (typcha .eq. 'C') then
                                    val(poscmp, nrsp, nuanom(tymast, nrpg), jaux) = &
                                        dimag(zc(adsv-1+kaux))
                                end if
                            else if (fieldScalar .eq. 4) then
                                if (typcha .eq. 'R') then
                                    val(poscmp, nrsp, nuanom(tymast, nrpg), jaux) = &
                                        abs(zr(adsv-1+kaux))
                                else if (typcha .eq. 'C') then
                                    val(poscmp, nrsp, nuanom(tymast, nrpg), jaux) = &
                                        abs(zc(adsv-1+kaux))
                                end if
                            else if (fieldScalar .eq. 5) then
                                if (typcha .eq. 'R') then
                                    val(poscmp, nrsp, nuanom(tymast, nrpg), jaux) = &
                                        atan2(0.d0, zr(adsv-1+kaux))*r8rddg()
                                else if (typcha .eq. 'C') then
                                    val(poscmp, nrsp, nuanom(tymast, nrpg), jaux) = &
                                        atan2(dimag(zc(adsv-1+kaux)), dble(zc(adsv-1+kaux))) &
                                        *r8rddg()
                                end if
                            end if
                        end if
                    end do
                end do
            else
                do iaux = ideb, ifin
                    ima = profas(iaux)
                    jaux = jaux+1
                    do nrpg = 1, nbpg
                        do nrsp = 1, nbsp
                            call cesexi('C', adsd, adsl, ima, nrpg, nrsp, nrcmpr, kaux)
                            if ((kaux .gt. 0)) then
                                if (fieldScalar .eq. 1) then
                                    val(poscmp, nrsp, nrpg, jaux) = zr(adsv-1+kaux)
                                else if (fieldScalar .eq. 2) then
                                    if (typcha .eq. 'R') then
                                        val(poscmp, nrsp, nrpg, jaux) = zr(adsv-1+kaux)
                                    else if (typcha .eq. 'C') then
                                        val(poscmp, nrsp, nrpg, jaux) = dble(zc(adsv-1+kaux))
                                    end if
                                else if (fieldScalar .eq. 3) then
                                    if (typcha .eq. 'R') then
                                        val(poscmp, nrsp, nrpg, jaux) = 0.d0
                                    else if (typcha .eq. 'C') then
                                        val(poscmp, nrsp, nrpg, jaux) = dimag(zc(adsv-1+kaux))
                                    end if
                                else if (fieldScalar .eq. 4) then
                                    if (typcha .eq. 'R') then
                                        val(poscmp, nrsp, nrpg, jaux) = abs(zr(adsv-1+kaux))
                                    else if (typcha .eq. 'C') then
                                        val(poscmp, nrsp, nrpg, jaux) = abs(zc(adsv-1+kaux))
                                    end if
                                else if (fieldScalar .eq. 5) then
                                    if (typcha .eq. 'R') then
                                        val(poscmp, nrsp, nrpg, jaux) = &
                                            atan2(0.d0, zr(adsv-1+kaux))*r8rddg()
                                    else if (typcha .eq. 'C') then
                                        val(poscmp, nrsp, nrpg, jaux) = &
                                            atan2(dimag(zc(adsv-1+kaux)), dble(zc(adsv-1+kaux))) &
                                            *r8rddg()
                                    end if
                                end if
                            end if
                        end do
                    end do
                end do
            end if
        end do
    end if
!
    if (niv .gt. 1) then
        call cpu_time(end)
        write (ifm, *) '<IRCMVA> FIN DE CREATION DES VALEURS EN ', end-start, " SEC"
    end if
!
end subroutine
