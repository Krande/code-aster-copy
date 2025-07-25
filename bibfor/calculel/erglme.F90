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
subroutine erglme(jceld, iavale, option, iord, ligrel, &
                  longt, nbgr, resuco, resuc1)
    implicit none
#include "jeveux.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/ltcrsd.h"
#include "asterfort/ltnotb.h"
#include "asterfort/nbelem.h"
#include "asterfort/rsadpa.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
!
    integer(kind=8) :: jceld, iavale, iord, longt, nbgr
    character(len=*) :: resuco
    character(len=19) :: resuc1, ligrel
    character(len=*) :: option
! person_in_charge: josselin.delmas at edf.fr
! =====================================================================
! ERREUR GLOBALE AU MAILLAGE - MECANIQUE
! **     **                    **
! =====================================================================
!
!     BUT:
!         CALCUL DES ESTIMATEURS GLOBAUX POUR LA MECANIQUE
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   JCELD  : ADRESSE DU DESCRIPTEUR DU CHAMP LOCAL
! IN   IAVALE : ADRESSE DES CHAMPS LOCAUX
! IN   OPTION :    'ERZ1_ELEM' OU 'ERZ2_ELEM'
!               OU 'QIZ1_ELEM' OU 'QIZ2_ELEM'
!               OU 'ERME_ELEM' OU 'QIRE_ELEM'
! IN   IORD   : NUMERO D'ORDRE
! IN   LIGREL : NOM DU LIGREL
! IN   LONGT  : NOMBRE DE CHAMPS LOCAUX
! IN   NBGR   : NOMBRE DE GRELS
! IN   RESUCO : NOM DU CONCEPT ENTRANT
! IN   RESUC1 : NOM DU CONCEPT RESULTAT DE LA COMMANDE CALC_ERREUR
!
!      SORTIE :
!-------------
!
! ......................................................................
!
!
!
    integer(kind=8) :: nbpar
    parameter(nbpar=5)
!
    integer(kind=8) :: ifi, nel
    integer(kind=8) :: mode, k, j, iret, iad, idecgr
    integer(kind=8) :: ladpa
!
    real(kind=8) :: err0, nors
    real(kind=8) :: listr(3), nu0, nuvo, nusa, nuno, termvo, termsa, termno
!
    character(len=3) :: typpar(nbpar)
    character(len=8) :: k8bid
    character(len=16) :: nompar(nbpar)
    character(len=19) :: nomt19
!
    complex(kind=8) :: cbid
!
! ----------------------------------------------------------------------
!
    call jemarq()
    cbid = (0.d0, 0.d0)
!
    ifi = iunifi('RESULTAT')
!
    nu0 = 0.d0
    err0 = 0.d0
    nors = 0.d0
    termvo = 0.d0
    termsa = 0.d0
    termno = 0.d0
    do 10, j = 1, nbgr
        mode = zi(jceld-1+zi(jceld-1+4+j)+2)
        if (mode .eq. 0) goto 10
        nel = nbelem(ligrel, j)
        idecgr = zi(jceld-1+zi(jceld-1+4+j)+8)
        do 20, k = 1, nel
            iad = iavale-1+idecgr+(k-1)*longt
            err0 = err0+zr(iad+1-1)**2
            nors = nors+zr(iad+3-1)**2
            termvo = termvo+zr(iad+4-1)**2
            termno = termno+zr(iad+6-1)**2
            termsa = termsa+zr(iad+8-1)**2
20          continue
10      end do
        if (option .eq. 'QIRE_ELEM') then
            err0 = sqrt(err0)
            termvo = sqrt(termvo)
            termsa = sqrt(termsa)
            termno = sqrt(termno)
        else
            if ((err0+nors) .gt. 0.d0) then
                nu0 = 100.d0*sqrt(err0/(err0+nors))
            else
                nu0 = 0.d0
            end if
            err0 = sqrt(err0)
            if ((termvo+nors) .gt. 0.d0) then
                nuvo = 100.d0*sqrt(termvo/(termvo+nors))
            else
                nuvo = 0.d0
            end if
            termvo = sqrt(termvo)
            if ((termsa+nors) .gt. 0.d0) then
                nusa = 100.d0*sqrt(termsa/(termsa+nors))
            else
                nusa = 0.d0
            end if
            termsa = sqrt(termsa)
            if ((termno+nors) .gt. 0.d0) then
                nuno = 100.d0*sqrt(termno/(termno+nors))
            else
                nuno = 0.d0
            end if
            termno = sqrt(termno)
            nors = sqrt(nors)
        end if
!
! ON ECRIT LES TERMES DANS LES PARAMETRES DE LA SD RESULTAT
!
        if (option .eq. 'ERME_ELEM') then
            call rsadpa(resuc1, 'E', 1, 'ERREUR_ERRE', iord, &
                        0, sjv=ladpa, styp=k8bid)
        else if (option .eq. 'ERZ1_ELEM') then
            call rsadpa(resuc1, 'E', 1, 'ERREUR_ERZ1', iord, &
                        0, sjv=ladpa, styp=k8bid)
        else if (option .eq. 'ERZ2_ELEM') then
            call rsadpa(resuc1, 'E', 1, 'ERREUR_ERZ2', iord, &
                        0, sjv=ladpa, styp=k8bid)
        else if (option .eq. 'QIRE_ELEM') then
            call rsadpa(resuc1, 'E', 1, 'ERREUR_QIRE', iord, &
                        0, sjv=ladpa, styp=k8bid)
        else if (option .eq. 'QIZ1_ELEM') then
            call rsadpa(resuc1, 'E', 1, 'ERREUR_QIZ1', iord, &
                        0, sjv=ladpa, styp=k8bid)
        else if (option .eq. 'QIZ2_ELEM') then
            call rsadpa(resuc1, 'E', 1, 'ERREUR_QIZ2', iord, &
                        0, sjv=ladpa, styp=k8bid)
        end if
        zr(ladpa) = err0
!
! ON CREE LA TABLE DES RESULTATS
!
        nompar(1) = 'NUME_ORDR'
        nompar(2) = 'OPTION'
        nompar(3) = 'ERRE_RELA'
        nompar(4) = 'ERRE_ABSO'
        nompar(5) = 'NORM_SIGM'
        typpar(1) = 'I'
        typpar(2) = 'K16'
        typpar(3) = 'R'
        typpar(4) = 'R'
        typpar(5) = 'R'
!
        call jeexin(resuc1//'.LTNT', iret)
        if (iret .eq. 0) call ltcrsd(resuc1, 'G')
!
        nomt19 = ' '
        call ltnotb(resuc1, 'ESTI_GLOB', nomt19)
        call jeexin(nomt19//'.TBBA', iret)
        if (iret .eq. 0) then
            call tbcrsd(nomt19, 'G')
            call tbajpa(nomt19, nbpar, nompar, typpar)
        end if
!
        listr(1) = nu0
        listr(2) = err0
        listr(3) = nors
!
        call tbajli(nomt19, nbpar, nompar, [iord], listr, &
                    [cbid], option, 0)
!
        write (ifi, *) ' '
        write (ifi, *) '***************************************************'
        if (option .eq. 'ERZ1_ELEM') then
            write (ifi, *) ' MECANIQUE : ESTIMATEUR D''ERREUR ZZ1'
        else if (option .eq. 'ERZ2_ELEM') then
            write (ifi, *) ' MECANIQUE : ESTIMATEUR D''ERREUR ZZ2'
        else if (option .eq. 'ERME_ELEM') then
            write (ifi, *) ' MECANIQUE : ESTIMATEUR D''ERREUR EN RESIDU ' &
                //'EXPLICITE'
        else if (option .eq. 'QIRE_ELEM') then
            write (ifi, *) ' MECANIQUE : ESTIMATEUR D''ERREUR EN QUANTITE ' &
                //'D''INTERET'
            write (ifi, *) '             - METHODE RESIDU EXPLICITE -'
        else if (option .eq. 'QIZ1_ELEM') then
            write (ifi, *) ' MECANIQUE : ESTIMATEUR D''ERREUR EN QUANTITE ' &
                //'D''INTERET'
            write (ifi, *) '             - METHODE ZZ1 -'
        else if (option .eq. 'QIZ2_ELEM') then
            write (ifi, *) ' MECANIQUE : ESTIMATEUR D''ERREUR EN QUANTITE ' &
                //'D''INTERET'
            write (ifi, *) '             - METHODE ZZ2 -'
        end if
        write (ifi, *) '***************************************************'
!
        write (ifi, 111) 'SD RESULTAT: ', resuco
        write (ifi, 110) ' NUMERO D''ORDRE: ', iord
!
        if (option .eq. 'QIRE_ELEM') then
            write (ifi, *) 'ERREUR            ABSOLUE'
            write (ifi, 108) ' TOTALE          ', err0
            write (ifi, 108) ' TERME VOLUMIQUE ', termvo
            write (ifi, 108) ' TERME SAUT      ', termsa
            write (ifi, 108) ' TERME NORMAL    ', termno
        else if ((option .eq. 'ERZ1_ELEM') .or. (option .eq. 'ERZ2_ELEM')) &
            then
            write (ifi, *) 'ERREUR            ABSOLUE'//&
         &            '         RELATIVE   NORME DE SIGMA'
            write (ifi, 108) ' TOTALE          ', err0, nu0, nors
        else
            write (ifi, *) 'ERREUR            ABSOLUE'//&
         &            '         RELATIVE   NORME DE SIGMA'
            write (ifi, 108) ' TOTALE          ', err0, nu0, nors
            write (ifi, 108) ' TERME VOLUMIQUE ', termvo, nuvo
            write (ifi, 108) ' TERME SAUT      ', termsa, nusa
            write (ifi, 108) ' TERME NORMAL    ', termno, nuno
        end if
!
108     format(a17, 1p, d16.8, 1x, 0p, f7.3, ' %', 1x, 1p, d16.8)
110     format(1p, a17, i5)
111     format(a17, a8)
!
        call jedema()
        end subroutine
