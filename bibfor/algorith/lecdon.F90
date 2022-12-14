! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
! *   LOGICIEL CODE_ASTER - COUPLAGE ASTER/EDYOS - Copyright EDF 2009  *
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
subroutine lecdon(ficext, unitpa, prdeff)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/getvtx.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ulisop.h"
#include "asterfort/ulopen.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer :: unitpa
    aster_logical :: prdeff, ficext
! person_in_charge: nicolas.greffet at edf.fr
! **********************************************************************
! **********************************************************************
!  LECDON : FONCTION
!  -----------------
! CE SSP PERMET DE LIRE LES DONNEES RELATIVES AUX PALIERS CONTENUES DANS
!  LE FICHIER FIC_DON, SOIT:
!          - LE NOMBRE DE PALIERS
!          - (POUR CHAQUE PALIER) LE NOEUD ASTER D'APPLICATION ET
!          LE TYPE DE PALIER
!=======================================================================
!  REFERENCES BIBLIOGRAPHIQUES
!  ---------------------------
! ======================================================================
!  DEVELOPPEMENTS ET CORRECTIONS D'ANOMALIES
!  -----------------------------------------
!  DATE: 13/02/09   AUTEUR: P. VAUGRANTE    ANOMALIE: DEVELOPPEMENT
!  DATE:            AUTEUR:                 ANOMALIE:
!  DATE:            AUTEUR:                 ANOMALIE:
!  DATE:            AUTEUR:                 ANOMALIE:
! ======================================================================
!  VARIABLES UTILISEES
!  -------------------
!
!  ___________________________________________________________________
! !    NOM   !   TYPE      !                  ROLE                    !
! !__________!_____________!__________________________________________!
! !          !             !                                          !
! ! CARAC"I" !  CHARACTER  !  VARIABLE TAMPON PERMETTANT DE FAIRE     !
! !          !             !  CORRESPONDRE A UN NOMBRE SON EQUIVALENT !
! !          !             !  EN CARACERES                            !
! !          !             !                                          !
! ! UNITPA   !  ENTIER     !  UNITE DE LECTURE                        !
! !          !             !                                          !
! ! CTYPE    !  CHARACTER*6!  TYPE DU PALIER LU                       !
! !          !             !                                          !
! ! IPAL     !  ENTIER     !  INDICE DE BOUCLE SUR LES PALIERS        !
! !          !             !                                          !
! ! NUMPAL   !  ENTIER     !  NUMERO DU PALIER LU                     !
! !          !             !                                          !
! ! NOMPRG   !  CHARACTER  !  NOM DU SSP (POUR ECRITURE DANS ERRCOU)  !
! !          !             !                                          !
! ! PALMAX   !  ENTIER     !  NOMBRE MAXIMUM DE PALIERS (PARAMETER)   !
! !          !             !                                          !
! ! DIMNAS   !  CHARACTER  !  NOMBRE DE DDL POUR UN NOEUD             !
! !          !             !                                          !
! !__________!_____________!__________________________________________!
!
!
!
!
!
!
!  COMMON ZI (TYPE: INTEGER) (NOM = 'NPAL')
!  _____________________________________________________________________
! !              !             !                                       !
! ! NBPAL        !  ADR        !  NOMBRE DE PALIERS POUR L'ETUDE       !
! !              !             !                                       !
! !______________!_____________!_______________________________________!
!
!
!
!  COMMON ZK8 (TYPE: CHARACTER*8) (NOM = 'C_PAL')
!  _____________________________________________________________________
! !              !             !                                       !
! ! TYPPAL(IPAL) ! ADR+(IPAL-1)!  TYPE DU PALIER CONSIDERE             !
! !              !             !                                       !
! ! FINPAL(IPAL) !  ADR+PALMAX !  TERMINAISON POUR LE PALIER CONSIDERE !
! !              !  +(IPAL-1)  !  PALIER N??I => _I                     !
! !              !             !                                       !
! ! CNPAL(IPAL)  ! ADR+2*PALMAX!  NOM DU NOEUD ASTER POUR LE PALIER    !
! !              !  +(IPAL-1)  !  CONSIDERE                            !
! !______________!_____________!_______________________________________!
!
!
!
!
!
!
!     VARIABLES INTERNES
!     ==================
    integer :: ifm, niv
    character(len=8) :: nomprg
    parameter(nomprg='LECDON')
!
    character(len=1) :: carac1
    character(len=2) :: carac2
!
    integer :: ipal, numpal, ngr, n2
    integer :: nbpal, iret
    character(len=6) :: ctype
!
    integer :: palmax
    parameter (palmax=20)
    character(len=6) :: typpal(palmax)
    character(len=3) :: finpal(palmax)
    character(len=8) :: cnpal(palmax), cnod
!
    integer :: zcpal, znpal
    character(len=16) :: k16nom
    character(len=24) :: cpal, npal
!
    call jemarq()
    niv = 0
    call infdbg('YACS_EDYOS', ifm, niv)
!
!     ASSIGNATION DES NOMS POUR LES ADRESSES DANS LES COMMON ASTER
!     ------------------------------------------------------------
    cpal='C_PAL'
    npal='N_PAL'
!
!     RESERVATION MEMOIRE POUR LES "COMMON"  ASTER
!     --------------------------------------------
    call jeexin(cpal, iret)
    if (iret .eq. 0) then
        call wkvect(cpal, 'G V K8', (3*palmax), zcpal)
    else
        call jeveuo(cpal, 'E', zcpal)
    endif
    call jeexin(npal, iret)
    if (iret .eq. 0) then
        call wkvect(npal, 'G V I', (1+palmax), znpal)
    else
        call jeveuo(npal, 'E', znpal)
    endif
!
    if (ficext) then
!
!     LECTURE DU FICHIER FIC_DON
!     --------------------------
        if (niv .ge. 2) write(ifm, *)'ASTEREDYOS: ', nomprg,&
                        ' DEBUT LECTURE DU FICHIER FIC_DON UNITE LOGIQUE ', unitpa
        k16nom ='                '
        if (ulisop ( unitpa, k16nom ) .eq. 0) then
            call ulopen(unitpa, ' ', ' ', 'NEW', 'O')
        endif
        read(unitpa,*)nbpal
        if (niv .ge. 2) write(ifm, *)'ASTEREDYOS: ', nomprg, ' ON A LU NBPAL =', nbpal
        if (nbpal .gt. palmax) then
            call utmess('F', 'EDYOS_43')
        endif
!
!     REMPLISSAGE "COMMON" ASTER POUR LE NOMBRE DE PALIERS
!     ----------------------------------------------------
        zi(znpal)=nbpal
!
!     BOUCLE DE LECTURE SUR LES PALIERS
!     ---------------------------------
        do ipal = 1, nbpal
            cnod='        '
            read(unitpa,*)numpal,cnod,ctype
            typpal(numpal)=ctype
            cnpal(ipal)=cnod
!
            if (ipal .lt. 10) then
                write(carac1,'(I1)')ipal
                finpal(ipal)='_'//carac1
            else
                write(carac2,'(I2)')ipal
                finpal(ipal)='_'//carac2
            endif
!
!
!   REMPLISSAGE "COMMON" ASTER POUR LES PALIERS (TYPE,TERMINAISON,NOEUD)
!   --------------------------------------------------------------------
            zk8(zcpal+(ipal-1))=ctype
            zk8(zcpal+palmax+(ipal-1))=finpal(ipal)
            zk8(zcpal+(2*palmax)+(ipal-1))=cnpal(ipal)
            if (niv .ge. 2) write(ifm, * )'ASTEREDYOS : LECDON : CTYPE - FINPAL - CNPAL=', ctype,&
                            ' -- ', finpal(ipal), ' -- ', cnpal(ipal)
!
!   REMPLISSAGE "COMMON" ASTER POUR LES NUMEROS DES NOEUDS DES PALIERS
!   ------------------------------------------------------------------
            zi(znpal+1+(ipal-1))=ipal
        end do
!
!     FIN DE BOUCLE DE LECTURE SUR LES PALIERS
!     ----------------------------------------
        prdeff = .true.
!
!     ECRITURE DES VARIABLES LUES
!     ---------------------------
        if (niv .ge. 2) then
            write(ifm,*)'ASTEREDYOS: ',nomprg,&
     &              ' - FIN LECTURE DU FICHIER FIC_DON '
            write(ifm,*)'ASTEREDYOS: ',nomprg,' - NOMBRE DE PALIERS: ',nbpal
            do ipal = 1, nbpal
                write(ifm,*)'ASTEREDYOS PALIER :',ipal,' TYPE :',&
                typpal(ipal), ' NOEUD ASTER : ',cnpal(ipal)
            end do
        endif
!
        if (ulisop ( unitpa, k16nom ) .ne. 0) call ulopen(-unitpa, ' ', ' ', 'NEW', 'O')
!
    else
!
!   LECTURE DES INFOS DEPUIS LE FICHIER DE COMMANDE
!   --------------------------------------------------------------------
        call getfac('PALIER_EDYOS', nbpal)
        if (nbpal .gt. palmax) then
            call utmess('F', 'EDYOS_43')
        endif
        zi(znpal)=nbpal
        do ipal = 1, nbpal
            call getvtx('PALIER_EDYOS', 'GROUP_NO', iocc=ipal, nbval=0, nbret=n2)
            if (abs(n2) .eq. 0) then
                call getvtx('PALIER_EDYOS', 'NOEUD', iocc=ipal, nbval=0, nbret=n2)
                if (abs(n2) .eq. 0) then
                    call utmess('F', 'EDYOS_49')
                else
                    ngr = -n2
                    ngr = 1
                    call getvtx('PALIER_EDYOS', 'NOEUD', iocc=ipal, nbval=ngr, vect=cnpal(ipal),&
                                nbret=n2)
                endif
            else
                ngr = -n2
                ngr = 1
                call getvtx('PALIER_EDYOS', 'GROUP_NO', iocc=ipal, nbval=ngr, vect=cnpal(ipal),&
                            nbret=n2)
            endif
            call getvtx('PALIER_EDYOS', 'TYPE_EDYOS', iocc=ipal, nbval=ngr, vect=ctype,&
                        nbret=n2)
            zk8(zcpal+(2*palmax)+(ipal-1))=cnpal(ipal)
            zk8(zcpal+(ipal-1))=ctype
        end do
!
        niv = 3
!
        do ipal = 1, nbpal
!
            if (ipal .lt. 10) then
                write(carac1,'(I1)')ipal
                finpal(ipal)='_'//carac1
            else
                write(carac2,'(I2)')ipal
                finpal(ipal)='_'//carac2
            endif
            zk8(zcpal+palmax+(ipal-1))=finpal(ipal)
            if (niv .ge. 2) write(ifm, * )'ASTEREDYOS : LECDON : CTYPE - FINPAL - CNPAL=',&
                            zk8(zcpal+(ipal-1)), ' -- ', finpal(ipal), ' -- ',&
                            zk8(zcpal+(2*palmax)+(ipal-1))
            zi(znpal+1+(ipal-1))=ipal
        end do
        if (niv .ge. 2) then
            write(ifm,*)'ASTEREDYOS: ',nomprg,&
     &              ' - FIN LECTURE ARGUMENTS PALIERS '
            write(ifm,*)'ASTEREDYOS: ',nomprg,' - NOMBRE PALIERS: ',&
            nbpal
        endif
    endif
!
!
    call jedema()
!
end subroutine
