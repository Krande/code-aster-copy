! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine ornorm(noma, listma, nbmail, reorie, norien,&
                  command)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/indiis.h"
#include "asterfort/infniv.h"
#include "asterfort/iorim1.h"
#include "asterfort/iorim2.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmavo.h"
#include "asterfort/utmess.h"
!
    integer :: listma(*), nbmail, norien
    aster_logical :: reorie
    character(len=8) :: noma
    character(len=24), optional :: command
!.======================================================================
!
!   ORNORM  --  LE BUT EST QUE TOUTES LES MAILLES DE LA LISTE SOIENT
!               ORIENTEES COMME LA PREMIERE MAILLE DE LA LISTE.
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NOMA           IN    K8      NOM DU MAILLAGE
!    LISTMA         IN    I       LISTE DES MAILLES A REORIENTER
!    NBMAIL         IN    I       NB DE MAILLES DE LA LISTE
!    NORIEN        VAR            NOMBRE DE MAILLES REORIENTEES
!.========================= DEBUT DES DECLARATIONS ====================
! -----  VARIABLES LOCALES
    integer :: nutyma, iliste
    integer :: ima, numail, numa, norieg, lliste
    integer :: im1, im2, ico, ibid(1)
    integer :: p1, p2, ifm, niv, p3, p4
    integer :: jdesm1, jdesm2
    integer :: nbmavo, indi, im3, nconex, zero
    aster_logical :: dime1, dime2
    character(len=1) :: lect
    character(len=2) :: kdim
    character(len=8) :: typel, nomail
    character(len=24) :: mailma, nomavo
    character(len=24) :: valk(2), cmd
    integer, pointer :: ori1(:) => null()
    integer, pointer :: ori2(:) => null()
    integer, pointer :: ori3(:) => null()
    integer, pointer :: ori4(:) => null()
    character(len=8), pointer :: ori5(:) => null()
    integer, pointer :: typmail(:) => null()
!
#define pasori(ima) ori1(ima).eq.0
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    cmd = ' '
    if (present(command)) cmd = command
!
    call jemarq()
    if (nbmail .eq. 0) goto 999
!
    call infniv(ifm, niv)
!
    zero = 0
    mailma = noma//'.NOMMAI'
    lect = 'L'
    if (reorie) lect = 'E'
!
! --- VECTEUR DU TYPE DES MAILLES DU MAILLAGE :
!     ---------------------------------------
    call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
!
! --- APPEL A LA CONNECTIVITE :
!     -----------------------
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', p2)
    call jeveuo(noma//'.CONNEX', lect, p1)
!
!     ALLOCATIONS :
!     -----------
    AS_ALLOCATE(vi=ori1, size=nbmail)
    AS_ALLOCATE(vi=ori2, size=nbmail)
    AS_ALLOCATE(vi=ori3, size=nbmail)
    AS_ALLOCATE(vi=ori4, size=nbmail)
    AS_ALLOCATE(vk8=ori5, size=nbmail)
!
! --- VERIFICATION DU TYPE DES MAILLES
! --- (ON DOIT AVOIR DES MAILLES DE PEAU) :
!     -----------------------------------
    dime1 = .false.
    dime2 = .false.
    do ima = 1, nbmail
        ori1(ima) = 0
        numa = listma(ima)
        ori3(ima) = zi(p2+numa)-zi(p2-1+numa)
        ori4(ima) = zi(p2+numa-1)
        jdesm1 = zi(p2+numa-1)
!
! ---   TYPE DE LA MAILLE COURANTE :
!       --------------------------
        nutyma = typmail(numa)
        call jenuno(jexnum('&CATA.TM.NOMTM', nutyma), typel)
        ori5(ima) = typel
!
        if (typel(1:4) .eq. 'QUAD') then
            dime2 = .true.
        else if (typel(1:4).eq.'TRIA') then
            dime2 = .true.
        else if (typel(1:3).eq.'SEG') then
            dime1 = .true.
        else
            call jenuno(jexnum(mailma, numa), nomail)
            valk(1) = nomail
            valk(2) = typel
            call utmess('F', 'MODELISA5_94', nk=2, valk=valk)
        endif
        if (dime1 .and. dime2) then
            call utmess('F', 'MODELISA5_98')
        endif
    end do
!
    if (dime2 .and. cmd(1:10) .eq. 'ORIE_LIGNE') then
        call utmess('F', 'MODELISA5_92')
    endif
!
!
! --- RECUPERATION DES MAILLES VOISINES DU GROUP_MA :
!     ---------------------------------------------
    kdim ='  '
    if (dime1) kdim ='1D'
    if (dime2) kdim ='2D'
    nomavo = '&&ORNORM.MAILLE_VOISINE '
    call utmavo(noma, kdim, listma, nbmail, 'V',&
                nomavo, zero, ibid)
    call jeveuo(jexatr(nomavo, 'LONCUM'), 'L', p4)
    call jeveuo(nomavo, 'L', p3)
!
    norieg = 0
!
! --- LA BOUCLE 100 DEFINIT LES CONNEXES
!
    nconex = 0
    do ima = 1, nbmail
        numail = listma(ima)
! ----- SI LA MAILLE N'EST PAS ORIENTEE ON L'ORIENTE
        if (pasori(ima)) then
            if (niv .eq. 2) then
                call jenuno(jexnum(mailma, numail), nomail)
                write (ifm,*) 'LA MAILLE ',nomail,&
     &                    ' SERT A ORIENTER UN NOUVEAU GROUPE CONNEXE'
            endif
            nconex = nconex + 1
            if (nconex .gt. 1) then
                if (cmd .ne. ' ') then
                    call utmess('F', 'MODELISA6_2')
                else
                    call utmess('F', 'MODELISA5_99')
                endif
            endif
            ori1(ima) = 1
            lliste = 0
            iliste = 0
            ori2(lliste+1) = ima
!
! ------- ON ORIENTE TOUTES LES MAILLES DU CONNEXE
!
200         continue
!
            im1 = ori2(iliste+1)
            jdesm1 = ori4(im1)
! ------- ON ESSAYE D'ORIENTER LES MAILLES VOISINES
            nbmavo = zi(p4+im1)-zi(p4-1+im1)
            do im3 = 1, nbmavo
                indi = zi(p3+zi(p4+im1-1)-1+im3-1)
                im2 = indiis ( listma, indi, 1, nbmail )
                if (im2 .eq. 0) goto 210
                numail = listma(im2)
                if (pasori(im2)) then
                    jdesm2 = ori4(im2)
!             VERIFICATION DE LA CONNEXITE ET REORIENTATION EVENTUELLE
                    if (dime1) ico = iorim1 ( zi(p1+jdesm1-1), zi(p1+jdesm2-1), reorie)
                    if (dime2) ico = iorim2 (&
                                     zi(p1+jdesm1-1), ori3(im1), zi(p1+jdesm2-1), ori3(im2),&
                                     reorie&
                                     )
!             SI MAILLES CONNEXES
                    if (ico .ne. 0) then
                        ori1(im2) = 1
                        lliste = lliste + 1
                        ori2(lliste+1) = im2
                        if (reorie .and. niv .eq. 2) then
                            call jenuno(jexnum(mailma, numail), nomail)
                            if (ico .lt. 0) then
                                write (ifm,*) 'LA MAILLE ',nomail,' A ETE REORIENTEE'
                            else
                                write (ifm,*) 'LA MAILLE ',nomail,' EST ORIENTEE'
                            endif
                        endif
                    endif
!
!             SI ORIENTATIONS CONTRAIRES
                    if (ico .lt. 0) norieg = norieg + 1
!
                endif
210             continue
            end do
            iliste = iliste + 1
            if (iliste .le. lliste) goto 200
        endif
    end do
!
    norien = norien + norieg
!
    AS_DEALLOCATE(vi=ori1)
    AS_DEALLOCATE(vi=ori2)
    AS_DEALLOCATE(vi=ori3)
    AS_DEALLOCATE(vi=ori4)
    AS_DEALLOCATE(vk8=ori5)
    call jedetr(nomavo)
!
999 continue
    call jedema()
end subroutine
