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
subroutine mmdelm(mdnode, xadj, adjncy, dhead, dforw, &
                  dbakw, qsize, llist, marker, maxint, &
                  tag, parent)
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
!
!--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDELM
!  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
!   C.ROSE AJ    DE PARENT + MODIFICA STRUCTURE : IF THEN ELSE ET DO WHI
!***************************************************************
!***************************************************************
!**     MMDELM ..... MULTIPLE MINIMUM DEGREE ELIMINATION     ***
!***************************************************************
!***************************************************************
!
!     PURPOSE - THIS ROUTINE ELIMINATES THE NODE MDNODE OF
!        MINIMUM DEGREE FROM THE ADJACENCY STRUCTURE, WHICH
!        IS STORED IN THE QUOTIENT GRAPH FORMAT.  IT ALSO
!        TRANSFORMS THE QUOTIENT GRAPH REPRESENTATION OF THE
!        ELIMINATION GRAPH.
!
!     INPUT PARAMETERS -
!        MDNODE - NODE OF MINIMUM DEGREE.
!        MAXINT - ESTIMATE OF MAXIMUM REPRESENTABLE (SHORT)
!                 INTEGER.
!        TAG    - TAG VALUE.
!
!     UPDATED PARAMETERS -
!        (XADJ,ADJNCY) - UPDATED ADJACENCY STRUCTURE.
!        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
!        QSIZE  - SIZE OF SUPERNODE.
!        MARKER - MARKER VECTOR.
!        LLIST  - TEMPORARY LINKED LIST OF ELIMINATED NABORS.
!
!***************************************************************
!
    integer(kind=8) :: adjncy(*), dbakw(*), dforw(*)
    integer(kind=8) :: llist(*), marker(*), qsize(*), dhead(*)
    integer(kind=8) :: xadj(*), parent(*)
    integer(kind=8) :: elmnt, i, istop, istrt, j, jstop, jstrt, link, maxint, mdnode
    integer(kind=8) :: nabor, node, npv, nqnbrs, nxnode, pvnode, rlmt, rloc, rnode, tag
    integer(kind=8) :: xqnbr
!
!***************************************************************
!
!        -----------------------------------------------
!        FIND REACHABLE SET AND PLACE IN DATA STRUCTURE.
!        -----------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    marker(mdnode) = tag
    istrt = xadj(mdnode)
    istop = xadj(mdnode+1)-1
!        -------------------------------------------------------
!        ELMNT POINTS TO THE BEGINNING OF THE LIST OF ELIMINATED
!        NABORS OF MDNODE, AND RLOC GIVES THE STORAGE LOCATION
!        FOR THE NEXT REACHABLE NODE.
!        -------------------------------------------------------
    elmnt = 0
    rloc = istrt
    rlmt = istop
    do i = istrt, istop
        nabor = adjncy(i)
        if (nabor .eq. 0) goto 120
        if (marker(nabor) .lt. tag) then
            marker(nabor) = tag
            if (dforw(nabor) .ge. 0) then
                adjncy(rloc) = nabor
                rloc = rloc+1
            else
                llist(nabor) = elmnt
                elmnt = nabor
                parent(nabor) = mdnode
            end if
        end if
    end do
120 continue
!            -----------------------------------------------------
!            MERGE WITH REACHABLE NODES FROM GENERALIZED ELEMENTS.
!            -----------------------------------------------------
!      DO WHILE (ELMNT.GT.0)
130 continue
    if (elmnt .gt. 0) then
        adjncy(rlmt) = -elmnt
        link = elmnt
140     continue
        jstrt = xadj(link)
        jstop = xadj(link+1)-1
        do j = jstrt, jstop
            node = adjncy(j)
            link = -node
            if (node .lt. 0) goto 140
            if (node .eq. 0) goto 170
            if (marker(node) .lt. tag .and. dforw(node) .ge. 0) then
                marker(node) = tag
!                            ---------------------------------
!                            USE STORAGE FROM ELIMINATED NODES
!                            IF NECESSARY.
!                            ---------------------------------
!            DO WHILE (RLOC.GE.RLMT)
150             continue
                if (rloc .ge. rlmt) then
                    link = -adjncy(rlmt)
                    rloc = xadj(link)
                    rlmt = xadj(link+1)-1
                    goto 150
! FIN DO WHILE
                end if
                adjncy(rloc) = node
                rloc = rloc+1
            end if
        end do
170     continue
        elmnt = llist(elmnt)
        goto 130
! FIN DO WHILE
    end if
    if (rloc .le. rlmt) adjncy(rloc) = 0
!        --------------------------------------------------------
!        FOR EACH NODE IN THE REACHABLE SET, DO THE FOLLOWING ...
!        --------------------------------------------------------
    link = mdnode
180 continue
    istrt = xadj(link)
    istop = xadj(link+1)-1
    do i = istrt, istop
        rnode = adjncy(i)
        link = -rnode
        if (rnode .eq. 0) goto 220
        if (rnode .lt. 0) goto 180
!                --------------------------------------------
!                IF RNODE IS IN THE DEGREE LIST STRUCTURE ...
!                --------------------------------------------
        pvnode = dbakw(rnode)
        if (pvnode .ne. 0 .and. pvnode .ne. (-maxint)) then
!                    -------------------------------------
!                    THEN REMOVE RNODE FROM THE STRUCTURE.
!                    -------------------------------------
            nxnode = dforw(rnode)
            if (nxnode .gt. 0) dbakw(nxnode) = pvnode
            if (pvnode .gt. 0) dforw(pvnode) = nxnode
            npv = -pvnode
            if (pvnode .lt. 0) dhead(npv) = nxnode
        end if
!                ----------------------------------------
!                PURGE INACTIVE QUOTIENT NABORS OF RNODE.
!                ----------------------------------------
        jstrt = xadj(rnode)
        jstop = xadj(rnode+1)-1
        xqnbr = jstrt
        do j = jstrt, jstop
            nabor = adjncy(j)
            if (nabor .eq. 0) goto 200
            if (marker(nabor) .lt. tag) then
                adjncy(xqnbr) = nabor
                xqnbr = xqnbr+1
            end if
        end do
200     continue
!                ----------------------------------------
!                IF NO ACTIVE NABOR AFTER THE PURGING ...
!                ----------------------------------------
        nqnbrs = xqnbr-jstrt
        if (nqnbrs .le. 0) then
!                    -----------------------------
!                    THEN MERGE RNODE WITH MDNODE.
!                    -----------------------------
            qsize(mdnode) = qsize(mdnode)+qsize(rnode)
            qsize(rnode) = 0
            marker(rnode) = maxint
            dforw(rnode) = -mdnode
            dbakw(rnode) = -maxint
        else
!                --------------------------------------
!                ELSE FLAG RNODE FOR DEGREE UPDATE, AND
!                ADD MDNODE AS A NABOR OF RNODE.
!                --------------------------------------
            dforw(rnode) = nqnbrs+1
            dbakw(rnode) = 0
            adjncy(xqnbr) = mdnode
            xqnbr = xqnbr+1
            if (xqnbr .le. jstop) adjncy(xqnbr) = 0
        end if
!
    end do
220 continue
end subroutine
