#ifndef MPICONTAINERUTILITIES_H_
#define MPICONTAINERUTILITIES_H_

/**
 * @file MPIContainerUtilities.h
 * @brief Fichier entete contenant des utilitaires de manipulation de containers STL en parallèle
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

/* person_in_charge: nicolas.sellenet at edf.fr */
#include <set>
#include <numeric>
#include <algorithm>

#ifdef ASTER_HAVE_MPI

#include "astercxx.h"

#include "aster_mpi.h"
#include "ParallelUtilities/MPIInfos.h"
#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/JeveuxVector.h"

class MPIContainerUtilities {
  private:
    int _nbProcs;
    int _rank;
    aster_comm_t *_commWorld;

    template<typename T>
      struct dependent_false : std::false_type {};
    template<typename T> static MPI_Datatype mpi_type() {
      static_assert(dependent_false<T>::value, "Unknown MPI type");
      throw std::runtime_error("Unknown MPI type");
      return MPI_CHAR;
    }

  public:
    MPIContainerUtilities();

    template < int length >
    std::vector< JeveuxString< length > >
    gatheringVectorsOnAllProcs( std::vector< JeveuxString< length > > &toGather ) const {
        typedef JeveuxString< length > JeveuxChar;
        VectorInt sizes( _nbProcs, 0 );
        VectorInt sizes2( _nbProcs, 0 );
        sizes[0] = toGather.size();
        aster_mpi_allgather( sizes.data(), 1, MPI_INTEGER, sizes2.data(), 1, MPI_INTEGER,
                             _commWorld );
        int sum = 0;
        for ( auto taille : sizes2 )
            sum += taille;

        std::vector< JeveuxChar > toReturn;
        for ( int rank = 0; rank < _nbProcs; ++rank ) {
            JeveuxChar *retour = new JeveuxChar[sizes2[rank]];

            if ( rank == _rank )
                for ( int position = 0; position < sizes2[rank]; ++position )
                    retour[position] = toGather[position];

            aster_mpi_bcast( retour, length * sizes2[rank], MPI_CHAR, rank, _commWorld );
            for ( int position = 0; position < sizes2[rank]; ++position )
                toReturn.push_back( JeveuxChar( retour[position] ) );
            delete[] retour;
        }

        return toReturn;
    };

    /// Gather values from all processes. Same data count from each
    /// process (wrapper for MPI_Allgather)
    template<typename T>
    void all_gather(const std::vector<T>& in_values,
                             std::vector<T>& out_values);

    /// Gather values from each process (variable count per process)
    template<typename T>
    void all_gatherv(const std::vector<T>& in_values,
                             std::vector<std::vector<T>>& out_values);

    /// Gather values from each process (variable count per process)
    template<typename T>
    void all_gatherv(const std::vector<T>& in_values,
                              std::vector<T>& out_values);

    /// Gather values from each process (variable count per process)
    template<typename T>
    void all_gatherv(const std::vector<T>& in_values,
                              JeveuxVector<T>& out_values);

    /// Gather values, one primitive from each process (MPI_Allgather)
    template<typename T>
    void all_gather(const T in_value, std::vector<T>& out_values);

    /// Gather values, one primitive from each process (MPI_Allgather).
    /// Specialization for std::string
    void all_gatherv(const std::string& in_values,
                           VectorString& out_values);
    /// Gather values, one primitive from each process (MPI_Allgather).
    /// Specialization for VectorString
    void all_gather(const VectorString& in_values,
                           VectorString& out_values);

    /// AllReduce values, one value from each process (MPI_AllReduce).
    template<typename T>
    void all_reduce(const T in_value, T& out_value, MPI_Op op);
};

//---------------------------------------------------------------------------
template<> inline MPI_Datatype MPIContainerUtilities::mpi_type<float>() { return MPI_FLOAT; }
template<> inline MPI_Datatype MPIContainerUtilities::mpi_type<double>() { return MPI_DOUBLE; }
template<> inline MPI_Datatype MPIContainerUtilities::mpi_type<short int>() { return MPI_SHORT; }
template<> inline MPI_Datatype MPIContainerUtilities::mpi_type<int>() { return MPI_INT; }
template<> inline MPI_Datatype MPIContainerUtilities::mpi_type<long int>() { return MPI_LONG; }
template<> inline MPI_Datatype MPIContainerUtilities::mpi_type<unsigned int>()
                                                                    { return MPI_UNSIGNED; }
template<> inline MPI_Datatype MPIContainerUtilities::mpi_type<unsigned long int>()
                                                                    { return MPI_UNSIGNED_LONG; }
template<> inline MPI_Datatype MPIContainerUtilities::mpi_type<long long>()
                                                                    { return MPI_LONG_LONG; }
template<> inline MPI_Datatype MPIContainerUtilities::mpi_type<bool>()
                                                                    { return MPI_LOGICAL; }
//---------------------------------------------------------------------------
inline void MPIContainerUtilities::all_gatherv(const std::string& in_values,
                        VectorString& out_values) {
    const std::size_t comm_size = _nbProcs;

    // Get data size on each process
    VectorInt pcounts(comm_size);
    int local_size = in_values.size();
    MPI_Allgather(&local_size, 1, MPI_INT, pcounts.data(), 1, MPI_INT, _commWorld->id);

    // Build offsets
    VectorInt offsets(comm_size + 1, 0);
    for (std::size_t i = 1; i <= comm_size; ++i)
    offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather
    const std::size_t n = std::accumulate(pcounts.begin(), pcounts.end(), 0);
    std::vector<char> _out(n);
    MPI_Allgatherv(const_cast<char*>(in_values.data()), in_values.size(),
                MPI_CHAR, _out.data(), pcounts.data(), offsets.data(),
                MPI_CHAR, _commWorld->id);

    // Rebuild
    out_values.resize(comm_size);
    for (std::size_t p = 0; p < comm_size; ++p)
    {
    out_values[p] = std::string(_out.begin() + offsets[p],
                                _out.begin() + offsets[p + 1]);
    }
}
//---------------------------------------------------------------------------
inline void MPIContainerUtilities::all_gather(const VectorString& in_values,
                        VectorString& out_values) {
    const std::size_t comm_size = _nbProcs;

    // Gather
    std::set<std::string> stringSet;
    for (auto str: in_values) {
        VectorString tmp(comm_size);
        all_gatherv(str, tmp);
        std::copy(tmp.begin(), tmp.end(), std::inserter(stringSet, stringSet.end()));
    }

    // Rebuild
    out_values.resize(stringSet.size());
    std::copy(stringSet.begin(), stringSet.end(), out_values.begin());
    std::sort(out_values.begin(), out_values.end());
}
//-------------------------------------------------------------------------
template<typename T>
    void MPIContainerUtilities::all_gather(const std::vector<T>& in_values,
                    std::vector<T>& out_values) {
    out_values.resize(in_values.size()*_nbProcs);
    MPI_Allgather(const_cast<T*>(in_values.data()), in_values.size(),
                mpi_type<T>(),
                out_values.data(), in_values.size(), mpi_type<T>(),
                _commWorld->id);
}
//---------------------------------------------------------------------------
template<typename T>
    void MPIContainerUtilities::all_gatherv(const std::vector<T>& in_values,
                    std::vector<std::vector<T>>& out_values) {

    const std::size_t comm_size = _nbProcs;

    // Get data size on each process
    VectorInt pcounts;
    const int local_size = in_values.size();
    MPIContainerUtilities::all_gather(local_size, pcounts);
    assert(pcounts.size() == comm_size);

    // Build offsets
    VectorInt offsets(comm_size + 1, 0);
    for (std::size_t i = 1; i <= comm_size; ++i)
    offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather data
    const std::size_t n = std::accumulate(pcounts.begin(), pcounts.end(), 0);
    std::vector<T> recvbuf(n);
    MPI_Allgatherv(const_cast<T*>(in_values.data()), in_values.size(),
                mpi_type<T>(),
                recvbuf.data(), pcounts.data(), offsets.data(),
                mpi_type<T>(), _commWorld->id);

    // Repack data
    out_values.resize(comm_size);
    for (std::size_t p = 0; p < comm_size; ++p)
    {
    out_values[p].resize(pcounts[p]);
    for (int i = 0; i < pcounts[p]; ++i)
        out_values[p][i] = recvbuf[offsets[p] + i];
    }
}
//---------------------------------------------------------------------------
template<typename T>
    void MPIContainerUtilities::all_gatherv(const std::vector<T>& in_values,
                    std::vector<T>& out_values) {

    const std::size_t comm_size = _nbProcs;

    // Get data size on each process
    VectorInt pcounts;
    const int local_size = in_values.size();
    MPIContainerUtilities::all_gather(local_size, pcounts);
    assert(pcounts.size() == comm_size);

    // Build offsets
    VectorInt offsets(comm_size + 1, 0);
    for (std::size_t i = 1; i <= comm_size; ++i)
    offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather data
    const std::size_t n = std::accumulate(pcounts.begin(), pcounts.end(), 0);
    out_values.resize(n);
    MPI_Allgatherv(const_cast<T*>(in_values.data()), in_values.size(),
                mpi_type<T>(),
                out_values.data(), pcounts.data(), offsets.data(),
                mpi_type<T>(), _commWorld->id);
}
//---------------------------------------------------------------------------
template<typename T>
    void MPIContainerUtilities::all_gatherv(const std::vector<T>& in_values,
                    JeveuxVector<T>& out_values) {

    const std::size_t comm_size = _nbProcs;

    // Get data size on each process
    VectorInt pcounts;
    const int local_size = in_values.size();
    MPIContainerUtilities::all_gather(local_size, pcounts);
    assert(pcounts.size() == comm_size);

    // Build offsets
    VectorInt offsets(comm_size + 1, 0);
    for (std::size_t i = 1; i <= comm_size; ++i)
    offsets[i] = offsets[i - 1] + pcounts[i - 1];

    // Gather data
    const std::size_t n = std::accumulate(pcounts.begin(), pcounts.end(), 0);

    if( out_values->size() < n)
    {
        if( out_values->isAllocated() )
            out_values->deallocate();
        out_values->allocate(n);
    }
    out_values->updateValuePointer();

    MPI_Allgatherv(const_cast<T*>(in_values.data()), in_values.size(),
                mpi_type<T>(),
                out_values->getDataPtr(), pcounts.data(), offsets.data(),
                mpi_type<T>(), _commWorld->id);

}
//---------------------------------------------------------------------------
template<typename T>
    void MPIContainerUtilities::all_gather(const T in_value,
                                           std::vector<T>& out_values) {
    out_values.resize(_nbProcs);
    MPI_Allgather(const_cast<T*>(&in_value), 1, mpi_type<T>(),
                out_values.data(), 1, mpi_type<T>(), _commWorld->id);
}
//---------------------------------------------------------------------------
template<typename T>
    void MPIContainerUtilities::all_reduce(const T in_value,
                                          T& out_value, MPI_Op op) {
    MPI_Allreduce(const_cast<T*>(&in_value), const_cast<T*>(&out_value),
                1, mpi_type<T>(), op, _commWorld->id);
}
//---------------------------------------------------------------------------


#endif /* MPICONTAINERUTILITIES_H_ */

#endif /* ASTER_HAVE_MPI */
