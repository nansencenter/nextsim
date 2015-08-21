/* \file issmmpi.h
 * \brief: header file that defines all the mpi wrappers that ISSM requires. The goal is to control
 * which MPI layer we are using at compile time: the standard mpi, the autodiff mpi or no mpi at all.
 */

#ifndef _ISSM_MPI_H_
#define _ISSM_MPI_H_

#include <cstddef>
#include <cassert>
//#include "../../shared/Numerics/types.h"

#ifdef _HAVE_MPI_
	/*Include header files: {{{*/
	#ifdef _HAVE_AMPI_
		#include <ampi/ampi.h>
	#else
		#include <mpi.h>
	#endif
	/*}}}*/
	/*MPI defines: *{{{*/
	// types
	typedef MPI_Comm             ISSM_MPI_Comm;
	typedef MPI_Datatype         ISSM_MPI_Datatype;
	typedef MPI_Op               ISSM_MPI_Op;
	typedef MPI_Status           ISSM_MPI_Status;

	// data types
	#define ISSM_MPI_CHAR       MPI_CHAR
	#ifdef _HAVE_AMPI_
		#define ISSM_MPI_DOUBLE    AMPI_ADOUBLE // corresponds to IssmDouble
	#else
		#define ISSM_MPI_DOUBLE    MPI_DOUBLE  // corresponds to IssmDouble
	#endif
	#define ISSM_MPI_PDOUBLE    MPI_DOUBLE  // corresponds to IssmPDouble
	#define ISSM_MPI_INT        MPI_INT

	// operations
	#define ISSM_MPI_MAX        MPI_MAX
	#define ISSM_MPI_MIN        MPI_MIN
	#define ISSM_MPI_PROD       MPI_PROD
	#define ISSM_MPI_SUM        MPI_SUM

	// others
	#define ISSM_MPI_COMM_WORLD    MPI_COMM_WORLD
	#define ISSM_MPI_STATUS_IGNORE MPI_STATUS_IGNORE
	#define ISSM_MPI_ANY_TAG       MPI_ANY_TAG
	#define ISSM_MPI_ANY_SOURCE    MPI_ANY_SOURCE

    /*other include files: */
	#include "./commops/commops.h"
	/*}}}*/
#else
	/*Our ISSM MPI defines: {{{*/
	// types
	typedef int                  ISSM_MPI_Comm;
	typedef int                  ISSM_MPI_Datatype;
	typedef int                  ISSM_MPI_Op;
	typedef int                  ISSM_MPI_Status;

	// data types
	#define ISSM_MPI_CHAR       1
	#define ISSM_MPI_DOUBLE     2
	#define ISSM_MPI_PDOUBLE    3
	#define ISSM_MPI_INT        4

	// operations
	#define ISSM_MPI_MAX        1
	#define ISSM_MPI_MIN        2
	#define ISSM_MPI_PROD       3
	#define ISSM_MPI_SUM        4

	// others
	#define ISSM_MPI_COMM_WORLD    1
	extern ISSM_MPI_Status ourIssmMPIStatusIgnore;
	#define ISSM_MPI_STATUS_IGNORE &ourIssmMPIStatusIgnore
	#define ISSM_MPI_ANY_TAG       2
	#define ISSM_MPI_ANY_SOURCE    3
	/*}}}*/
# endif

template <class T> ISSM_MPI_Datatype TypeToMPIType(){assert(false);};
template <> inline ISSM_MPI_Datatype TypeToMPIType<double>(){return ISSM_MPI_DOUBLE;};
template <> inline ISSM_MPI_Datatype TypeToMPIType<int>(){return ISSM_MPI_INT;};
template <> inline ISSM_MPI_Datatype TypeToMPIType<char>(){return ISSM_MPI_CHAR;};
template <class T> int ISSM_MPI_Bcast(T *buffer, int count,int root, ISSM_MPI_Comm comm){  /*{{{*/

	int rc=0;

	/*Get MPI type*/
	ISSM_MPI_Datatype datatype = TypeToMPIType<T>();

#ifdef _HAVE_MPI_
# ifdef _HAVE_AMPI_
	rc=AMPI_Bcast(buffer,
				count,
				datatype,
				root,
				comm);
# else
	rc=MPI_Bcast(buffer,
				count,
				datatype,
				root,
				comm);
# endif
#else
	// nothing to be done here
#endif
	return rc;
}/*}}}*/

/* interfaces  {{{*/
int ISSM_MPI_Allgather(void *sendbuf, int sendcount, ISSM_MPI_Datatype sendtype, void *recvbuf, int recvcount, ISSM_MPI_Datatype recvtype, ISSM_MPI_Comm comm);
int ISSM_MPI_Allgatherv(void *sendbuf, int sendcount, ISSM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts, int *displs, ISSM_MPI_Datatype recvtype, ISSM_MPI_Comm comm);
int ISSM_MPI_Allreduce(void *sendbuf, void *recvbuf, int count, ISSM_MPI_Datatype datatype, ISSM_MPI_Op op, ISSM_MPI_Comm comm);
int ISSM_MPI_Barrier(ISSM_MPI_Comm comm);
int ISSM_MPI_Bcast(void *buffer, int count, ISSM_MPI_Datatype datatype, int root, ISSM_MPI_Comm comm);
int ISSM_MPI_Comm_free(ISSM_MPI_Comm *comm);
int ISSM_MPI_Comm_rank(ISSM_MPI_Comm comm, int *rank);
int ISSM_MPI_Comm_size( ISSM_MPI_Comm comm, int *size);
int ISSM_MPI_Finalize(void);
int ISSM_MPI_Gather(void *sendbuf, int sendcnt, ISSM_MPI_Datatype sendtype, void *recvbuf, int recvcnt, ISSM_MPI_Datatype recvtype, int root, ISSM_MPI_Comm comm);
int ISSM_MPI_Gatherv(void *sendbuf, int sendcnt, ISSM_MPI_Datatype sendtype, void *recvbuf, int *recvcnts, int *displs, ISSM_MPI_Datatype recvtype, int root, ISSM_MPI_Comm comm);
int ISSM_MPI_Init(int *argc, char ***argv);
int ISSM_MPI_Recv(void *buf, int count, ISSM_MPI_Datatype datatype, int source, int tag, ISSM_MPI_Comm comm, ISSM_MPI_Status *status);
int ISSM_MPI_Reduce(void *sendbuf, void *recvbuf, int count, ISSM_MPI_Datatype datatype, ISSM_MPI_Op op, int root, ISSM_MPI_Comm comm);
int ISSM_MPI_Scatter(void *sendbuf, int sendcnt, ISSM_MPI_Datatype sendtype, void *recvbuf, int recvcnt, ISSM_MPI_Datatype recvtype, int root, ISSM_MPI_Comm comm);
int ISSM_MPI_Scatterv(void *sendbuf, int *sendcnts, int *displs, ISSM_MPI_Datatype sendtype, void *recvbuf, int recvcnt, ISSM_MPI_Datatype recvtype, int root, ISSM_MPI_Comm comm);
int ISSM_MPI_Send(void *buf, int count, ISSM_MPI_Datatype datatype, int dest, int tag, ISSM_MPI_Comm comm);
double ISSM_MPI_Wtime(void);

// special for Adol-C locations when buffers are allocated with new
// this could end up in the xNew template specialized for adoubles
// so as to not litter the code with it.
void ISSM_MPI_ContiguousInAdolc(size_t aSize);
/*}}}*/

#endif  //#ifndef _ISSM_MPI_H_
