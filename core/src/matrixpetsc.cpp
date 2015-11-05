/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   matrixpetsc.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Jul  6 16:39:32 2015
 */

#include <vectorpetsc.hpp>
#include <matrixpetsc.hpp>

/**
 * @class MatrixPetsc
 * @brief Wrapper for petsc matrices
 *
 * @see
 *
 */

namespace Nextsim
{

MatrixPetsc::MatrixPetsc( Communicator const& comm )
	:
    M_comm( comm ),
    M_is_initialized( false )
{
#if 0
    int ierr = 0;

    ierr = MatCreate(M_comm,&M_mat);
    CHKERRABORT( comm, ierr );

    ierr = MatSetType(M_mat,MATSEQAIJ);
    CHKERRABORT( comm, ierr );

    M_is_initialized = false;
#endif
}

MatrixPetsc::~MatrixPetsc()
{
    this->clear();
}

Mat
MatrixPetsc::mat() const
{
    ASSERT(M_is_initialized, "null petsc matrix");

    return M_mat;
}

Mat&
MatrixPetsc::mat()
{
    ASSERT(M_is_initialized, "null petsc matrix");

    return M_mat;
}

void
MatrixPetsc::init( const size_type m, const size_type n, const size_type nnz )
{
	if (m==0 || n==0)
		return;

    {
        if (M_is_initialized)
            this->clear();
    }

    int ierr = 0;
    int nrow = static_cast<int> (m);
    int ncol = static_cast<int> (n);
    int n_nz  = static_cast<int> (nnz);

    // ierr = PetscInitialize(0,(char ***)"", PETSC_NULL, PETSC_NULL);
    // CHKERRABORT( M_comm, ierr );

    ierr = MatCreateSeqAIJ ( M_comm, nrow, ncol, n_nz, PETSC_NULL, &M_mat );
    CHKERRABORT( M_comm, ierr );


    ierr = MatSetFromOptions ( M_mat );
    CHKERRABORT( M_comm, ierr );

    ierr = MatSetOption( M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE );
    MatSetOption( M_mat,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE );
    CHKERRABORT( M_comm, ierr );

    // generates an error for new matrix entry
    ierr = MatSetOption ( M_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE );
    CHKERRABORT( M_comm, ierr );

    M_is_initialized = true;

    this->zero();
}

void
MatrixPetsc::zero()
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int ierr = 0;
    PetscBool is_assembled;
    MatAssembled( M_mat, &is_assembled );

    //std::cout<<"Assemble= "<< is_assembled <<"\n";
    if ( is_assembled )
    {
        ierr = MatZeroEntries( M_mat );
        CHKERRABORT( M_comm,ierr );
    }
}

void
MatrixPetsc::close() const
{
    int ierr = 0;

    ierr = MatAssemblyBegin(M_mat,MAT_FINAL_ASSEMBLY);
    CHKERRABORT( M_comm, ierr );

    ierr = MatAssemblyEnd(M_mat,MAT_FINAL_ASSEMBLY);
    CHKERRABORT( M_comm, ierr );
}

void
MatrixPetsc::setValue(size_type const& i, size_type const& j, value_type const& value)
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int ierr=0, i_val=i, j_val=j;

    PetscScalar petsc_value = static_cast<PetscScalar>( value );
    ierr = MatSetValues( M_mat, 1, &i_val, 1, &j_val,
                         &petsc_value, INSERT_VALUES );
    CHKERRABORT( M_comm,ierr );
}

void
MatrixPetsc::addValue(size_type const& i, size_type const& j, value_type const& value)
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int ierr=0, i_val=i, j_val=j;

    PetscScalar petsc_value = static_cast<PetscScalar>( value );
    ierr = MatSetValues( M_mat, 1, &i_val, 1, &j_val,
                         &petsc_value, ADD_VALUES );
    CHKERRABORT( M_comm,ierr );
}

void
MatrixPetsc::setMatrix(int* rows, int nrows,
                       int* cols, int ncols,
                       value_type* data)
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int ierr=0;

    ierr = MatSetValues( M_mat,
                         nrows, /*( int* )*/ rows,
                         ncols, /*( int* )*/ cols,
                         /*( PetscScalar* )*/ data,
                         INSERT_VALUES );

    CHKERRABORT( M_comm,ierr );
}

void
MatrixPetsc::addMatrix(int* rows, int nrows,
                       int* cols, int ncols,
                       value_type* data)
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int ierr=0;

    ierr = MatSetValues( M_mat,
                         nrows, /*( int* )*/ rows,
                         ncols, /*( int* )*/ cols,
                         /*( PetscScalar* )*/ data,
                         ADD_VALUES );

    CHKERRABORT( M_comm,ierr );
}

typename MatrixPetsc::value_type
MatrixPetsc::l1Norm() const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int ierr=0;
    //value_type petsc_value;
    value_type value;

    this->close();

    ierr = MatNorm( M_mat, NORM_INFINITY, &value );
    CHKERRABORT( M_comm,ierr );

    //value = static_cast<real_type>( petsc_value );
    return value;
}

typename MatrixPetsc::value_type
MatrixPetsc::linftyNorm() const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int ierr=0;
    //double petsc_value;
    value_type value;

    this->close();

    ierr = MatNorm( M_mat, NORM_INFINITY, &value );
    CHKERRABORT( M_comm,ierr );

    //value = static_cast<real_type>( petsc_value );
    return value;
}

typename MatrixPetsc::value_type
MatrixPetsc::operator () (size_type const& i, size_type const& j) const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    const PetscScalar *petsc_row;
    const PetscInt    *petsc_cols;

    value_type value=0.;

    int
        ierr=0,
        ncols=0,
        i_val=static_cast<int>( i ),
        j_val=static_cast<int>( j );

    ierr = MatGetRow( M_mat, i_val, &ncols, &petsc_cols, &petsc_row );
    CHKERRABORT( M_comm,ierr );

    // Perform a binary search to find the contiguous index in
    // petsc_cols (resp. petsc_row) corresponding to global index j_val
    std::pair<const int*, const int*> p =
        std::equal_range ( &petsc_cols[0], &petsc_cols[0] + ncols, j_val );

    // Found an entry for j_val
    if ( p.first != p.second )
    {
        // The entry in the contiguous row corresponding
        // to the j_val column of interest
        const int j = std::distance ( const_cast<int*>( &petsc_cols[0] ),
                                      const_cast<int*>( p.first ) );

        assert ( j < ncols );
        assert ( petsc_cols[j] == j_val );

        value = static_cast<value_type> ( petsc_row[j] );

        ierr  = MatRestoreRow( M_mat, i_val,
                               &ncols, &petsc_cols, &petsc_row );
        CHKERRABORT( this->comm(),ierr );

        return value;
    }

    return 0;
}

MatrixPetsc&
MatrixPetsc::operator = (MatrixPetsc const& M)
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    M.close();
    int ierr= MatDuplicate( M.mat(),MAT_COPY_VALUES,&M_mat );
    CHKERRABORT( M_comm,ierr );

    return *this;
}

typename MatrixPetsc::size_type
MatrixPetsc::size1() const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int petsc_m=0, petsc_n=0, ierr=0;

    ierr = MatGetSize ( M_mat, &petsc_m, &petsc_n );
    CHKERRABORT( M_comm,ierr );
    return static_cast<size_type>( petsc_m );

}

typename MatrixPetsc::size_type
MatrixPetsc::size2() const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int petsc_m=0, petsc_n=0, ierr=0;

    ierr = MatGetSize ( M_mat, &petsc_m, &petsc_n );
    CHKERRABORT( this->comm(),ierr );
    return static_cast<size_type>( petsc_n );

}

void
MatrixPetsc::printScreen() const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int ierr = 0;
    ierr = MatView(M_mat,PETSC_VIEWER_STDOUT_WORLD);
    CHKERRABORT( M_comm, ierr );
}

void
MatrixPetsc::printMatlab(std::string const& filename) const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    // if (filename.substr(filename.size()-2, filename.size()) != ".m")
    //     filename += ".m";

    int ierr = 0;

    //ASSERT(ierr == 0, "MatrixPetsc not properly initialized");

    //std::cout<<"PATH= "<< boost::filesystem::path(filename).stem().string() <<"\n";
    //std::cout<<"CHAR= " << filename.substr(filename.size()-2, filename.size()) <<"\n";

    std::string matfilename = Environment::nextsimDir().string() + "/matlab/" + filename;

    fs::path path(matfilename);
    if ( !fs::exists(path) )
        fs::create_directories(path.parent_path());


    //PetscObjectSetName((PetscObject)M_mat,boost::filesystem::path("out_"+matfilename).stem().string().c_str());
    PetscObjectSetName((PetscObject)M_mat,("var_" + path.stem().string()).c_str());

    PetscViewer petsc_viewer;

    ierr = PetscViewerCreate ( M_comm, &petsc_viewer );
    CHKERRABORT( M_comm, ierr );

    if ( matfilename != "NULL" )
    {
        ierr = PetscViewerASCIIOpen( M_comm,
                                     matfilename.c_str(),
                                     &petsc_viewer );
        CHKERRABORT( M_comm,ierr );

        ierr = PetscViewerSetFormat ( petsc_viewer,
                                      PETSC_VIEWER_ASCII_MATLAB );
        CHKERRABORT( M_comm,ierr );

        ierr = MatView ( M_mat, petsc_viewer );
        CHKERRABORT( M_comm,ierr );
    }

    ierr = PetscViewerDestroy( &petsc_viewer );
    CHKERRABORT( M_comm,ierr );
}


void
MatrixPetsc::printBinary(std::string const& filename) const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int ierr = 0;

    PetscObjectSetName((PetscObject)M_mat,boost::filesystem::path("out_"+filename).stem().string().c_str());
    PetscViewer petsc_viewer;

    ierr = PetscViewerCreate ( M_comm, &petsc_viewer );
    CHKERRABORT( M_comm, ierr );

    if ( filename != "NULL" )
    {
        ierr = PetscViewerBinaryOpen(M_comm,
                                     filename.c_str(),
                                     FILE_MODE_WRITE,
                                     &petsc_viewer);
        CHKERRABORT( M_comm,ierr );

        ierr = MatView ( M_mat, petsc_viewer );
        CHKERRABORT( M_comm,ierr );
    }

    ierr = PetscViewerDestroy( &petsc_viewer );
    CHKERRABORT( M_comm,ierr );
}

void
MatrixPetsc::loadBinary(std::string const& filename)
{
    this->clear();

    int ierr = 0;

    PetscObjectSetName((PetscObject)M_mat,boost::filesystem::path("in_"+filename).stem().string().c_str());
    PetscViewer petsc_viewer;

    ierr = PetscViewerCreate ( M_comm, &petsc_viewer );
    CHKERRABORT( M_comm, ierr );

    if ( filename != "NULL" )
    {
        ierr = PetscViewerBinaryOpen(M_comm,
                                     filename.c_str(),
                                     FILE_MODE_READ,
                                     &petsc_viewer);
        CHKERRABORT( M_comm,ierr );

        ierr = MatLoad ( M_mat, petsc_viewer );
        CHKERRABORT( M_comm,ierr );

        M_is_initialized = true;
    }

    ierr = PetscViewerDestroy( &petsc_viewer );
    CHKERRABORT( M_comm,ierr );
}

void
MatrixPetsc::clear()
{
    int ierr = 0;

    PetscBool pinit;
    PetscInitialized( &pinit );

    if (pinit && M_is_initialized)
    {
        ierr = MatDestroy(&M_mat);
        CHKERRABORT( M_comm, ierr );

        M_is_initialized = false;
    }
}

bool
MatrixPetsc::closed() const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    int ierr=0;
    PetscBool assembled;

    ierr = MatAssembled( M_mat, &assembled );
    CHKERRABORT( M_comm,ierr );

    return ( assembled == PETSC_TRUE );
}

bool
MatrixPetsc::isSymmetric() const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    PetscBool b;
    MatIsSymmetric( M_mat, 1e-13, &b );

    return b;
}

void
MatrixPetsc::diagonal( VectorPetsc& out ) const
{
    int ierr = MatGetDiagonal( M_mat, out.vec() );
    CHKERRABORT( M_comm,ierr );
}

void
MatrixPetsc::on(std::vector<int> const& flags, VectorPetsc& rhs)
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");

    ASSERT(this->size2() == rhs.size(), "invalid right-hand side");

    this->close();
    // apply homogeneous dirichlet boundary conditions
    VectorPetsc values;
    values.init(rhs.size());
    //values.zero();

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    MatSetOption( M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE );
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE );
#else
    MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS );
#endif

    MatZeroRowsColumns(M_mat, flags.size(), flags.data(), 1.0, values.vec(), rhs.vec() );

#if 0
    // first step: elimination of dirichlet rows
    for (int flag : flags)
    {
        const PetscInt* petsc_columns;
        const PetscScalar* petsc_data;
        int ncols;
        int ierr = MatGetRow( M_mat, flag, &ncols, &petsc_columns,PETSC_NULL);
        CHKERRABORT( M_comm, ierr );

        std::vector<double> data(ncols,0.);
        for (int jj=0; jj<ncols; ++jj)
        {
            if (petsc_columns[jj]==flag)
                data[jj] = 1.;
        }

        // this->close();

        this->setMatrix(&flag, 1,
                        (int *)petsc_columns, ncols, &data[0]);

        // this->close();

        rhs(flag) = 0.;
    }

    // in-place matrix transposition
    int ierr = MatTranspose( M_mat, MAT_INITIAL_MATRIX, &M_mat );
    CHKERRABORT( M_comm, ierr );

    // dirichlet boundary application should be symmetric
    // second step: elimination of dirichlet columns
    for (int flag : flags)
    {
        const PetscInt* petsc_columns;
        const PetscScalar* petsc_data;
        int ncols;
        int ierr = MatGetRow( M_mat, flag, &ncols, &petsc_columns,PETSC_NULL);
        //int ierr = MatGetRow( M_mat, flag, &ncols, &petsc_columns,&petsc_data);
        CHKERRABORT( M_comm, ierr );

        std::vector<double> data(ncols,0.);
        for (int jj=0; jj<ncols; ++jj)
        {
            if (petsc_columns[jj]==flag)
                data[jj] = 1.;
        }

#if 0
        if (tt == 0)
        {
            std::cout<<"flag = "<< flag <<"\n";
            for (int jj=0; jj<ncols; ++jj)
                std::cout<<"petsc_data["<< jj << "]= "<< petsc_data[jj] <<"\n";

            for (int jj=0; jj<ncols; ++jj)
                std::cout<<"petsc_columns["<< jj << "]= "<< petsc_columns[jj] <<"\n";
        }
        ++tt;
#endif

        this->close();

        this->setMatrix(&flag, 1,
                            (int *)petsc_columns, ncols, &data[0]);

        this->close();
    }
#endif
}

typename MatrixPetsc::value_type
MatrixPetsc::energy(VectorPetsc& u) const
{
    ASSERT(M_is_initialized, "MatrixPetsc not properly initialized");
    ASSERT(this->size2() == u.size(), "invalid right-hand side");

    this->close();

    PetscScalar e;
    VectorPetsc v;
    v.init(u.size());
    MatMult( M_mat, u.vec(), v.vec() );
    VecDot( u.vec(), v.vec(), &e );

    return e;
}

} // Nextsim
