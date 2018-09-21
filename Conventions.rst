neXtSIM conventions
===================

Git branching and merging
-------------------------

For the code we adopt the following system for branching and merging:

1. **master** branch: numbered releases of the code. Never edited. Merged from *develop* and *hot fix* branches (see notes on workflow below). Long living.
2. **develop** branch: rather stable version of code under development. Never edited. Merged from topic specific issue branches. Long living.
3. **issue<NNN>_<short-heading>** - issue specific branches (NNN = issue number). Main working area. Short living. Branched from, and merged back into develop.
4. **hotfix<NNN>_<short-heading>** - branches that are specific to a hotfix issue. Hotfixes are bugfixes on master that should be fixed as soon as possible.

.. note::

   1. Never edit code in the master or develop branch. Always make a new branch for your edits.
   2. A new branch should be very specific to only one problem. It should be short living.
   3. Commit often.
   4. Branch often.
   5. Branch only from master or from develop.
   6. Create pull requests for your branches and **always** assign a reviewer to merge and delete the branch, and close the issue. Always include head developer in pull requests, but don't rely on only him to review your edits.

How to report and handle new issues (bugs, improvements, new features, etc.)
----------------------------------------------------------------------------

If you discover a bug check that no one else has reported the same issue on https://github.com/nansencenter/nextsim/issues. If the bug has not been reported create an issue and assign someone to fix it (possibly yourself). Please notify people in person if you're assigning issues to them.

If you would like to suggest an improvement or a new feature check that no one else has made a similar request on https://github.com/nansencenter/nextsim/issues. If this is not the case create a new issue, assigning or mentioning anyone you think could be affected or interested by your suggestion, always including the head developer.

If you have been assigned an issue on https://github.com/nansencenter/nextsim/issues address it using the following steps. For issues requiring a hotfix (bugs on the master branch that should be fixed as soon as possible):

        1. Branching off from **master**, create an issue branch on your local system named **hotfix<NNN>_<short-heading>** where NNN is the issue number from GitHub. This will be the main (short living) working area.
        2. Write the necessary code.
                   a. Make sure you test your modifications well. 
                   b. Feel free to commit and push your issue branch often
        3. Once the issue is fixed merge the **master** branch back into your issue branch and resolve any conflicts.
        4. Create a pull request on GitHub to merge the issue branch back into **master**. Always include at least one reviewer who will then merge and delete the issue branch.
        5. Merge the **master** branch into **develop** in your local repository, resolve conflicts, test, and push the updated **develop** branch.
        6. Close the issue.

For issues not requiring a hotfix (less urgent bug-fixes and feature requests):

        1. Branching off from **develop**, create an issue branch on your local system named **issue<NNN>_<short-heading>** where NNN is the issue number from GitHub. This will be the main (short living) working area.
        2. Write the necessary code.
                   a. Make sure you test your modifications well. 
                   b. Feel free to commit and push your issue branch often
        3. Once the issue is fixed merge the **develop** branch back into your issue branch and resolve any conflicts.
        4. Create a pull request on GitHub to merge the issue branch back into **develop**. Always include at least one reviewer who will then merge and delete the issue branch, and close the issue.

Code conventions
-------------------

.. note:: This section requires substantially more work

        * neXtSIM is written using ISO C++11
        * Indentation is done using a set of four (4) spaces (no tabs)
        * All array opperations should be done using std::vectors - not C-style arrays
        * The use of C-style pointers, new, and delete is strongly discouraged
        * Names of (class) global variables start with M_, except for diagnostic variables which start with D_. The names of prognostic global variables should start with P_, but this is not yet implemented.
        * Variables should be named using the underscore-as-space convention (such_as_this), while functions should be named using the capital-instead-of-space convention (suchAsThis).
        * Names and values of physical constants recide in ``model/constants.hpp``
        * Runtime modifiable model parameters are set in ``model/options.cpp`` and can be set using option files.


Commenting the code; conventions
--------------------------------

	Doxygen is used to generate a documentation of NeXtSIM, which can be found here: file:///Users/verans/Developer/nextsim/doc/html/index.html        
	When developing the code, re-writing or creating new functions, please comment it as follow.

	* Create a header for the function, structured as below:
	  // ---------------------------------------------------
	  //! The function does (...)
	  //! * Notes, particularities, etc.
	  //! Called by the <functionWhichCallsThatParticularFunction>() function.

	  Precisions:
	  - Any comment preceded by "//!" will be treated as a description of a function/variable in Doxygen and will appear in the documentation.
	  - The first line of the header should describe briefly what the function does.
	  - The following lines should add any valuable information or comment on the function. "*" produces a bullet point in Doxygen. 
	  - The last line should indicate the function() or functions() that directly call the function in question. 
            The () should be included at the end of the function's name. This automatically creates a link to this function in Doxygen, enabling the user to trace back the origin and order of call for the functions.
	  - The header should be placed just before the function declaration, i.e., the type and function name and opening bracket, { .

	* Place logical and relevant comments within the function, especially if the function is long, and calls many functions itself. 
          Enumerate the different logical steps taken within the function using "//! - 1) ...description...". 
	  This creates en indented, numbered list in the Doxygen documentation that indicates that function's structure. 
	  It is better to avoid placing "//!" in front of comments that express a doubt or a "to do" type of idea, as those will interfere with the documentation of the current version of the code. 

	* After the closing bracket, recall the name of the function as follow:
          //<functionName>
	  This will not appear in the Doxygen documentation, but eases the reading when scrolling through the code, i.e. the finiteelement.cpp file. 

	* Leave 2 blank lines between each function declaration.

	
	Below is an example of a documented function:

	//------------------------------------------------------------------------------------------------------
	//! Initializes constants, dataset descriptions, the time, mesh, variables, forcings, bathymetry, moorings and drifters.
	//! * Also outputs restarts for debugging.
	//! Called by the run() function.
	void
	FiniteElement::init()
	{

    	//! - 1) Initializes everything that doesn't depend on the mesh (constants, dataset descriptions and time) using the initOptAndParam() function,

    		M_comm.barrier();

    		pcpt = 0;
    		mesh_adapt_step=0;
    		had_remeshed=false;

    		this->initOptAndParam();
    		M_current_time = time_init /*+ pcpt*time_step/(24*3600.0)*/;

    	//! - 2) Initializes the mesh using the initMesh() function,
   
		this->initMesh();

    		if (M_rank==0)
    		{
		LOG(INFO) << "-----------------------Simulation started on "<< Nextsim::current_time_local() <<"\n";
       	 	LOG(INFO) <<"TIMESTEP= "<< time_step <<"\n";
        	LOG(INFO) <<"DURATION= "<< duration <<"\n";
    		}

    		// We need to set the scale_coeff et al after initialising the mesh - this was previously done in initConstants
    		// The mean resolution of the small_arctic_10km mesh is 7446.71 m. Using 74.5 gives scale_coef = 0.100022, for that mesh
    		boost::mpi::broadcast(M_comm, M_res_root_mesh, 0);


		(.........)

 	}//init
        
