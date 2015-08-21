/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, BamgQuadtree.h)*/
#ifndef _BAMGQUADTREE_H
#define _BAMGQUADTREE_H

#include "./DataSet.h"

namespace bamg {

	const int  MaxDepth = 30;
	const long MaxISize = ( 1L << MaxDepth);  // = 2^30 : 010000000000..000 (bitwise operation)

	class BamgVertex;

	class BamgQuadtree{

		private:

			/*A quadtree box contains a maximum of 4 vertices. 4 other quadtree boxes are
			 * created if a fifth vertex is added to the same box. A Quadtree box is therefore
			 * composed of EITHER:
			 * - up to 4 vertices
			 * - 4 "sub" quadtree boxes*/
			class BamgQuadtreeBox: public Object{ 
				public:
					int nbitems; // number of current vertices in the box
					union{
						BamgQuadtreeBox* b[4];
						BamgVertex*  v[4];
					};
					/*Object functions*/
					void    Echo()       {_error_("not implemented yet"); };
					void    DeepEcho()   {_error_("not implemented yet"); };
					int     Id()         {_error_("not implemented yet"); };
					int     ObjectEnum() {_error_("not implemented yet"); };
					Object *copy()       {_error_("not implemented yet"); };
					void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ _error_("not implemented yet!"); };
			};

			/*BamgQuadtree private Fields*/
			DataSet* boxcontainer;

		public:

			/*BamgQuadtree public Fields*/
			BamgQuadtreeBox* root;
			long         NbBamgQuadtreeBox;
			long         NbVertices;

			BamgQuadtree();
			BamgQuadtree(Mesh *t,long nbv=-1);
			~BamgQuadtree();

			BamgVertex      *NearestVertex(Icoor1 i,Icoor1 j);
			BamgQuadtreeBox *NewBamgQuadtreeBox(void);
			BamgVertex      *ToClose(BamgVertex &,double ,Icoor1,Icoor1);
			void             Add(BamgVertex &w);
	};
}
#endif
