#include <limits.h>
#include <string.h>
#include <stdlib.h>

#include "./include.h"
#include "./DataSet.h"

namespace bamg {

	/*MACROS {{{*/
	/* 
	 * 
	 *    J    j
	 *    ^    ^
	 *    |    | +--------+--------+
	 *    |    | |        |        |
	 * 1X |    | |   2    |   3    |
	 *    |    | |        |        |
	 *    |    | +--------+--------+
	 *    |    | |        |        |
	 * 0X |    | |   0    |   1    |
	 *    |    | |        |        |
	 *    |    | +--------+--------+
	 *    |    +-----------------------> i
	 *    |         
	 *    |----------------------------> I
	 *              X0        X1  
	 *
	 * box 0 -> I=0 J=0 IJ=00  = 0
	 * box 1 -> I=1 J=0 IJ=01  = 1
	 * box 2 -> I=0 J=1 IJ=10  = 2
	 * box 3 -> I=1 J=1 IJ=11  = 3
	 */
#define INTER_SEG(a,b,x,y) (((y) > (a)) && ((x) <(b)))
#define ABS(i) ((i)<0 ?-(i) :(i))
#define MAX1(i,j) ((i)>(j) ?(i) :(j))
#define NORM(i1,j1,i2,j2) MAX1(ABS((i1)-(j1)),ABS((i2)-(j2)))

	//IJ(i,j,l) returns the box number of i and j with respect to l
	//if !j&l and !i&l -> 0 (box zero: lower left )
	//if !j&l and  i&l -> 1 (box one:  lower right)
	//if  j&l and !i&l -> 2 (box two:  upper left )
	//if  j&l and  i&l -> 3 (box three:upper right)
#define IJ(i,j,l)  ((j&l) ? ((i&l) ? 3:2 ) :((i&l) ? 1:0 ))

	//I_IJ(k,l) returns l if first  bit of k is 1, else 0
#define I_IJ(k,l)  ((k&1) ? l:0)
	//J_IJ(k,l) returns l if second bit of k is 1, else 0
#define J_IJ(k,l)  ((k&2) ? l:0)
	/*}}}*/
	/*DOCUMENTATION What is a BamgQuadtree? {{{
	 * A Quadtree is a very simple way to group vertices according
	 * to their locations. A square that holds all the points of the mesh
	 * (or the geometry) is divided into 4 boxes. As soon as one box
	 * hold more than 4 vertices, it is divided into 4 new boxes, etc...
	 * There cannot be more than MAXDEEP (=30) subdivision.
	 * This process is like a Dichotomy in dimension 2
	 *
	 *  + - -  -    - -    -    - - + -   - + - + - + - -     - - +
	 *  |                           |       |   | X |             |
	 *                                      + - + - +
	 *  |                           |       |   |   |             |
	 *                              + -   - + - + - +             +
	 *  |                           |       |       |             |
	 *                         
	 *  |                           |       |       |             |
	 *  + - -  -    - -    -    - - + -   - + -   - + - -     - - +
	 *  |                           |               |             |
	 *                         
	 *  |                           |               |             |
	 *                         
	 *  |                           |               |             |
	 *  |                           |               |             |
	 *  + - -  -    - -    -    - - + -   -   -   - + - -     - - +
	 *  |                           |                             |
	 *                         
	 *  |                           |                             |
	 *                         
	 *  |                           |                             |
	 *                         
	 *  |                           |                             |
	 *  |                           |                             |
	 *  |                           |                             |
	 *  |                           |                             |
	 *  |                           |                             |
	 *  + - -  -    - -    -    - - + -   -   -   -   - -     - - +
	 *
	 * The coordinate system used in a quadtree are integers to avoid
	 * round-off errors. The vertex in the lower left box has the coordinates
	 * (0 0) 
	 * The upper right vertex has the follwing coordinates:
	 * 2^30 -1           2^30 -1        in decimal
	 * 0 1 1 1 .... 1    0 1 1 1 .... 1 in binary
	 *  \--   29  --/     \--   29  --/
	 * Using binaries is therefore very easy to locate a vertex in a box:
	 * we just need to look at the bits from the left to the right (See ::Add)
	 }}}*/

	/*Constructors/Destructors*/
	BamgQuadtree::BamgQuadtree(){/*{{{*/

		/*Number of boxes and vertices*/
		NbBamgQuadtreeBox=0;
		NbVertices=0;

		/*Create container*/
		boxcontainer=new DataSet();

		/*Create Root, pointer toward the main box*/
		root=NewBamgQuadtreeBox();

		}
	/*}}}*/
	BamgQuadtree::BamgQuadtree(Mesh * t,long nbv){ /*{{{*/

		/*Number of boxes and vertices*/
		NbBamgQuadtreeBox=0;
		NbVertices=0;

		/*Create container*/
		boxcontainer=new DataSet();

		/*Create Root, pointer toward the main box*/
		root=NewBamgQuadtreeBox();

		/*Check Sizes*/
		_assert_(MaxISize>MaxICoor);

		/*Add all vertices of the mesh*/
		if (nbv==-1) nbv=t->nbv;
		for (int i=0;i<nbv;i++) Add(t->vertices[i]);

	}
	/*}}}*/
	BamgQuadtree::~BamgQuadtree() {/*{{{*/
		delete boxcontainer;
		root=NULL;
	}
	/*}}}*/

	/*Methods*/
	void  BamgQuadtree::Add(BamgVertex &w){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, BamgQuadtree.cpp/Add)*/
		BamgQuadtreeBox** pb=NULL;
		BamgQuadtreeBox*  b=NULL;

		/*Get integer coodinate of current point w*/
		long i=w.i.x, j=w.i.y;

		/*Initialize level*/
		long level=MaxISize;

		/*Get inital box (the largest)*/
		pb = &root;

		/*Find the smallest box where w is located*/
		while((b=*pb) && (b->nbitems<0)){ 

			//shift b->nbitems by -1
			b->nbitems--;

			//shifted righ by one bit: level=00000010 -> 00000001
			level >>= 1;

			//Get next subbox according to the bit value (level)
			pb = &b->b[IJ(i,j,level)];
		}

		/*OK, we have found b, a Subbox holding vertices (might be full)
		  check that the vertex is not already in the box*/
		if (b){      
			if (b->nbitems > 3 &&  b->v[3] == &w) return;
			if (b->nbitems > 2 &&  b->v[2] == &w) return;
			if (b->nbitems > 1 &&  b->v[1] == &w) return;
			if (b->nbitems > 0 &&  b->v[0] == &w) return;
		}

		/*check that l is not 0 (this should not happen as MaxDepth = 30)*/
		_assert_(level>0);

		/*Now, try to add the vertex, if the subbox is full (nbitems=4), we have to divide it
		  in 4 new subboxes*/
		while ((b= *pb) && (b->nbitems == 4)){ // the BamgQuadtreeBox is full

			/*Copy the 4 vertices in the current BamgQuadtreebox*/
			BamgVertex* v4[4];
			v4[0]= b->v[0];
			v4[1]= b->v[1];
			v4[2]= b->v[2];
			v4[3]= b->v[3];

			/*set nbitems as negative 
			 * (box full -> holds 4 pointers toward subboxes and not 4 vertices)*/
			b->nbitems = -b->nbitems;

			/*Initialize the 4 pointers toward the 4 subboxes*/
			b->b[0]=b->b[1]=b->b[2]=b->b[3]=NULL;

			/*level = 0010000 -> 0001000*/
			level >>= 1;

			/*Put the four vertices in the new boxes*/
			for (int k=0;k<4;k++){

				int          ij;
				/*bb is the new "sub"box of b where v4[k] is located*/
				BamgQuadtreeBox *bb = b->b[ij=IJ(v4[k]->i.x,v4[k]->i.y,level)];

				// alloc the BamgQuadtreeBox
				if (!bb) bb=b->b[ij]=NewBamgQuadtreeBox(); 

				/*Copy the current vertex*/
				bb->v[bb->nbitems++] = v4[k];
			}

			/*Get the subbox where w (i,j) is located*/
			pb = &b->b[IJ(i,j,level)];
		}

		/*alloc the BamgQuadtreeBox if necessary*/
		if (!(b=*pb)) b=*pb= NewBamgQuadtreeBox();

		/*Add w*/
		b->v[b->nbitems++]=&w;

		//Increase NbVertices by one (we have one new vertex)
		NbVertices++;
	}
	/*}}}*/
	BamgVertex*  BamgQuadtree::NearestVertex(Icoor1 i,Icoor1 j) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, BamgQuadtree.cpp/NearestVertex)*/

		/*Intermediaries*/
		BamgQuadtreeBox *pb[MaxDepth];
		int          pi[MaxDepth];
		Icoor1       ii[MaxDepth];
		Icoor1       jj[MaxDepth];
		int          level;
		long         n0;
		BamgQuadtreeBox *b;
		long         h0;
		long         h = MaxISize;
		long         hb= MaxISize;
		Icoor1       i0=0,j0=0;

		/*initial output as NULL (no vertex found)*/
		BamgVertex*  nearest_v=NULL;

		/*Project w coordinates (i,j) onto [0,MaxISize-1] x [0,MaxISize-1] -> (iplus,jplus)*/
		Icoor1 iplus( i<MaxISize ? (i<0?0:i) : MaxISize-1);
		Icoor1 jplus( j<MaxISize ? (j<0?0:j) : MaxISize-1);

		/*Get initial Quadtree box (largest)*/
		b = root;

		/*if the tree is empty, return NULL pointer*/
		if (!root->nbitems) return nearest_v; 

		/*else, find the smallest non-empty BamgQuadtreeBox containing  the point (i,j)*/
		while((n0=b->nbitems)<0){

			Icoor1       hb2 = hb >> 1;             //size of the current box
			int          k   = IJ(iplus,jplus,hb2); //box number (0,1,2 or 3)
			BamgQuadtreeBox *b0  = b->b[k];             //pointer toward current box

			/* break if NULL box or empty (Keep previous box b)*/
			if (( b0 == NULL) || (b0->nbitems == 0)) break;

			/*Get next Quadtree box*/
			b=b0;	
			i0 += I_IJ(k,hb2); // i orign of BamgQuadtreeBox (macro)
			j0 += J_IJ(k,hb2); // j orign of BamgQuadtreeBox 
			hb = hb2;          // size of the box (in Int)
		}

		/*The box b, is the smallest box containing the point (i,j) and
		 * has the following properties:
		 * - n0: number of items (>0 if vertices, else boxes)
		 * - hb: box size (int)
		 * - i0: x coordinate of the lower left corner
		 * - j0: y coordinate of the lower left corner*/

		/* if the current subbox is holding vertices, we are almost done*/
		if (n0>0){  
			//loop over the vertices of the box and find the closest vertex
			for(int k=0;k<n0;k++){

				/*get integer coordinates of current vertex*/
				I2 i2=b->v[k]->i;

				/*Compute norm with w*/
				h0=NORM(iplus,i2.x,jplus,i2.y);

				/*is it smaller than previous value*/
				if (h0<h){
					h = h0;
					nearest_v = b->v[k];
				}
			}
			/*return closest vertex*/
			return nearest_v;
		}

		/* general case: the current box is empty, we have to go backwards
			and find the closest not-empty box and find the closest vertex*/

		/*Initialize search variables*/
		pb[0]=b;                             //pointer toward the box b
		pi[0]=b->nbitems>0?(int)b->nbitems:4;//number of boxes in b
		ii[0]=i0;                            //i coordinate of the box lowest left corner
		jj[0]=j0;                            //j coordinate of the box lowest left corner

		/*initialize h: smallest box size, containing a vertex close to w*/
		h=hb;

		/*Main loop*/
		level=0;
		do {

			/*get current box*/
			b= pb[level];

			/*Loop over the items in current box (if not empty!)*/
			while (pi[level]){

				/*We are looping now over the items of b. k is the current index (in [0 3])*/
				pi[level]--;
				int k=pi[level];

				/*if the current subbox is holding vertices (b->nbitems<0 is subboxes)*/
				if (b->nbitems>0){
					I2 i2=b->v[k]->i;
					h0 = NORM(iplus,i2.x,jplus,i2.y);
					if (h0<h){
						h=h0;
						nearest_v=b->v[k];
					}
				}
				/*else: current box b is pointing toward 4 boxes
				 * test sub-box k and go deeper into the tree if it is non empty
				 * and contains the point w modulo a size h that is either the size of the smallest
				 * non empty box containing w, or the closest point to w (so far) */
				else{
					BamgQuadtreeBox* b0=b;

					/*if the next box exists:*/
					if((b=b->b[k])){

						/*Get size (hb) and coordinates of the current sub-box lowest left corner*/
						hb>>=1;
						Icoor1 iii = ii[level]+I_IJ(k,hb);
						Icoor1 jjj = jj[level]+J_IJ(k,hb);

						/*if the current point (iplus,jplus) is in b (modulo h), this box is good:
						 * it is holding vertices that are close to w */
						if (INTER_SEG(iii,iii+hb,iplus-h,iplus+h) && INTER_SEG(jjj,jjj+hb,jplus-h,jplus+h)){
							level++;
							pb[level]= b;
							pi[level]= b->nbitems>0 ?(int)  b->nbitems : 4  ;
							ii[level]= iii;
							jj[level]= jjj;
						}

						//else go backwards
						else{
							//shifted righ by one bit: hb=001000000 -> 01000000
							b=b0;
							hb<<=1;
						}
					}
					else{
						/*Current box is NULL, go to next subbox of b (k=k-1)*/
						b=b0;
					}
				}
			}

			/*We have found a vertex, now, let's try the other boxes of the previous level
			 * in case there is a vertex closest to w that has not yet been tested*/
			hb <<= 1;
		} while (level--);

		/*return nearest_v, nearest vertex*/
		return nearest_v;

	}
	/*}}}*/
	BamgQuadtree::BamgQuadtreeBox* BamgQuadtree::NewBamgQuadtreeBox(void){/*{{{*/

		/*Output*/
		BamgQuadtreeBox* newbox=NULL;

		/*Create and initialize a new box*/
		newbox=new BamgQuadtreeBox;
		newbox->nbitems=0;
		newbox->b[0]=NULL;
		newbox->b[1]=NULL;
		newbox->b[2]=NULL;
		newbox->b[3]=NULL;

		/*Add root to the container*/
		boxcontainer->AddObject(newbox);

		/*Increase counter*/
		NbBamgQuadtreeBox++;

		/*currentbox now points toward next quadtree box*/
		return newbox;
	}/*}}}*/
	BamgVertex*   BamgQuadtree::ToClose(BamgVertex & v,double seuil,Icoor1 hx,Icoor1 hy){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, BamgQuadtree.cpp/ToClose)*/

		const Icoor1 i=v.i.x;
		const Icoor1 j=v.i.y;
		const R2 X(v.r);
		const Metric  Mx(v.m);

		BamgQuadtreeBox * pb[ MaxDepth ];
		int  pi[ MaxDepth  ];
		Icoor1 ii[  MaxDepth ], jj [ MaxDepth];
		int l=0; // level
		BamgQuadtreeBox * b;
		Icoor1 hb =  MaxISize;
		Icoor1 i0=0,j0=0;

		//  BamgVertex *vn=0;

		if (!root->nbitems)
		 return 0; // empty tree 

		// general case -----
		pb[0]=root;
		pi[0]=root->nbitems>0 ?(int)  root->nbitems : 4  ;
		ii[0]=i0;
		jj[0]=j0;
		do {    
			b= pb[l];
			while (pi[l]--){ 	      
				int k = pi[l];

				if (b->nbitems>0){ // BamgVertex BamgQuadtreeBox none empty
					I2 i2 =  b->v[k]->i;
					if ( ABS(i-i2.x) <hx && ABS(j-i2.y) <hy )
					  {
						R2 XY(X,b->v[k]->r);
						if(LengthInterpole(Mx(XY), b->v[k]->m(XY)) < seuil){
							return b->v[k]; 
						}
					  }
				}
				else{ // Pointer BamgQuadtreeBox 
					BamgQuadtreeBox *b0=b;
					if ((b=b->b[k])){
						hb >>=1 ; // div by 2
						long iii = ii[l]+I_IJ(k,hb);
						long jjj = jj[l]+J_IJ(k,hb);

						if  (INTER_SEG(iii,iii+hb,i-hx,i+hx) && INTER_SEG(jjj,jjj+hb,j-hy,j+hy)){
							pb[++l]=  b;
							pi[l]= b->nbitems>0 ?(int)  b->nbitems : 4  ;
							ii[l]= iii;
							jj[l]= jjj;

						}
						else{
							b=b0;
							hb <<=1 ;
						}
					}
					else{
						b=b0;
					}
				}
			}
			hb <<= 1; // mul by 2 
		} while (l--);

		return 0;
	}
	/*}}}*/
}
