/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   connectivity.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Oct 15 13:45:33 2018
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cmath>
#include <string>


#ifndef __MeshConnectivity_HPP
#define __MeshConnectivity_HPP 1

namespace Nextsim
{


int columnCompare(std::vector<int> const& col, int m, int n, int i, int j);
void columnSort(std::vector<int>& col, int m, int n);
void columnSwap(std::vector<int>& col, int m, int n, int icol1, int icol2);
void sortHeapExternal(int n, int *indx, int *i, int *j, int isgn);
void elementConnectivity(std::vector<int> const& triangle_nodes, std::vector<int>& triangle_neighbors, int triangle_num);

int
columnCompare(std::vector<int> const& col, int m, int n, int i, int j)
{
	if (i == j)
	{
		return 0;
	}

	int k = 1;

	while (k <= m)
	{
		if (col[k-1+(i-1)*m] < col[k-1+(j-1)*m])
		{
			return (-1);
		}
		else if (col[k-1+(j-1)*m] < col[k-1+(i-1)*m])
		{
			return 1;
		}

		k = k + 1;
	}

	return 0;
}

void
columnSort(std::vector<int>& col, int m, int n)
{
	int i = 0;
	int indx = 0;
	int isgn = 0;
	int j = 0;

	for ( ; ; )
	{
		sortHeapExternal(n, &indx, &i, &j, isgn);

		if (0 < indx)
		{
			columnSwap(col, m, n, i, j);
		}

		else if (indx < 0)
		{
			isgn = columnCompare(col, m, n, i, j);
		}
		else if (indx == 0)
		{
			break;
		}
	}
}

void
columnSwap(std::vector<int>& col, int m, int n, int icol1, int icol2)
{
	if (icol1 == icol2)
	{
		return;
	}

	int offset = 1;
	int t;

	for (int i = 0; i < m; i++)
	{
		t = col[i+(icol1-offset)*m];
		col[i+(icol1-offset)*m] = col[i+(icol2-offset)*m];
		col[i+(icol2-offset)*m] = t;
	}
}


void
sortHeapExternal(int n, int *indx, int *i, int *j, int isgn )
{
	static int i_save = 0;
	static int j_save = 0;
	static int k = 0;
	static int k1 = 0;
	static int n1 = 0;

	//  INDX = 0: This is the first call.

	if (*indx == 0)
	{

		i_save = 0;
		j_save = 0;
		k = n / 2;
		k1 = k;
		n1 = n;
	}

	//  INDX < 0: The user is returning the results of a comparison.

	else if (*indx < 0)
	{
		if (*indx == -2)
		{
			if (isgn < 0)
			{
				i_save = i_save + 1;
			}

			j_save = k1;
			k1 = i_save;
			*indx = -1;
			*i = i_save;
			*j = j_save;
			return;
		}

		if (0 < isgn)
		{
			*indx = 2;
			*i = i_save;
			*j = j_save;
			return;
		}

		if (k <= 1)
		{
			if (n1 == 1)
			{
				i_save = 0;
				j_save = 0;
				*indx = 0;
			}
			else
			{
				i_save = n1;
				j_save = 1;
				n1 = n1 - 1;
				*indx = 1;
			}
			*i = i_save;
			*j = j_save;
			return;
		}
		k = k - 1;
		k1 = k;
	}

	//  0 < INDX: the user was asked to make an interchange.

	else if (*indx == 1)
	{
		k1 = k;
	}

	for ( ; ; )
	{

		i_save = 2 * k1;

		if (i_save == n1)
		{
			j_save = k1;
			k1 = i_save;
			*indx = -1;
			*i = i_save;
			*j = j_save;
			return;
		}
		else if (i_save <= n1)
		{
			j_save = i_save + 1;
			*indx = -2;
			*i = i_save;
			*j = j_save;
			return;
		}

		if (k <= 1)
		{
			break;
		}

		k = k - 1;
		k1 = k;
	}

	if (n1 == 1)
	{
		i_save = 0;
		j_save = 0;
		*indx = 0;
		*i = i_save;
		*j = j_save;
	}
	else
	{
		i_save = n1;
		j_save = 1;
		n1 = n1 - 1;
		*indx = 1;
		*i = i_save;
		*j = j_save;
	}

}

void
elementConnectivity(std::vector<int> const& triangle_nodes, std::vector<int>& triangle_neighbors, int triangle_num)
{
	triangle_neighbors.resize(3*triangle_num);
	std::vector<int> col(4*(3*triangle_num));
	int triangle_order = 3;
	int i;
	int icol;
    int j;
    int k;
    int side1;
    int side2;
    int tri1;
    int tri2;

    for (int tri = 0; tri < triangle_num; tri++)
    {
        i = triangle_nodes[0+tri*triangle_order];
        j = triangle_nodes[1+tri*triangle_order];
        k = triangle_nodes[2+tri*triangle_order];

        if (i < j)
        {
            col[0+(3*tri+0)*4] = i;
            col[1+(3*tri+0)*4] = j;
            col[2+(3*tri+0)*4] = 3;
            col[3+(3*tri+0)*4] = tri + 1;
        }
        else
        {
            col[0+(3*tri+0)*4] = j;
            col[1+(3*tri+0)*4] = i;
            col[2+(3*tri+0)*4] = 3;
            col[3+(3*tri+0)*4] = tri + 1;
        }

        if (j < k)
        {
            col[0+(3*tri+1)*4] = j;
            col[1+(3*tri+1)*4] = k;
            col[2+(3*tri+1)*4] = 1;
            col[3+(3*tri+1)*4] = tri + 1;
        }
        else
        {
            col[0+(3*tri+1)*4] = k;
            col[1+(3*tri+1)*4] = j;
            col[2+(3*tri+1)*4] = 1;
            col[3+(3*tri+1)*4] = tri + 1;
        }

        if (k < i)
        {
            col[0+(3*tri+2)*4] = k;
            col[1+(3*tri+2)*4] = i;
            col[2+(3*tri+2)*4] = 2;
            col[3+(3*tri+2)*4] = tri + 1;
        }
        else
        {
            col[0+(3*tri+2)*4] = i;
            col[1+(3*tri+2)*4] = k;
            col[2+(3*tri+2)*4] = 2;
            col[3+(3*tri+2)*4] = tri + 1;
        }
    }

    columnSort(col, 4, 3*triangle_num);

    for (int j = 0; j < triangle_num; j++)
    {
	    for (i = 0; i < 3; i++)
	    {
		    triangle_neighbors[i+j*3] = -1;
	    }
    }

    icol = 1;

    for ( ; ; )
    {
	    if (3 * triangle_num <= icol)
	    {
		    break;
	    }

	    if (col[0+(icol-1)*4] != col[0+icol*4] ||
	         col[1+(icol-1)*4] != col[1+icol*4])
	    {
		    icol = icol + 1;
		    continue;
	    }

	    side1 = col[2+(icol-1)*4];
	    tri1 =  col[3+(icol-1)*4];
	    side2 = col[2+ icol   *4];
	    tri2 =  col[3+ icol   *4];

	    triangle_neighbors[side1-1+(tri1-1)*3] = tri2;
	    triangle_neighbors[side2-1+(tri2-1)*3] = tri1;

	    icol = icol + 2;
    }
}


} // Nextsim
#endif
