/*
	To store a variable indexd by phonon scattering
	* 
	Isak Kilen @ 2018
*/

#ifndef __INDEXSET2D_SCREENING_H_INCLUDED__
#define __INDEXSET2D_SCREENING_H_INCLUDED__


#include <iostream>
#include "classSortWeight2.h"

using namespace std;


class indexSet2d_screening
{
	public:
		indexSet2d_screening()
		{
			weight_e= NULL; 
			weight_h= NULL; 
			i_q 	= NULL; 
			i_kq 	= NULL; 
			b_q 	= NULL; 
			b_kq 	= NULL; 
			index_length_row = 0;
			index_length_col = NULL;
			
			indexSort = NULL;
		}
		
		~indexSet2d_screening()
		{
			if (weight_e != NULL)
			{
				for(int i = 0; i < index_length_row; i++)
				{
					delete [] weight_e[i];
					delete [] weight_h[i];
					delete [] i_kq[i];
					delete [] i_q[i];

					for(int j = 0; j < index_length_col[i]; j++)
					{
						delete [] b_q[i][j];
						delete [] b_kq[i][j];
					}

					delete [] b_q[i];
					delete [] b_kq[i];
				}
				delete [] weight_e;
				delete [] weight_h;
				delete [] i_kq;
				delete [] i_q;
				delete [] b_q;
				delete [] b_kq;
				delete [] index_length_col;

				delete indexSort;
			}
		}

		indexSet2d_screening(const indexSet2d_screening &obj)
		{
			index_length_row = obj.index_length_row;

			if (obj.weight_e!= NULL)
			{
				index_length_col = new int[index_length_row];
				for(int i = 0; i < index_length_row; i++)
				{
					index_length_col[i] = obj.index_length_col[i];
				}

				weight_e= new double*[index_length_row];
				weight_h= new double*[index_length_row];
				i_q 	= new int*[index_length_row]; 
				i_kq 	= new int*[index_length_row]; 
				b_q 	= new double**[index_length_row]; 
				b_kq 	= new double**[index_length_row]; 
				
				for(int i = 0; i < index_length_row; i++)
				{
					weight_e[i]	= new double[index_length_col[i]];
					weight_h[i]	= new double[index_length_col[i]];
					i_q[i]		= new int[index_length_col[i]];
					i_kq[i]		= new int[index_length_col[i]];
					b_q[i]		= new double*[index_length_col[i]];
					b_kq[i]		= new double*[index_length_col[i]];
					
					for(int j = 0; j < index_length_col[i]; j++)
					{
						weight_e[i][j] 	= obj.weight_e[i][j];
						weight_h[i][j] 	= obj.weight_h[i][j];
						i_q[i][j]	= obj.i_q[i][j];
						i_kq[i][j]	= obj.i_kq[i][j];

						b_q[i][j] = new double[3];
						b_kq[i][j] = new double[3];
						for(int k = 0; k < 3; k++)
						{
							b_q[i][j][k]	= obj.b_q[i][j][k];
							b_kq[i][j][k]	= obj.b_kq[i][j][k];
						}
					}
				}

				indexSort = new sortWeight2(*obj.indexSort);
			}
		}
		
		indexSet2d_screening(int N, int *array_M)
		{
			index_length_row = N;
			index_length_col = new int[index_length_row];
			for(int i = 0; i < index_length_row; i++)
			{
				index_length_col[i] = array_M[i];
			}
			
			weight_e= new double*[index_length_row];
			weight_h= new double*[index_length_row];
			i_q 	= new int*[index_length_row]; 
			i_kq 	= new int*[index_length_row]; 
			b_q 	= new double**[index_length_row]; 
			b_kq 	= new double**[index_length_row]; 

			indexSort = new sortWeight2();
			
			for(int i = 0; i < index_length_row; i++)
			{
				weight_e[i]	= new double[index_length_col[i]];
				weight_h[i]	= new double[index_length_col[i]];
				i_q[i]		= new int[index_length_col[i]];
				i_kq[i]		= new int[index_length_col[i]];
				b_q[i]		= new double*[index_length_col[i]];
				b_kq[i]		= new double*[index_length_col[i]];
				
				for(int j = 0; j < index_length_col[i]; j++)
				{
					weight_e[i][j] 	= 0;
					weight_h[i][j] 	= 0;
					i_q[i][j]	= 0;
					i_kq[i][j]	= 0;
					b_q[i][j]	= new double[3];
					b_kq[i][j]	= new double[3];
					for(int k = 0; k < 3; k++)
					{
						b_q[i][j][k]	= 0.0;
						b_kq[i][j][k]	= 0.0;
					}
				}
			}
		}
		
		void updateWeight(int k, int index, double Fe, double Fh, int q, double beta0_q, double beta1_q, double beta2_q, int kq, double beta0_kq, double beta1_kq, double beta2_kq)
		{
			weight_e[k][index] 	= Fe;
			weight_h[k][index] 	= Fh;
			i_q[k][index]		= q;
			i_kq[k][index]		= kq;
			b_q[k][index][0] 	= beta0_q;
			b_kq[k][index][0] 	= beta0_kq;
			b_q[k][index][1] 	= beta1_q;
			b_kq[k][index][1] 	= beta1_kq;
			b_q[k][index][2] 	= beta2_q;
			b_kq[k][index][2] 	= beta2_kq;
		}

		void sortIndices()
		{
			// Initialize temporary arrays
			double **weight_dummy	= new double*[index_length_row];
			int **i_dummy 		= new int*[index_length_row]; 
			double ***b_dummy 	= new double**[index_length_row]; 

			for(int i = 0; i < index_length_row; i++)
			{
				weight_dummy[i]	= new double[index_length_col[i]];
				i_dummy[i]	= new int[index_length_col[i]];
				b_dummy[i]	= new double*[index_length_col[i]];
				
				for(int j = 0; j < index_length_col[i]; j++)
				{
					weight_dummy[i][j]	= 0;
					i_dummy[i][j]		= 0;
					b_dummy[i][j]		= new double[3];
					for(int k = 0; k < 3; k++)
					{
						b_dummy[i][j][k]= 0.0;
					}
				}
			}
			for(int i = 0; i < index_length_row; i++)
			{
				// Sort indices, does not change arrays
				indexSort->initSort(i_q[i], i_kq[i], index_length_col[i]);
				int *index_list = indexSort->runSort();

				// Sort Q vectors
				for(int j = 0; j < index_length_col[i]; j++)
				{
					weight_dummy[i][j]	= weight_e[i][index_list[j]];
					i_dummy[i][j]		= i_q[i][index_list[j]];
					for(int k = 0; k < 3; k++)
					{
						b_dummy[i][j][k]= b_q[i][index_list[j]][k];
					}
				}
				for(int j = 0; j < index_length_col[i]; j++)
				{
					weight_e[i][j]	= weight_dummy[i][j];
					i_q[i][j]	= i_dummy[i][j];
					for(int k = 0; k < 3; k++)
					{
						b_q[i][j][k]= b_dummy[i][j][k];
					}
				}
				// Sort K-Q vectors
				for(int j = 0; j < index_length_col[i]; j++)
				{
					weight_dummy[i][j]	= weight_h[i][index_list[j]];
					i_dummy[i][j]		= i_kq[i][index_list[j]];
					for(int k = 0; k < 3; k++)
					{
						b_dummy[i][j][k]= b_kq[i][index_list[j]][k];
					}
				}
				for(int j = 0; j < index_length_col[i]; j++)
				{
					weight_h[i][j]	= weight_dummy[i][j];
					i_kq[i][j]	= i_dummy[i][j];
					for(int k = 0; k < 3; k++)
					{
						b_kq[i][j][k]= b_dummy[i][j][k];
					}
				}
			}

			// Delete dummy arrays
			for(int i = 0; i < index_length_row; i++)
			{
				delete [] weight_dummy[i];
				delete [] i_dummy[i];

				for(int j = 0; j < index_length_col[i]; j++)
				{
					delete [] b_dummy[i][j];
				}

				delete [] b_dummy[i];
			}
			delete [] weight_dummy;
			delete [] i_dummy;
			delete [] b_dummy;
			
		}

		
		inline double F_e(int i, int index)
		{
			return weight_e[i][index];
		}
		inline double F_h(int i, int index)
		{
			return weight_h[i][index];
		}
		
		inline int I_kq(int i, int index)
		{
			return i_kq[i][index];
		}

		inline int I_q(int i, int index)
		{
			return i_q[i][index];
		}
		
		inline void beta_kq(int i, int index, double *b0, double *b1, double *b2)
		{
			*b0 = b_kq[i][index][0];
			*b1 = b_kq[i][index][1];
			*b2 = b_kq[i][index][2];
		}

		inline void beta_q(int i, int index, double *b0, double *b1, double *b2)
		{
			*b0 = b_q[i][index][0];
			*b1 = b_q[i][index][1];
			*b2 = b_q[i][index][2];
		}
		
		inline int get_num_cols(int i)
		{
			return index_length_col[i];
		}
	
	private:
	
		double **weight_e; // Primary variable
		double **weight_h; // Primary variable
		int **i_kq; // k-q or k+q
		int **i_q; // k-q or k+q
		double ***b_kq; // Interpolation constant
		double ***b_q; // Interpolation constant
		int index_length_row; // Number of row elements
		int *index_length_col; // Number of col elements
		
		sortWeight2 *indexSort;
};

#endif







