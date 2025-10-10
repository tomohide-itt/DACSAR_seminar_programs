#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

// main関数
int main()
{
	int dim = 2;
	int num_node = 4;
	int num_elem = 1;
	int npe = 4; // node per element
	int gp_order = 3; // 積分次数

	// 節点座標
	std::vector<std::vector<double>> xy( num_node+1, std::vector<double>(dim, 0.0 ) );
	xy[1][0] = 0.0;
	xy[1][1] = 0.0;
	xy[2][0] = 1.0;
	xy[2][1] = 0.0;
	xy[3][0] = 1.0;
	xy[3][1] = 1.0;
	xy[4][0] = 0.0;
	xy[4][1] = 1.0;

	// 要素の情報
	std::vector<std::vector<int>> nod( num_elem+1, std::vector<int>( npe ) );
	nod[1][0] = 1;
	nod[1][1] = 2;
	nod[1][2] = 3;
	nod[1][3] = 4;

	// Ｄマトリクスをつくる
	double E  = 1.0e+6; //ヤング率
	double nu = 0.333;	//ポアソン比
	double K = E/3.0/(1.0-2.0*nu); //体積弾性係数
	double G = E/2.0/(1.0+nu); //せん断弾性係数
	std::vector<double> Dmat( 16, 0.0 );
	// D = K(1x1) +2G(I - 1/3(1x1))
	std::vector<double> delxdel( 16, 0.0 );
	delxdel[ 0] = 1.0;
	delxdel[ 1] = 1.0;
	delxdel[ 2] = 1.0;
	delxdel[ 4] = 1.0;
	delxdel[ 5] = 1.0;
	delxdel[ 6] = 1.0;
	delxdel[ 8] = 1.0;
	delxdel[ 9] = 1.0;
	delxdel[10] = 1.0;
	//
	std::vector<double> idn( 16, 0.0 );
	idn[ 0] = 1.0;
	idn[ 5] = 1.0;
	idn[10] = 1.0;
	idn[15] = 0.5;
	//
	for( int i=0; i<16; i++ )
	{
		Dmat[i] = (K-2.0/3.0*G)*delxdel[i] + 2.0*G*idn[i];
	}
	//+++
	std::cout << "Dmat:" << std::endl;
	std::cout << std::scientific << std::setprecision(5);
	for( int i=0; i<4; i++ )
	{
		for( int j=0; j<4; j++ )
		{
			std::cout << std::setw(15) << Dmat[i*4+j];
		}
		std::cout << std::endl;
	}
	//---

	// 積分点の位置
	std::vector<double> gp_pos( gp_order );
	gp_pos[0] = -sqrt( 3.0/5.0 );
	gp_pos[1] =  0.0;
	gp_pos[2] =  sqrt( 3.0/5.0 );

	// 積分点の重み
	std::vector<double> gp_wei( gp_order );
	gp_wei[0] = 5.0/9.0;
	gp_wei[1] = 8.0/9.0;
	gp_wei[2] = 5.0/9.0;

	int num_gp = pow( gp_order, dim );

	for( int ng=0; ng<num_gp; ng++ )
	{
		double xi, eta;
		if( ng==0 ){ xi=gp_pos[0]; eta=gp_pos[0]; }
		if( ng==1 ){ xi=gp_pos[0]; eta=gp_pos[1]; }
		if( ng==2 ){ xi=gp_pos[0]; eta=gp_pos[2]; }
		if( ng==3 ){ xi=gp_pos[1]; eta=gp_pos[0]; }
		if( ng==4 ){ xi=gp_pos[1]; eta=gp_pos[1]; }
		if( ng==5 ){ xi=gp_pos[1]; eta=gp_pos[2]; }
		if( ng==6 ){ xi=gp_pos[2]; eta=gp_pos[0]; }
		if( ng==7 ){ xi=gp_pos[2]; eta=gp_pos[1]; }
		if( ng==8 ){ xi=gp_pos[2]; eta=gp_pos[2]; }
		// dN/dr (4x2)
		std::vector<std::vector<double>> dNdr( num_node, std::vector<double>(dim,0.0) );
		dNdr[0][0] = -0.25*( 1.0-eta );
		dNdr[1][0] =  0.25*( 1.0-eta );
		dNdr[2][0] =  0.25*( 1.0+eta );
		dNdr[3][0] = -0.25*( 1.0+eta );
		dNdr[0][1] = -0.25*( 1.0-xi );
		dNdr[1][1] = -0.25*( 1.0+xi );
		dNdr[2][1] =  0.25*( 1.0+xi );
		dNdr[3][1] =  0.25*( 1.0+xi );
		//+++
		std::cout << "dNdr[ng=" << ng << "]:" << std::endl;
		for( int i=0; i<dim; i++ )
		{
			for( int j=0; j<num_node; j++ )
			{
				std::cout << std::setw(15) << dNdr[j][i];
			}
			std::cout << std::endl;
		}
		//---
		// J = dNdrT xy (2x2)
		std::vector<std::vector<double>> J( dim, std::vector<double>(dim, 0.0) );
		for( int i=0; i<dim; i++ )
		{
			for( int j=0; j<dim; j++ )
			{
				J[i][j] = 0.0;
				for( int k=0; k<4; k++ )
				{
					J[i][j] += dNdr[k][i] * xy[k+1][j];
				}
			}
		}
		//+++
		std::cout << "J[ng=" << ng << "]:" << std::endl;
		for( int i=0; i<dim; i++ )
		{
			for( int j=0; j<dim; j++ )
			{
				std::cout << std::setw(15) << J[i][j];
			}
			std::cout << std::endl;
		}
		//---
		// detJ
		double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
		// Jinv (2x2)
		std::vector<std::vector<double>> Jinv( dim, std::vector<double>(dim,0.0) );
		Jinv[0][0] =  J[1][1]/detJ;
		Jinv[0][1] = -J[1][0]/detJ;
		Jinv[1][0] = -J[0][1]/detJ;
		Jinv[1][1] =  J[0][0]/detJ;
		// derivN (4x2)
		std::vector<std::vector<double>> derivN( num_node, std::vector<double>(dim, 0.0) );
		for( int i=0; i<num_node; i++ )
		{
			for( int j=0; j<dim; j++ )
			{
				derivN[i][j] = 0.0;
				for( int k=0; k<dim; k++ )
				{
					derivN[i][j] += Jinv[j][k] * dNdr[i][k];
				}
			}
		}
		//+++
		std::cout << "derivN[ng=" << ng << "]:" << std::endl;
		for( int i=0; i<dim; i++ )
		{
			for( int j=0; j<num_node; j++ )
			{
				std::cout << std::setw(15) << derivN[j][i];
			}
			std::cout << std::endl;
		}
		//---

	}
		

	

	return 0;
}


