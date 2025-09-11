#include <iostream>

// 関数の前方宣言
int sum_int( int a, int b );

// main関数
int main()
{
	int a = sum_int( 3, 5 );
	std::cout << a << std::endl;
	return 0;
}

// 二つの整数を足して結果を返す関数
int sum_int( int a, int b )
{
	return a + b;
}


