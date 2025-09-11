#include <iostream>

// 関数の前方宣言
void show_sum_int( int a, int b );

// main関数
int main()
{
	show_sum_int( 3, 5 );
	return 0;
}

// 二つの整数を足して表示する関数
void show_sum_int( int a, int b )
{
	std::cout << a << " + " << b << " = " << a+b << std::endl;
}


