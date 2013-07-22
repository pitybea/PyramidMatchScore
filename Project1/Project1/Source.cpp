#include <stdlib.h>
#include <stdio.h>


class sth
{
public:
	int myi;
	void operator ++ ()
	{
		myi++;
	};

	void operator + (int b)
	{
		myi+=b;
	};

	void operator () (int b,int c)
	{
			myi+=b;
			myi+=c;
	}
};

int main()
{
	sth a;

	a.myi=15;
	a(5,6);
	printf("%d",a.myi);

	getchar();
	return 0;
}