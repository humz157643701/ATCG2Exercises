#include <iostream>
#include <Window.h>
int main(void)
{
	try
	{	
		Window wnd;
		wnd.run();
	}
	catch (std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknown exception\n";
	}
	int dummy = 0;
}