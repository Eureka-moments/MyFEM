#include "Test.h"
#include <QtWidgets/QApplication>
//#include "global.h"

int main(int argc, char* argv[])
{
    
    QApplication a(argc, argv);
    Test w;
    w.show();
	int rvalue = a.exec();
    w.close_click();
    return rvalue;
}