/*
 * mon.cpp
 *
 *  Created on: Oct 23, 2018
 *      Author: xuyihan
 */
#include "mon.h"
#include <iostream>
#include <vector>
using namespace std;

void Mon::show_details() {
	wait();
	cout << "input1: " << a << endl;
	cout << "input2: " << b << endl;
	cout << "result: " << f << endl;
	sc_assert(f == a * b);
	wait();
	cout << "input1: " << a << endl;
	cout << "input2: " << b << endl;
	cout << "result: " << f << endl;
	sc_assert(f == a * b);
	wait();
	cout << "input1: " << a << endl;
	cout << "input2: " << b << endl;
	cout << "result: " << f << endl;
	sc_assert(f == a * b);
	wait();
	cout << "input1: " << a << endl;
	cout << "input2: " << b << endl;
	cout << "result: " << f << endl;
	sc_assert(f == a * b);
	wait();
	cout << "input1: " << a << endl;
	cout << "input2: " << b << endl;
	cout << "result: " << f << endl;
	sc_assert(f == a * b);
	wait();
	cout << "input1: " << a << endl;
	cout << "input2: " << b << endl;
	cout << "result: " << f << endl;
	sc_assert(f == a * b);
	wait();
	cout << "input1: " << a << endl;
	cout << "input2: " << b << endl;
	cout << "result: " << f << endl;
	sc_assert(f == a * b);
	sc_stop;
}



