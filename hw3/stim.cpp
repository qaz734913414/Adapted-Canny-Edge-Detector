/*
 * Stim.cpp
 *
 *  Created on: Oct 23, 2018
 *      Author: xuyihan
 */
#include "stim.h"
#include <iostream>
#include <vector>

void Stim::stimulus() {
	wait();
	/*
	A = 100;
	B = 200;
	wait();
	A = -10;
	B = 23;
	wait();
	A = 25;
	B = -3;
	wait();*/
	a = 1;
	b = 6;
	wait();
	a = 2;
	b = 6;
	wait();
	a = 3;
	b = 6;
	wait();
	a = 4;
	b = 6;
	wait();
	a = 5;
	b = 6;
	wait();
	a = 6;
	b = 6;
	wait();
	a = 7;
	b = 6;
	wait();
	sc_stop();
}

