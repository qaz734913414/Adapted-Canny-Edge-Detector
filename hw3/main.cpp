/*
 * main.cpp
 *
 *  Created on: Oct 23, 2018
 *      Author: xuyihan
 */
#include "systemc.h"
#include "top.h"
int sc_main(int argc, char* argv[]) {
	Top top("top");
	sc_start();
	return 0;
}

