/*
 * mon.h
 *
 *  Created on: Oct 23, 2018
 *      Author: xuyihan
 */

#ifndef MON_H_
#define MON_H_

#include "systemc.h"
SC_MODULE(Mon)
{
	sc_in<bool> Clk;
	sc_in<int> a;
	sc_in<int> b;
	sc_in<int> f;

	void show_details();
	SC_CTOR(Mon)
	{
		SC_THREAD(show_details);
		sensitive << Clk.neg();
	}
};

#endif /* MON_H_ */
