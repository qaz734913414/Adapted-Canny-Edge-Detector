/*
 * Stim.h
 *
 *  Created on: Oct 23, 2018
 *      Author: xuyihan
 */

#ifndef STIM_H_
#define STIM_H_
#include "systemc.h"
SC_MODULE(Stim)
{
	sc_in<bool> Clk;
	sc_out<int> a;
	sc_out<int> b;

	void stimulus();
	SC_CTOR(Stim)
	{
		SC_THREAD(stimulus);
		sensitive << Clk.pos();
	}
};

#endif /* STIM_H_ */
