/*
 * Top.h
 *
 *  Created on: Oct 23, 2018
 *      Author: xuyihan
 */

#ifndef TOP_H_
#define TOP_H_

#include "systemc.h"
#include "mult.h"
#include "mon.h"
#include "stim.h"

SC_MODULE(Top)
{
	sc_signal<int> asig, bsig, fsig;
	sc_clock testclk;
	Stim stim1;
	Mult uut;
	Mon mon1;

	SC_CTOR(Top)
	: testclk("testclk", 10, SC_NS),
	stim1("stim1"),
	uut("uut"),
	mon1("mon1")
	{
		stim1.a(asig);
		stim1.b(bsig);
		stim1.Clk(testclk);
		uut.a(asig);
		uut.b(bsig);
		uut.f(fsig);
		mon1.a.bind(asig);
		mon1.b.bind(bsig);
		mon1.f.bind(fsig);
		mon1.Clk.bind(testclk);
	}
};
#endif /* TOP_H_ */
