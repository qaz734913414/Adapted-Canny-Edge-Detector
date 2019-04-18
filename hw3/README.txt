Description: 
	The Stim module sense the clock signal, when clock is high, it generates two inputs a and b, and sends them to the mult and mon module. The mult module sense the change of a and b, if change happened, the action function compute the multiply of a and b, and outputs the result f to mon module. The mon module will sense the falling edge of the clock, when the mult module finish multiplying, the show_details function will print the result and the corresponding input

Problem I met:
	In the stimulate function, I forgot to add a wait() function before sc_stop(), which result in an absence of execution of the last set of input. Also, I did not change the uppercase variables to match the rest of variables.