# Turbo Codes Simulator
I wrote this little simulator during a project for the class of Information Theory and Codes during my graduate studies. The code is written in C99 and it is self-contained.
The main purpose was to demonstrate the incredible performances of turbo-codes (near Shannon Capacity) instead of obtaining a fast encoding\decoding.

## Important Notes
This simulator has been written in two weeks so it is not well documented nor optimized in terms of efficiency and the code might seem a bit dirty. Anyway, there are plenty of possibilities to extend code in a very clean manner.

## Files description
<b>./conv_encoder.h</b><br>
&emsp;This library implements a SINGLE convolutional encoder with or without feedback. 

<b>./shiftRegister.h</b><br>
%emsp;Internally used by conv_encoder to provide an implementation of a a linear shift register block.

<b>./graph.h</b><br>
&emsp;Implements the trellis structures to describe the convolutional code. It also includes a BCJR implementation (In earlier version a training code for the Viterbi Algorithm was also implemented)

<b>./graph_turbo.h</b><br>
&emsp;Implements the BCJR algorithm specifically modified to work with turbo decoder. It is assumed that the SNR is exactly known.

<b>./main.c</b><br>
&emsp;Main file with some routines to test performances. This were used to build the Eb/N0 vs. BER plot and to test performance in sending binary files. Be aware that this tests could require several minutes in order to complete. Catching BER <1E-5 could be very challenging!

<b>./max_lookup.h</b><br>
&emsp;Lookup table for MAX*

<b>./mt/*</b><br>
&emsp;Marsenne Twister used to generate random numbers in order to perform a better monte-carlo simulation. (Needed to catch BER < 1E-5)

<b>./scramblers/*</b><br>
&emsp;CCSDS scramblers

