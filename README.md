# Do-it-yourself EEG Code

Written in Processing, this code aims to give you really easy-to-use functions
that help process the data from an EEG.

# EEG Basic

This is a simple visualizer that shows the raw input waveform,
the FFT of that raw input, a longer-time sample buffer with a zoom in of the
lower frequencies, and an exponentially decaying average of the EEG spectrum
in certain critical EEG bands.  It assumes the data is being fed in from the
microphone slot on the computer.