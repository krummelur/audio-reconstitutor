# audio-reconstitutor
An dft based audio effect.

Takes an input file (.wav, mono) splits into chunks and performs dft on it, then generates signals for the most prominent frequencies for each chunk and writes it to an output file.  
the result will change depending on the length of audio to analyze for each chunk, the chunk length, the chunk blending, and the signal shape.
