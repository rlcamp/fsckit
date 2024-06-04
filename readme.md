# fsckit/unfsckit

This is a work-in-progress, proof-of-concept digital modulation loosely inspired by LoRa[1], which has some properties desirable for acoustic comms over unreliable channels. Care will need to be taken to ensure that we do not violate the LoRa patent if using anything similar to this for anything important. There is lots of relevant stuff in the academic literature describing how to time-align LoRa. I have not read the LoRa patent itself.

## Concepts

### Frequency-shifted chirp keying

The basic concept is identical to multiple-frequency-shift-keying (MFSK) with an additional linearly-increasing-with-wraparound component to the transmitted waveform. Each message starts with a preamble consisting of 2 or more upsweeps followed by 2 downsweeps, allowing the receiver to detect the beginning of a message and simulateously disambiguate timing and carrier frequency offset[2].

Adding the chirped component to MFSK modulation also ensures that all symbols equally sweep the entire passband, providing resilience against non-flat responses in the transmitter, propagation channel, and/or receiver.

### MFSK detection

Once thus synced, the receiver can simply de-chirp the received signal by multiplication with a conjugate of the expected unmodulated chirp sequence. The output of this de-chirping is simple MFSK modulation, which can be optimally[1] decoded using block-wise fast Fourier transforms on a critically sampled subband - that is, each symbol consists of S complex samples, over which a length-S complex-to-complex FFT is performed. The loudest of the resulting S bins is the index of the encoded data symbol.

Suboptimal demodulation in the high-SNR limit can be performed without an FFT by doing simple FM demodulation of the dechirped waveform within each symbol period.

### Spreading factor

This modulation falls under the category of "spread spectrum" techniques in that there is a knob that controls spectral redundancy. Since each symbol encodes one of S unique values (and therefore, N = log2(S) bits) and has a duration linearly proportional to S, the redundancy factor is S / log2(S), or equivalently, 2^N / N. Values of N of interest to us range from 4 to 6 (compare to LoRa which uses values of N from 7 to 12 inclusive, in 125-500 kHz of bandwidth around 1 GHz).

### Todo

- Better false alarm mitigation. This will require both better resistance to incorrectly leaving the detecting-a-preamble state, as well as a way to quickly return to it rather than demodulating a long data packet that does not exist. A checksum on the packet length might be most surefire way to achieve this, at the expense of adding overhead to every packet

- Improve preamble detection at low SNR

- Determine the minimum necessary amount of bandpass filtering in actual practice - the eight-pole Butterworth bandpass filter downstream of the high pass filter may be overkill

- Identify and eliminate slow libc calls. The code has a number of remainderf() and lrintf() calls within the hot loop which can be a bottleneck or not, depending on how well they are implemented and optimized for the finite-math-only case within various libc's we care about.

- Figure out what is covered by the patent and make sure we are in the clear

- Add some forward error correction. LoRa uses one of several very simple Hamming codes. We can probably do better (and may be forced to in order to not step on the patent) but bear in mind this needs to work on a microcontroller with almost no available memory. Techniques which rely on a Fourier transform as the primitive would be convenient, as the demodulator already uses one.

### References

[1] L. Vangelista, "Frequency shift chirp modulation: The LoRa modulation," IEEE Signal Processing Letters, vol. 24, no. 12, pp. 1818–1821, Dec 2017.

[2] C. Bernier, F. Dehmas and N. Deparis, "Low Complexity LoRa Frame Synchronization for Ultra-Low Power Software-Defined Radios," in IEEE Transactions on Communications, vol. 68, no. 5, pp. 3140-3152, May 2020, doi: 10.1109/TCOMM.2020.2974464 https://cea.hal.science/cea-02280910/file/TCOM.pdf
