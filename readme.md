# Chirp FSK modulator/demodulator

This is a proof-of-concept digital modulation loosely inspired by LoRa[1], which has some properties desirable for acoustic communications over unreliable channels.

## Concepts

### Frequency-shifted chirp keying

The basic concept is identical to multiple-frequency-shift-keying (MFSK) with an additional linearly-increasing-with-wraparound component to the transmitted waveform. Each message starts with a preamble consisting of 3 or more upsweeps followed by 2 downsweeps, allowing the receiver to detect the beginning of a message and simulateously disambiguate timing and carrier frequency offset[2].

Adding the chirped component to MFSK modulation also ensures that all symbols equally sweep the entire passband, providing resilience against non-flat responses in the transmitter, propagation channel, and/or receiver.

### MFSK detection

Once synced, the receiver can simply de-chirp the received signal by multiplying by the conjugate of the expected unmodulated chirp. The output of this de-chirping is simple MFSK modulation, which can be optimally[1] decoded using block-wise fast Fourier transforms on a critically sampled subband (that is, each symbol consists of S complex samples, over which a length-S complex-to-complex FFT is performed). The loudest of the resulting S bins is the index of the encoded data symbol. Gray coding is used to ensure that single-shift errors result in only single-bit errors.

If the next-downstream logic is a soft-decision forward error correcting decoder, the individual bits can be soft-decided as follows. The dechirping and FFT are performed as before. For each output bit, a power-if-one and a power-if-zero accumulator are initialized. The magnitude squared of the FFT bins of all shifts are inspected. For each bin whose Gray-coded index is set, that bin's magnitude squared is added to the power-if-one accumulator, otherwise it is added to the power-if-zero accumulator. After summing the power in all bins into one or the other of these accumulators, a soft decision between +1.0 and -1.0 can be obtained by dividing the difference of these accumulators by their sum. The details of the normalizations of this method have not yet been proven to be optimal.

### Spreading factor

This modulation falls under the category of "spread spectrum" techniques in that there is a knob that controls spectral redundancy. Since each symbol encodes one of S unique values (and therefore, N = log2(S) bits) and has a duration linearly proportional to S, the redundancy factor is S / log2(S), or equivalently, 2^N / N. Values of N of interest to us range from 3 to 6 (compare to LoRa which uses values of N from 7 to 12 inclusive, in 125-500 kHz of bandwidth around 1 GHz).

### Forward error correction

The current implementation uses a Hamming 7,4 code with a controllable amount of interleaving between it and the chirp layer. Without interleaving, the Hamming layer would not actually add any value vs simply using a larger spreading factor in the chirp layer. The interleaving factor should be matched to the expected channel conditions (that is, the effective time between each of the 7 bits should be greater than the expected duration of an interference burst).

After deinterleaving, the receiver does a dumb simple exhaustive naive maximum-likelihood brute force search for the most likely uncorrupted code word using soft bit decisions as inputs. We can probably do better, but bear in mind this needs to work on a microcontroller with almost no available memory.

### Todo

- Pin down the required bandpass filter parameters, there are three knobs here that interact

- Identify and eliminate slow libc calls. The code has a number of remainderf() and lrintf() calls within the hot loop which can be a bottleneck or not, depending on how well they are implemented and optimized for the finite-math-only case within various libc's we care about.

- Better forward error correction. We're currently brute-force soft decoding a Hamming 7,4 layer in exponential time. This could be replaced with a Hadamard 8,4 layer which can be soft-decoded for a lot less compute effort

### License

ISC license.

### References

[1] L. Vangelista, "Frequency shift chirp modulation: The LoRa modulation," IEEE Signal Processing Letters, vol. 24, no. 12, pp. 1818–1821, Dec 2017.

[2] C. Bernier, F. Dehmas and N. Deparis, "Low Complexity LoRa Frame Synchronization for Ultra-Low Power Software-Defined Radios," in IEEE Transactions on Communications, vol. 68, no. 5, pp. 3140-3152, May 2020, doi: 10.1109/TCOMM.2020.2974464 https://cea.hal.science/cea-02280910/file/TCOM.pdf

## Notes for maintainers

### `dprintf`

Debug text is emitted to `stderr` via `dprintf(2, ...)` instead of via `fprintf(stderr, ...)` to allow embedded usage without incurring any `malloc`/`free` activity within the loop.
