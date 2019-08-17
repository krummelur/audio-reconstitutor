package se.fredriksonsound.audioreconstitutor;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;

import java.io.*;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import edu.illinois.cs.cs125.lib.wavfile.WavFile;
import edu.illinois.cs.cs125.lib.wavfile.WavFileException;


public class main {
    static double sampleRate;

    //What algorithm to used to generate the signal
    static final AudioSignal.SignalGenerator generatorAlgorithm = Math::sin;
    //How many signals should be generated for each chunk
    static final int NUM_VOICES = 512;
    //How long to fade in and out each chunk
    static final int FADE_IN_OUT_SAMPLES = 2048 * 12;
    //Should chunks overlap or not, enabling overlap also gives each signal a random phase offset to minimize interference.
    static boolean SHOULD_OVERLAP_CHUNKS = true;
    //How many samples of the input should be analyzed for each chunk (starts from first sample of chunk)
    static final int FFTSIZE = 1024 * 32;
    //How many samples should be in each chunk
    static final int SAMPLES_PER_CHUNK = FFTSIZE / 2;

    public static void main(String[] args) {
        AudioFileHandling af = new AudioFileHandling();
        AudioSignal as = new AudioSignal();
        double[] inputSignal = AudioFileHandling.loadAudioFile("sample.wav");
        int totalOutputLength = (inputSignal.length / FFTSIZE) * FFTSIZE;
        double[] outputSignal = new double[totalOutputLength];
        DoubleFFT_1D fft_1d = new DoubleFFT_1D(FFTSIZE);


        double currentSignalChunk[] = new double[FFTSIZE];
        int offset = 0;
        while (offset + FFTSIZE < totalOutputLength) {

            System.arraycopy(inputSignal, offset, currentSignalChunk, 0, FFTSIZE);

            for (int j = 0; j < currentSignalChunk.length; j++) {
                currentSignalChunk[j] = currentSignalChunk[j] * (0.54f - 0.46f * Math.cos(2 * Math.PI * j / (FFTSIZE - 1)));
            }

            fft_1d.realForward(currentSignalChunk);

            double consolidatedOutput[] = new double[currentSignalChunk.length / 2];
            for (int i = 0; i < currentSignalChunk.length / 2; i += 1) {
                consolidatedOutput[i] = Math.sqrt(currentSignalChunk[i * 2] * currentSignalChunk[i * 2] + currentSignalChunk[i * 2 + 1] * currentSignalChunk[i * 2 + 1]);
            }

            TreeMap<Double, Double> peakMap = new TreeMap<>();
            int adjacentSamples = 5;
            for (int i = adjacentSamples; (i < consolidatedOutput.length - adjacentSamples) && (peakMap.size() < NUM_VOICES); i++) {
                boolean isPeak = true;
                for (int j = 1; j < adjacentSamples; j++) {
                    isPeak &= !(consolidatedOutput[i - j] >= consolidatedOutput[i]);
                    isPeak &= !(consolidatedOutput[i + j] >= consolidatedOutput[i]);
                }
                //Will not catch cases where there is a plateau
                if (isPeak) {
                    peakMap.put(consolidatedOutput[i], sampleRate * i / (currentSignalChunk.length));
                }
            }

            double[] currentChunkOutput = new double[SAMPLES_PER_CHUNK + FADE_IN_OUT_SAMPLES * (SHOULD_OVERLAP_CHUNKS ? 1 : 0)];
            if (offset != 0 && SHOULD_OVERLAP_CHUNKS) {
                System.arraycopy(outputSignal, offset, currentChunkOutput, 0, FADE_IN_OUT_SAMPLES);
            }

            double maxAmplitude = peakMap.lastKey();
            if (!peakMap.isEmpty()) {
                for (Map.Entry<Double, Double> entry : peakMap.entrySet()) {
                    double[] newTone = as.generateAudioWithFrequency(
                            entry.getValue(),
                            entry.getKey() / maxAmplitude / 12,
                            currentChunkOutput.length / sampleRate,
                            sampleRate,
                            FADE_IN_OUT_SAMPLES,
                            generatorAlgorithm);
                    for (int i = 0; i < currentChunkOutput.length; i++) {
                        currentChunkOutput[i] += newTone[i];
                    }
                }
                System.arraycopy(currentChunkOutput, 0, outputSignal, offset, Math.min(currentChunkOutput.length, outputSignal.length - offset));
                offset += SAMPLES_PER_CHUNK;
                System.out.println(String.format("%.2f", (double) offset / totalOutputLength * 100) + " %");
            }
        }
        af.saveAudioFile("output.wav", outputSignal);
        System.exit(0);
    }


}

class AudioFileHandling {
    public static double[] loadAudioFile(String inputFilePath) {
        WavFile wf = null;
        try {
            wf = WavFile.openWavFile(new File(inputFilePath));
            main.sampleRate = wf.getSampleRate();
            int length = (int) wf.getNumFrames();
            double[] audioSamples = new double[length];
            wf.readFrames(audioSamples, length);
            return audioSamples;
        } catch (Exception e) {
            throw new RuntimeException(e);
        } finally {
            if (wf != null) {
                try {
                    wf.close();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
    }

    public static void saveAudioFile(String outputFilePath, double[] outputSignal) {
        WavFile wf = null;
        try {
            wf = WavFile.newWavFile(new File(outputFilePath), 1, (long) outputSignal.length+1, 16, (long) main.sampleRate);
            wf.writeFrames(outputSignal, outputSignal.length);
        } catch (IOException e) {
            throw new RuntimeException(e);
        } catch (WavFileException e) {
            throw new RuntimeException(e);
        }
        finally {
            if (wf != null) {
                try {
                    wf.close();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
    }
}


class AudioSignal {
    Random random = new Random();
    public static final double PI = Math.PI;

    double[] generateAudioWithFrequency(double frequency, double amplitude, double lengthSeconds, double sampleRate,int fadeInOutSamples, AudioSignal.SignalGenerator gen) {
        //double fadeSinDivisor = (FADE_IN_OUT_SAMPLES)/(Math.PI);
        double fadeSinDivisor = (fadeInOutSamples * 2) / (Math.PI);
        double phase = main.SHOULD_OVERLAP_CHUNKS ? random.nextDouble() * Math.PI * 2 : 0;
        double[] signal = new double[(int) (sampleRate * lengthSeconds)];

        for (int i = 0; i < signal.length; i++) {
            signal[i] = gen.generate(phase + ((Math.PI / sampleRate) * (double) i * frequency * 2)) * amplitude;
            if (i < fadeInOutSamples) {
                //signal[i] = signal[i] * (1+Math.sin((double)i/fadeSinDivisor-Math.PI/2))/2;
                signal[i] = signal[i] * (Math.sin((double) i / fadeSinDivisor));
            }
            if (i > signal.length - fadeInOutSamples) {
                //signal[i] = signal[i] * (1+Math.sin(-1.0*((double)(i-signal.length)/fadeSinDivisor+Math.PI/2)))/2;
                signal[i] = signal[i] * (Math.sin(-1.0 * ((double) (i - signal.length) / fadeSinDivisor)));
            }
        }
        return signal;
    }

    public static double square(double t) {
        return (Math.floor(Math.sin(t)) * 2 + 1);
    }

    public static double tri(double t) {
        return -1 * Math.abs(((t + PI / 2) % (2 * PI) - PI) / (PI / 2)) + 1;
    }

    public static double saw(double t) {
        return ((t / PI - 1) % (2) - 1);
    }

    interface SignalGenerator {
        double generate(double t);
    }
}

