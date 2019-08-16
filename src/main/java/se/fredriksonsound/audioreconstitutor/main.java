package se.fredriksonsound.audioreconstitutor;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import edu.illinois.cs.cs125.lib.wavfile.WavFile;
import edu.illinois.cs.cs125.lib.wavfile.WavFileException;


public class main {
    static double sampleRate = 44100.0;
    static int numVoices = 500;
    static int fadeInOutSamples = 2048;

    public static void main(String[] args) {
        AudioStuff as = new AudioStuff();
        double[]signal = loadAudioFile("sample.wav");
        int FFTSIZE = 1024*32;
        int samplesPerChunk = FFTSIZE/8;
        int totalOutputLength = (signal.length / FFTSIZE) * FFTSIZE;
        double[] outputSignal = new double[totalOutputLength];
        DoubleFFT_1D fft_1d = new DoubleFFT_1D(FFTSIZE);
        int offset = 0;
        double currentSignalChunk[] = new double[FFTSIZE];


        while (offset + FFTSIZE < totalOutputLength) {

            System.arraycopy(signal, offset, currentSignalChunk, 0, FFTSIZE);

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
            for (int i = adjacentSamples; i < consolidatedOutput.length - adjacentSamples && peakMap.size() < numVoices; i++) {
                boolean isPeak = true;
                for (int j = 1; j < adjacentSamples; j++) {
                    if (consolidatedOutput[i - j] >= consolidatedOutput[i])
                        isPeak = false;
                    if (consolidatedOutput[i + j] >= consolidatedOutput[i])
                        isPeak = false;
                }
                //Will not catch edge cases where there is a plateau
                if (isPeak && (consolidatedOutput[i] > 30 || true)) {
                    peakMap.put(consolidatedOutput[i], sampleRate * i / (currentSignalChunk.length));
                }
            }
            //peakMap.forEach();

            double[] currentChunkOutput = new double[samplesPerChunk + fadeInOutSamples];
            if(offset != 0) {
                System.arraycopy(outputSignal, offset, currentChunkOutput, 0, fadeInOutSamples);
            }
            System.out.println( String.format("%.2f", (double)offset / totalOutputLength * 100) + " %");
            double maxAmplitude = peakMap.lastKey();
            if (!peakMap.isEmpty()) {
                int currentEntry = 0;
                for (Map.Entry<Double, Double> entry : peakMap.entrySet()) {
                    double[] newTone = as.createSampleWithFrequency(
                            entry.getValue(),
                            sampleRate,
                            entry.getKey() / maxAmplitude / 12,
                            currentChunkOutput.length / sampleRate,
                            fadeInOutSamples,
                            t -> as.saw(t));
                    for (int i = 0; i < currentChunkOutput.length; i++) {
                        currentChunkOutput[i] += newTone[i];
                    }
                }
                System.arraycopy(currentChunkOutput, 0, outputSignal, offset, Math.min(currentChunkOutput.length, outputSignal.length-offset));
                offset += samplesPerChunk;
            }
        }
        WavFile wf = null;
        try {
            wf = WavFile.newWavFile(new File("output.wav"), 1, (long) outputSignal.length, 16, (long) sampleRate);
            wf.writeFrames(outputSignal, outputSignal.length);
            wf.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        } catch (WavFileException e) {
            throw new RuntimeException(e);
        }

        System.exit(0);
    }

    public static double[] loadAudioFile(String inputFilePath) {
        WavFile wf = null;
        try {
            wf = WavFile.openWavFile(new File(inputFilePath));
            sampleRate = wf.getSampleRate();
            int length = (int) wf.getNumFrames();
            double[] audioSamples = new double[length];
            System.out.println("VALIDBITS: " + wf.getValidBits());
            wf.readFrames(audioSamples, length);
            wf.close();
            return audioSamples;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}

class AudioStuff {
    public static final double PI = Math.PI;
    Random random = new Random();

    double[] createSampleWithFrequency(double frequency, double sampleRate, double amplitude, double lengthSeconds, int fadeInOutSamples, SignalGenerator gen) {
        //frequency = 220;
        //double fadeSinDivisor = (fadeInOutSamples)/(Math.PI);
        double fadeSinDivisor = (fadeInOutSamples*2)/(Math.PI);
        double phase = random.nextDouble()*Math.PI*2;
        //phase = 0;
        double[] signal = new double[(int) (sampleRate * lengthSeconds)];
        for (int i = 0; i < signal.length; i++) {

            signal[i] = gen.generate(phase + ((Math.PI / sampleRate) * (double)i * frequency * 2)) * amplitude;
            //signal[i] = 0.2;
            if(i < fadeInOutSamples) {
                //signal[i] = signal[i] * (1+Math.sin((double)i/fadeSinDivisor-Math.PI/2))/2;
                signal[i] = signal[i] * (Math.sin((double)i/fadeSinDivisor));
            }
            if(i > signal.length - fadeInOutSamples) {
                //signal[i] = signal[i] * (1+Math.sin(-1.0*((double)(i-signal.length)/fadeSinDivisor+Math.PI/2)))/2;
                signal[i] = signal[i] * (Math.sin(-1.0*((double)(i-signal.length)/fadeSinDivisor)));
            }
        }
        return signal;
    }

    public double square(double t) {
        return (Math.floor(Math.sin(t)) * 2 + 1);
    }

    public double tri(double t) {
        return -1*Math.abs( ((t+PI/2)%(2*PI)-PI)/(PI/2))+1;
    }

    public double saw(double t) {
        return ((t/PI-1)%(2)-1);
    }

    interface SignalGenerator {
        double generate(double t);
    }
}

