# Frequency-estimation-in-audio-signal-using-MATLAB

MATLAB script for detecting speech formants using three AR methods: Yule–Walker, Levinson–Durbin, and Burg. Supports file input or live recording, applies pre-emphasis, estimates AR models, plots spectra, and marks up to seven formants. Useful for speech and DSP analysis.

## Features

- **Two input modes** – load an existing audio file (WAV, MP3, FLAC, OGG) or capture audio directly from a microphone.
- **Pre-emphasis filter** – compensates for spectral tilt inherent in speech signals.
- **Three AR estimation methods**:
  - *Yule–Walker* – solves the autocorrelation normal equations.
  - *Levinson–Durbin* – efficient recursive solution of the same equations.
  - *Burg* – minimises forward and backward prediction errors, yielding stable models with improved frequency resolution.
- **Formant detection** – identifies up to 7 spectral peaks in the 80–4500 Hz range (covers F1–F4 for typical speech).
- **Side-by-side plots** – a three-panel figure displays the normalised AR spectrum for each method with formant frequencies annotated.

## Requirements

| Dependency | Version |
|---|---|
| MATLAB | R2018b or newer |
| Signal Processing Toolbox | any compatible version |

## Usage

1. Open MATLAB and navigate to the directory containing `frequency_estimation.m`.
2. Run the function:

   ```matlab
   frequency_estimation()
   ```

3. When prompted, enter:
   - `1` to load an audio file via a file-browser dialog, **or**
   - `2` to record audio from the default microphone (you will also be asked for the recording duration in seconds).

4. A figure window will appear showing the AR spectrum and detected formants for all three methods.

## Algorithm Overview

```
Audio input
    │
    ▼
Mono conversion + normalisation
    │
    ▼
Pre-emphasis filter  H(z) = 1 − 0.97 z⁻¹
    │
    ▼
AR coefficient estimation (order p = 16)
    ├── Yule–Walker  (aryule)
    ├── Levinson–Durbin  (levinson)
    └── Burg  (arburg)
    │
    ▼
All-pole frequency response  |H(e^jω)|  (freqz, N = 2048 points)
    │
    ▼
Peak picking in 80–4500 Hz  →  formant frequencies
    │
    ▼
Annotated spectrum plots (3 subplots)
```

## File Structure

```
frequency_estimation.m   Main script and helper function (findFormants)
README.md                This file
```
