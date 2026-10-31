#include <Wire.h>
#include "MAX30105.h"
#include <Adafruit_GFX.h>
#include <Adafruit_SSD1306.h>
#include <math.h>

#define SDA_PIN    21
#define SCL_PIN    22
#define BUTTON_PIN 4

MAX30105 sensor;
Adafruit_SSD1306 oled(128, 64, &Wire, -1);

// ---------------- Sensor settings ----------------
const byte AMP_RED = 60;           // 0-255. Tune so red/IR DC are ~50k-200k
const byte AMP_IR  = 60;
const byte SAMPLE_AVG  = 4;
const byte LED_MODE    = 2;        // 2 = red + IR only
const int  SAMPLE_RATE = 400;      // 400 Hz / 4 averaging = 100 Hz per channel
const int  PULSE_WIDTH = 411;
const int  ADC_RANGE   = 16384;

const int   FS      = 100;
const float DT_MS   = 1000.0f / FS;
const int   WINDOW_SEC = 30;                 // keep <= 45 (RAM)
const int   N          = FS * WINDOW_SEC;
const int   SETTLE_SAMPLES = 2 * FS;         // ignore first 2 s after finger placed
const int   EDGE = FS;                       // trim 1 s at each end after filtering
const long  FINGER_THRESH = 50000;
const int   MAX_BEATS = 120;

// ---------------- Quality-gate thresholds ----------------
const float DC_MIN = 30000, DC_MAX = 250000;
const float PI_MIN_PCT   = 0.05f;   // IR perfusion index, %
const float SNR_MIN      = 3.0f;    // pulse RMS / high-frequency noise RMS
const float CORR_MIN     = 0.85f;   // red vs IR waveform similarity
const float DRIFT_MAX_PCT = 5.0f;   // baseline wander (pressure/motion)
const int   NN_MIN       = 15;      // minimum valid beat-to-beat intervals
const float KEEP_FRAC_MIN = 0.80f;  // fraction of intervals that are plausible
const int   WAVE_MIN     = 8;       // minimum beats usable for morphology
const int   BA_MIN       = 5;

enum { F_DC = 1, F_PERF = 2, F_SNR = 4, F_CORR = 8, F_DRIFT = 16,
       F_BEATS = 32, F_INTERVALS = 64, F_WAVE = 128 };

// ---------------- Hb model (EDIT AFTER CALIBRATION) ----------------
// Hb = baseline(sex) + sum_k WEIGHT[k] * (feature_k - FEAT_MEAN[k]) / FEAT_SD[k]
// Until you fit these from volunteer data, keep MODEL_CALIBRATED = false.
enum { SEX_UNSET = 0, SEX_MALE = 1, SEX_FEMALE = 2 };
const float BASE_HB_MALE   = 14.5f;   // replace with mean lab Hb of YOUR male volunteers
const float BASE_HB_FEMALE = 12.5f;   // replace with mean lab Hb of YOUR female volunteers
const float ANEMIA_CUTOFF_MALE   = 13.0f;   // WHO adult, g/dL
const float ANEMIA_CUTOFF_FEMALE = 12.0f;
const bool  MODEL_CALIBRATED = false;
const float MODEL_MAE = 0.0f;               // cross-validated error (g/dL) from your Python fit

const int K = 12;
const char *FEAT_NAMES[K] = { "dc_ratio", "R_red", "R_ir", "R_ratio", "pi_pct", "hr_bpm",
                              "rmssd_ms", "rise_ms", "width50_ms", "area_norm", "ba_ratio", "delay_ms" };
float FEAT_MEAN[K] = {0,0,0,0,0,0,0,0,0,0,0,0};
float FEAT_SD[K]   = {1,1,1,1,1,1,1,1,1,1,1,1};
float WEIGHT[K]    = {0,0,0,0,0,0,0,0,0,0,0,0};

const bool PRETTY_PRINT = true;
const bool LOG_FAILED_ROWS = false;

// ---------------- Buffers ----------------
uint32_t bufRed[N], bufIR[N];
float fRed[N], fIr[N], scratch[N];
float gPos[MAX_BEATS], gRaw[MAX_BEATS];
bool  gValid[MAX_BEATS];
float gRise[MAX_BEATS], gW50[MAX_BEATS], gArea[MAX_BEATS], gBA[MAX_BEATS];

// ---------------- State ----------------
int      count = 0, settle = 0;
bool     needRemoval = false, waitShown = false;
uint32_t holdUntil = 0, lastButton = 0;
int      subjectSex = SEX_UNSET;
String   subjectId = "none";
float    refHb = NAN;
String   cmdBuf;

struct Result {
  uint16_t failMask;
  bool  ok;
  float redDc, irDc, dcRatio, redAc, irAc, Rred, Rir, Rratio, spo2, pi;
  float hr, sdnn, rmssd;
  int   nBeats, nNN;
  float rise, width50, areaNorm, ba, delayMs;
  float snr, corr, drift;
  float prior, hb;
};

// =====================================================================
// Function declared
// =====================================================================
float medianOf(const float *v, int n) {
  if (n <= 0) return NAN;
  float t[MAX_BEATS];
  if (n > MAX_BEATS) n = MAX_BEATS;
  for (int i = 0; i < n; i++) t[i] = v[i];
  for (int i = 1; i < n; i++) {
    float key = t[i]; int j = i - 1;
    while (j >= 0 && t[j] > key) { t[j + 1] = t[j]; j--; }
    t[j + 1] = key;
  }
  return (n & 1) ? t[n / 2] : 0.5f * (t[n / 2 - 1] + t[n / 2]);
}

char sexChar() { return subjectSex == SEX_MALE ? 'M' : (subjectSex == SEX_FEMALE ? 'F' : '?'); }

// Detrend (1 s centred moving average), invert (more blood = less light, so
// systolic peaks point up), smooth (5-pt). Returns pulse RMS (the "AC" part).
float prepareSignal(const uint32_t *x, float *f, float dc, float &driftPct, float &hfRms) {
  const int HALF = FS / 2;
  float bMin = 1e12f, bMax = -1e12f;
  for (int i = 0; i < N; i++) {
    int lo = i - HALF < 0 ? 0 : i - HALF;
    int hi = i + HALF > N - 1 ? N - 1 : i + HALF;
    double s = 0;
    for (int k = lo; k <= hi; k++) s += x[k];
    float base = (float)(s / (hi - lo + 1));
    scratch[i] = -((float)x[i] - base);
    if (i >= EDGE && i < N - EDGE) {
      if (base < bMin) bMin = base;
      if (base > bMax) bMax = base;
    }
  }
  for (int i = 0; i < N; i++) {
    int lo = i - 2 < 0 ? 0 : i - 2;
    int hi = i + 2 > N - 1 ? N - 1 : i + 2;
    float s = 0;
    for (int k = lo; k <= hi; k++) s += scratch[k];
    f[i] = s / (hi - lo + 1);
  }
  double sq = 0, hf = 0;
  for (int i = EDGE; i < N - EDGE; i++) {
    sq += (double)f[i] * f[i];
    float m = 0;
    for (int k = -2; k <= 2; k++) m += (float)x[i + k];
    m /= 5.0f;
    float d = (float)x[i] - m;
    hf += (double)d * d;
  }
  int n = N - 2 * EDGE;
  hfRms = (float)sqrt(hf / n);
  driftPct = 100.0f * (bMax - bMin) / dc;
  return (float)sqrt(sq / n);
}

// Dominant pulse period (in samples) by autocorrelation over lags 0.3-1.5 s (40-200 bpm).
// The full-beat period correlates better with itself than the half-beat spacing between
// the systolic peak and the smaller diastolic bump, so this does not lock onto the bump.
int estimatePeriod(const float *f) {
  const int LMIN = 30, LMAX = 150;
  double best = -1e30;
  int bestLag = LMIN;
  for (int lag = LMIN; lag <= LMAX; lag++) {
    double s = 0;
    for (int i = EDGE; i < N - EDGE - LMAX; i++) s += (double)f[i] * f[i + lag];
    if (s > best) { best = s; bestLag = lag; }
  }
  return bestLag;
}

// Peak detection on the filtered IR pulse. Window and refractory time scale with the
// estimated period; parabolic interpolation gives sub-sample peak positions.
int findPeaks(const float *f, float rms, float *pos, int period) {
  int W = (int)(0.35f * period);
  if (W < 10) W = 10;
  const int MIN_DIST = (int)(0.6f * period);
  float thr = 0.5f * rms;
  int n = 0, lastIdx = -1000;
  for (int i = EDGE; i < N - EDGE; i++) {
    float v = f[i];
    if (v < thr) continue;
    bool isMax = true;
    for (int k = 1; k <= W; k++) {
      if (f[i - k] > v || f[i + k] >= v) { isMax = false; break; }
    }
    if (!isMax) continue;
    float y0 = f[i - 1], y1 = f[i], y2 = f[i + 1];
    float den = y0 - 2 * y1 + y2;
    float delta = (den != 0) ? 0.5f * (y0 - y2) / den : 0;
    if (delta > 1) delta = 1;
    if (delta < -1) delta = -1;
    float p = i + delta;
    if (i - lastIdx < MIN_DIST) {
      if (n > 0 && v > f[lastIdx]) { pos[n - 1] = p; lastIdx = i; }
      continue;
    }
    if (n < MAX_BEATS) { pos[n++] = p; lastIdx = i; }
  }
  // Prominence filter: drop peaks much smaller than the typical peak (diastolic bumps, noise)
  if (n >= 4) {
    float h[MAX_BEATS];
    for (int k = 0; k < n; k++) h[k] = f[(int)roundf(pos[k])];
    float mh = medianOf(h, n);
    int m = 0;
    for (int k = 0; k < n; k++) if (h[k] >= 0.5f * mh) pos[m++] = pos[k];
    n = m;
  }
  return n;
}

// Red-vs-IR delay by cross-correlation. Positive = IR waveform lags red.
float redIrDelay(const float *fr, const float *fi, float &corrPeak) {
  const int L = 10;
  double c[2 * L + 1], er = 0, ei = 0;
  for (int i = EDGE; i < N - EDGE; i++) { er += (double)fr[i] * fr[i]; ei += (double)fi[i] * fi[i]; }
  int best = 0;
  for (int lag = -L; lag <= L; lag++) {
    double s = 0;
    for (int i = EDGE; i < N - EDGE; i++) s += (double)fr[i] * fi[i + lag];
    c[lag + L] = s;
    if (s > c[best]) best = lag + L;
  }
  corrPeak = (float)(c[best] / sqrt(er * ei + 1e-9));
  float delta = 0;
  if (best > 0 && best < 2 * L) {
    double den = c[best - 1] - 2 * c[best] + c[best + 1];
    if (den != 0) delta = (float)(0.5 * (c[best - 1] - c[best + 1]) / den);
  }
  return ((best - L) + delta) * DT_MS;
}

float computeHb(const float *feat, float &prior) {
  prior = (subjectSex == SEX_MALE) ? BASE_HB_MALE : BASE_HB_FEMALE;
  float hb = prior;
  for (int k = 0; k < K; k++) hb += WEIGHT[k] * (feat[k] - FEAT_MEAN[k]) / FEAT_SD[k];
  if (hb < 5) hb = 5;
  if (hb > 20) hb = 20;
  return hb;
}

// =====================================================================
// Full analysis of one 30 s window
// =====================================================================
Result analyze() {
  Result r;
  memset(&r, 0, sizeof(r));
  r.hr = r.sdnn = r.rmssd = r.rise = r.width50 = r.areaNorm = r.ba = NAN;
  r.prior = r.hb = NAN;

  double sr = 0, si = 0;
  for (int i = 0; i < N; i++) { sr += bufRed[i]; si += bufIR[i]; }
  r.redDc = (float)(sr / N);
  r.irDc  = (float)(si / N);
  r.dcRatio = r.redDc / r.irDc;
  if (r.redDc < DC_MIN || r.redDc > DC_MAX || r.irDc < DC_MIN || r.irDc > DC_MAX) r.failMask |= F_DC;

  float driftR, hfR, hfI;
  r.redAc = prepareSignal(bufRed, fRed, r.redDc, driftR, hfR);
  r.irAc  = prepareSignal(bufIR, fIr, r.irDc, r.drift, hfI);
  r.Rred = r.redAc / r.redDc;
  r.Rir  = r.irAc / r.irDc;
  r.Rratio = r.Rred / r.Rir;
  r.pi = 100.0f * r.Rir;
  r.snr = r.irAc / (hfI + 1e-6f);

  // SpO2: standard empirical curve for R = (ACred/DCred)/(ACir/DCir). Uncalibrated estimate.
  float s = -45.060f * r.Rratio * r.Rratio + 30.354f * r.Rratio + 94.845f;
  r.spo2 = s > 100 ? 100 : (s < 70 ? 70 : s);

  if (r.pi < PI_MIN_PCT) r.failMask |= F_PERF;
  if (r.snr < SNR_MIN) r.failMask |= F_SNR;
  if (r.drift > DRIFT_MAX_PCT || driftR > DRIFT_MAX_PCT) r.failMask |= F_DRIFT;

  r.delayMs = redIrDelay(fRed, fIr, r.corr);
  if (r.corr < CORR_MIN) r.failMask |= F_CORR;

  // ---- Beats, heart rate, HRV ----
  int period = estimatePeriod(fIr);
  float hrAuto = 60.0f * FS / period;
  int nPeaks = findPeaks(fIr, r.irAc, gPos, period);
  r.nBeats = nPeaks;
  if (nPeaks < 4) {
    r.failMask |= F_BEATS | F_INTERVALS | F_WAVE;
  } else {
    int nRaw = nPeaks - 1;
    for (int i = 0; i < nRaw; i++) gRaw[i] = (gPos[i + 1] - gPos[i]) * DT_MS;
    float med = medianOf(gRaw, nRaw);
    double sum = 0;
    int nNN = 0;
    for (int i = 0; i < nRaw; i++) {
      gValid[i] = gRaw[i] >= 300 && gRaw[i] <= 1500 && fabsf(gRaw[i] - med) <= 0.25f * med;
      if (gValid[i]) { sum += gRaw[i]; nNN++; }
    }
    r.nNN = nNN;
    if (nNN >= 3) {
      float mean = (float)(sum / nNN);
      double sq = 0;
      for (int i = 0; i < nRaw; i++) if (gValid[i]) sq += (gRaw[i] - mean) * (gRaw[i] - mean);
      r.sdnn = (float)sqrt(sq / (nNN - 1));
      r.hr = 60000.0f / mean;
      // Cross-check against the autocorrelation rate: a mismatch means peaks were mis-detected
      if (fabsf(r.hr - hrAuto) > 0.2f * hrAuto) r.failMask |= F_INTERVALS;
      double ds = 0;
      int nd = 0;
      for (int i = 0; i < nRaw - 1; i++)
        if (gValid[i] && gValid[i + 1]) { float d = gRaw[i + 1] - gRaw[i]; ds += d * d; nd++; }
      r.rmssd = (nd >= 10) ? (float)sqrt(ds / nd) : NAN;
    }
    if (nNN < NN_MIN || nNN < KEEP_FRAC_MIN * nRaw || isnan(r.rmssd)) r.failMask |= F_INTERVALS;
    if (nPeaks < NN_MIN + 1) r.failMask |= F_BEATS;

    // ---- Waveform morphology on the IR pulse ----
    int nM = 0, nB = 0;
    const int S = 3;
    for (int k = 1; k < nPeaks; k++) {
      if (!gValid[k - 1]) continue;
      int pk = (int)roundf(gPos[k]);
      int back = (int)(0.7f * (gPos[k] - gPos[k - 1]));
      int start = pk - back;
      if (start < S + 1) start = S + 1;
      int foot = start;
      for (int i = start; i < pk; i++) if (fIr[i] < fIr[foot]) foot = i;
      float amp = fIr[pk] - fIr[foot];
      if (amp <= 0) continue;
      float half = fIr[foot] + 0.5f * amp;

      float tUp = NAN, tDn = NAN;
      for (int i = foot + 1; i <= pk; i++)
        if (fIr[i] >= half) { tUp = (i - 1) + (half - fIr[i - 1]) / (fIr[i] - fIr[i - 1]); break; }
      int lim = pk + 80 > N - 2 ? N - 2 : pk + 80;
      for (int i = pk + 1; i <= lim; i++)
        if (fIr[i] <= half) { tDn = (i - 1) + (fIr[i - 1] - half) / (fIr[i - 1] - fIr[i]); break; }
      if (isnan(tUp) || isnan(tDn)) continue;

      int Lb = (int)roundf(gPos[k] - gPos[k - 1]);
      int endI = foot + Lb;
      if (endI > N - 1) continue;
      double area = 0;
      for (int i = foot; i <= endI; i++) { float v = fIr[i] - fIr[foot]; if (v > 0) area += v; }

      if (nM < MAX_BEATS) {
        gRise[nM] = (gPos[k] - foot) * DT_MS;
        gW50[nM]  = (tDn - tUp) * DT_MS;
        gArea[nM] = (float)(area / (amp * (Lb + 1)));
        nM++;
      }

      // Second derivative (acceleration plethysmogram): b/a ratio
      int aIdx = foot;
      float aV = -1e30f;
      for (int i = foot; i <= pk; i++) {
        float d2 = fIr[i + S] - 2 * fIr[i] + fIr[i - S];
        if (d2 > aV) { aV = d2; aIdx = i; }
      }
      int bEnd = aIdx + 25 > N - 1 - S ? N - 1 - S : aIdx + 25;
      float bV = 1e30f;
      for (int i = aIdx; i <= bEnd; i++) {
        float d2 = fIr[i + S] - 2 * fIr[i] + fIr[i - S];
        if (d2 < bV) bV = d2;
      }
      if (aV > 0 && bV < 0 && nB < MAX_BEATS) gBA[nB++] = bV / aV;
    }
    if (nM >= WAVE_MIN) {
      r.rise = medianOf(gRise, nM);
      r.width50 = medianOf(gW50, nM);
      r.areaNorm = medianOf(gArea, nM);
    } else r.failMask |= F_WAVE;
    if (nB >= BA_MIN) r.ba = medianOf(gBA, nB);
    else r.failMask |= F_WAVE;
  }

  r.ok = (r.failMask == 0);

  // ---- Hb estimate (sex-specific baseline + optical corrections) ----
  if (subjectSex != SEX_UNSET) {
    float feat[K] = { r.dcRatio, r.Rred, r.Rir, r.Rratio, r.pi, r.hr,
                      r.rmssd, r.rise, r.width50, r.areaNorm, r.ba, r.delayMs };
    bool finite = true;
    for (int k = 0; k < K; k++) if (isnan(feat[k])) finite = false;
    float prior;
    float hb = computeHb(feat, prior);
    r.prior = prior;
    if (MODEL_CALIBRATED && r.ok && finite) r.hb = hb;
  }
  return r;
}
