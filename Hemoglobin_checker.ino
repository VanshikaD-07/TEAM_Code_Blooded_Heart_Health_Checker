#include <Wire.h>
#include "MAX30105.h"
#include <Adafruit_GFX.h>
#include <Adafruit_SSD1306.h>
#include <math.h>

#define SDA_PIN    21
#define SCL_PIN    22
#define BUTTON_PIN 4
#define CSV_PATH   "/hb_data.csv"

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

// =====================================================================
// Output
// =====================================================================

void printFailReasons(uint16_t m) {
  if (m & F_DC)        Serial.println("  - DC level out of range: adjust AMP_RED/AMP_IR or finger contact");
  if (m & F_PERF)      Serial.println("  - Pulse too weak (perfusion index low): warm finger, lighter pressure");
  if (m & F_SNR)       Serial.println("  - Too noisy: block room light, keep finger still");
  if (m & F_CORR)      Serial.println("  - Red and IR waveforms disagree: motion or poor contact");
  if (m & F_DRIFT)     Serial.println("  - Baseline drifting: pressure changed or finger moved");
  if (m & F_BEATS)     Serial.println("  - Too few heartbeats detected");
  if (m & F_INTERVALS) Serial.println("  - Beat intervals irregular or too few valid for HRV");
  if (m & F_WAVE)      Serial.println("  - Not enough clean beats for waveform features");
}

const char *CSV_HEADER =
    "id,sex,ref_hb,ts_ms,quality_ok,fail_mask,hr_bpm,rmssd_ms,sdnn_ms,spo2,pi_pct,"
    "red_dc,ir_dc,dc_ratio,R_red,R_ir,R_ratio,rise_ms,width50_ms,area_norm,ba_ratio,"
    "delay_ms,snr,corr,drift_pct,n_beats,n_nn,hb_est";

bool fsReady = false;

void addNum(String &s, float v, int dec, bool last = false) {
  if (!isnan(v)) s += String(v, dec);
  if (!last) s += ',';
}

String csvRow(const Result &r) {
  String s = subjectId;
  s += ',';
  s += sexChar(); s += ',';
  addNum(s, refHb, 1);
  s += String(millis()); s += ',';
  s += String(r.ok ? 1 : 0); s += ',';
  s += String(r.failMask); s += ',';
  addNum(s, r.hr, 1); addNum(s, r.rmssd, 1); addNum(s, r.sdnn, 1); addNum(s, r.spo2, 1); addNum(s, r.pi, 3);
  addNum(s, r.redDc, 0); addNum(s, r.irDc, 0); addNum(s, r.dcRatio, 4);
  addNum(s, r.Rred, 5); addNum(s, r.Rir, 5); addNum(s, r.Rratio, 4);
  addNum(s, r.rise, 1); addNum(s, r.width50, 1); addNum(s, r.areaNorm, 4); addNum(s, r.ba, 4);
  addNum(s, r.delayMs, 1); addNum(s, r.snr, 1); addNum(s, r.corr, 3); addNum(s, r.drift, 2);
  s += String(r.nBeats); s += ',';
  s += String(r.nNN); s += ',';
  addNum(s, r.hb, 2, true);
  return s;
}

// Append one row to the ESP32's flash file (header is written when the file is created)
void saveRow(const String &row) {
  if (!fsReady) { Serial.println("(flash not available: row not saved on board)"); return; }
  bool isNew = !LittleFS.exists(CSV_PATH);
  File f = LittleFS.open(CSV_PATH, FILE_APPEND);
  if (!f) { Serial.println("(flash write failed)"); return; }
  if (isNew) f.println(CSV_HEADER);
  f.println(row);
  f.close();
  Serial.println("(saved to ESP32 flash: type 'dump' to print the file)");
}

// Print the whole saved file to the Serial Monitor (or to log_serial.py, which saves it on your laptop)
void dumpFile() {
  if (!fsReady || !LittleFS.exists(CSV_PATH)) { Serial.println("No saved data yet."); return; }
  File f = LittleFS.open(CSV_PATH, FILE_READ);
  Serial.println("----- BEGIN hb_data.csv -----");
  while (f.available()) Serial.write(f.read());
  Serial.println("----- END hb_data.csv -----");
  f.close();
}

void clearFile() {
  if (fsReady && LittleFS.exists(CSV_PATH)) LittleFS.remove(CSV_PATH);
  Serial.println("Saved data cleared.");
}

void spaceInfo() {
  if (!fsReady) { Serial.println("Flash not available."); return; }
  Serial.printf("Flash file system: %u of %u bytes used\n",
                (unsigned)LittleFS.usedBytes(), (unsigned)LittleFS.totalBytes());
}


void printCsvHeader() { Serial.println(CSV_HEADER); }

void printReport(const Result &r) {
  if (PRETTY_PRINT) {
    Serial.println();
    Serial.println("========== MEASUREMENT COMPLETE ==========");
    Serial.printf("Subject: %s   Sex: %c   Reference Hb: ", subjectId.c_str(), sexChar());
    if (isnan(refHb)) Serial.println("not set"); else Serial.printf("%.1f g/dL\n", refHb);
    Serial.printf("QUALITY: %s\n", r.ok ? "PASS" : "FAIL - retry");
    if (!r.ok) printFailReasons(r.failMask);

    Serial.println("\n-- Heart (from IR pulse peaks) --");
    Serial.printf("Heart rate            : %.1f bpm\n", r.hr);
    Serial.printf("HRV RMSSD             : %.1f ms  (beat-to-beat variability, 30 s window)\n", r.rmssd);
    Serial.printf("HRV SDNN              : %.1f ms  (overall variability)\n", r.sdnn);
    Serial.printf("Beats / valid intervals: %d / %d\n", r.nBeats, r.nNN);

    Serial.println("\n-- Oxygen and perfusion --");
    Serial.printf("SpO2 (approx, uncal.) : %.1f %%\n", r.spo2);
    Serial.printf("Perfusion index (IR)  : %.3f %%  (pulse strength)\n", r.pi);

    Serial.println("\n-- Light levels --");
    Serial.printf("Red DC / IR DC        : %.0f / %.0f counts\n", r.redDc, r.irDc);
    Serial.printf("DC ratio (red/IR)     : %.4f\n", r.dcRatio);
    Serial.printf("AC/DC red, IR         : %.5f , %.5f\n", r.Rred, r.Rir);
    Serial.printf("AC/DC ratio of ratios : %.4f  (R used for SpO2)\n", r.Rratio);

    Serial.println("\n-- Pulse waveform (median over beats) --");
    Serial.printf("Rise time             : %.1f ms  (foot to peak)\n", r.rise);
    Serial.printf("Width at 50%% height   : %.1f ms\n", r.width50);
    Serial.printf("Area index            : %.4f  (pulse fullness, 0-1)\n", r.areaNorm);
    Serial.printf("b/a ratio             : %.4f  (2nd-derivative wave ratio)\n", r.ba);

    Serial.println("\n-- Red vs IR timing --");
    Serial.printf("IR delay vs red       : %.1f ms  (experimental; 10 ms sampling steps)\n", r.delayMs);

    Serial.println("\n-- Signal quality --");
    Serial.printf("Pulse SNR index       : %.1f   Red/IR correlation: %.3f   Baseline drift: %.2f %%\n",
                  r.snr, r.corr, r.drift);

    Serial.println("\n-- Hemoglobin --");
    if (subjectSex == SEX_UNSET) {
      Serial.println("Set sex first (type m or f) to get the sex-specific baseline.");
    } else if (!MODEL_CALIBRATED) {
      Serial.printf("NOT CALIBRATED. Sex-mean prior for %c is %.1f g/dL: a population average,\n"
                    "not a measurement. Fit the model on volunteer data, then set WEIGHT[],\n"
                    "FEAT_MEAN[], FEAT_SD[], BASE_HB_* and MODEL_CALIBRATED = true.\n",
                    sexChar(), r.prior);
    } else if (isnan(r.hb)) {
      Serial.println("No estimate: quality check failed or a feature is missing.");
    } else {
      float cut = (subjectSex == SEX_MALE) ? ANEMIA_CUTOFF_MALE : ANEMIA_CUTOFF_FEMALE;
      Serial.printf("Hb estimate: %.1f g/dL (+/- %.1f)   %s (cutoff %.1f)\n", r.hb, MODEL_MAE,
                    r.hb < cut ? "LOW-HB RISK" : "within normal range", cut);
    }
    Serial.println("Screening aid only, not a diagnosis.");
    Serial.println("\nCSV row:");
  }
  if (r.ok || LOG_FAILED_ROWS) { String row = csvRow(r); Serial.println(row); saveRow(row); }
  else if (PRETTY_PRINT) Serial.println("(row not logged: quality failed)");
}

void showWait() {
  oled.clearDisplay();
  oled.setTextSize(1);
  oled.setTextColor(SSD1306_WHITE);
  oled.setCursor(0, 0);
  oled.println("Place finger on");
  oled.println("sensor, hold still");
  oled.printf("\nSex: %c   ID: %s\n", sexChar(), subjectId.c_str());
  if (subjectSex == SEX_UNSET) oled.println("Type m/f or press btn");
  oled.display();
  waitShown = true;
}

void showProgress(int sec) {
  oled.clearDisplay();
  oled.setCursor(0, 0);
  oled.println("Measuring...");
  oled.printf("Hold still %d/%d s\n", sec, WINDOW_SEC);
  oled.display();
  waitShown = false;
}

void showResult(const Result &r) {
  oled.clearDisplay();
  oled.setCursor(0, 0);
  if (!r.ok) {
    oled.println("QUALITY: RETRY");
    if (r.failMask & F_DC)        oled.println("Check LED/contact");
    if (r.failMask & F_PERF)      oled.println("Pulse too weak");
    if (r.failMask & (F_SNR | F_CORR)) oled.println("Noisy signal");
    if (r.failMask & F_DRIFT)     oled.println("Finger moved");
    if (r.failMask & (F_BEATS | F_INTERVALS | F_WAVE)) oled.println("Too few clean beats");
  } else {
    oled.printf("HR %.0f  SpO2 %.0f%%\n", r.hr, r.spo2);
    oled.printf("RMSSD %.0f ms\n", r.rmssd);
    oled.printf("PI %.2f%%  Sex %c\n", r.pi, sexChar());
    if (subjectSex == SEX_UNSET) oled.println("Set sex: m/f");
    else if (!MODEL_CALIBRATED) oled.println("Hb: needs calibration");
    else if (isnan(r.hb)) oled.println("Hb: n/a");
    else {
      float cut = (subjectSex == SEX_MALE) ? ANEMIA_CUTOFF_MALE : ANEMIA_CUTOFF_FEMALE;
      oled.printf("Hb %.1f g/dL\n%s\n", r.hb, r.hb < cut ? "LOW-HB RISK" : "Normal range");
    }
  }
  oled.println("Remove finger");
  oled.display();
  waitShown = false;
}

// =====================================================================
// Input handling
// =====================================================================
void setSex(int s) {
  subjectSex = s;
  Serial.printf("Sex set to %c\n", sexChar());
  showWait();
}

void processCmd(String s) {
  s.trim();
  if (s.length() == 0) return;
  char c = tolower(s[0]);
  if (c == 'm') setSex(SEX_MALE);
  else if (c == 'f') setSex(SEX_FEMALE);
  else if (c == 'r') { refHb = s.substring(1).toFloat(); Serial.printf("Reference Hb set to %.1f g/dL\n", refHb); }
  else if (c == 'i') { subjectId = s.substring(1); subjectId.trim(); refHb = NAN;
                       Serial.printf("Subject ID: %s (reference Hb cleared)\n", subjectId.c_str()); showWait(); }
  else if (c == 'x') { refHb = NAN; Serial.println("Reference Hb cleared"); }
  else Serial.println("Commands: m | f | r<Hb> | i<ID> | x");
}

void handleSerial() {
  while (Serial.available()) {
    char c = Serial.read();
    if (c == '\n' || c == '\r') { if (cmdBuf.length()) { processCmd(cmdBuf); cmdBuf = ""; } }
    else if (cmdBuf.length() < 24) cmdBuf += c;
  }
}

void handleButton() {
  if (digitalRead(BUTTON_PIN) == LOW && millis() - lastButton > 300) {
    lastButton = millis();
    setSex(subjectSex == SEX_MALE ? SEX_FEMALE : SEX_MALE);
  }
}

// =====================================================================
void setup() {
  Serial.begin(115200);
  pinMode(BUTTON_PIN, INPUT_PULLUP);
  Wire.begin(SDA_PIN, SCL_PIN);
  Wire.setClock(400000);          // use 100000 if you see I2C errors

  if (!oled.begin(SSD1306_SWITCHCAPVCC, 0x3C)) Serial.println("OLED not found");

  if (!sensor.begin(Wire, I2C_SPEED_FAST)) {
    Serial.println("MAX3010x not found - check wiring");
    oled.clearDisplay(); oled.setCursor(0, 0); oled.println("Sensor not found"); oled.display();
    while (1) delay(1000);
  }
  sensor.setup(60, SAMPLE_AVG, LED_MODE, SAMPLE_RATE, PULSE_WIDTH, ADC_RANGE);
  sensor.setPulseAmplitudeRed(AMP_RED);
  sensor.setPulseAmplitudeIR(AMP_IR);

  fsReady = LittleFS.begin(true);   // true = format on first use
  if (!fsReady) Serial.println("LittleFS mount failed: rows will only print to serial");
  Serial.println("Hb screening prototype ready. Commands: m | f | r<Hb> | i<ID> | x | dump | clear | space");
  printCsvHeader();
  showWait();
}

void loop() {
  handleSerial();
  handleButton();
  sensor.check();

  while (sensor.available()) {
    uint32_t r  = sensor.getFIFORed();
    uint32_t ir = sensor.getFIFOIR();
    sensor.nextSample();

    if (ir < FINGER_THRESH) {            // no finger
      if (needRemoval) needRemoval = false;
      if (count > 0 || settle > 0) { count = 0; settle = 0; waitShown = false; }
      if (!waitShown && millis() > holdUntil) showWait();
      continue;
    }
    if (needRemoval) continue;           // result is showing; wait for finger removal
    if (settle < SETTLE_SAMPLES) { settle++; continue; }

    bufRed[count] = r; bufIR[count] = ir; count++;
    if (count % FS == 0 && count < N) showProgress(count / FS);

    if (count >= N) {
      showProgress(WINDOW_SEC);
      Result res = analyze();
      printReport(res);
      showResult(res);
      count = 0; settle = 0;
      needRemoval = true;
      holdUntil = millis() + 30000;       // keep result on screen 8 s after removal
      sensor.clearFIFO();
      break;
    }
  }
}
